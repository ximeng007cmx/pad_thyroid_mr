# -*- coding: utf-8 -*-
"""
Step 3.1  结局 SNP 提取 + Harmonise（甲减/TSH/FT4 → PAD）

输入
  work/instruments_{hypo,tsh,ft4}.tsv          Step 2 工具变量
  data/outcome/pad_GCST90018890.h.tsv.gz       PAD 结局（Sakaue 2021，rsid 在末列，无 n 列）

关键口径（Step 2 实测结论）
  1. 对齐一律用 rsID；位置仅作兜底（用 SNP 字段=面板坐标=正确 GRCh38，与 PAD 精确相等）。
  2. PAD 列序：chromosome|base_pair_location|effect_allele|other_allele|beta|standard_error|
     effect_allele_frequency|p_value|variant_id|hm_coordinate_conversion|hm_code|rsid（无 n 列）。
  3. 链方向（2026-10-07 实测）：PAD 与 1000G 正向参考一致；暴露 GCST90572789/90/91 等位标注
     混合链（variant_id/hm_code 证实）。互补对齐后 EAF 差中位 0.02-0.03 与同链组一致 →
     数据自洽，harmonise 时对非回文位点做链互补翻转（COMP 分支），无歧义。

输出
  work/pad_sub.tsv                        结局中命中的全部 SNP 行（一次性提取，复用）
  step3/checks/harmonise_<trait>.tsv      逐 SNP 决策台账
  step3/results/harmonised_<trait>.tsv    harmonise 后的 MR 输入（Step 3.2+ 用）
"""
import os, sys, csv, gzip, math, time, datetime

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
PAD  = os.path.join(ROOT, "data", "outcome", "pad_GCST90018890.h.tsv.gz")
WORK = os.path.join(ROOT, "work")
S3   = os.path.join(ROOT, "step3")
CHK  = os.path.join(S3, "checks")
RES  = os.path.join(S3, "results")
LOG  = os.path.join(S3, "logs")
for d in (CHK, RES, LOG):
    os.makedirs(d, exist_ok=True)
STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF  = os.path.join(LOG, "s3_harmonise_%s.log" % STAMP)

TRAITS = ["hypo", "tsh", "ft4"]

def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line)
    with open(LOGF, "a", encoding="utf-8") as f:
        f.write(line + "\n")

# ------------------------------------------------ 1) 读工具变量
inst = {}
for t in TRAITS:
    rows = list(csv.DictReader(open(os.path.join(WORK, "instruments_%s.tsv" % t),
                                    encoding="utf-8"), delimiter="\t"))
    for r in rows:
        r["F"] = float(r["F"])
        r["beta"] = float(r["beta"]); r["se"] = float(r["se"])
        r["eaf"] = float(r["EAF"]);  r["p"] = float(r["pval"])
    inst[t] = rows
    log("%s: 工具 %d 个" % (t, len(rows)))

# rsID -> (trait, row)：rs_id 优先，rsid 兜底
rs_map = {}
for t in TRAITS:
    for r in inst[t]:
        for key in (r["rs_id_gwas"], r["rsid_gwas"]):
            if key and key not in rs_map:
                rs_map[key] = (t, r)
# 位置兜底：SNP 字段 = 面板坐标 = 正确 GRCh38
pos_map = {}
for t in TRAITS:
    for r in inst[t]:
        pos_map.setdefault(r["SNP"], (t, r))

log("rsID 索引 %d 条；位置索引 %d 条" % (len(rs_map), len(pos_map)))

# ------------------------------------------------ 2) 扫描 PAD，提取命中行
sub_path = os.path.join(WORK, "pad_sub.tsv")
if os.path.exists(sub_path):
    log("pad_sub.tsv 已存在，跳过扫描（如需重提请删除该文件）")
else:
    t0 = time.time()
    hit_rs = hit_pos = 0
    n = 0
    out = open(sub_path + ".tmp", "w", encoding="utf-8", newline="")
    w = csv.writer(out, delimiter="\t")
    w.writerow(["rsid", "chr", "pos", "ea", "oa", "beta", "se", "eaf", "pval"])
    with gzip.open(PAD, "rt", encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        hdr = next(r)
        assert hdr[0] == "chromosome" and hdr[-1] == "rsid", hdr
        i_chr, i_pos, i_ea, i_oa = 0, 1, 2, 3
        i_beta, i_se, i_eaf, i_p = 4, 5, 6, 7
        i_rs = len(hdr) - 1
        for row in r:
            n += 1
            rsid = row[i_rs]
            got = None
            if rsid in rs_map:
                got = "rs"
            else:
                key = "%s:%s" % (row[i_chr], row[i_pos])
                if key in pos_map:
                    got = "pos"
            if got is None:
                continue
            try:
                w.writerow([rsid, row[i_chr], row[i_pos], row[i_ea].upper(), row[i_oa].upper(),
                            row[i_beta], row[i_se], row[i_eaf], row[i_p]])
            except ValueError:
                continue
            if got == "rs":
                hit_rs += 1
            else:
                hit_pos += 1
            if n % 4000000 == 0:
                log("  扫描 %dM 行 ..." % (n // 1000000))
    out.close()
    os.replace(sub_path + ".tmp", sub_path)
    log("PAD 扫描完成：%d 行，%.0f s；rsID 命中 %d，位置兜底命中 %d"
        % (n, time.time() - t0, hit_rs, hit_pos))

# 读回
pad_rows = {}
for d in csv.DictReader(open(sub_path, encoding="utf-8"), delimiter="\t"):
    if d["rsid"]:
        pad_rows.setdefault(d["rsid"], d)
    pad_rows.setdefault("POS:%s:%s" % (d["chr"], d["pos"]), d)
log("pad_sub.tsv 共 %d 个唯一 rsID" % len([k for k in pad_rows if not k.startswith("POS:")]))

# ------------------------------------------------ 3) Harmonise
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}   # 链互补（2026-10-07 实测：暴露文件混合链，PAD 正向）
HCOL = ["SNP", "chr", "pos_panel", "pos_gwas", "rs_id", "rsid_gwas", "match_route",
        "A1", "A2", "beta_X", "se_X", "eaf_X", "p_X", "F_X",
        "beta_Y", "se_Y", "eaf_Y", "p_Y",
        "harmonised", "action", "eaf_diff", "beta_XY_ratio", "se_XY_ratio", "note"]

def harmonise(t):
    rep = open(os.path.join(CHK, "harmonise_%s.tsv" % t), "w", encoding="utf-8", newline="")
    wr = csv.DictWriter(rep, fieldnames=HCOL, delimiter="\t", extrasaction="ignore")
    wr.writeheader()
    ok_rows = []
    stat = {"total": 0, "rs_id": 0, "rsid_gwas": 0, "pos": 0, "no_outcome": 0,
            "allele_mismatch": 0, "palindromic_ambiguous": 0, "se_bad": 0,
            "kept": 0, "eaf_discordant": 0, "flipped": 0, "complement": 0,
            "eaf_gap_015": 0}
    for r in inst[t]:
        stat["total"] += 1
        rec = {k: "" for k in HCOL}
        rec.update({"SNP": r["SNP"], "chr": r["chr"], "pos_panel": r["SNP"].split(":")[1],
                    "pos_gwas": r["pos"], "rs_id": r["rs_id_gwas"], "rsid_gwas": r["rsid_gwas"],
                    "beta_X": r["beta"], "se_X": r["se"], "eaf_X": r["eaf"], "p_X": r["p"],
                    "F_X": r["F"], "A1": r["EA"], "A2": r["OA"]})
        # --- 找结局
        y = None
        for key, route in ((r["rs_id_gwas"], "rs_id"), (r["rsid_gwas"], "rsid_gwas")):
            if key and key in pad_rows and not key.startswith("POS:"):
                y, rec["match_route"] = pad_rows[key], route
                stat[route] += 1
                break
        if y is None:
            y0 = pad_rows.get("POS:" + r["SNP"])
            if y0 is not None:
                y, rec["match_route"] = y0, "pos"
                stat["pos"] += 1
        if y is None:
            stat["no_outcome"] += 1
            rec.update({"harmonised": "FALSE", "action": "no_outcome_snp"})
            wr.writerow(rec); continue
        rec.update({"beta_Y": y["beta"], "se_Y": y["se"], "eaf_Y": y["eaf"], "p_Y": y["pval"]})
        try:
            bY, seY, eafY = float(y["beta"]), float(y["se"]), float(y["eaf"])
        except ValueError:
            stat["se_bad"] += 1
            rec.update({"harmonised": "FALSE", "action": "outcome_numeric_bad"})
            wr.writerow(rec); continue
        if not (seY > 0) or not math.isfinite(bY):
            stat["se_bad"] += 1
            rec.update({"harmonised": "FALSE", "action": "outcome_se_nonpositive"})
            wr.writerow(rec); continue
        e1, o1 = r["EA"], r["OA"]            # 暴露 EA/OA
        e2, o2 = y["ea"], y["oa"]            # 结局 EA/OA
        strand_tag = ""
        if {e1, o1} != {e2, o2}:
            # 链互补分支：结局翻到暴露的链（非回文时无歧义）
            # 依据 2026-10-07 验证：暴露文件等位标注混合链（variant_id/hm_code 证实），
            # 互补对齐后 eaf 差中位 0.02-0.03，与同链组一致 → 数据自洽，可安全翻转
            if e1 in COMP and o1 in COMP and {COMP[e1], COMP[o1]} == {e2, o2}:
                e2, o2 = COMP[e2], COMP[o2]
                strand_tag = "complement_"
            else:
                stat["allele_mismatch"] += 1
                rec.update({"harmonised": "FALSE", "action": "allele_mismatch",
                            "note": "X=%s/%s Y=%s/%s" % (e1, o1, y["ea"], y["oa"])})
                wr.writerow(rec); continue
        if strand_tag:
            stat["complement"] += 1
        action = strand_tag + "straight"
        if e2 == e1:
            pass                              # 已同向
        else:                                 # 效应等位相反 → 翻转结局
            bY, eafY = -bY, 1.0 - eafY
            action = strand_tag + "flipped"
            stat["flipped"] += 1
        pal = (e1 + o1) in ("AT", "TA", "CG", "GC")
        diff = abs(r["eaf"] - eafY)
        rec["eaf_diff"] = round(diff, 4)
        if pal and diff > 0.1:
            stat["palindromic_ambiguous"] += 1
            rec.update({"harmonised": "FALSE", "action": "palindromic_ambiguous"})
            wr.writerow(rec); continue
        if not pal and diff > 0.2:
            stat["eaf_discordant"] += 1
            rec.update({"harmonised": "FALSE", "action": "eaf_discordant"})
            wr.writerow(rec); continue
        seR = seY / abs(r["beta"])
        bR = bY / r["beta"]
        if diff > 0.15:
            stat["eaf_gap_015"] += 1
            rec["note"] = "eaf_gap>0.15 QC留意"
        rec.update({"harmonised": "TRUE", "action": action, "beta_XY_ratio": bR, "se_XY_ratio": seR})
        wr.writerow(rec)
        ok_rows.append({"SNP": r["SNP"], "rs_id": r["rs_id_gwas"], "chr": r["chr"],
                        "pos": r["SNP"].split(":")[1], "A1": e1, "A2": o1,
                        "beta_X": r["beta"], "se_X": r["se"], "eaf_X": r["eaf"], "p_X": r["p"],
                        "F_X": r["F"], "beta_Y": bY, "se_Y": seY, "eaf_Y": eafY, "p_Y": y["pval"],
                        "beta_XY": bR, "se_XY": seR, "palindromic": pal, "action": action})
        stat["kept"] += 1
    rep.close()
    outp = os.path.join(RES, "harmonised_%s.tsv" % t)
    with open(outp, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(ok_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(ok_rows)
    log("%s: 候选 %d → 保留 %d | 无结局 %d | 等位基因不符 %d | 链互补 %d | 回文歧义 %d | EAF 不符 %d | 翻转 %d | EAF差>0.15 %d"
        % (t, stat["total"], stat["kept"], stat["no_outcome"], stat["allele_mismatch"],
           stat["complement"], stat["palindromic_ambiguous"], stat["eaf_discordant"],
           stat["flipped"], stat["eaf_gap_015"]))
    return stat

allstat = {}
for t in TRAITS:
    allstat[t] = harmonise(t)

with open(os.path.join(RES, "harmonise_summary.json"), "w", encoding="utf-8") as f:
    import json
    json.dump(allstat, f, ensure_ascii=False, indent=2)
log("完成。日志: %s" % LOGF)
