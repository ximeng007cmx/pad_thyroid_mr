# -*- coding: utf-8 -*-
"""
Step 4  FinnGen R13 I9_PAD 复制分析（甲减/TSH/FT4 → PAD）

输入
  work/instruments_{hypo,tsh,ft4}.tsv                     Step 2 工具变量
  data/replication/finngen_R13_I9_PAD.gz                  FinnGen R13 I9_PAD
     列：#chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
     效应等位 = alt；GRCh38（已验证 chr1:13668 与 PAD/Sakaue 文件一致）
  FinnGen R13 manifest：I9_PAD = 22,244 例 / 460,490 对照（用于 n_eff）

⚠ 重叠声明（见 deliverables/Step4_样本重叠影响评估.md）
  甲减暴露 GWAS 含 FinnGen DF10 → 甲减→FinnGen PAD 存在样本重叠，**其结果不能作为复制证据**，
  仅作参考（输出中标 OVERLAP_INVALID）。TSH/FT4 暴露不含芬兰 → 为独立复制（VALID）。

输出
  work/finngen_pad_sub.tsv                        命中的 FinnGen 行（一次性提取，复用）
  step4/checks/harmonise_finngen_<t>.tsv          逐 SNP 决策台账
  step4/results/harmonised_finngen_<t>.tsv        MR 输入
  step4/results/mr_finngen_<t>.json               全部统计量
  step4/results/mr_finngen_summary.tsv            三性状 × 方法汇总（含重叠标注）
"""
import os, sys, csv, gzip, math, time, datetime, json

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
sys.path.insert(0, os.path.join(ROOT, "step3", "scripts"))
from mrlib import (p_norm, p_chi2, ivw, q_stat, ivw_random_dl, egger,
                   weighted_median, weighted_mode, ivw_bootstrap_se,
                   presso_equiv, p_binom_two)

FIN   = os.path.join(ROOT, "data", "replication", "finngen_R13_I9_PAD.gz")
WORK  = os.path.join(ROOT, "work")
S4    = os.path.join(ROOT, "step4")
CHK   = os.path.join(S4, "checks"); RES = os.path.join(S4, "results"); LOG = os.path.join(S4, "logs")
for d in (CHK, RES, LOG): os.makedirs(d, exist_ok=True)
STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF  = os.path.join(LOG, "s4_finngen_%s.log" % STAMP)

TRAITS = ["hypo", "tsh", "ft4"]
TRAIT_CN = {"hypo": "甲减", "tsh": "TSH", "ft4": "FT4"}
# 重叠标注：甲减暴露含 FinnGen DF10 → FinnGen PAD 复制无效
OVERLAP = {"hypo": "OVERLAP_INVALID(暴露含FinnGen DF10)", "tsh": "VALID", "ft4": "VALID"}
# FinnGen R13 I9_PAD
FG_CASE, FG_CTRL = 22244, 460490
N_EFF_FG = 4.0 * FG_CASE * FG_CTRL / (FG_CASE + FG_CTRL)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def log(m):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), m)
    print(line)
    with open(LOGF, "a", encoding="utf-8") as f: f.write(line + "\n")

log("FinnGen R13 I9_PAD: %d 例 / %d 对照，n_eff=%.0f" % (FG_CASE, FG_CTRL, N_EFF_FG))
log("重叠标注: hypo=%s | tsh=%s | ft4=%s" % (OVERLAP["hypo"], OVERLAP["tsh"], OVERLAP["ft4"]))

# ---------------- 1) 读工具变量
inst = {}
for t in TRAITS:
    rows = list(csv.DictReader(open(os.path.join(WORK, "instruments_%s.tsv" % t),
                                    encoding="utf-8"), delimiter="\t"))
    for r in rows:
        r["F"] = float(r["F"]); r["beta"] = float(r["beta"]); r["se"] = float(r["se"])
        r["eaf"] = float(r["EAF"]); r["p"] = float(r["pval"])
        r["R2"] = float(r["R2"]); r["n_eff"] = float(r["n_eff"])
    inst[t] = rows
    log("%s: 工具 %d 个" % (t, len(rows)))

rs_map, pos_map = {}, {}
for t in TRAITS:
    for r in inst[t]:
        for k in (r["rs_id_gwas"], r["rsid_gwas"]):
            if k and k not in rs_map: rs_map[k] = True
        pos_map.setdefault("%s:%s" % (r["chr"], r["SNP"].split(":")[1]), True)
log("rsID 索引 %d；位置索引 %d" % (len(rs_map), len(pos_map)))

# ---------------- 2) 扫描 FinnGen，提取命中行
sub = os.path.join(WORK, "finngen_pad_sub.tsv")
if os.path.exists(sub):
    log("finngen_pad_sub.tsv 已存在，跳过扫描")
else:
    t0 = time.time(); n = hit_rs = hit_pos = 0
    out = open(sub + ".tmp", "w", encoding="utf-8", newline="")
    w = csv.writer(out, delimiter="\t")
    w.writerow(["rsid", "chr", "pos", "alt", "ref", "beta", "se", "af_alt", "pval", "gene"])
    with gzip.open(FIN, "rt", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        assert hdr[0].lstrip("#") == "chrom" and hdr[4] == "rsids", hdr
        for row in rd:
            n += 1
            rsids = [x for x in row[4].split(",") if x]
            key = None
            for x in rsids:
                if x in rs_map: key = x; break
            route = "rs"
            if key is None:
                pk = "%s:%s" % (row[0], row[1])
                if pk in pos_map: route = "pos"; key = pk
            if key is None: continue
            w.writerow([key, row[0], row[1], row[3].upper(), row[2].upper(),
                        row[8], row[9], row[10], row[6], row[5]])
            if route == "rs": hit_rs += 1
            else: hit_pos += 1
            if n % 4000000 == 0: log("  扫描 %dM 行 ..." % (n // 1000000))
    out.close(); os.replace(sub + ".tmp", sub)
    log("FinnGen 扫描完成：%d 行，%.0f s；rsID 命中 %d，位置兜底 %d" % (n, time.time() - t0, hit_rs, hit_pos))

pad = {}
for d in csv.DictReader(open(sub, encoding="utf-8"), delimiter="\t"):
    if not d["rsid"].startswith("POS:"): pad.setdefault(d["rsid"], d)
    pad.setdefault("POS:" + d["chr"] + ":" + d["pos"], d)
log("命中唯一 rsID %d 个" % len([k for k in pad if not k.startswith("POS:")]))

# ---------------- 3) Harmonise（与 Step 3 同一套判定规则）
HCOL = ["SNP", "chr", "pos_panel", "rs_id", "rsid_gwas", "match_route",
        "A1", "A2", "beta_X", "se_X", "eaf_X", "F_X",
        "beta_Y", "se_Y", "eaf_Y", "p_Y", "gene_Y",
        "harmonised", "action", "eaf_diff", "beta_XY", "se_XY", "note"]

def harmonise(t):
    rep = open(os.path.join(CHK, "harmonise_finngen_%s.tsv" % t), "w", encoding="utf-8", newline="")
    wr = csv.DictWriter(rep, fieldnames=HCOL, delimiter="\t", extrasaction="ignore"); wr.writeheader()
    ok = []
    st = dict(total=0, kept=0, no_outcome=0, allele_mismatch=0, complement=0, flipped=0,
              palindromic_ambiguous=0, eaf_discordant=0, eaf_gap=0)
    for r in inst[t]:
        st["total"] += 1
        rec = {k: "" for k in HCOL}
        rec.update({"SNP": r["SNP"], "chr": r["chr"], "pos_panel": r["SNP"].split(":")[1],
                    "rs_id": r["rs_id_gwas"], "rsid_gwas": r["rsid_gwas"],
                    "beta_X": r["beta"], "se_X": r["se"], "eaf_X": r["eaf"], "F_X": r["F"],
                    "A1": r["EA"], "A2": r["OA"]})
        y = None
        for k, route in ((r["rs_id_gwas"], "rs_id"), (r["rsid_gwas"], "rsid_gwas")):
            if k and k in pad and not k.startswith("POS:"):
                y, rec["match_route"] = pad[k], route; break
        if y is None:
            y0 = pad.get("POS:" + r["SNP"])
            if y0 is not None: y, rec["match_route"] = y0, "pos"
        if y is None:
            st["no_outcome"] += 1
            rec.update({"harmonised": "FALSE", "action": "no_outcome_snp"}); wr.writerow(rec); continue
        rec.update({"beta_Y": y["beta"], "se_Y": y["se"], "eaf_Y": y["af_alt"],
                    "p_Y": y["pval"], "gene_Y": y["gene"]})
        try:
            bY, seY, eafY = float(y["beta"]), float(y["se"]), float(y["af_alt"])
        except ValueError:
            rec.update({"harmonised": "FALSE", "action": "outcome_numeric_bad"}); wr.writerow(rec); continue
        if not (seY > 0) or not math.isfinite(bY):
            rec.update({"harmonised": "FALSE", "action": "outcome_se_nonpositive"}); wr.writerow(rec); continue
        e1, o1 = r["EA"], r["OA"]; e2, o2 = y["alt"], y["ref"]
        tag = ""
        if {e1, o1} != {e2, o2}:
            if e1 in COMP and o1 in COMP and {COMP[e1], COMP[o1]} == {e2, o2}:
                e2, o2 = COMP[e2], COMP[o2]; tag = "complement_"
            else:
                st["allele_mismatch"] += 1
                rec.update({"harmonised": "FALSE", "action": "allele_mismatch",
                            "note": "X=%s/%s Y=%s/%s" % (e1, o1, y["alt"], y["ref"])})
                wr.writerow(rec); continue
        if tag: st["complement"] += 1
        action = tag + "straight"
        if e2 != e1:
            bY, eafY = -bY, 1.0 - eafY; action = tag + "flipped"; st["flipped"] += 1
        pal = (e1 + o1) in ("AT", "TA", "CG", "GC")
        diff = abs(r["eaf"] - eafY)
        rec["eaf_diff"] = round(diff, 4)
        if pal and diff > 0.1:
            st["palindromic_ambiguous"] += 1
            rec.update({"harmonised": "FALSE", "action": "palindromic_ambiguous"}); wr.writerow(rec); continue
        if not pal and diff > 0.2:
            st["eaf_discordant"] += 1
            rec.update({"harmonised": "FALSE", "action": "eaf_discordant"}); wr.writerow(rec); continue
        bR, seR = bY / r["beta"], seY / abs(r["beta"])
        if diff > 0.15:
            st["eaf_gap"] += 1; rec["note"] = "eaf_gap>0.15"
        rec.update({"harmonised": "TRUE", "action": action, "beta_XY": bR, "se_XY": seR})
        wr.writerow(rec)
        ok.append({"SNP": r["SNP"], "rs_id": r["rs_id_gwas"], "chr": r["chr"],
                   "A1": e1, "A2": o1, "beta_X": r["beta"], "se_X": r["se"], "eaf_X": r["eaf"],
                   "F_X": r["F"], "R2_X": r["R2"], "n_eff_X": r["n_eff"],
                   "beta_Y": bY, "se_Y": seY, "eaf_Y": eafY, "p_Y": y["pval"],
                   "beta_XY": bR, "se_XY": seR, "gene_Y": y["gene"], "action": action})
        st["kept"] += 1
    rep.close()
    with open(os.path.join(RES, "harmonised_finngen_%s.tsv" % t), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(ok[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(ok)
    log("%s: 候选 %d → 保留 %d | 无结局 %d | 等位不符 %d | 链互补 %d | 回文歧义 %d | EAF不符 %d | 翻转 %d | EAF差>0.15 %d"
        % (t, st["total"], st["kept"], st["no_outcome"], st["allele_mismatch"], st["complement"],
           st["palindromic_ambiguous"], st["eaf_discordant"], st["flipped"], st["eaf_gap"]))
    return ok, st

# ---------------- 4) MR
def direction(rows):
    r2x = sum(r["R2_X"] for r in rows); k = len(rows)
    r2y = sum((r["beta_Y"]/r["se_Y"])**2 / ((r["beta_Y"]/r["se_Y"])**2 + N_EFF_FG) for r in rows)
    nX = rows[0]["n_eff_X"]
    return {"R2_X": r2x, "R2_Y": r2y, "nR2_X": nX*r2x, "nR2_Y": N_EFF_FG*r2y,
            "direction": "forward (X->Y)" if nX*r2x > N_EFF_FG*r2y else "reverse?"}

summary = []
for t in TRAITS:
    rows, st = harmonise(t)
    bX = [r["beta_X"] for r in rows]; sX = [r["se_X"] for r in rows]
    bY = [r["beta_Y"] for r in rows]; sY = [r["se_Y"] for r in rows]
    bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
    k = len(rows)
    out = {"trait": t, "trait_cn": TRAIT_CN[t], "nsnp": k, "overlap_status": OVERLAP[t],
           "outcome": "FinnGen R13 I9_PAD (22,244 cases / 460,490 controls)", "harmonise": st}
    bF, seF = ivw(bR, sR)
    out["ivw_fixed"] = {"beta": bF, "se": seF, "p": p_norm(bF/seF)}
    rnd = ivw_random_dl(bR, sR); out["ivw_random"] = rnd
    eg = egger(bX, bY, sY); out["mr_egger"] = eg
    bm = weighted_median(bR, sR); sem = ivw_bootstrap_se(bR, sR, "median")
    out["weighted_median"] = {"beta": bm, "se": sem, "p": p_norm(bm/sem)}
    bmo = weighted_mode(bR, sR); semo = ivw_bootstrap_se(bR, sR, "mode", nb=500)
    out["weighted_mode"] = {"beta": bmo, "se": semo, "p": p_norm(bmo/semo)}
    out["presso_equiv"] = presso_equiv(bR, sR)
    out["directionality"] = direction(rows)
    npos = sum(1 for x in bR if x > 0)
    out["sign_test"] = {"n_positive": npos, "n_negative": k-npos,
                        "p_binomial_two_sided": p_binom_two(min(npos, k-npos), k)}
    # 留一法
    loo = []
    for i in range(k):
        keep = [j for j in range(k) if j != i]
        r2 = ivw_random_dl([bR[j] for j in keep], [sR[j] for j in keep])
        loo.append({"SNP": rows[i]["SNP"], "nsnp": k-1, "beta": r2["beta"], "se": r2["se"], "p": r2["p"]})
    with open(os.path.join(RES, "loo_finngen_%s.tsv" % t), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["SNP","nsnp","beta","se","p"], delimiter="\t"); w.writeheader(); w.writerows(loo)
    with open(os.path.join(RES, "mr_finngen_%s.json" % t), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    log("%s [%s]: IVW固定 β=%.4f (p=%.3g) | IVW随机 β=%.4f (p=%.3g) | Egger β=%.4f 截距p=%.3g | Q_p=%.3g"
        % (t, OVERLAP[t], bF, out["ivw_fixed"]["p"], rnd["beta"], rnd["p"],
           eg["beta"] if eg else float("nan"), eg["intercept_pval"] if eg else float("nan"), rnd["Q_pval"]))
    def row(method, d, note=""):
        summary.append({"trait": t, "trait_cn": TRAIT_CN[t], "overlap": OVERLAP[t], "method": method,
                        "nsnp": d.get("nsnp", k), "beta": round(d["beta"], 6), "se": round(d["se"], 6),
                        "or": round(math.exp(d["beta"]), 6), "p": "%.4g" % d["p"], "note": note})
    row("IVW(固定效应)", out["ivw_fixed"])
    row("IVW(随机效应 DL)", rnd, "Q=%s df=%d Q_p=%.3g tau2=%.3g" % ("%.2f" % rnd["Q"], rnd["Q_df"], rnd["Q_pval"], rnd["tau2"]))
    if eg: row("MR-Egger", eg, "intercept=%s p=%.3g" % (round(eg["intercept"],6), eg["intercept_pval"]))
    row("加权中位数", out["weighted_median"], "bootstrap se 2000")
    row("加权众数", out["weighted_mode"], "bootstrap se 500")
    pe = out["presso_equiv"]
    dn = ""
    if pe["distortion"]:
        dc = pe["distortion"]; dn = "离群 %d, 校正beta=%s p=%.3g, distortion p=%.3g" % (
            dc["n_outliers"], round(dc["beta_corrected"],6), dc["p_corrected"], dc["distortion_p"])
    summary.append({"trait": t, "trait_cn": TRAIT_CN[t], "overlap": OVERLAP[t],
                    "method": "MR-PRESSO等价", "nsnp": k, "beta": "", "se": "", "or": "",
                    "p": "%.4g" % pe["global_Q_pval"], "note": dn})
    dj = out["directionality"]
    summary.append({"trait": t, "trait_cn": TRAIT_CN[t], "overlap": OVERLAP[t],
                    "method": "方向性(R²对比)", "nsnp": k, "beta": "", "se": "", "or": "", "p": "",
                    "note": "R2X=%.4f R2Y=%.5f nR2X=%.0f nR2Y=%.0f → %s" % (dj["R2_X"], dj["R2_Y"], dj["nR2_X"], dj["nR2_Y"], dj["direction"])})
    stt = out["sign_test"]
    summary.append({"trait": t, "trait_cn": TRAIT_CN[t], "overlap": OVERLAP[t], "method": "符号一致性",
                    "nsnp": k, "beta": "", "se": "", "or": "", "p": "%.4g" % stt["p_binomial_two_sided"],
                    "note": "正 %d / 负 %d" % (stt["n_positive"], stt["n_negative"])})

with open(os.path.join(RES, "mr_finngen_summary.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["trait","trait_cn","overlap","method","nsnp","beta","se","or","p","note"], delimiter="\t")
    w.writeheader(); w.writerows(summary)
log("完成。汇总: %s" % os.path.join(RES, "mr_finngen_summary.tsv"))
