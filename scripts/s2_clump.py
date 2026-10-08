# -*- coding: utf-8 -*-
"""
Step 2 主流程：LD clumping → 工具变量清单 → F / R² 统计

流程
  A. 校验面板坐标版本（hg38 抽查 rs7412 / rs429358 等已知位点）
  B. 把面板 .bim 的 SNP ID 统一改写成 "chr:pos"（原文件备份为 *.rsid.orig）
     —— 这样只用位置匹配，彻底绕开 rsID 版本差异与重名问题
  C. 由三套候选位点生成 PLINK clump 输入（chr:pos + p）
  D. 调 PLINK 1.9 --clump：r² < 0.001，窗口 10 000 kb
  E. 解析索引 SNP，回填 beta/se/eaf/n，算 z、F、R²、整体 F
  F. 输出 TSV / XLSX / Markdown

依赖：PLINK 1.9（tools/plink.exe）、openpyxl、标准库。**不需要 numpy / pandas。**

用法：
  python s2_clump.py --check-build     # 只抽查坐标版本
  python s2_clump.py                   # 全流程
"""
import os, sys, re, csv, json, math, shutil, subprocess, argparse, time

ROOT   = r"D:/WorkSpace_Study/pad_thyroid_mr"
PANEL  = os.path.join(ROOT, "ref", "1kg_eur_hg38", "allchr.EUR.biallelicsnps.cbb")
PLINK  = os.path.join(ROOT, "tools", "plink.exe")
WORK   = os.path.join(ROOT, "work")
DELIV  = os.path.join(ROOT, "deliverables")
os.makedirs(WORK, exist_ok=True)
os.makedirs(DELIV, exist_ok=True)

TRAITS = [
    # tag, 中文名, 类型(quant/binary), 样本信息, 候选文件
    ("hypo", "甲状腺功能减退症（Hypothyroidism）", "binary",  dict(ncase=113393, nctrl=1065268)),
    ("tsh",  "促甲状腺激素（TSH）",               "quant",   dict(n=482873)),
    ("ft4",  "游离甲状腺素（FT4 / Thyroxine）",    "quant",   dict(n=191449)),
]

# 面板坐标版本抽查：hg38 下的期望位置
# rs7412 / rs429358 = APOE（hg38 经典坐标）；rs1801133 = MTHFR，
# 已用 Ensembl REST 核实其 GRCh38 位置为 1:11796321（11856378 是 GRCh37 位置，勿混淆）
BUILD_PROBES = {
    "rs7412":   ("19", 44908822),   # APOE
    "rs429358": ("19", 44908684),   # APOE
    "rs1801133":("1",  11796321),   # MTHFR（Ensembl REST 核实，GRCh38）
}


# ---------------------------------------------------------------- 工具
def read_candidates(tag):
    """读 work/cand_<tag>.tsv -> list[dict]，并用 chr:pos 建索引。"""
    path = os.path.join(WORK, "cand_%s.tsv" % tag)
    rows, by_pos = [], {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            d["key"] = "%s:%s" % (d["chr"], d["pos"])
            for k in ("beta", "se", "eaf", "p", "n", "z", "F"):
                d[k] = float(d[k])
            rows.append(d)
            by_pos[d["key"]] = d
    return rows, by_pos


def check_build():
    """抽查面板坐标版本。rsID 从 .rsid.orig（或 ID 映射表）里找——
    主 .bim 的 ID 已改写成 chr:pos，不再含 rsID。"""
    bim = PANEL + ".bim"
    bak = bim + ".rsid.orig"
    idmap = os.path.join(WORK, "panel_id_map.tsv")
    src = bak if os.path.exists(bak) else bim
    found = {}
    if os.path.exists(idmap):
        with open(idmap, "r", encoding="utf-8") as f:
            for d in csv.DictReader(f, delimiter="\t"):
                if d["rs_id"] in BUILD_PROBES:
                    found[d["rs_id"]] = (d["chr"], int(d["pos"]))
    elif os.path.exists(src):
        with open(src, "r", encoding="utf-8") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    p = line.split()
                if len(p) < 6:
                    continue
                if p[1] in BUILD_PROBES:
                    found[p[1]] = (p[0], int(p[3]))
    else:
        print("[!] 找不到 %s 或 %s" % (bim, src))
        return False
    print("=== 面板坐标版本抽查 ===")
    ok = True
    for rsid, (ec, ep) in BUILD_PROBES.items():
        if rsid not in found:
            print("  %-10s 面板中未找到" % rsid)
            ok = False
            continue
        gc, gp = found[rsid]
        hit = (str(gc) == ec and gp == ep)
        print("  %-10s 面板=%s:%d   期望hg38=%s:%d   %s"
              % (rsid, gc, gp, ec, ep, "✅" if hit else "❌"))
        ok = ok and hit
    print("结论：", "面板确认为 hg38 ✅" if ok else "坐标版本存疑，需人工确认 ❌")
    return ok


def rewrite_bim_ids():
    """把 .bim 的 SNP ID 改成 chr:pos，原文件备份。
    同时落盘 work/panel_id_map.tsv（idx/rsid/chr/pos/a1/a2/id），供 rsID 匹配用。"""
    bim = PANEL + ".bim"
    bak = bim + ".rsid.orig"
    idmap = os.path.join(WORK, "panel_id_map.tsv")
    if os.path.exists(bak) and os.path.exists(idmap):
        print("[bim] 已改写为 chr:pos（原文件在 %s，ID 映射表已存在）" % os.path.basename(bak))
        return True
    if not os.path.exists(bak):
        shutil.copy2(bim, bak)          # 此时 bim 仍带 rsID
    # 一律从 bak（含 rsID 的原始版本）重建
    n = dup = 0
    seen_pos = {}
    rows = []
    with open(bak, "r", encoding="utf-8") as fi:
        for idx, line in enumerate(fi):
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                p = line.split()
                if len(p) < 6:
                    continue
            chrom, rsid, pos, a1, a2 = p[0], p[1], p[3], p[4], p[5]
            key = (chrom, pos)
            if key in seen_pos:
                dup += 1
                sid = "%s:%s:dup%d" % (chrom, pos, dup)
            else:
                seen_pos[key] = True
                sid = "%s:%s" % (chrom, pos)
            rows.append((idx, rsid, chrom, pos, a1, a2, sid))
            n += 1
    with open(bim, "w", encoding="utf-8", newline="\n") as fo:
        fo.write("\n".join("\t".join((r[2], r[6], "0", r[3], r[4], r[5])) for r in rows) + "\n")
    with open(idmap, "w", encoding="utf-8", newline="\n") as fm:
        fm.write("idx\trs_id\tchr\tpos\ta1\ta2\tid\n")
        for r in rows:
            fm.write("\t".join(str(x) for x in r) + "\n")
    print("[bim] 改写完成：%d 个变异（%d 个同位置重复已加后缀），ID 映射表 -> %s" % (n, dup, idmap))
    return True


def build_panel_maps():
    """返回 rs2id（rsID -> 面板 SNP ID）与 pos2ids（chr:pos -> [面板 SNP ID]）。"""
    idmap = os.path.join(WORK, "panel_id_map.tsv")
    rs2id, pos2ids = {}, {}
    with open(idmap, "r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            rs2id.setdefault(d["rs_id"], d["id"])
            pos2ids.setdefault("%s:%s" % (d["chr"], d["pos"]), []).append(d["id"])
    return rs2id, pos2ids


def match_candidates(tag, by_pos, rs2id, pos2ids):
    """把候选位点对齐到面板。
    主路由 = rsID（rs_id 优先，其次 rsid）；兜底 = 位置匹配 pos-1
    （已三方验证暴露文件坐标整体 +1），且要求该位置在面板唯一。
    同一面板位点被多个候选命中时保留 p 最小者。
    返回 (matched: 面板SNP ID -> 候选行, stats)。"""
    matched = {}
    stats = {"rs_id": 0, "rsid": 0, "pos_shift": 0, "no_match": 0, "dup_panel": 0}
    cands = sorted(by_pos.values(), key=lambda d: d["p"])   # p 升序 → 先到先得
    for d in cands:
        pid, route = None, None
        for key, rt in ((d["rs_id"], "rs_id"), (d["rsid"], "rsid")):
            if key and key in rs2id:
                pid, route = rs2id[key], rt
                break
        if pid is None:
            k = "%s:%d" % (d["chr"], int(d["pos"]) - 1)
            ids = pos2ids.get(k)
            if ids and len(ids) == 1:
                pid, route = ids[0], "pos_shift"
            elif ids:
                stats["dup_panel"] += 1
        if pid is None:
            stats["no_match"] += 1
            continue
        if pid in matched:
            stats["dup_panel"] += 1
            continue
        stats[route] += 1
        d["match_route"] = route
        matched[pid] = d
    return matched, stats


def make_clump_input(tag, matched):
    path = os.path.join(WORK, "clump_in_%s.txt" % tag)
    keys = sorted(matched.keys(), key=lambda k: (int(k.split(":")[0]), int(k.split(":")[1])))
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("SNP\tP\n")
        for k in keys:
            f.write("%s\t%.6g\n" % (k, matched[k]["p"]))
    return path


def run_plink(tag, clump_in):
    out = os.path.join(WORK, "clump_%s" % tag)
    cmd = [PLINK, "--bfile", PANEL, "--clump", clump_in,
           "--clump-p1", "5e-8", "--clump-p2", "1",
           "--clump-r2", "0.001", "--clump-kb", "10000",
           "--clump-snp-field", "SNP", "--clump-field", "P",
           "--allow-no-sex", "--out", out]
    print("[plink] %s" % " ".join(cmd[1:]))
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(PLINK))
    print(p.stdout[-2500:])
    if p.returncode != 0:
        print(p.stderr[-1500:])
        raise RuntimeError("plink 退出码 %s" % p.returncode)
    print("[plink] 用时 %.1f s" % (time.time() - t0))
    return out + ".clumped"


def parse_clumped(path):
    """返回 (索引SNP列表, {索引SNP: [被剔除SNP]})"""
    index, members = [], {}
    if not os.path.exists(path):
        return index, members
    with open(path, "r", encoding="utf-8") as f:
        hdr = None
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                p = line.split()
            if len(p) < 3:
                continue
            if p[0].startswith("CHR") or p[0] == "CHR":
                hdr = p
                continue
            snp = p[2]
            index.append(snp)
            sp2 = p[-1] if len(p) >= 12 else ""
            mem = []
            if sp2 and sp2 != "NONE":
                for item in sp2.split(","):
                    m = re.match(r"^(.+?)\(([\d.eE+-]+)\)$", item.strip())
                    if m:
                        mem.append((m.group(1), float(m.group(2))))
            members[snp] = mem
    return index, members


# ---------------------------------------------------------------- 统计
def r2_from_z(z, neff):
    """TwoSampleMR get_r_from_pn 的代数等价式：R² = z² / (z² + n - 2)。
    二值性状用 n_eff = 4·n_case·n_ctrl/(n_case+n_ctrl) 代入，属观测尺度近似（已在报告中标注）。"""
    return (z * z) / (z * z + neff - 2.0)


def build_instruments(tag, meta, matched, index, members):
    rows = []
    for snp in index:
        d = matched.get(snp)
        if d is None:          # PLINK 报告了索引 SNP 但不在我们的映射里（不应发生）
            continue
        z = d["z"]
        F = z * z
        neff = (4.0 * meta["ncase"] * meta["nctrl"] / (meta["ncase"] + meta["nctrl"])
                if tag == "hypo" else meta["n"])
        r2 = r2_from_z(z, neff)
        maf = min(d["eaf"], 1 - d["eaf"])
        rows.append({
            "SNP": snp, "chr": d["chr"], "pos": d["pos"],
            "EA": d["ea"], "OA": d["oa"],
            "EAF": d["eaf"], "MAF": maf,
            "beta": d["beta"], "se": d["se"], "pval": d["p"],
            "z": z, "F": F, "R2": r2,
            "n": d["n"], "n_eff": neff,
            "rsid_gwas": d["rsid"], "rs_id_gwas": d["rs_id"],
            "match_route": d.get("match_route", ""),
            "n_clumped_out": len(members.get(snp, [])),
            "clumped_out": ";".join("%s(%.3f)" % t for t in members.get(snp, [])[:40]),
            "palindromic": (d["ea"] + d["oa"]) in ("AT", "TA", "CG", "GC"),
            "low_maf": maf < 0.01,
        })
    rows.sort(key=lambda r: r["pval"])
    return rows


def summarize(tag, meta, rows):
    k = len(rows)
    r2sum = sum(r["R2"] for r in rows)
    neff = rows[0]["n_eff"] if rows else 0
    Fs = [r["F"] for r in rows]
    meanF = sum(Fs) / k if k else 0
    overallF = (r2sum / (1 - r2sum)) * ((neff - k - 1) / k) if k and r2sum < 1 else float("nan")
    return {
        "trait": tag, "n_instruments": k,
        "r2_sum": r2sum, "r2_pct": r2sum * 100,
        "mean_F": meanF, "median_F": sorted(Fs)[k // 2] if k else 0,
        "min_F": min(Fs) if k else 0,
        "n_F_lt_10": sum(1 for x in Fs if x < 10),
        "overall_F": overallF,
        "n_palindromic": sum(1 for r in rows if r["palindromic"]),
        "n_low_maf": sum(1 for r in rows if r["low_maf"]),
        "neff": neff,
    }


# ---------------------------------------------------------------- 输出
COLS = ["SNP", "chr", "pos", "EA", "OA", "EAF", "MAF", "beta", "se", "pval",
        "z", "F", "R2", "n", "n_eff", "rsid_gwas", "rs_id_gwas", "match_route",
        "n_clumped_out", "palindromic", "low_maf", "clumped_out"]


def write_tsv(tag, rows):
    path = os.path.join(WORK, "instruments_%s.tsv" % tag)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def safe_sheet(name):
    """Excel 表名不允许 [ ] : * ? / \\ 且 ≤31 字符"""
    for ch in r'[]:*?/\\':
        name = name.replace(ch, "-")
    return (name[:28] or "sheet")


def write_xlsx(allrows, sums):
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    hdr_fill = PatternFill("solid", fgColor="1F4E79")
    hdr_font = Font(color="FFFFFF", bold=True)

    # 汇总表
    ws = wb.active
    ws.title = "汇总"
    ws.append(["性状", "工具数", "累计 R²(%)", "平均 F", "中位 F", "最小 F",
               "F<10 个数", "整体 F", "回文位点", "MAF<0.01", "有效样本量"])
    for c in ws[1]:
        c.fill, c.font = hdr_fill, hdr_font
    for s in sums:
        name = dict((t[0], t[1]) for t in TRAITS)[s["trait"]]
        ws.append([name, s["n_instruments"], round(s["r2_pct"], 3),
                   round(s["mean_F"], 1), round(s["median_F"], 1), round(s["min_F"], 1),
                   s["n_F_lt_10"], round(s["overall_F"], 1) if s["overall_F"] == s["overall_F"] else "NA",
                   s["n_palindromic"], s["n_low_maf"], int(s["neff"])])
    ws.freeze_panes = "A2"
    for i, w in enumerate([30, 10, 12, 10, 10, 10, 10, 12, 10, 10, 14], 1):
        ws.column_dimensions[get_column_letter(i)].width = w

    for tag, rows in allrows.items():
        name = dict((t[0], t[1]) for t in TRAITS)[tag]
        ws = wb.create_sheet(safe_sheet(name))
        ws.append(COLS)
        for c in ws[1]:
            c.fill, c.font = hdr_fill, hdr_font
            c.alignment = Alignment(horizontal="center")
        for r in rows:
            ws.append([r.get(c) for c in COLS])
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        widths = {"SNP": 16, "clumped_out": 46, "rsid_gwas": 16, "rs_id_gwas": 16}
        for i, c in enumerate(COLS, 1):
            ws.column_dimensions[get_column_letter(i)].width = widths.get(c, 11)
        for row in ws.iter_rows(min_row=2):
            if row[18].value is True:      # palindromic
                for c in row:
                    c.fill = PatternFill("solid", fgColor="FFF2CC")
            if row[19].value is True:      # low maf
                for c in row:
                    c.fill = PatternFill("solid", fgColor="FCE4D6")
    out = os.path.join(DELIV, "Step2_工具变量清单.xlsx")
    wb.save(out)
    return out


def write_md(allrows, sums):
    name = dict((t[0], t[1]) for t in TRAITS)
    L = ["# Step 2 工具变量筛选结果", "",
         "> 面板：1000G Phase 3 hg38 EUR（Zenodo 19068199，294 例；CC-BY-4.0）",
         "> clumping：PLINK v1.90b7，r² < 0.001，窗口 10 000 kb",
         "> 候选阈值：p < 5×10⁻⁸（常染色体、双等位 SNP）", "",
         "## 汇总", "",
         "| 性状 | 工具数 | 累计 R²(%) | 平均 F | 中位 F | 最小 F | F<10 | 整体 F | 回文 | MAF<0.01 | n_eff |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for s in sums:
        L.append("| %s | %d | %.3f | %.1f | %.1f | %.1f | %d | %s | %d | %d | %d |" % (
            name[s["trait"]], s["n_instruments"], s["r2_pct"], s["mean_F"], s["median_F"],
            s["min_F"], s["n_F_lt_10"],
            ("%.1f" % s["overall_F"]) if s["overall_F"] == s["overall_F"] else "NA",
            s["n_palindromic"], s["n_low_maf"], int(s["neff"])))
    L += ["", "## 口径说明", "",
          "- **F = z² = (beta/se)²**，逐位点精确值，无分布假设。MR 惯例要求 F > 10。",
          "- **R²** 用 TwoSampleMR `get_r_from_pn` 的代数等价式 `R² = z²/(z² + n − 2)`；",
          "  二值性状（甲减）以 `n_eff = 4·n_case·n_ctrl/(n_case+n_ctrl)` 代入，属**观测尺度近似**，",
          "  用于相对比较与整体 F 计算，不等同于责任量表（liability scale）R²。",
          "- **整体 F** = (R²/(1−R²))·((n_eff − k − 1)/k)，k 为工具数。",
          "- 回文位点（A/T、C/G）与 MAF<0.01 位点**暂未剔除**，已在清单中标记，留待 Step 3 协调时处理。",
          "- 样本重叠提示：暴露（Nat Genet 2025）与结局 GCST90018890（Sakaue 2021）均含 UKB，Step 3 需评估。",
          ""]
    for tag, rows in allrows.items():
        L += ["## %s —— 前 15 个工具（共 %d）" % (name[tag], len(rows)), "",
              "| SNP(chr:pos) | EA/OA | EAF | beta | se | p | z | F | R² | 剔除数 |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for r in rows[:15]:
            L.append("| %s | %s/%s | %.4f | %.4f | %.4f | %.2e | %.2f | %.1f | %.5f | %d |" % (
                r["SNP"], r["EA"], r["OA"], r["EAF"], r["beta"], r["se"], r["pval"],
                r["z"], r["F"], r["R2"], r["n_clumped_out"]))
        L.append("")
    out = os.path.join(DELIV, "Step2_工具变量清单.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-build", action="store_true")
    ap.add_argument("--skip-build-check", action="store_true")
    a = ap.parse_args()

    if a.check_build:
        sys.exit(0 if check_build() else 1)

    if not os.path.exists(PANEL + ".bed") or not os.path.exists(PANEL + ".bim"):
        print("[!] 面板尚未下载完整，请先跑 s2_download_panel.py")
        sys.exit(2)
    if not os.path.exists(PLINK):
        print("[!] 找不到 PLINK：%s" % PLINK)
        sys.exit(2)

    if not a.skip_build_check:
        check_build()
    rewrite_bim_ids()
    rs2id, pos2ids = build_panel_maps()
    print("[面板] rsID 映射 %d 条，位置 %d 个" % (len(rs2id), len(pos2ids)))

    allrows, sums, mstats_all = {}, [], {}
    for tag, cn, kind, meta in TRAITS:
        print("\n===== %s =====" % cn)
        cand, by_pos = read_candidates(tag)
        print("[候选] %d 个（p<5e-8，常染色体双等位 SNP）" % len(cand))
        matched, mstats = match_candidates(tag, by_pos, rs2id, pos2ids)
        mstats_all[tag] = mstats
        print("[对齐] rs_id=%d  rsid=%d  位置兜底=%d  无匹配=%d  面板重复=%d"
              % (mstats["rs_id"], mstats["rsid"], mstats["pos_shift"],
                 mstats["no_match"], mstats["dup_panel"]))
        cin = make_clump_input(tag, matched)
        clumped = run_plink(tag, cin)
        index, members = parse_clumped(clumped)
        print("[clump] 独立工具 %d 个" % len(index))
        rows = build_instruments(tag, meta, matched, index, members)
        s = summarize(tag, meta, rows)
        print("[统计] 平均F=%.1f 最小F=%.1f 累计R²=%.2f%% 整体F=%.1f"
              % (s["mean_F"], s["min_F"], s["r2_pct"], s["overall_F"]))
        write_tsv(tag, rows)
        allrows[tag], _ = rows, sums.append(s)

    x = write_xlsx(allrows, sums)
    m = write_md(allrows, sums)
    with open(os.path.join(WORK, "step2_summary.json"), "w", encoding="utf-8") as f:
        json.dump({"summary": sums, "match_stats": mstats_all}, f, ensure_ascii=False, indent=2)
    print("\n交付：\n  %s\n  %s" % (x, m))


if __name__ == "__main__":
    main()
