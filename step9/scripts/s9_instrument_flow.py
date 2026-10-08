# -*- coding: utf-8 -*-
"""P2-10a: 复核工具筛选各阶段计数（对应 STROBE-MR 10a 流程图）。
直接复用 s2_clump.py 的匹配逻辑（import），不重写、不改数。
输出：step9/results/instrument_selection_flow.tsv
"""
import os, sys, io, csv, importlib.util, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"

# 动态载入 s2_clump.py（复用 read_candidates / build_panel_maps / match_candidates）
spec = importlib.util.spec_from_file_location("s2", os.path.join(ROOT, "scripts", "s2_clump.py"))
s2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s2)

rs2id, pos2ids = s2.build_panel_maps()

DIAG = {  # 来自 work/cand_<t>.diag.txt（实测）
    "hypo": {"scanned": 7344354, "p5e8": 16247, "autosomal": 16247, "biallelic": 12782, "non_snp_excluded": 3465},
    "tsh":  {"scanned": 6151997, "p5e8": 17246, "autosomal": 17246, "biallelic": 13731, "non_snp_excluded": 3515},
    "ft4":  {"scanned": 5901316, "p5e8": 1690,  "autosomal": 1690,  "biallelic": 1474,  "non_snp_excluded": 216},
}
HARM = {  # 来自 step3/results/harmonise_summary.json
    "hypo": {"total": 222, "kept": 221, "eaf_discordant": 1, "no_outcome": 0},
    "tsh":  {"total": 258, "kept": 256, "eaf_discordant": 0, "no_outcome": 2},
    "ft4":  {"total": 65,  "kept": 65,  "eaf_discordant": 0, "no_outcome": 0},
}
NAME = {"hypo":"甲减 (hypothyroidism)", "tsh":"TSH", "ft4":"FT4"}

rows = []
for tag in ("hypo", "tsh", "ft4"):
    cand, by_pos = s2.read_candidates(tag)
    n_rows = len(cand)
    n_keys = len(by_pos)                      # chr:pos 去重后的候选数
    dup_by_pos = n_rows - n_keys               # 因 chr:pos 重复被覆盖
    matched, mstats = s2.match_candidates(tag, by_pos, rs2id, pos2ids)
    n_matched = len(matched)
    d = DIAG[tag]; h = HARM[tag]
    rows.append({
        "trait": tag, "trait_cn": NAME[tag],
        "scanned": d["scanned"],
        "p_lt_5e8": d["p5e8"],
        "biallelic_snp": d["biallelic"],
        "excluded_non_snp": d["non_snp_excluded"],
        "cand_rows": n_rows, "cand_unique_pos": n_keys, "cand_dup_by_pos": dup_by_pos,
        "mapped_to_panel": n_matched,
        "route_rs_id": mstats["rs_id"], "route_rsid": mstats["rsid"], "route_pos_fallback": mstats["pos_shift"],
        "excluded_unmapped": d["biallelic"] - n_matched,
        "excluded_unmapped_nomatch": mstats["no_match"],
        "excluded_unmapped_dup": (d["biallelic"] - n_matched) - mstats["no_match"],
        "after_clumping": h["total"],
        "excluded_ld": n_matched - h["total"],
        "after_harmonise": h["kept"],
        "excluded_harmonise": h["total"] - h["kept"],
        "excl_harmonise_eaf": h["eaf_discordant"],
        "excl_harmonise_missing_outcome": h["no_outcome"],
    })

out = os.path.join(ROOT, "step9", "results", "instrument_selection_flow.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
cols = list(rows[0].keys())
with open(out, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
    for r in rows: w.writerow(r)

print(f"saved: {out}\n")
for r in rows:
    print("="*72)
    print(r["trait_cn"])
    print(f"  扫描总变异            {r['scanned']:,}")
    print(f"  p < 5e-8              {r['p_lt_5e8']:,}")
    print(f"  常染色体双等位 SNP     {r['biallelic_snp']:,}   (排除非 SNP/indel: {r['excluded_non_snp']:,})")
    print(f"    候选表实际行数       {r['cand_rows']:,}   (chr:pos 重复覆盖: {r['cand_dup_by_pos']:,})")
    print(f"  对齐参考面板          {r['mapped_to_panel']:,}   (rs_id {r['route_rs_id']:,} + rsid {r['route_rsid']} + 位置兜底 {r['route_pos_fallback']})")
    print(f"    排除未对齐           {r['excluded_unmapped']:,}   (面板无匹配/位置歧义 {r['excluded_unmapped_nomatch']:,} + 重复命中 {r['excluded_unmapped_dup']})")
    print(f"  LD clumping 后独立工具 {r['after_clumping']:,}   (LD 归并排除 {r['excluded_ld']:,})")
    print(f"  harmonise 后保留       {r['after_harmonise']:,}   (排除 {r['excluded_harmonise']}: EAF 不符 {r['excl_harmonise_eaf']} / 结局缺失 {r['excl_harmonise_missing_outcome']})")
