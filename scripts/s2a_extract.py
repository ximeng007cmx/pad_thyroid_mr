# -*- coding: utf-8 -*-
"""
Step 2.1  从三套暴露 GWAS 的 harmonised 文件中提取候选工具变量（p < 5e-8）
纯标准库实现（gzip + csv），不依赖 pandas。

输入：
  data/exposure/hypo_GCST90572791.h.tsv.gz   （Hypothyroidism，二值）
  data/exposure/tsh_GCST90572789.h.tsv.gz    （TSH，连续）
  data/exposure/ft4_GCST90572790.h.tsv.gz    （FT4，连续）
输出：
  work/cand_<trait>.tsv      全部 p<5e-8 记录
  work/cand_<trait>.diag.txt 诊断报告
"""
import gzip, os, sys, csv, time, math
from collections import Counter

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
DATA = os.path.join(ROOT, "data", "exposure")
WORK = os.path.join(ROOT, "work")
os.makedirs(WORK, exist_ok=True)

PVAL_THRESH = 5e-8

TRAITS = [
    ("hypo", "hypo_GCST90572791.h.tsv.gz"),
    ("tsh",  "tsh_GCST90572789.h.tsv.gz"),
    ("ft4",  "ft4_GCST90572790.h.tsv.gz"),
]

# GWAS-SSF 列序（实测）
COL = {"chromosome": 0, "base_pair_location": 1, "effect_allele": 2, "other_allele": 3,
       "beta": 4, "standard_error": 5, "effect_allele_frequency": 6, "p_value": 7,
       "rsid": 8, "rs_id": 9, "n": 10, "hm_coordinate_conversion": 11, "hm_code": 12,
       "variant_id": 13}

HDR = ["chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "rsid", "rs_id", "n",
       "hm_conv", "hm_code", "variant_id", "z", "F"]


def run(tag, fname):
    path = os.path.join(DATA, fname)
    out = os.path.join(WORK, f"cand_{tag}.tsv")
    diag = os.path.join(WORK, f"cand_{tag}.diag.txt")

    t0 = time.time()
    stats = Counter()
    n_rows = 0
    ratio_rsid_used = Counter()      # rsid / rs_id 哪一列有值
    samples_rsid = []
    samples_rs_id = []

    fout = open(out, "w", newline="", encoding="utf-8")
    w = csv.writer(fout, delimiter="\t")
    w.writerow(HDR)

    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        assert header[COL["p_value"]] == "p_value", header
        for row in r:
            n_rows += 1
            if len(row) < 14:
                stats["短行"] += 1
                continue
            try:
                p = float(row[COL["p_value"]])
            except ValueError:
                stats["p值不可解析"] += 1
                continue
            if not (p < PVAL_THRESH):
                continue
            stats["p<5e-8"] += 1

            chrom = row[COL["chromosome"]]
            if chrom not in {str(i) for i in range(1, 23)}:
                stats["非常染色体"] += 1
                continue
            stats["常染色体"] += 1

            ea = row[COL["effect_allele"]].upper()
            oa = row[COL["other_allele"]].upper()
            if len(ea) != 1 or len(oa) != 1 or ea not in "ACGT" or oa not in "ACGT":
                stats["非SNP(含 indel/多碱基)"] += 1
                continue
            if ea == oa:
                stats["等位基因相同"] += 1
                continue
            stats["双等位SNP"] += 1

            rsid, rs_id = row[COL["rsid"]], row[COL["rs_id"]]
            if rsid and rsid != ".":
                ratio_rsid_used["rsid"] += 1
                if len(samples_rsid) < 3:
                    samples_rsid.append(rsid)
            elif rs_id and rs_id != ".":
                ratio_rsid_used["rs_id"] += 1
                if len(samples_rs_id) < 3:
                    samples_rs_id.append(rs_id)
            else:
                ratio_rsid_used["均缺失"] += 1

            try:
                beta = float(row[COL["beta"]])
                se = float(row[COL["standard_error"]])
                eaf = float(row[COL["effect_allele_frequency"]])
                n = float(row[COL["n"]])
            except ValueError:
                stats["数值字段不可解析"] += 1
                continue
            if se <= 0 or not math.isfinite(beta) or not math.isfinite(se):
                stats["se非正/非有限"] += 1
                continue
            stats["完整记录"] += 1

            z = beta / se
            F = z * z
            w.writerow([chrom, row[COL["base_pair_location"]], ea, oa,
                        f"{beta:.6g}", f"{se:.6g}", f"{eaf:.6g}", f"{p:.6g}",
                        rsid, rs_id, f"{n:.0f}", row[COL["hm_coordinate_conversion"]],
                        row[COL["hm_code"]], row[COL["variant_id"]],
                        f"{z:.6g}", f"{F:.6g}"])
    fout.close()

    with open(diag, "w", encoding="utf-8") as d:
        d.write(f"# 候选位点提取诊断 — {tag} ({fname})\n")
        d.write(f"扫描总行数            : {n_rows:,}\n")
        d.write(f"阈值                  : p < {PVAL_THRESH:g}\n")
        d.write(f"耗时                  : {time.time()-t0:.1f} s\n\n")
        for k, v in stats.items():
            d.write(f"{k:<22}: {v:,}\n")
        d.write(f"\nrsID 列取值分布       : {dict(ratio_rsid_used)}\n")
        d.write(f"rsid 样例             : {samples_rsid}\n")
        d.write(f"rs_id 样例            : {samples_rs_id}\n")

    print(f"[{tag}] 总行 {n_rows:,} | p<5e-8 {stats['p<5e-8']:,} | "
          f"常染色体 {stats['常染色体']:,} | 双等位SNP {stats['双等位SNP']:,} | "
          f"完整记录 {stats['完整记录']:,} | {time.time()-t0:.1f}s")
    return stats["完整记录"]


if __name__ == "__main__":
    print("=== Step 2.1 候选位点提取 ===")
    for tag, fn in TRAITS:
        run(tag, fn)
    print("=== 完成 ===")
