# -*- coding: utf-8 -*-
"""
Step 8 前置验证：共定位候选区域的**变异覆盖度**检查（决定 coloc 可行性）。

为什么先做这一步
  coloc/coloc.susie 需要区域内「全部变异」的完整汇总统计（不只是显著位点）。
  甲状腺文件（GWAS Catalog 新格式，Rand 2025）位点密度偏低，因此先量化：
    ① 区域内各自有多少变异；
    ② 两两可配对（坐标可对齐）的变异有多少；
    ③ 主效 lead SNP 本身是否在文件里。

坐标口径（沿用 Step 6 实测结论）
  · 面板/instruments 的 SNP ID 位置  = 真实 GRCh38 1-based
  · 甲状腺文件：SNP 位置 = 真实 + 1，indel 位置 = 真实 + 0
  · 血管结局文件（SSF）：hm_pos = 真实
  计数时：甲状腺按 file_pos ∈ [start, end+1]，结局按 file_pos ∈ [start, end]；
  求交集时甲状腺用 {p-1, p} 双候选映射回真实坐标。

用法：
  python step8/scripts/s8_precheck.py                # 默认 top10 位点
  python step8/scripts/s8_precheck.py --top 15 --win 1000000
"""
import os
import io
import sys
import csv
import gzip
import json
import time
import argparse
import datetime

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
WORK = os.path.join(ROOT, "work")
RES = os.path.join(ROOT, "step8", "results")
LOGD = os.path.join(ROOT, "step8", "logs")
os.makedirs(RES, exist_ok=True)
os.makedirs(LOGD, exist_ok=True)

LOG = None


def log(m):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(line, flush=True)
    if LOG:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")


# ---------------------------------------------------------------- 数据集定义
THY_COLS = dict(chr="chromosome", pos="base_pair_location", ea="effect_allele",
                oa="other_allele", beta="beta", se="standard_error",
                eaf="effect_allele_frequency", p="p_value", rsid="rsid", rs_id="rs_id")
SSF_COLS = dict(chr="hm_chrom", pos="hm_pos", rsid="hm_rsid", ea="hm_effect_allele",
                oa="hm_other_allele", beta="hm_beta", se="standard_error",
                eaf="hm_effect_allele_frequency", p="p_value")
MVP_COLS = dict(chr="Chromosome", pos="Position", rsid="SNP_ID", ea="EA",
                oa="Allele2", beta="Effect", se="SE", eaf="EAF", p="PValue")
SAK_COLS = dict(chr="chromosome", pos="base_pair_location", rsid="rsid",
                ea="effect_allele", oa="other_allele", beta="beta",
                se="standard_error", eaf="effect_allele_frequency", p="p_value")

DATASETS = [
    # tag, 显示名, 路径, 列定义, 是否有双 rsid 列, 坐标偏移候选(甲状腺=0/1)
    ("thy_hypo", "甲状腺·甲减 (Rand2025)", "data/exposure/hypo_GCST90572791.h.tsv.gz", THY_COLS, True, (0, 1)),
    ("thy_tsh", "甲状腺·TSH (Rand2025)", "data/exposure/tsh_GCST90572789.h.tsv.gz", THY_COLS, True, (0, 1)),
    ("thy_ft4", "甲状腺·FT4 (Rand2025)", "data/exposure/ft4_GCST90572790.h.tsv.gz", THY_COLS, True, (0, 1)),
    ("cad_nouk", "CAD (Nikpay2015; GCST003116)", "data/panvascular/cad_GCST003116.h.tsv.gz", SSF_COLS, False, (0,)),
    ("is_sv", "IS-SV (MEGASTROKE; GCST005841)", "data/panvascular/is_sv_GCST005841.h.tsv.gz", SSF_COLS, False, (0,)),
    ("is_laa", "IS-LAA (MEGASTROKE; GCST005840)", "data/panvascular/is_laa_GCST005840.h.tsv.gz", SSF_COLS, False, (0,)),
    ("mvp_meta", "MVP PAD (trans-ancestry)", "raw/mvp_pad/MVP.te.PAD.dbGAP.txt.gz", MVP_COLS, False, (0,)),
]


def _open(path):
    p = os.path.join(ROOT, path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def norm_chr(c):
    c = (c or "").strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c.upper()


# ---------------------------------------------------------------- 1 选候选位点
def pick_candidates(top, gap):
    p = os.path.join(WORK, "instruments_hypo.tsv")
    rows = []
    with io.open(p, encoding="utf-8") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            snp = d["SNP"]
            c, pos = snp.split(":")[0], int(snp.split(":")[1])
            rows.append(dict(snp=snp, chr=c, pos=pos, pval=float(d["pval"]),
                             ea=d["EA"], oa=d["OA"], eaf=d["EAF"],
                             rsid=d.get("rsid_gwas", ""), rs_id=d.get("rs_id_gwas", "")))
    rows.sort(key=lambda r: r["pval"])
    picked = []
    for r in rows:
        if any(r["chr"] == q["chr"] and abs(r["pos"] - q["pos"]) < gap for q in picked):
            continue
        r["mhc"] = (r["chr"] == "6" and 25_000_000 <= r["pos"] <= 35_000_000)
        picked.append(r)
        if len(picked) >= top:
            break
    return picked


# ---------------------------------------------------------------- 2 单文件扫描
def scan(tag, path, cols, dual_rsid, offsets, regions):
    """regions: list of dict(id, chr, start, end)；返回 {rid: dict(n, pos_set, rs)}"""
    out = {r["id"]: dict(n=0, pos=set(), rs=set()) for r in regions}
    # 建立 (chr, pos) -> rid 的快速索引（按 chr 分组，用区间判断）
    by_chr = {}
    for r in regions:
        by_chr.setdefault(r["chr"], []).append(r)
    n = 0
    with _open(path) as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {k: hdr.index(v) for k, v in cols.items() if v in hdr}
        ci, pi = ix.get("chr"), ix.get("pos")
        if ci is None or pi is None:
            log("    [警告] %s 缺 chr/pos 列，跳过" % tag)
            return out
        rs_cols = [ix[k] for k in ("rsid", "rs_id") if k in ix] if dual_rsid else \
                  ([ix["rsid"]] if "rsid" in ix else [])
        for row in rd:
            n += 1
            try:
                c = norm_chr(row[ci])
                pos = int(row[pi])
            except (TypeError, ValueError, IndexError):
                continue
            cands = by_chr.get(c)
            if not cands:
                continue
            for r in cands:
                lo, hi = r["start"], r["end"]
                if any(lo <= pos + o <= hi for o in offsets):
                    d = out[r["id"]]
                    d["n"] += 1
                    d["pos"].add(pos)
                    for rc in rs_cols:
                        v = row[rc]
                        if v and v != "NA":
                            d["rs"].add(v)
                    break
    log("    %-28s 扫 %d 行，命中候选区 %d 个变异" % (
        tag, n, sum(out[r["id"]]["n"] for r in regions)))
    return out


def main():
    global LOG
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--win", type=int, default=1_000_000)
    ap.add_argument("--rescan", action="store_true", help="忽略缓存，强制重新扫描")
    a = ap.parse_args()
    LOG = os.path.join(LOGD, "precheck_%s.log" % datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

    log("=" * 78)
    log("Step 8 前置验证：共定位候选区变异覆盖度（窗口 ±%d bp）" % a.win)
    log("=" * 78)

    cand = pick_candidates(a.top, gap=2_000_000)
    log("[1] 候选位点 %d 个（甲减工具按 p 排序、间隔 ≥2 Mb 去重）" % len(cand))
    regions = []
    for i, r in enumerate(cand):
        regions.append(dict(id="R%02d" % (i + 1), chr=r["chr"],
                            start=r["pos"] - a.win, end=r["pos"] + a.win, lead=r))
        log("   R%02d  chr%s:%d  p=%.3g  %s%s" % (
            i + 1, r["chr"], r["pos"], r["pval"], r["rs_id"],
            "  [MHC]" if r["mhc"] else ""))

    RAW = os.path.join(RES, "precheck_raw.json")
    if os.path.exists(RAW) and not a.rescan:
        log("[2] 载入扫描缓存（跳过重扫；加 --rescan 可强制重扫）")
        with io.open(RAW, encoding="utf-8") as f:
            raw = json.load(f)
        res = {}
        for tag, _, _, _, _, _ in DATASETS:
            res[tag] = {rid: dict(n=d["n"], pos=set(d["pos"]), rs=set(d["rs"]))
                        for rid, d in raw.get(tag, {}).items()}
    else:
        log("[2] 扫描各数据集（单遍，统计区域内变异数）")
        res = {}
        for tag, name, path, cols, dual, offs in DATASETS:
            res[tag] = scan(tag, path, cols, dual, offs, regions)
        with io.open(RAW, "w", encoding="utf-8") as f:
            json.dump({tag: {rid: dict(n=d["n"], pos=sorted(d["pos"]), rs=sorted(d["rs"]))
                             for rid, d in res[tag].items()}
                       for tag, _, _, _, _, _ in DATASETS}, f)
        log("    扫描缓存 → %s" % os.path.basename(RAW))

    # ---------------------------------------------------------- 3 汇总
    log("[3] 汇总")
    primary_out = ["cad_nouk", "mvp_meta", "is_sv", "is_laa"]
    rows_out = []
    for r in regions:
        rid = r["id"]
        lead = r["lead"]
        row = dict(id=rid, chr=r["chr"], pos=lead["pos"], rs=lead["rs_id"],
                   rs_alt=lead["rsid"], pval=lead["pval"], mhc=lead["mhc"])
        for tag, _, _, _, _, _ in DATASETS:
            d = res[tag][rid]
            row["n_" + tag] = d["n"]
            row["lead_" + tag] = bool(lead["rs_id"] in d["rs"] or lead["rsid"] in d["rs"])
        # thy hypo 的「真实坐标」集合（SNP=+1 需回退 1）
        thy = set()
        for p in res["thy_hypo"][rid]["pos"]:
            thy.add(p); thy.add(p - 1)
        thy_rs = res["thy_hypo"][rid]["rs"]
        for otag in primary_out:
            opos = res[otag][rid]["pos"]
            row["ov_" + otag] = len(thy & opos) if thy and opos else 0
            # rsID 交集（不受装配版本影响的稳健指标）
            row["ovrs_" + otag] = len(thy_rs & res[otag][rid]["rs"])
        rows_out.append(row)

    # 落盘
    jp = os.path.join(RES, "precheck_coverage.json")
    with open(jp, "w", encoding="utf-8") as f:
        json.dump(dict(window=a.win, candidates=[
            dict(id=r["id"], chr=r["chr"], pos=rows_out[i]["pos"], rs=rows_out[i]["rs"],
                 pval=rows_out[i]["pval"], mhc=rows_out[i]["mhc"],
                 **{k: v for k, v in rows_out[i].items()
                    if k.startswith(("n_", "lead_", "ov_"))})
            for i, r in enumerate(regions)]), f, ensure_ascii=False, indent=2)
    tp = os.path.join(RES, "precheck_coverage.tsv")
    with open(tp, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows_out)

    # 控制台明细
    for row in rows_out:
        log("")
        log("── %s  chr%s:%d (%s)  p=%.3g%s" % (
            row["id"], row["chr"], row["pos"], row["rs"], row["pval"],
            "  [MHC，共定位通常排除]" if row["mhc"] else ""))
        for tag, name, _, _, _, _ in DATASETS:
            log("     %-28s 变异 %6d  lead%s" % (
                name, row["n_" + tag], "✓" if row["lead_" + tag] else "✗"))
        log("     与甲减可配对（位置）：CAD %d | MVP %d | IS-SV %d | IS-LAA %d" % (
            row["ov_cad_nouk"], row["ov_mvp_meta"], row["ov_is_sv"], row["ov_is_laa"]))
        log("     与甲减可配对（rsID）：CAD %d | MVP %d | IS-SV %d | IS-LAA %d" % (
            row["ovrs_cad_nouk"], row["ovrs_mvp_meta"], row["ovrs_is_sv"], row["ovrs_is_laa"]))

    # 总体统计
    log("")
    log("── 总体（%d 个候选区，含/不含 MHC）" % len(rows_out))
    for tag, name, _, _, _, _ in DATASETS:
        vals = [row["n_" + tag] for row in rows_out]
        leads = sum(1 for row in rows_out if row["lead_" + tag])
        log("   %-28s 每区中位 %6d（min %d / max %d），lead 命中 %d/%d" % (
            name, sorted(vals)[len(vals) // 2], min(vals), max(vals), leads, len(vals)))
    log("")
    log("结果 → %s" % os.path.basename(tp))
    log("=" * 78)
    log("完成")


if __name__ == "__main__":
    main()
