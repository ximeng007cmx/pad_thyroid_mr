# -*- coding: utf-8 -*-
"""
Step 8 共定位主程序：coloc.abf（纯 Python，见 coloc_abf.py）。

设计（与前置验证结论一致）
  · 只做**甲减**主效位点，窗口 ±1 Mb（10 个位点，去重后间隔 ≥2 Mb）
  · 结局只用**无样本重叠**队列：CAD-Nikpay、MVP-PAD、IS-SV、IS-LAA
  · 变体集合 = 两文件**共同存在**、等位基因同向/互补、MAF ≥ 0.01
  · 后验敏感性：sd.prior ∈ {0.15, 0.2}，p12 ∈ {1e-6, 1e-5, 1e-4}

坐标口径（实测）
  · 甲状腺（Rand2025）：SNP 位置 = 真实 GRCh38 + 1，indel = 真实
  · CAD-Nikpay / IS-SV / IS-LAA：hm_pos = 真实 GRCh38
  · ⚠ MVP-PAD：**GRCh37(hg19)** → 关闭位置匹配，只走 rsID

性能：每个文件**单遍扫描**收集全部候选区（共 5 次扫描），不做 区域×文件 次扫描。

用法：python step8/scripts/s8_coloc.py
"""
import os
import sys
import csv
import gzip
import json
import time
import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from coloc_abf import coloc_abf  # noqa: E402

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
RES = os.path.join(ROOT, "step8", "results")
LOGD = os.path.join(ROOT, "step8", "logs")
os.makedirs(RES, exist_ok=True)
os.makedirs(LOGD, exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
WIN = 1_000_000
MAF_MIN = 0.01

THY_REL = "data/exposure/hypo_GCST90572791.h.tsv.gz"
THY_COLS = dict(chr="chromosome", pos="base_pair_location", ea="effect_allele",
                oa="other_allele", beta="beta", se="standard_error",
                eaf="effect_allele_frequency", p="p_value", rsid="rsid", rs_id="rs_id")
SSF_COLS = dict(chr="hm_chrom", pos="hm_pos", ea="hm_effect_allele", oa="hm_other_allele",
                beta="hm_beta", se="standard_error", eaf="hm_effect_allele_frequency",
                p="p_value", rsid="hm_rsid")
MVP_COLS = dict(chr="Chromosome", pos="Position", ea="EA", oa="Allele2",
                beta="Effect", se="SE", eaf="EAF", p="PValue", rsid="SNP_ID")

OUTCOMES = [
    # key, 显示名, 路径, 列, 是否允许位置匹配, 备注
    ("cad_nouk", "CAD (Nikpay2015; GCST003116)",
     "data/panvascular/cad_GCST003116.h.tsv.gz", SSF_COLS, True, "GRCh38"),
    ("mvp_meta", "PAD (MVP trans-ancestry meta)",
     "raw/mvp_pad/MVP.te.PAD.dbGAP.txt.gz", MVP_COLS, False, "GRCh37 → 仅 rsID"),
    ("is_sv", "IS-SV (MEGASTROKE; GCST005841)",
     "data/panvascular/is_sv_GCST005841.h.tsv.gz", SSF_COLS, True, "GRCh38"),
    ("is_laa", "IS-LAA (MEGASTROKE; GCST005840)",
     "data/panvascular/is_laa_GCST005840.h.tsv.gz", SSF_COLS, True, "GRCh38"),
]

LOG = None


def log(m):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(line, flush=True)
    if LOG:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")


def _open(rel):
    p = os.path.join(ROOT, rel)
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.endswith(".gz") \
        else open(p, "r", encoding="utf-8", errors="replace")


def nchr(c):
    c = (c or "").strip()
    return c[3:].upper() if c.lower().startswith("chr") else c.upper()


def scan_by_rs(rel, cols, rs2rid):
    """按 rsID 收集变异（不依赖坐标窗口）——用于 GRCh37 文件（MVP）"""
    out = {}
    want = set(rs2rid)
    n = hit = 0
    with _open(rel) as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {k: hdr.index(v) for k, v in cols.items() if v in hdr}
        for row in rd:
            n += 1
            s = row[ix["rsid"]]
            if s not in want:
                continue
            rid = rs2rid[s]
            try:
                b = float(row[ix["beta"]]); se = float(row[ix["se"]])
            except (TypeError, ValueError):
                continue
            if not (se > 0) or not np.isfinite(b):
                continue
            try:
                eaf = float(row[ix["eaf"]])
            except (TypeError, ValueError):
                eaf = None
            out.setdefault(rid, []).append(
                dict(pos=int(row[ix["pos"]]), beta=b, se=se, eaf=eaf,
                     ea=(row[ix["ea"]] or "").upper(), oa=(row[ix["oa"]] or "").upper(), rs=s))
            hit += 1
    log("    %-40s 扫 %d 行 → 按 rsID 命中 %d 个变异" % (os.path.basename(rel), n, hit))
    return out


def scan_all_regions(rel, cols, regions, dual_rsid=False):
    """单遍扫描，返回 {rid: [rec, ...]}"""
    out = {r["id"]: [] for r in regions}
    by_chr = {}
    for r in regions:
        by_chr.setdefault(r["chr"], []).append(r)
    n = 0
    with _open(rel) as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {k: hdr.index(v) for k, v in cols.items() if v in hdr}
        ci, pi = ix["chr"], ix["pos"]
        for row in rd:
            n += 1
            c = nchr(row[ci])
            cands = by_chr.get(c)
            if not cands:
                continue
            try:
                pos = int(row[pi])
            except (TypeError, ValueError):
                continue
            tgt = None
            for r in cands:
                if r["start"] <= pos <= r["end"]:
                    tgt = r["id"]
                    break
            if tgt is None:
                continue
            try:
                b = float(row[ix["beta"]]); s = float(row[ix["se"]])
            except (TypeError, ValueError):
                continue
            if not (s > 0) or not np.isfinite(b):
                continue
            try:
                eaf = float(row[ix["eaf"]])
            except (TypeError, ValueError):
                eaf = None
            rec = dict(pos=pos, beta=b, se=s, eaf=eaf,
                       ea=(row[ix["ea"]] or "").upper(), oa=(row[ix["oa"]] or "").upper())
            if "rsid" in ix and row[ix["rsid"]] not in ("", "NA"):
                rec["rs"] = row[ix["rsid"]]
            if dual_rsid and "rs_id" in ix and row[ix["rs_id"]] not in ("", "NA"):
                rec["rs_id"] = row[ix["rs_id"]]
            out[tgt].append(rec)
    log("    %-40s 扫 %d 行 → 候选区合计 %d 变异" % (
        os.path.basename(rel), n, sum(len(v) for v in out.values())))
    return out


def build_index(thy):
    by_rs, by_pos = {}, {}
    for r in thy:
        for k in ("rs_id", "rs"):
            if r.get(k):
                by_rs.setdefault(r[k], r)
        by_pos.setdefault(r["pos"] - 1, r)   # SNP：文件位 = 真实+1
        by_pos.setdefault(r["pos"], r)       # indel：文件位 = 真实
    return by_rs, by_pos


def align(thy, out):
    te, to = thy["ea"], thy["oa"]
    oe, oo = out["ea"], out["oa"]
    if {oe, oo} == {te, to}:
        flip = (oe != te)
    elif te in COMP and to in COMP and {COMP[oe], COMP[oo]} == {te, to}:
        flip = (COMP[oe] != te)
    else:
        return None
    f = thy["eaf"]
    if f is not None:
        f = min(f, 1.0 - f)
        if f < MAF_MIN:
            return None
    return (thy["beta"], thy["se"], (-out["beta"] if flip else out["beta"]), out["se"])


def main():
    global LOG
    LOG = os.path.join(LOGD, "coloc_%s.log" % datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    log("=" * 78)
    log("Step 8 共定位 coloc.abf：甲减（Rand2025）× 无重叠血管结局")
    log("=" * 78)

    pj = os.path.join(RES, "precheck_coverage.json")
    with open(pj, "r", encoding="utf-8") as f:
        pc = json.load(f)
    cands = pc["candidates"]
    regions = [dict(id=c["id"], chr=c["chr"], start=c["pos"] - WIN,
                    end=c["pos"] + WIN, lead=c) for c in cands]
    log("[1] 候选区 %d 个（窗口 ±%d bp）" % (len(regions), WIN))

    log("[2] 单遍扫描：甲状腺（甲减）")
    thy_all = scan_all_regions(THY_REL, THY_COLS, regions, dual_rsid=True)
    # 甲状腺各区的 rsID 集合 → 供 GRCh37 文件按 rsID 取变异
    rs2rid = {}
    for rid, recs in thy_all.items():
        for r in recs:
            for k in ("rs_id", "rs"):
                if r.get(k):
                    rs2rid.setdefault(r[k], rid)
    log("   甲状腺候选区 rsID 合计 %d 个（用于 GRCh37 文件匹配）" % len(rs2rid))
    out_all = {}
    for key, name, rel, cols, allow_pos, note in OUTCOMES:
        log("[2] 扫描：%s（%s）" % (name, note))
        if allow_pos:
            out_all[key] = scan_all_regions(rel, cols, regions)
        else:
            out_all[key] = scan_by_rs(rel, cols, rs2rid)

    log("[3] 逐区逐结局做 coloc.abf")
    rows = []
    for reg in regions:
        rid, c, pos, lead = reg["id"], reg["chr"], reg["lead"]["pos"], reg["lead"]
        log("")
        log("── %s chr%s:%d (%s)%s" % (rid, c, pos, lead.get("rs", ""),
                                      "  [MHC]" if lead.get("mhc") else ""))
        thy = thy_all[rid]
        by_rs, by_pos = build_index(thy)
        log("   甲状腺区变异 %d" % len(thy))
        for key, name, rel, cols, allow_pos, note in OUTCOMES:
            outv = out_all[key].get(rid, [])
            pairs, n_rs, n_pos = [], 0, 0
            for r in outv:
                t = by_rs.get(r.get("rs", ""))
                if t is not None:
                    n_rs += 1
                elif allow_pos:
                    t = by_pos.get(r["pos"])
                    if t is not None:
                        n_pos += 1
                if t is None:
                    continue
                a = align(t, r)
                if a is None:
                    continue
                pairs.append(a)
            if not pairs:
                log("   %-34s 区变异 %6d → 可配对 0（跳过）" % (name, len(outv)))
                rows.append(dict(region=rid, chr=c, pos=pos, rs=lead.get("rs", ""),
                                 mhc=int(bool(lead.get("mhc"))), outcome=key,
                                 outcome_label=name, assembly=note, n=0,
                                 n_out_region=len(outv), m_rs=0, m_pos=0,
                                 verdict="无可配对变异"))
                continue
            d15 = coloc_abf(pairs, 1e-4, 1e-4, 1e-5, 0.15)
            d20 = coloc_abf(pairs, 1e-4, 1e-4, 1e-5, 0.20)
            s6 = coloc_abf(pairs, 1e-4, 1e-4, 1e-6, 0.15)
            s4 = coloc_abf(pairs, 1e-4, 1e-4, 1e-4, 0.15)
            log("   %-34s 区变异 %6d → 可配对 %4d (rs %4d/pos %3d) | "
                "H4=%.3f H3=%.3f H1=%.3f H2=%.3f H0=%.3f | %s" % (
                    name, len(outv), len(pairs), n_rs, n_pos,
                    d15["H4"], d15["H3"], d15["H1"], d15["H2"], d15["H0"], d15["verdict"]))
            rows.append(dict(
                region=rid, chr=c, pos=pos, rs=lead.get("rs", ""),
                mhc=int(bool(lead.get("mhc"))), outcome=key, outcome_label=name,
                assembly=note, n=len(pairs), n_out_region=len(outv),
                m_rs=n_rs, m_pos=n_pos,
                H0=d15["H0"], H1=d15["H1"], H2=d15["H2"], H3=d15["H3"], H4=d15["H4"],
                verdict=d15["verdict"], H4_sd020=d20["H4"],
                H4_p12_1e6=s6["H4"], H4_p12_1e4=s4["H4"]))

    keys = sorted({k for r in rows for k in r})
    with open(os.path.join(RES, "coloc_results.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    log("")
    log("汇总 → coloc_results.tsv（%d 行）" % len(rows))
    # 结论摘要
    h4 = [r for r in rows if r.get("H4") and r["H4"] >= 0.8]
    h3 = [r for r in rows if r.get("H3") and r["H3"] >= 0.8]
    log("PP.H4 ≥ 0.8：%d 条；PP.H3 ≥ 0.8：%d 条（共 %d 条有效比较）" % (
        len(h4), len(h3), sum(1 for r in rows if r.get("n"))))
    log("=" * 78)
    log("完成")


if __name__ == "__main__":
    main()
