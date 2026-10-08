# -*- coding: utf-8 -*-
"""
Step 6 反向 MR：以动脉粥样硬化性血管结局为「暴露」、甲状腺功能为「结局」，
检验是否存在反向因果（PAD / CAD -> 甲减 / TSH / FT4）。

流程（每个暴露独立走一遍）
  1. 扫描暴露 GWAS，收集 p<5e-8 候选（rsID 优先，FinnGen 的 rsids 取逗号前第一个）
  2. rsID -> 面板 SNP ID（chr:pos，GRCh38）：流式扫描 work/panel_id_map.tsv，只留候选
  3. PLINK 1.90 --clump（r²<0.001，窗口 10 000 kb）取独立工具
  4. 扫描三个甲状腺功能 GWAS，提取索引工具的效应（rsID 优先，位置兜底按 +1 偏移）
  5. 逐性状 harmonise + MR（IVW 固定/随机 DL、MR-Egger、加权中位数/众数、
     Cochran Q、MR-PRESSO 等价、LOO、方向性 R² 对比）

与正向分析的口径完全一致（工具 p<5e-8、clump r²<0.001/10 Mb；
等位基因集合匹配 → 链互补翻转 → 回文/EAF 判别）。

用法：
  python step6/scripts/s6_rev_pipeline.py --exposure all
  python step6/scripts/s6_rev_pipeline.py --exposure cad_nikpay
"""
import os
import sys
import csv
import gzip
import json
import math
import time
import argparse
import datetime
import subprocess

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
sys.path.insert(0, os.path.join(ROOT, "step3", "scripts"))
from mrlib import (p_norm, ivw, ivw_random_dl, egger, weighted_median,
                   weighted_mode, ivw_bootstrap_se, presso_equiv, p_binom_two)

PLINK = os.path.join(ROOT, "tools", "plink.exe")
PANEL = os.path.join(ROOT, "ref", "1kg_eur_hg38", "allchr.EUR.biallelicsnps.cbb")
PANEL_MAP = os.path.join(ROOT, "work", "panel_id_map.tsv")
WORK = os.path.join(ROOT, "work")
RES = os.path.join(ROOT, "step6", "results")
LOGD = os.path.join(ROOT, "step6", "logs")
os.makedirs(RES, exist_ok=True)
os.makedirs(LOGD, exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
TRAITS = ["hypo", "tsh", "ft4"]
TRAIT_CN = {"hypo": "甲减", "tsh": "TSH", "ft4": "FT4"}
P_THRESH = 5e-8

# GWAS Catalog SSF harmonised 标准列
SSF = dict(rsid="hm_rsid", chr="hm_chrom", pos="hm_pos", ea="hm_effect_allele",
           oa="hm_other_allele", beta="hm_beta", se="standard_error",
           eaf="hm_effect_allele_frequency", p="p_value")

# ---------------- 暴露 = 血管结局
EXPOSURES = {
    "pad_sakaue": dict(
        path=os.path.join(ROOT, "data", "outcome", "pad_GCST90018890.h.tsv.gz"),
        gz=True, sep="\t",
        cols=dict(rsid="rsid", chr="chromosome", pos="base_pair_location",
                  ea="effect_allele", oa="other_allele", beta="beta",
                  se="standard_error", eaf="effect_allele_frequency", p="p_value"),
        label="PAD (Sakaue 2021; GCST90018890)", short="PAD-Sakaue",
        n_case=7114, n_control=475964, note="暴露队列含 UKB"),
    "pad_finngen": dict(
        path=os.path.join(ROOT, "data", "replication", "finngen_R13_I9_PAD.gz"),
        gz=True, sep="\t",
        cols=dict(rsid="rsids", chr="#chrom", pos="pos", ea="alt", oa="ref",
                  beta="beta", se="sebeta", eaf="af_alt", p="pval"),
        label="PAD (FinnGen R13 I9_PAD)", short="PAD-FinnGen",
        n_case=22244, n_control=460490, note="暴露队列含 FinnGen Freeze 10"),
    "cad_nikpay": dict(
        path=os.path.join(ROOT, "data", "panvascular", "cad_GCST003116.h.tsv.gz"),
        gz=True, sep="\t", cols=SSF,
        label="CAD (Nikpay 2015; GCST003116)", short="CAD-Nikpay",
        n_case=60801, n_control=123504, note="与暴露无样本重叠"),
}

# ---------------- 结局 = 甲状腺功能
OUTCOMES = {
    "hypo": dict(path=os.path.join(ROOT, "data", "exposure", "hypo_GCST90572791.h.tsv.gz"),
                 label="Hypothyroidism (Rand 2025)"),
    "tsh": dict(path=os.path.join(ROOT, "data", "exposure", "tsh_GCST90572789.h.tsv.gz"),
                label="TSH (Rand 2025)"),
    "ft4": dict(path=os.path.join(ROOT, "data", "exposure", "ft4_GCST90572790.h.tsv.gz"),
                label="FT4 (Rand 2025)"),
}
THY = dict(rsid="rsid", rs_id="rs_id", chr="chromosome", pos="base_pair_location",
           ea="effect_allele", oa="other_allele", beta="beta",
           se="standard_error", eaf="effect_allele_frequency", p="p_value")

# 各命中路由的优先级（越小越可信）
ROUTE_RANK = {"rs_id": 0, "rsid": 1, "pos+1": 2, "pos+0": 3}


def _alleles_ok(e1, o1, e2, o2):
    """暴露(X)等位基因 vs 结局(Y)等位基因：返回 same / complement / None"""
    if {e1, o1} == {e2, o2}:
        return "same"
    if e1 in COMP and o1 in COMP and {COMP[e1], COMP[o1]} == {e2, o2}:
        return "complement"
    return None


def _open(path, gz):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if gz else \
           open(path, "r", encoding="utf-8", errors="replace")


LOG = None


def log(m):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(line, flush=True)
    if LOG:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")


# ---------------------------------------------------------------- 1 扫描暴露
def scan_exposure(spec):
    """扫描暴露 GWAS，返回 p<5e-8 候选 list[dict]"""
    c = spec["cols"]
    cand = []
    n_total = 0
    with _open(spec["path"], spec["gz"]) as f:
        rd = csv.reader(f, delimiter=spec["sep"])
        hdr = next(rd)
        ix = {}
        for k, name in c.items():
            if name not in hdr:
                raise KeyError("暴露文件缺少列 %r（实际表头：%s）" % (name, hdr[:20]))
            ix[k] = hdr.index(name)
        for row in rd:
            n_total += 1
            try:
                p = float(row[ix["p"]])
            except (ValueError, IndexError):
                continue
            if not (p < P_THRESH):
                continue
            # 部分 harmonised 文件存在 beta/se 为 NA 的行（无法 harmonise 的位点）→ 跳过
            try:
                beta = float(row[ix["beta"]]); se = float(row[ix["se"]])
            except (ValueError, IndexError):
                continue
            if not (se > 0) or not math.isfinite(beta):
                continue
            rs = (row[ix["rsid"]] or "").split(",")[0].strip()
            cand.append(dict(rsid=rs,
                             chr=row[ix["chr"]], pos=row[ix["pos"]],
                             ea=(row[ix["ea"]] or "").upper(), oa=(row[ix["oa"]] or "").upper(),
                             beta=beta, se=se,
                             eaf=row[ix["eaf"]], p=p))
    # 同 rsID 只留 p 最小
    best = {}
    for d in cand:
        if d["rsid"] and d["rsid"] != "NA":
            if d["rsid"] not in best or d["p"] < best[d["rsid"]]["p"]:
                best[d["rsid"]] = d
    log("  扫描 %d 行 → p<5e-8 命中 %d 行，去重后 %d 个 rsID" % (n_total, len(cand), len(best)))
    return list(best.values())


# ---------------------------------------------------------------- 2 rsID -> 面板ID
def map_to_panel(cand):
    """流式扫描 panel_id_map.tsv，为候选 rsID 找面板 SNP ID（chr:pos）"""
    want = set(d["rsid"] for d in cand if d.get("rsid"))
    got = {}
    with open(PANEL_MAP, "r", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for d in rd:
            if d["rs_id"] in want:
                got.setdefault(d["rs_id"], d)
                if len(got) >= len(want):
                    break
    hit = 0
    for d in cand:
        m = got.get(d["rsid"])
        if m:
            d["panel_id"] = m["id"]
            d["panel_chr"] = m["chr"]
            d["panel_pos"] = int(m["pos"])
            hit += 1
        else:
            d["panel_id"] = None
    log("  面板映射：%d / %d 命中（%.1f%%）" % (hit, len(cand), 100.0 * hit / max(1, len(cand))))
    return cand


# ---------------------------------------------------------------- 3 clump
def run_clump(key, cand):
    rows = [d for d in cand if d.get("panel_id")]
    if not rows:
        return []
    inp = os.path.join(WORK, "rev_clump_in_%s.txt" % key)
    with open(inp, "w", encoding="utf-8", newline="\n") as f:
        f.write("SNP\tP\n")
        for d in sorted(rows, key=lambda x: x["p"]):
            f.write("%s\t%.6g\n" % (d["panel_id"], d["p"]))
    out = os.path.join(WORK, "rev_clump_%s" % key)
    cmd = [PLINK, "--bfile", PANEL, "--clump", inp,
           "--clump-p1", "5e-8", "--clump-p2", "1",
           "--clump-r2", "0.001", "--clump-kb", "10000",
           "--clump-snp-field", "SNP", "--clump-field", "P",
           "--allow-no-sex", "--out", out]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(PLINK))
    if p.returncode != 0:
        log("  [plink] 失败：%s" % (p.stderr or "")[-400:])
        raise RuntimeError("plink rc=%s" % p.returncode)
    idx = []
    clumped = out + ".clumped"
    if os.path.exists(clumped):
        with open(clumped, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split()
                if len(parts) < 3 or parts[0] == "CHR":
                    continue
                idx.append(parts[2])          # 第 3 列 = 索引 SNP
    log("  clump：%d 个候选 → %d 个独立工具" % (len(rows), len(idx)))
    return idx


# ---------------------------------------------------------------- 4 扫描结局（甲状腺）
def scan_outcome(tag, tools):
    """单遍流式扫描甲状腺 GWAS，多路由命中工具。

    坐标规律（本步实测，chr1 n≈35 万，rsID 精确配对 + Ensembl 三方核验）：
      · 面板 .bim 位置 = 真实 GRCh38 1-based 坐标（VCF POS 直写，非 0-based）
      · 甲状腺文件 SNP   ：甲状腺pos = 真实 GRCh38 + 1   （87.0%）
      · 甲状腺文件 indel ：甲状腺pos = 真实 GRCh38       （99.9%）
    故位置路由依次尝试  panel_pos+1  与  panel_pos。
    rsID 路由用甲状腺文件的 rs_id / rsid 两列，且必须等位基因同向或链互补才接受。
    """
    spec = OUTCOMES[tag]
    by_rs, by_pos = {}, {}
    for pid, d in tools.items():
        for col in ("rs_id", "rsid"):
            k = d.get(col)
            if k and k != "NA":
                by_rs.setdefault(k, pid)
        p = int(d["panel_pos"])
        by_pos.setdefault((d["panel_chr"], p + 1), pid)   # SNP 情形
        by_pos.setdefault((d["panel_chr"], p), pid)       # indel 情形

    cands = {}          # pid -> (rank, route, strand, rec)
    n = 0
    with _open(spec["path"], True) as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {k: hdr.index(v) for k, v in THY.items() if v in hdr}
        for row in rd:
            n += 1
            c = row[ix["chr"]]
            p = row[ix["pos"]]
            e = (row[ix["ea"]] or "").upper()
            o = (row[ix["oa"]] or "").upper()
            hits = []
            for col in ("rs_id", "rsid"):
                k = row[ix[col]] if col in ix else ""
                if k and k != "NA" and k in by_rs:
                    hits.append((col, by_rs[k]))
            if (c, p) in by_pos:
                pv = int(p)
                hits.append(("pos+1" if pv == int(tools[by_pos[(c, p)]]["panel_pos"]) + 1
                             else "pos+0", by_pos[(c, p)]))
            for route, pid in hits:
                d = tools[pid]
                st = _alleles_ok(d["ea"], d["oa"], e, o)
                if st is None:
                    continue
                rank = ROUTE_RANK[route]
                if pid not in cands or rank < cands[pid][0]:
                    cands[pid] = (rank, route, st, dict(
                        rsid=row[ix["rsid"]] if "rsid" in ix else "",
                        chr=c, pos=p, ea=e, oa=o, beta=row[ix["beta"]],
                        se=row[ix["se"]], eaf=row[ix["eaf"]], p=row[ix["p"]]))
    out = {}
    for pid, (rank, route, st, rec) in cands.items():
        rec["route"] = route
        rec["strand"] = st
        out[pid] = rec
    rts = {}
    for r in out.values():
        rts[r["route"]] = rts.get(r["route"], 0) + 1
    log("  结局 %s：扫描 %d 行 → 命中工具 %d / %d（%s）"
        % (TRAIT_CN[tag], n, len(out), len(tools),
           "、".join("%s×%d" % (k, v) for k, v in sorted(rts.items())) or "无"))
    return out


# ---------------------------------------------------------------- 5 harmonise + MR
def harmonise(tools, yrows):
    """暴露 = 血管结局工具, 结局 = 甲状腺。返回 (rows, stats)"""
    st = dict(total=0, kept=0, no_outcome=0, allele_mismatch=0, complement=0,
              flipped=0, eaf_missing=0, palindromic_ambiguous=0, eaf_discordant=0,
              palindromic_no_eaf=0)
    ok = []
    for pid, d in tools.items():
        st["total"] += 1
        y = yrows.get(pid)
        if y is None:
            st["no_outcome"] += 1
            continue
        try:
            bY, seY = float(y["beta"]), float(y["se"])
        except (TypeError, ValueError):
            st["no_outcome"] += 1
            continue
        if not (seY > 0) or not math.isfinite(bY) or not math.isfinite(d["beta"]) or not (d["se"] > 0):
            st["no_outcome"] += 1
            continue
        try:
            eafY = float(y["eaf"]); have_eaf = True
        except (TypeError, ValueError):
            eafY = None; have_eaf = False; st["eaf_missing"] += 1
        e1, o1 = d["ea"], d["oa"]          # 暴露（血管）等位基因
        e2, o2 = y["ea"], y["oa"]          # 结局（甲状腺）等位基因
        acted = "straight"
        if {e1, o1} != {e2, o2}:
            if e1 in COMP and o1 in COMP and {COMP[e1], COMP[o1]} == {e2, o2}:
                e2, o2 = COMP[e2], COMP[o2]; st["complement"] += 1; acted = "complement"
            else:
                st["allele_mismatch"] += 1
                continue
        if e2 != e1:                        # 把结局对齐到暴露的效应等位基因
            bY = -bY
            if have_eaf:
                eafY = 1.0 - eafY
            st["flipped"] += 1
            acted = "flip" if acted == "straight" else acted + "_flip"
        pal = (e1 + o1) in ("AT", "TA", "CG", "GC")
        if have_eaf:
            try:
                eafX = float(d["eaf"])
            except (TypeError, ValueError):
                eafX, have_eaf = None, False
        if have_eaf:
            diff = abs(eafX - eafY)
            if pal and diff > 0.1:
                st["palindromic_ambiguous"] += 1; continue
            if not pal and diff > 0.2:
                st["eaf_discordant"] += 1; continue
        elif pal:
            st["palindromic_no_eaf"] += 1; continue
        ok.append(dict(SNP=d["panel_id"], rs=pid, chr=d["panel_chr"], pos=d["panel_pos"],
                       A1=e1, A2=o1,
                       beta_X=d["beta"], se_X=d["se"],
                       eaf_X=(eafX if eafX is not None else ""),
                       beta_Y=bY, se_Y=seY, eaf_Y=(eafY if eafY is not None else ""),
                       p_Y=y["p"], beta_XY=bY / d["beta"], se_XY=seY / abs(d["beta"]),
                       action=acted, route=y.get("route", ""), strand=y.get("strand", ""),
                       pos_Y=y.get("pos", "")))
        st["kept"] += 1
    return ok, st


def _safe(fn, *a, **kw):
    """工具数量过少时部分方法会退化/除零，统一兜底为 None（结果里标 "NA"）。"""
    try:
        return fn(*a, **kw)
    except Exception:
        return None


def _pv(b, se, se_key="se"):
    """安全的 p 值：beta/se 缺失或 se<=0 时返回 None"""
    try:
        if b is None or se is None or not math.isfinite(b) or not math.isfinite(se) or se <= 0:
            return None
        return p_norm(b / se)
    except Exception:
        return None


def f4(x):
    """%.4f 的安全版本"""
    try:
        return "NA" if x is None else "%.4f" % float(x)
    except (TypeError, ValueError):
        return "NA"


def g3(x):
    """%.3g 的安全版本"""
    try:
        return "NA" if x is None else "%.3g" % float(x)
    except (TypeError, ValueError):
        return "NA"


def mr_report(rows, n_case, n_eff_out, label_x, trait, overlap_note=None):
    bX = [r["beta_X"] for r in rows]; sX = [r["se_X"] for r in rows]
    bY = [r["beta_Y"] for r in rows]; sY = [r["se_Y"] for r in rows]
    bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
    k = len(rows)
    iv = _safe(ivw, bR, sR) or (None, None)
    bF, seF = iv
    rnd = _safe(ivw_random_dl, bR, sR)
    eg = _safe(egger, bX, bY, sY) if k >= 3 else None
    bm = _safe(weighted_median, bR, sR) if k >= 3 else None
    sem = _safe(ivw_bootstrap_se, bR, sR, "median") if k >= 3 else None
    bmo = _safe(weighted_mode, bR, sR) if k >= 4 else None
    semo = _safe(ivw_bootstrap_se, bR, sR, "mode", nb=500) if k >= 4 else None
    npos = sum(1 for x in bR if x > 0)
    r2x = sum((r["beta_X"] / r["se_X"]) ** 2 / ((r["beta_X"] / r["se_X"]) ** 2 + n_case) for r in rows)
    r2y = sum((r["beta_Y"] / r["se_Y"]) ** 2 / ((r["beta_Y"] / r["se_Y"]) ** 2 + n_eff_out) for r in rows)
    return {"label": label_x, "trait": trait, "nsnp": k, "n_case": n_case, "n_eff_out": n_eff_out,
            "overlap_note": overlap_note,
            "ivw_fixed": {"beta": bF, "se": seF, "p": _pv(bF, seF)},
            "ivw_random": rnd, "mr_egger": eg,
            "weighted_median": {"beta": bm, "se": sem, "p": _pv(bm, sem)},
            "weighted_mode": {"beta": bmo, "se": semo, "p": _pv(bmo, semo)},
            "presso_equiv": _safe(presso_equiv, bR, sR),
            "sign_test": {"n_positive": npos, "n_negative": k - npos,
                          "p": p_binom_two(min(npos, k - npos), k)},
            "directionality": {"R2_X": r2x, "R2_Y": r2y,
                               "nR2_X": n_case * r2x, "nR2_Y": n_eff_out * r2y,
                               "direction": "forward (X->Y)" if n_case * r2x > n_eff_out * r2y else "reverse?"}}


def write_loo(key, t, rows):
    bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
    k = len(rows); loo = []
    if k >= 2:
        for i in range(k):
            keep = [j for j in range(k) if j != i]
            r2 = _safe(ivw_random_dl, [bR[j] for j in keep], [sR[j] for j in keep])
            if r2:
                loo.append({"SNP": rows[i]["SNP"], "nsnp": k - 1,
                            "beta": r2["beta"], "se": r2["se"], "p": r2["p"]})
    p = os.path.join(RES, "revloo_%s_%s.tsv" % (key, t))
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["SNP", "nsnp", "beta", "se", "p"], delimiter="\t")
        w.writeheader(); w.writerows(loo)
    return p


HCOL = ["SNP", "rs", "chr", "pos", "A1", "A2", "beta_X", "se_X", "eaf_X",
        "beta_Y", "se_Y", "eaf_Y", "p_Y", "beta_XY", "se_XY", "action",
        "route", "strand", "pos_Y"]


def main():
    global LOG
    ap = argparse.ArgumentParser()
    ap.add_argument("--exposure", default="all", help="all 或逗号分隔的 key")
    ap.add_argument("--skip-clump", action="store_true")
    a = ap.parse_args()
    LOG = os.path.join(LOGD, "rev_%s.log" % datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

    keys = list(EXPOSURES) if a.exposure == "all" else a.exposure.split(",")
    n_eff_out = {t: 4.0 * (113393 * 1065268) / (113393 + 1065268) if t == "hypo"
                 else (482873 if t == "tsh" else 191449) for t in TRAITS}
    # 注意：结局（甲状腺）的 n_eff 用其自身样本量；hypo 用 4·nc·nctrl/(nc+nctrl)
    for key in keys:
        spec = EXPOSURES[key]
        log("=" * 78)
        log("反向 MR：暴露 = %s（%s）" % (spec["label"], spec["note"]))
        cand = scan_exposure(spec)
        if not cand:
            log("  无候选，跳过"); continue
        cand = map_to_panel(cand)
        cpath = os.path.join(WORK, "rev_cand_%s.tsv" % key)
        with open(cpath, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["rsid", "chr", "pos", "ea", "oa", "beta",
                                              "se", "eaf", "p", "panel_id", "panel_chr", "panel_pos"],
                               delimiter="\t", extrasaction="ignore")
            w.writeheader(); w.writerows(cand)
        if a.skip_clump:
            idx = [d["panel_id"] for d in cand if d.get("panel_id")]
        else:
            idx = run_clump(key, cand)
        if not idx:
            log("  无独立工具，跳过"); continue
        tools = {d["panel_id"]: d for d in cand if d.get("panel_id") in set(idx)}
        log("  工具集合 %d 个" % len(tools))

        n_case = spec["n_case"]
        summary = []
        for t in TRAITS:
            if t == "hypo":
                neff_y = 4.0 * 113393 * 1065268 / (113393 + 1065268)
            elif t == "tsh":
                neff_y = 482873.0
            else:
                neff_y = 191449.0
            yrows = scan_outcome(t, tools)
            rows, st = harmonise(tools, yrows)
            if not rows:
                log("  %s × %s：无可用工具（结局命中 %d/%d）"
                    % (spec["short"], TRAIT_CN[t], len(yrows), len(tools)))
                with open(os.path.join(RES, "rev_%s_%s.json" % (key, t)), "w",
                          encoding="utf-8") as f:
                    json.dump({"label": spec["label"], "trait": t, "nsnp": 0,
                               "overlap_note": spec["note"],
                               "outcome_label": OUTCOMES[t]["label"],
                               "coverage": {"n_tools": len(tools), "n_outcome_hit": len(yrows),
                                            "n_harmonised": 0},
                               "harmonise": st}, f, ensure_ascii=False, indent=2)
                summary.append(dict(exposure=key, trait=t, nsnp=0, beta="", se="",
                                    p="", p_random="", egger_int_p="", Q_p=""))
                continue
            rep = mr_report(rows, n_case, neff_y, spec["label"], t, spec["note"])
            rep["harmonise"] = st
            rts = {}
            for y in yrows.values():
                rts[y["route"]] = rts.get(y["route"], 0) + 1
            rep["coverage"] = {"n_tools": len(tools), "n_outcome_hit": len(yrows),
                               "hit_pct": round(100.0 * len(yrows) / max(1, len(tools)), 1),
                               "by_route": rts,
                               "n_harmonised": len(rows)}
            rep["outcome_label"] = OUTCOMES[t]["label"]
            with open(os.path.join(RES, "rev_%s_%s.json" % (key, t)), "w", encoding="utf-8") as f:
                json.dump(rep, f, ensure_ascii=False, indent=2)
            with open(os.path.join(RES, "rev_harmonised_%s_%s.tsv" % (key, t)), "w",
                      encoding="utf-8", newline="") as f:
                w = csv.DictWriter(f, fieldnames=HCOL, delimiter="\t", extrasaction="ignore")
                w.writeheader(); w.writerows(rows)
            write_loo(key, t, rows)
            _rnd = rep["ivw_random"] or {}
            _eg = rep["mr_egger"] or {}
            log("  %s → %s: %d SNP | IVW固定 β=%s (p=%s) | 随机 p=%s | Egger截距 p=%s | Q p=%s | 命中 %d/%d"
                % (spec["short"], TRAIT_CN[t], rep["nsnp"],
                   f4(rep["ivw_fixed"]["beta"]), g3(rep["ivw_fixed"]["p"]),
                   g3(_rnd.get("p")), g3(_eg.get("intercept_pval")), g3(_rnd.get("Q_pval")),
                   len(yrows), len(tools)))
            summary.append(dict(exposure=key, trait=t, nsnp=rep["nsnp"],
                                beta=rep["ivw_fixed"]["beta"], se=rep["ivw_fixed"]["se"],
                                p=rep["ivw_fixed"]["p"], p_random=_rnd.get("p"),
                                egger_int_p=_eg.get("intercept_pval"),
                                Q_p=_rnd.get("Q_pval")))
        if summary:
            sp = os.path.join(RES, "rev_%s_summary.tsv" % key)
            with open(sp, "w", encoding="utf-8", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(summary[0].keys()), delimiter="\t")
                w.writeheader(); w.writerows(summary)
            log("  汇总 → %s" % os.path.basename(sp))
    log("=" * 78)
    log("完成")


if __name__ == "__main__":
    main()
