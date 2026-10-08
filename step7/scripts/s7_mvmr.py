# -*- coding: utf-8 -*-
"""
Step 7 多变量 MR（MVMR）：把三个甲状腺性状（甲减 / TSH / FT4）同时放进模型，
拆解「甲减 → 动脉粥样硬化性血管疾病」的效应是否**独立于 TSH 与 FT4**。

方法（全部为两步法汇总数据 MVMR 的标准做法）
  1. 工具变量：三个性状单变量工具集的**并集**，再整体 clump 一次
     （r²<0.001、窗口 10 Mb、p1=5e-8），保证工具之间相互独立、无 LD。
  2. MV-IVW（Sanderson 2019）：
        β̂ = (Z'WZ)^{-1} Z'W y ,  W = diag(1/σ_Y²) ,
        Var(β̂) = (Z'WZ)^{-1}
  3. MV-Egger：加截距的 WLS，截距 p 作为整体水平多效性检验。
  4. 同质性/多效性：Q = Σ_j w_j (y_j − Z_jβ̂)² ，df = L − K。
  5. 条件 F 统计量（Sanderson 2021, Testing and correcting for weak and
     pleiotropic instruments in two-sample MVMR, bioRxiv 2020.04.02.021980，式 4 与式 7）：
        Q_{x_k}  = Σ_j (π̂_kj − δ̃_k′ Π̂_{−k,j})² / σ²_{x_k,j}
        σ²_{x_k,j} = σ²_{k,j} + Σ_{m≠k} δ̃²_{k,m} σ²_{m,j}
        F_TS,k   = Q_{x_k} / (L − (K−1))
     δ̃_k 由 π̂_k 对 Π̂_{−k} 的逆方差加权最小二乘得到（权重 1/σ²_{k,j}）。
    ⚠ 上式默认各暴露效应估计之间**协方差为 0**。本项目的三个甲状腺性状来自
      同一项 GWAS（Rand 2025），协方差不为 0 且无法从汇总数据估计，
      故条件 F 与 MV-Q 的零分布仅为近似 —— 结果中已显式标注该限制。

坐标口径（Step 6 实测结论，务必遵守）
  · 面板 .bim 位置 = 真实 GRCh38 1-based 坐标
  · 甲状腺文件（GWAS Catalog 新格式）：SNP 位置 = 真实坐标 + 1；indel = 真实坐标
  · 血管结局文件（SSF harmonised / Sakaue 格式）：位置 = 真实坐标
  因此位置路由对甲状腺文件注册 {P, P+1}，对结局文件只注册 {P}；并一律要求
  等位基因「同向或链互补」才接受，避免邻近位点误配。

用法：
  python step7/scripts/s7_mvmr.py                    # 全部结局
  python step7/scripts/s7_mvmr.py --outcomes sakaue,cad_nouk
"""
import os
import re
import sys
import csv
import gzip
import json
import math
import time
import argparse
import datetime
import subprocess

import numpy as np

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
sys.path.insert(0, os.path.join(ROOT, "step3", "scripts"))
from mrlib import p_norm, p_chi2, ivw, ivw_random_dl  # noqa: E402

PLINK = os.path.join(ROOT, "tools", "plink.exe")
PANEL = os.path.join(ROOT, "ref", "1kg_eur_hg38", "allchr.EUR.biallelicsnps.cbb")
WORK = os.path.join(ROOT, "work")
RES = os.path.join(ROOT, "step7", "results")
LOGD = os.path.join(ROOT, "step7", "logs")
os.makedirs(RES, exist_ok=True)
os.makedirs(LOGD, exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
TRAITS = ["hypo", "tsh", "ft4"]
TRAIT_CN = {"hypo": "甲减", "tsh": "TSH", "ft4": "FT4"}
P_THRESH = 5e-8

# ---------------------------------------------------------------- 甲状腺暴露文件
THY_COLS = dict(rsid="rsid", rs_id="rs_id", chr="chromosome", pos="base_pair_location",
                ea="effect_allele", oa="other_allele", beta="beta",
                se="standard_error", eaf="effect_allele_frequency", p="p_value")
THY = {
    "hypo": dict(path=os.path.join(ROOT, "data", "exposure", "hypo_GCST90572791.h.tsv.gz"),
                 label="Hypothyroidism (Rand 2025; GCST90572791)",
                 n_case=113393, n_control=1065268),
    "tsh":  dict(path=os.path.join(ROOT, "data", "exposure", "tsh_GCST90572789.h.tsv.gz"),
                 label="TSH (Rand 2025; GCST90572789)", n_case=None, n_control=None),
    "ft4":  dict(path=os.path.join(ROOT, "data", "exposure", "ft4_GCST90572790.h.tsv.gz"),
                 label="FT4 (Rand 2025; GCST90572790)", n_case=None, n_control=None),
}

# ---------------------------------------------------------------- 血管结局
SSF_COLS = dict(rsid="hm_rsid", chr="hm_chrom", pos="hm_pos",
                ea="hm_effect_allele", oa="hm_other_allele", beta="hm_beta",
                se="standard_error", eaf="hm_effect_allele_frequency", p="p_value")

OUTCOMES = {
    "sakaue": dict(
        path=os.path.join(ROOT, "data", "outcome", "pad_GCST90018890.h.tsv.gz"),
        cols=dict(rsid="rsid", chr="chromosome", pos="base_pair_location",
                  ea="effect_allele", oa="other_allele", beta="beta",
                  se="standard_error", eaf="effect_allele_frequency", p="p_value"),
        n_case=7114, n_control=475964, label="PAD (Sakaue2021; GCST90018890)",
        overlap="OVERLAP_UKB（暴露 GWAS 含 UKB）"),
    "cad_nouk": dict(
        path=os.path.join(ROOT, "data", "panvascular", "cad_GCST003116.h.tsv.gz"),
        cols=SSF_COLS, n_case=60801, n_control=123504,
        label="CAD (Nikpay2015 CARDIoGRAMplusC4D; GCST003116)", overlap=None),
    "cad_ukb": dict(
        path=os.path.join(ROOT, "data", "panvascular", "cad_GCST90013864.h.tsv.gz"),
        cols=SSF_COLS, n_case=29339, n_control=322724,
        label="CAD (Mbatchou2021 UKB; GCST90013864)",
        overlap="OVERLAP_UKB（暴露 GWAS 含 UKB）"),
    "mi_ukb": dict(
        path=os.path.join(ROOT, "data", "panvascular", "mi_GCST90038610.h.tsv.gz"),
        cols=SSF_COLS, n_case=11081, n_control=473517,
        label="MI (Donertas2021 UKB; GCST90038610)",
        overlap="OVERLAP_UKB（暴露 GWAS 含 UKB）"),
    "is_laa": dict(
        path=os.path.join(ROOT, "data", "panvascular", "is_laa_GCST005840.h.tsv.gz"),
        cols=SSF_COLS, n_case=4373, n_control=406111,
        label="IS-LAA (MEGASTROKE Malik2018; GCST005840)", overlap=None),
    "is_sv": dict(
        path=os.path.join(ROOT, "data", "panvascular", "is_sv_GCST005841.h.tsv.gz"),
        cols=SSF_COLS, n_case=5386, n_control=192662,
        label="IS-SV (MEGASTROKE Malik2018; GCST005841)", overlap=None),
    "mvp_meta": dict(
        path=os.path.join(ROOT, "raw", "mvp_pad", "MVP.te.PAD.dbGAP.txt.gz"),
        cols=dict(rsid="SNP_ID", chr="Chromosome", pos="Position",
                  ea="EA", oa="Allele2", beta="Effect", se="SE",
                  eaf="EAF", p="PValue"),
        n_case=31307, n_control=211753, label="MVP PAD (trans-ancestry meta)",
        overlap=None),
}

LOG = None


def log(m):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(line, flush=True)
    if LOG:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


_IDRE = re.compile(r"^([0-9XYMT]+):(\d+)")


def snp_chr_pos(snp):
    m = _IDRE.match(snp)
    return m.group(1), int(m.group(2))


# ---------------------------------------------------------------- 1 工具集
def load_univariable_instruments():
    """读 work/instruments_<t>.tsv（Step2 产出：各性状的独立工具）"""
    out = {}
    for t in TRAITS:
        p = os.path.join(WORK, "instruments_%s.tsv" % t)
        rows = {}
        with open(p, "r", encoding="utf-8") as f:
            for d in csv.DictReader(f, delimiter="\t"):
                rows[d["SNP"]] = dict(
                    SNP=d["SNP"], chr=d["chr"], pos=int(d["pos"]),
                    EA=d["EA"].upper(), OA=d["OA"].upper(),
                    beta=float(d["beta"]), se=float(d["se"]), p=float(d["pval"]),
                    rs_id=d.get("rs_id_gwas", ""), rsid=d.get("rsid_gwas", ""))
        out[t] = rows
        log("  %-4s 单变量工具 %d 个" % (TRAIT_CN[t], len(rows)))
    return out


def joint_clump(inst):
    """三性状工具并集 → 整体 clump"""
    best = {}
    for t in TRAITS:
        for snp, d in inst[t].items():
            if snp not in best or d["p"] < best[snp][0]:
                best[snp] = (d["p"], t)
    log("  并集去重后候选 %d 个" % len(best))
    inp = os.path.join(WORK, "mvmr_clump_in.txt")
    with open(inp, "w", encoding="utf-8", newline="\n") as f:
        f.write("SNP\tP\n")
        for snp, (p, _) in sorted(best.items(), key=lambda x: x[1][0]):
            f.write("%s\t%.6g\n" % (snp, max(min(float(p), 1.0), 1e-300)))
    out = os.path.join(WORK, "mvmr_clump")
    cmd = [PLINK, "--bfile", PANEL, "--clump", inp,
           "--clump-p1", "5e-8", "--clump-p2", "1",
           "--clump-r2", "0.001", "--clump-kb", "10000",
           "--clump-snp-field", "SNP", "--clump-field", "P",
           "--allow-no-sex", "--out", out]
    pr = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(PLINK))
    if pr.returncode != 0:
        raise RuntimeError("plink rc=%s %s" % (pr.returncode, (pr.stderr or "")[-400:]))
    idx = []
    with open(out + ".clumped", "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3 or parts[0] == "CHR":
                continue
            idx.append(parts[2])
    log("  联合 clump：%d 候选 → %d 个独立工具" % (len(best), len(idx)))
    # 记录每个工具来自哪些性状（p<5e-8 者）
    src = {}
    for snp in idx:
        src[snp] = [t for t in TRAITS if snp in inst[t]]
    return idx, src


# ---------------------------------------------------------------- 2 扫描
def _alleles_ok(e1, o1, e2, o2):
    if {e1, o1} == {e2, o2}:
        return "same"
    if e1 in COMP and o1 in COMP and {COMP[e1], COMP[o1]} == {e2, o2}:
        return "complement"
    return None


def scan_thyroid(tag, snps, ref_rs):
    """扫描甲状腺文件。ref_rs: {rsid: snp}（来自 instruments 表的 rs_id/rsid 列）"""
    spec = THY[tag]
    by_rs, by_pos = {}, {}
    for rsid, snp in ref_rs.items():
        if rsid and rsid != "NA":
            by_rs.setdefault(rsid, snp)
    for snp in snps:
        c, p = snp_chr_pos(snp)
        by_pos.setdefault((c, p), []).append(snp)      # indel 情形
        by_pos.setdefault((c, p + 1), []).append(snp)  # SNP 情形
    got = {}
    n = 0
    want = set(snps)
    with _open(spec["path"]) as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {k: hdr.index(v) for k, v in THY_COLS.items() if v in hdr}
        for row in rd:
            n += 1
            c = row[ix["chr"]]
            try:
                p = int(row[ix["pos"]])
            except (TypeError, ValueError):
                continue
            e = (row[ix["ea"]] or "").upper()
            o = (row[ix["oa"]] or "").upper()
            hits = []
            for col in ("rs_id", "rsid"):
                k = row[ix[col]] if col in ix else ""
                if k and k != "NA" and k in by_rs:
                    hits.append(("rs_id" if col == "rs_id" else "rsid", by_rs[k]))
            for snp in by_pos.get((c, p), []):
                pc, pp = snp_chr_pos(snp)
                hits.append(("pos+1" if p == pp + 1 else "pos+0", snp))
            for route, snp in hits:
                if snp in got or snp not in want:
                    continue
                got[snp] = dict(ea=e, oa=o, beta=row[ix["beta"]], se=row[ix["se"]],
                                eaf=row[ix["eaf"]] if "eaf" in ix else "NA",
                                p=row[ix["p"]], route=route, pos=p,
                                rsid=row[ix["rsid"]] if "rsid" in ix else "")
    log("    甲状腺 %-4s：扫 %d 行 → 命中 %d / %d" % (TRAIT_CN[tag], n, len(got), len(snps)))
    return got


def scan_outcome(key, snps, ref_rs=None):
    spec = OUTCOMES[key]
    by_pos = {}
    for snp in snps:
        c, p = snp_chr_pos(snp)
        by_pos.setdefault((c, p), []).append(snp)
    by_rs = {}
    for rsid, snp in (ref_rs or {}).items():
        if rsid and rsid != "NA" and snp in set(snps):
            by_rs.setdefault(rsid, snp)
    got = {}
    n = 0
    want = set(snps)
    with _open(spec["path"]) as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {k: hdr.index(v) for k, v in spec["cols"].items() if v in hdr}
        for row in rd:
            n += 1
            cc = row[ix["chr"]]
            try:
                pp = int(row[ix["pos"]])
            except (TypeError, ValueError):
                continue          # harmonise 未成功的行坐标为 NA
            hits = list(by_pos.get((cc, pp), []))
            k = row[ix["rsid"]] if "rsid" in ix else ""
            if k and k != "NA" and k in by_rs:
                hits.append(by_rs[k])
            if not hits:
                continue
            e = (row[ix["ea"]] or "").upper()
            o = (row[ix["oa"]] or "").upper()
            for snp in hits:
                if snp in got or snp not in want:
                    continue
                got[snp] = dict(ea=e, oa=o, beta=row[ix["beta"]], se=row[ix["se"]],
                                eaf=row[ix["eaf"]] if "eaf" in ix else "NA",
                                p=row[ix["p"]], route="pos+0")
    log("    结局 %-10s：扫 %d 行 → 命中 %d / %d" % (key, n, len(got), len(snps)))
    return got


# ---------------------------------------------------------------- 3 harmonise
def align(rec, ref_ea, ref_oa):
    """把 rec 的效应对齐到参考等位基因。返回 (beta, se, eaf) 或 None"""
    try:
        b = float(rec["beta"]); s = float(rec["se"])
    except (TypeError, ValueError):
        return None
    if not (s > 0) or not math.isfinite(b):
        return None
    e, o = rec["ea"], rec["oa"]
    if {e, o} == {ref_ea, ref_oa}:
        pass
    elif e in COMP and o in COMP and {COMP[e], COMP[o]} == {ref_ea, ref_oa}:
        e, o = COMP[e], COMP[o]
    else:
        return None
    flipped = (e != ref_ea)
    if flipped:
        b = -b
    try:
        f = float(rec["eaf"])
        if flipped:
            f = 1.0 - f
    except (TypeError, ValueError):
        f = None
    return b, s, f, ("flip" if flipped else "straight"), e, o


def build_ref_rs(inst):
    """gwas rsID -> 面板 SNP ID（三性状合并）"""
    ref_rs = {}
    for t in TRAITS:
        for snp, d in inst[t].items():
            for k in (d.get("rs_id"), d.get("rsid")):
                if k:
                    ref_rs.setdefault(k, snp)
    return ref_rs


def build_dataset(keys, inst, src, ref_rs):
    """返回 merged: {snp: dict(beta_K=[...], se_K=[...], ea, oa, ...)}"""
    thy = {t: scan_thyroid(t, keys, ref_rs) for t in TRAITS}
    # 参考等位基因优先取甲减；缺失则依次 TSH、FT4
    ref = {}
    for snp in keys:
        for t in TRAITS:
            r = thy[t].get(snp)
            if r:
                ref[snp] = (r["ea"], r["oa"])
                break
    merged = {}
    stat = dict(no_ref=0, drop_trait=0, keep=0)
    for snp in keys:
        if snp not in ref:
            stat["no_ref"] += 1
            continue
        re_, ro_ = ref[snp]
        bK, sK, eafK, acts = [], [], [], []
        ok = True
        for t in TRAITS:
            r = thy[t].get(snp)
            a = align(r, re_, ro_) if r else None
            if a is None:
                ok = False
                break
            bK.append(a[0]); sK.append(a[1]); eafK.append(a[2]); acts.append(a[3])
        if not ok:
            stat["drop_trait"] += 1
            continue
        stat["keep"] += 1
        merged[snp] = dict(SNP=snp, EA=re_, OA=ro_, beta_X=bK, se_X=sK, eaf_X=eafK,
                           action="|".join(acts),
                           src=";".join([TRAIT_CN[t] for t in TRAITS if snp in inst[t]]),
                           rs=thy["hypo"].get(snp, {}).get("rsid", ""),
                           pX=[thy[t].get(snp, {}).get("p", "NA") for t in TRAITS])
    log("    三性状齐全 %d / %d（缺参考等位基因 %d、缺少某性状 %d）"
        % (stat["keep"], len(keys), stat["no_ref"], stat["drop_trait"]))
    return merged, thy, stat


# ---------------------------------------------------------------- 4 MVMR 估计
def mvmr_point(Z, y, W):
    A = Z.T @ W @ Z
    V = np.linalg.inv(A)
    b = V @ (Z.T @ W @ y)
    return b, np.sqrt(np.diag(V))


def cond_f(Z, SeX):
    """Sanderson 2021 式 (4)(7)：两样本条件 F 统计量（协方差按 0 处理）"""
    L, K = Z.shape
    out = []
    for k in range(K):
        other = [m for m in range(K) if m != k]
        P = Z[:, other]
        zk = Z[:, k]
        sek = SeX[:, k]
        wk = np.diag(1.0 / sek ** 2)
        d = np.linalg.solve(P.T @ wk @ P, P.T @ wk @ zk)
        resid = zk - P @ d
        var = sek ** 2 + ((d ** 2)[None, :] * SeX[:, other] ** 2).sum(axis=1)
        Qk = float(np.sum(resid ** 2 / var))
        df = L - (K - 1)
        out.append(dict(exposure=TRAITS[k], Qx=Qk, df=df,
                        F=Qk / df if df > 0 else float("nan"),
                        p=p_chi2(Qk, df) if df > 0 else float("nan")))
    return out


def run_mvmr(merged, outcome_key, yrows):
    """构建完整矩阵并跑 MVMR"""
    snps, bK, sK, bY, sY, meta = [], [], [], [], [], []
    for snp, d in merged.items():
        r = yrows.get(snp)
        if not r:
            continue
        a = align(r, d["EA"], d["OA"])
        if a is None:
            continue
        try:
            bb = float(r["beta"]); ss = float(r["se"])
        except (TypeError, ValueError):
            continue
        if not (ss > 0) or not math.isfinite(bb):
            continue
        snps.append(snp)
        bK.append(d["beta_X"]); sK.append(d["se_X"])
        bY.append(a[0]); sY.append(a[1])
        meta.append(dict(SNP=snp, chr=snp_chr_pos(snp)[0], pos=snp_chr_pos(snp)[1],
                         EA=d["EA"], OA=d["OA"], src=d["src"], action=d["action"],
                         beta_obs=a[0], se_obs=a[1],
                         **{("beta_" + t): d["beta_X"][i] for i, t in enumerate(TRAITS)},
                         **{("se_" + t): d["se_X"][i] for i, t in enumerate(TRAITS)},
                         **{("p_" + t): d["pX"][i] for i, t in enumerate(TRAITS)}))
    L = len(snps)
    if L < 4:
        return None, {"nsnp": L, "note": "可用 SNP 不足 4 个，未估计"}
    Z = np.array(bK, dtype=float)
    SeX = np.array(sK, dtype=float)
    y = np.array(bY, dtype=float)
    Sy = np.array(sY, dtype=float)
    W = np.diag(1.0 / Sy ** 2)

    b, se = mvmr_point(Z, y, W)
    K = Z.shape[1]
    resid = y - Z @ b
    Q = float(np.sum(resid ** 2 / Sy ** 2))
    dfq = L - K
    mvmr = {}
    for i, t in enumerate(TRAITS):
        mvmr[t] = dict(beta=float(b[i]), se=float(se[i]), p=p_norm(b[i] / se[i]))
    # MV-Egger（加截距）
    D = np.column_stack([np.ones(L), Z])
    A = D.T @ W @ D
    bE = np.linalg.inv(A) @ (D.T @ W @ y)
    seE = np.sqrt(np.diag(np.linalg.inv(A)))
    egger = dict(intercept=float(bE[0]), intercept_se=float(seE[0]),
                 intercept_pval=p_norm(bE[0] / seE[0]),
                 beta={t: float(bE[i + 1]) for i, t in enumerate(TRAITS)},
                 se={t: float(seE[i + 1]) for i, t in enumerate(TRAITS)})
    cf = cond_f(Z, SeX)
    # 弱工具门槛：全部条件 F > 10
    ok = all(x["F"] > 10 for x in cf)
    return dict(meta=meta, nsnp=L, mvmr=mvmr, egger=egger, cond_f=cf,
                Q=Q, Q_df=dfq, Q_pval=p_chi2(Q, dfq) if dfq > 0 else None,
                weak=("通过（所有条件 F > 10）" if ok else "存在条件 F ≤ 10 的暴露，估计可能偏倚")), \
           {"nsnp": L}


def main():
    global LOG
    ap = argparse.ArgumentParser()
    ap.add_argument("--outcomes", default="all")
    a = ap.parse_args()
    LOG = os.path.join(LOGD, "mvmr_%s.log" % datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    log("=" * 78)
    log("Step 7 MVMR：甲减 / TSH / FT4 三暴露联合模型")
    log("=" * 78)

    log("[1] 读取单变量工具集")
    inst = load_univariable_instruments()
    log("[2] 联合 clump")
    keys, src = joint_clump(inst)
    if not keys:
        log("无工具，退出"); return
    log("[3] 扫描三个甲状腺文件")
    ref_rs = build_ref_rs(inst)
    merged, thy, stat = build_dataset(keys, inst, src, ref_rs)

    # 工具清单落盘
    ip = os.path.join(RES, "mvmr_instruments.tsv")
    with open(ip, "w", encoding="utf-8", newline="") as f:
        cols = ["SNP", "EA", "OA", "src"] + ["beta_" + t for t in TRAITS] + \
               ["se_" + t for t in TRAITS] + ["p_" + t for t in TRAITS]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for snp, d in merged.items():
            row = dict(SNP=snp, EA=d["EA"], OA=d["OA"], src=d["src"])
            for i, t in enumerate(TRAITS):
                row["beta_" + t] = d["beta_X"][i]
                row["se_" + t] = d["se_X"][i]
                row["p_" + t] = d["pX"][i]
            w.writerow(row)
    log("  工具清单 → %s" % os.path.basename(ip))

    outs = list(OUTCOMES) if a.outcomes == "all" else a.outcomes.split(",")
    summary = []
    for key in outs:
        spec = OUTCOMES[key]
        log("[4] 结局 %s（%s）" % (key, spec["label"]))
        yrows = scan_outcome(key, list(merged.keys()), ref_rs)
        res, info = run_mvmr(merged, key, yrows)
        if res is None:
            log("  跳过：%s" % info)
            continue
        univ = {}
        for t in TRAITS:
            jp = os.path.join(ROOT, "step4", "results", "mr_%s_%s.json" % (key, t))
            if os.path.exists(jp):
                with open(jp, "r", encoding="utf-8") as f:
                    j = json.load(f)
                univ[t] = {"beta": j["ivw_fixed"]["beta"], "se": j["ivw_fixed"]["se"],
                           "p": j["ivw_fixed"]["p"], "nsnp": j.get("nsnp")}
        payload = dict(outcome=key, label=spec["label"], n_case=spec["n_case"],
                       n_control=spec["n_control"], overlap=spec["overlap"],
                       nsnp=res["nsnp"], mvmr=res["mvmr"], mv_egger=res["egger"],
                       cond_f=res["cond_f"], Q=res["Q"], Q_df=res["Q_df"],
                       Q_pval=res["Q_pval"], weak_check=res["weak"],
                       univariable=univ,
                       limitation="三个暴露来自同一 GWAS，暴露效应估计间协方差不为零且无法从"
                                  "汇总数据估计；条件 F 与 MV-Q 的零分布仅为近似。")
        with open(os.path.join(RES, "mvmr_%s.json" % key), "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
        with open(os.path.join(RES, "mvmr_instruments_%s.tsv" % key), "w",
                  encoding="utf-8", newline="") as f:
            cols = list(res["meta"][0].keys())
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader(); w.writerows(res["meta"])
        log("   %d SNP | MV-IVW：%s" % (res["nsnp"], "；".join(
            "%s β=%.4f (p=%.3g)" % (TRAIT_CN[t], res["mvmr"][t]["beta"], res["mvmr"][t]["p"])
            for t in TRAITS)))
        log("   条件 F：%s | MV-Q=%.1f (df=%d, p=%.3g)"
            % ("、".join("%s=%.1f" % (TRAIT_CN[x["exposure"]], x["F"]) for x in res["cond_f"]),
               res["Q"], res["Q_df"], res["Q_pval"]))
        for t in TRAITS:
            summary.append(dict(outcome=key, exposure=t,
                                beta_mvmr=res["mvmr"][t]["beta"], se_mvmr=res["mvmr"][t]["se"],
                                p_mvmr=res["mvmr"][t]["p"],
                                beta_univ=(univ.get(t) or {}).get("beta"),
                                se_univ=(univ.get(t) or {}).get("se"),
                                p_univ=(univ.get(t) or {}).get("p"),
                                nsnp=res["nsnp"],
                                condF=[x["F"] for x in res["cond_f"] if x["exposure"] == t][0],
                                MVQ=res["Q"], MVQ_p=res["Q_pval"], overlap=spec["overlap"]))
    if summary:
        sp = os.path.join(RES, "mvmr_summary.tsv")
        with open(sp, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(summary)
        log("汇总 → %s" % os.path.basename(sp))
    log("=" * 78)
    log("完成")


if __name__ == "__main__":
    main()
