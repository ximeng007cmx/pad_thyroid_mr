# -*- coding: utf-8 -*-
"""
Step 9 附：样本重叠偏倚的「解析 + 蒙特卡洛模拟」——重叠比例 π ↔ 偏倚幅度曲线

================================================================================
一、问题
================================================================================
Step 5 §3 的配对设计观察到：同一暴露（甲减）在两个结局 GWAS 上给出的 β 相差
1.77×（CAD）/ 1.91×（PAD，MVP→Sakaue）/ 3.62×（PAD，MVP→FinnGen）。
此前文档把这些差异全部归因于「样本重叠」。本脚本用模拟定量回答：
**在项目实际的工具强度下，样本重叠最多能造成多大的偏倚？能否解释上述差异？**

================================================================================
二、理论（两样本 MR 的采样结构，无混杂结构假设）
================================================================================
暴露样本 S1（n1）、结局样本 S2（n2），重叠 S_o（n_o）。定义

    π ≡ n_o / n2        （结局样本中与暴露样本重叠的比例；π=1 即全重叠）

逐 SNP j 的估计误差 e_Xj = β̂_Xj − β_Xj、e_Yj = β̂_Yj − β_Yj。在
(a) 个体相互独立、(b) 基因型跨个体独立、(c) 抽样误差近似正态 三条假设下：

    Var(e_Xj) = se_Xj² = σ_X² / n1
    Var(e_Yj) = se_Yj² = σ_Y² / n2
    Cov(e_Xj, e_Yj) = n_o · Cov(X_i, Y_i) / (n1 · n2) = π · β_obs · se_Xj²

其中 β_obs ≡ Cov(X,Y)/σ_X² 是**观察性（含混杂）关联**。第三个等式是纯代数：
    n_o·Cov(X,Y)/(n1 n2) = (n_o/n2)·(Cov(X,Y)/n1) = π·β_obs·(σ_X²/n1) = π·β_obs·se_Xj²
**不依赖混杂的结构**（不需要「同一混杂因子同时进 X 与 Y」这类假设）。

由此得到 IVW 的渐近期望（β_Y = β·β_X，即先忽略多效性）：

        plim β̂ = β · (F_w + π·r) / (F_w + 1)            （精确，见 §四 MC）

其中
    F_w = Σ w_j β_Xj² / Σ w_j se_Xj²   （加权平均工具强度，本项目 ≈ 108–204）
    r   = β_obs / β                     （观察性关联相对于因果效应的倍数）

定义相对放大倍数 R(π) ≡ plim(π) / plim(0)，则

        **R(π) = 1 + π · r / F_w**                       （直线，过点 (0, 1)）

    ⇒ 反解：要产生 R 倍放大，需 π = (R−1)·F_w / r；
      在 π = 1（全重叠）时能买到的上限只有 R_max = 1 + r/F_w；
      若要 **任意** π ≤ 1 都解释 R，则需 r = (R−1)·F_w（即观察性关联必须是因果效应的这么多倍）。

★ 关于两种常见权重口径（**本次由 MC 纠正的一个推导错误，如实记录**）
  项目 s4_pipeline 用的是「比值均值型」IVW：对 R_j = β_Yj/β_Xj 取 w_j = 1/se_Rj²（se_Rj = se_Yj/|β̂_Xj|）的加权平均。
  另有一种「回归型」写法：β̂_Y ~ β̂_X 过原点、权重 1/se_Yj²。二者**代数恒等**：
      Σ (β̂_X²/se_Y²)(β̂_Y/β̂_X) / Σ (β̂_X²/se_Y²) = Σ (β̂_Xβ̂_Y/se_Y²) / Σ (β̂_X²/se_Y²)
  §四 MC 数值上亦完全重合（两列逐点相同）。若把 se_R 当作**固定**的已知量做一阶展开，
  会得到 plim = β·[1 + (1 − πr)/F_w] —— 符号与 MC 相反。原因：该权重含**随机的** β̂_X，
  不能当常数处理。本脚本保留两列输出以显式记录这一纠正，正式结论一律采用 MC 与
  R(π) = 1 + π·r/F_w（MC 最大偏差 0.004）。

**核心结论（见 §五）**：F_w ≈ 109 意味着重叠能造成的最大改变只有百分之几；
要让重叠解释 1.77×，需要 r = 0.774·F_w ≈ 84（观察性关联必须是因果效应的 84 倍）；
要让 3.62×，需要 r ≈ 287。而现实中甲减/TSH 与心血管结局的观察性 RR 只有 1.3–1.7 量级
（PMID 29978767、16828622）→ 相关系数量级为个位数。**重叠不可能解释这些差异。**

================================================================================
三、输入（只读，全部来自既有结果文件，不重跑 MR、不改数）
================================================================================
    step4/results/harmonised_<tag>_<t>.tsv   逐 SNP 的 beta_X / se_X / beta_Y / se_Y
    step3/results/mr_main_<t>.json           PAD 主分析（Sakaue 2021）
    step4/results/mr_<tag>_<t>.json          各复制层 IVW 点估计（观测到的 R 用）
    step9/scripts/mrlib 同口径：IVW = Σw(β_Y/β_X)/Σw，w = 1/se_R²，se_R = se_Y/|β_X|

================================================================================
四、输出
================================================================================
    step9/results/overlap_sim_curve.tsv     π 网格 × 配置 × r → 解析 R 与 MC R
    step9/results/overlap_sim_anchor.tsv    每个「无重叠 vs 有重叠」配对：观测 R、
                                            所需 π、所需 F
    step9/results/overlap_sim_summary.json  汇总（含 MC 与解析的一致性）
    step9/figures/overlap_bias_curve.png    π ↔ 相对偏倚曲线（含观测值对照）
    step9/figures/overlap_bias_mc.png       蒙特卡洛验证（点 + 区间 vs 解析线）
    step9/figures/overlap_bias_required.png 「要解释观测差异需要多弱/多大的重叠」
================================================================================
"""
import io
import os
import csv
import json
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
RES = os.path.join(ROOT, "step9", "results")
FIG = os.path.join(ROOT, "step9", "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.3, "savefig.bbox": "tight"})

RNG_SEED = 20261007
NREP = 10000
PIS = [i / 20.0 for i in range(0, 21)]          # π = 0.00 … 1.00，步长 0.05（解析曲线）
MC_PIS = [0.0, 0.25, 0.50, 0.75, 1.00]          # 蒙特卡洛只在 5 个点验证（省算力）
RS = [1.0, 2.0, 5.0, 10.0, 20.0]                # 观察性/因果 倍数 r 的扫描网格

TRAIT_CN = {"hypo": "Hypothyroidism", "tsh": "TSH", "ft4": "FT4"}


# ============================================================ 数据读取
def load_harm(tag, t, kind="tag"):
    """读取 harmonised 输入（列见 Step2/Step4 落盘约定）。
    kind == "main" 时主分析（PAD-Sakaue）的输入在 step3/results/harmonised_<t>.tsv。"""
    if kind == "main":
        p = os.path.join(ROOT, "step3", "results", "harmonised_%s.tsv" % t)
    else:
        p = os.path.join(ROOT, "step4", "results", "harmonised_%s_%s.tsv" % (tag, t))
    if not os.path.exists(p):
        return None
    rows = list(csv.DictReader(io.open(p, encoding="utf-8"), delimiter="\t"))
    if not rows:
        return None
    return {
        "bX": np.array([float(r["beta_X"]) for r in rows]),
        "sX": np.array([float(r["se_X"]) for r in rows]),
        "bY": np.array([float(r["beta_Y"]) for r in rows]),
        "sY": np.array([float(r["se_Y"]) for r in rows]),
        "n": len(rows),
    }


def mr_json(kind, tag, t):
    """读取既有 MR 结果 JSON（只取 IVW 固定效应点估计，用于观测到的 R）"""
    if kind == "main":
        p = os.path.join(ROOT, "step3", "results", "mr_main_%s.json" % t)
    elif kind == "finngen" and t == "tsh" and os.path.exists(
            os.path.join(ROOT, "step4", "results", "mr_finngen_r13_tsh.json")):
        p = os.path.join(ROOT, "step4", "results", "mr_finngen_r13_tsh.json")
    else:
        p = os.path.join(ROOT, "step4", "results", "mr_%s_%s.json" % (tag, t))
    if not os.path.exists(p):
        return None
    d = json.load(io.open(p, encoding="utf-8"))
    return d


# ============================================================ 工具强度
def f_weighted(w, bX, sX):
    """F_w = Σ w β_X² / Σ w se_X²"""
    return float(np.sum(w * bX ** 2) / np.sum(w * sX ** 2))


def ratio_ivw(w, bXh, bYh):
    """比值均值型 IVW：w_j = 1/se_Rj²，se_Rj = se_Yj/|β̂_Xj|（与项目 s4_pipeline 同口径）"""
    R = bYh / bXh
    return float(np.sum(w * R) / np.sum(w))


def reg_ivw(w2, bXh, bYh):
    """回归型 IVW：w_j = 1/se_Yj²，β̂ = Σwβ_Xβ_Y / Σwβ_X²（过原点加权回归）"""
    return float(np.sum(w2 * bXh * bYh) / np.sum(w2 * bXh ** 2))


# ============================================================ 解析式
# 唯一的解析解：R(π) = 1 + π·r/F_w  （IVW 两种权重口径代数等价，见模块 docstring）
def R_analytic(pi, r, F):
    """R(π) ≡ plim(π)/plim(0) = 1 + π·r/F_w"""
    return 1.0 + pi * r / F


def pi_required(R, r, F):
    """反解：给定 r 与 F，要产生 R 倍放大需要多大的重叠比例 π"""
    return (R - 1.0) * F / r


def r_required(R, F, pi=1.0):
    """反解：π 固定为 1 时，要产生 R 倍放大需要多大的 r = β_obs/β"""
    return (R - 1.0) * F / pi


def F_required(R, r, pi=1.0):
    """反解：给定 R/r/π，工具强度需要多弱"""
    return pi * r / (R - 1.0)


def R_max_at_full_overlap(r, F):
    """π = 1 时重叠能买到的最大放大倍数"""
    return 1.0 + r / F


# ============================================================ 蒙特卡洛
def simulate(bX, sX, sY, beta_true, beta_obs, pi, nrep=NREP, seed=RNG_SEED,
             Z=None, chunk=2500):
    """蒙特卡洛：从真实工具结构出发，按 (π, β_obs) 生成相关的抽样误差。

    e_Xj ~ N(0, se_Xj²)
    e_Yj = (c_j/se_Xj²)·e_Xj + N(0, se_Yj² − c_j²/se_Xj²)，c_j = π·β_obs·se_Xj²
    ⇒ Cov(e_Xj, e_Yj) = c_j（精确）

    返回 {ratio: 比值均值型 β̂ 分布, reg: 回归型 β̂ 分布, rho_max}。
    Z 可外部传入 (Z1, Z2) 以在多个 (r, π) 间复用同一批随机数 → 方差削减。
    """
    m = len(bX)
    c = pi * beta_obs * sX ** 2
    rho = c / (sX * sY)
    if np.any(np.abs(rho) >= 1.0):
        return None                                  # 超出生成模型的可行域
    slope = c / sX ** 2
    cvar = sY ** 2 - c ** 2 / sX ** 2
    if np.any(cvar <= 0):
        return None
    csd = np.sqrt(cvar)

    if Z is None:
        rng = np.random.default_rng(seed)
        Z1 = rng.standard_normal((nrep, m))
        Z2 = rng.standard_normal((nrep, m))
    else:
        Z1, Z2 = Z
    nrep = Z1.shape[0]

    w2 = 1.0 / sY ** 2
    bR = np.empty(nrep)
    bG = np.empty(nrep)
    for a in range(0, nrep, chunk):
        b = min(a + chunk, nrep)
        z1 = Z1[a:b]
        eX = z1 * sX[None, :]
        eY = slope[None, :] * eX + Z2[a:b] * csd[None, :]
        bXh = bX[None, :] + eX
        bYh = beta_true * bX[None, :] + eY
        w = bXh ** 2 / sY[None, :] ** 2              # 1/se_R²，与项目口径一致
        bR[a:b] = (w * (bYh / bXh)).sum(1) / w.sum(1)
        bG[a:b] = (w2[None, :] * bXh * bYh).sum(1) / (w2[None, :] * bXh ** 2).sum(1)
    return {"ratio": bR, "reg": bG, "rho_max": float(np.max(np.abs(rho)))}


# ============================================================ 配置表
# (key, 显示名, trait, 数据来源 kind/tag, 是否有重叠, 无重叠 β 锚点来源)
CONFIGS = [
    ("cad_nouk", "CAD (Nikpay 2015)", "hypo", "tag", "cad_nouk", False, "self"),
    ("cad_ukb", "CAD (UKB Mbatchou 2021)", "hypo", "tag", "cad_ukb", True, "cad_nouk"),
    ("mvp_meta", "PAD (MVP trans-ancestry)", "hypo", "tag", "mvp_meta", False, "self"),
    ("sakaue", "PAD (Sakaue 2021, main)", "hypo", "main", None, True, "mvp_meta"),
    ("finngen", "PAD (FinnGen R13)", "hypo", "tag", "finngen", True, "mvp_meta"),
    ("cad_nouk", "CAD (Nikpay 2015)", "tsh", "tag", "cad_nouk", False, "self"),
    ("cad_ukb", "CAD (UKB Mbatchou 2021)", "tsh", "tag", "cad_ukb", True, "cad_nouk"),
    ("cad_nouk", "CAD (Nikpay 2015)", "ft4", "tag", "cad_nouk", False, "self"),
    ("cad_ukb", "CAD (UKB Mbatchou 2021)", "ft4", "tag", "cad_ukb", True, "cad_nouk"),
]


def build_cfg_table():
    out = []
    for key, label, t, kind, tag, ov, ref in CONFIGS:
        h = load_harm(tag, t, kind)
        if h is None:
            continue
        d = mr_json(kind, tag, t)
        if d is None:
            continue
        w = 1.0 / h["sY"] ** 2
        F = f_weighted(w, h["bX"], h["sX"])
        out.append(dict(key=key, label=label, trait=t, tag=(tag or "sakaue"),
                        kind=kind, overlap=ov, ref=ref, harm=h, Fw=F,
                        beta=d["ivw_fixed"]["beta"], se=d["ivw_fixed"]["se"],
                        nsnp=len(h["bX"]),
                        F_mean=float(np.mean((h["bX"] / h["sX"]) ** 2)),
                        F_min=float(np.min((h["bX"] / h["sX"]) ** 2))))
    return out


# ============================================================ 主流程
def main():
    cfgs = build_cfg_table()
    if not cfgs:
        raise SystemExit("没有可用的配置")
    idx = {(c["key"], c["trait"]): c for c in cfgs}

    # ---------- 1. 观测到的「无重叠 vs 有重叠」配对（甲减为主要看）
    pairs = [
        dict(disease="CAD", base=("cad_nouk", "hypo"), over=("cad_ukb", "hypo"),
             note="CARDIoGRAMplusC4D → UKB（结局来源不同 + 暴露含 UKB）"),
        dict(disease="PAD", base=("mvp_meta", "hypo"), over=("sakaue", "hypo"),
             note="MVP（无重叠）→ Sakaue 2021（暴露含 UKB，部分重叠）"),
        dict(disease="PAD", base=("mvp_meta", "hypo"), over=("finngen", "hypo"),
             note="MVP（无重叠）→ FinnGen R13（暴露含 FinnGen DF10，重叠）"),
    ]

    anchor_rows = []
    for p in pairs:
        b = idx[p["base"]]
        o = idx[p["over"]]
        R_obs = o["beta"] / b["beta"]
        anchor_rows.append(dict(
            disease=p["disease"], base_label=b["label"], over_label=o["label"],
            beta_base=b["beta"], beta_over=o["beta"], R_obs=R_obs,
            Fw_base=b["Fw"], Fw_over=o["Fw"], note=p["note"]))
    # 旧的文档口径（Sakaue→FinnGen、Nikpay→UKB）保留作对照
    bc = idx[("cad_nouk", "hypo")]
    oc = idx[("cad_ukb", "hypo")]
    bs = idx[("sakaue", "hypo")]
    os_ = idx[("finngen", "hypo")]

    # ---------- 2. 解析曲线（细网格）+ 蒙特卡洛（粗网格验证）
    curve_rows = []
    mc_checks = []
    for ci, c in enumerate(cfgs):
        if c["trait"] != "hypo":
            continue
        btrue = c["beta"] if c["ref"] == "self" else idx[(c["ref"], c["trait"])]["beta"]
        h = c["harm"]
        rng = np.random.default_rng(RNG_SEED + 1000 * ci)
        Z = (rng.standard_normal((NREP, len(h["bX"]))),
             rng.standard_normal((NREP, len(h["bX"]))))
        b0 = simulate(h["bX"], h["sX"], h["sY"], btrue, 0.0, 0.0, Z=Z)   # π=0 基线
        m0r, m0g = float(np.mean(b0["ratio"])), float(np.mean(b0["reg"]))
        for r in RS:
            bobs = r * btrue
            for pi in PIS:
                Ra = R_analytic(pi, r, c["Fw"])
                row = dict(cfg=c["key"], label=c["label"], trait=c["trait"], Fw=c["Fw"],
                           r=r, pi=pi, R_analytic=Ra, R_mc=None, R_mc_ratio=None,
                           mc_lo=None, mc_hi=None, mc_rho_max=None, note="")
                if pi in MC_PIS:
                    sim = simulate(h["bX"], h["sX"], h["sY"], btrue, bobs, pi, Z=Z)
                    if sim is None:
                        row["note"] = "超出生成模型可行域（|rho|>=1）"
                    else:
                        row["R_mc"] = float(np.mean(sim["reg"]) / m0g)
                        row["R_mc_ratio"] = float(np.mean(sim["ratio"]) / m0r)
                        row["mc_lo"] = float(np.percentile(sim["reg"], 2.5) / m0g)
                        row["mc_hi"] = float(np.percentile(sim["reg"], 97.5) / m0g)
                        row["mc_rho_max"] = sim["rho_max"]
                        mc_checks.append(dict(
                            cfg=c["key"], r=r, pi=pi, Fw=c["Fw"],
                            R_analytic=Ra, R_mc=row["R_mc"],
                            R_mc_ratio=row["R_mc_ratio"], rho_max=sim["rho_max"]))
                curve_rows.append(row)

    # ---------- 3. 反解：观测到的 R 需要多少 π / 多大的 r / 多弱的 F
    inv_rows = []
    for a in anchor_rows:
        for r in RS:
            inv_rows.append(dict(
                disease=a["disease"], pair="%s -> %s" % (a["base_label"], a["over_label"]),
                R_obs=a["R_obs"], Fw=a["Fw_over"], r=r,
                pi_req=pi_required(a["R_obs"], r, a["Fw_over"]),
                r_req_at_pi1=r_required(a["R_obs"], a["Fw_over"], 1.0),
                F_req_at_pi1_r=r_required(a["R_obs"], a["Fw_over"], 1.0),
                F_req_at_pi1=F_required(a["R_obs"], r, 1.0),
                R_max_at_pi1=R_max_at_full_overlap(r, a["Fw_over"])))

    # ---------- 4. TSH→CAD 翻转的定量（无重叠下为负 → 重叠数据集成显著正向）
    # 注意：TSH 的 β_base < 0 而观察性关联 β_obs > 0，故不能用「β_obs = r·β（r>0）」的参数化；
    # 这里直接以 β_obs（TSH→CAD 的观察性 log-OR）为自变量。
    bs_, os_ = idx[("cad_nouk", "tsh")], idx[("cad_ukb", "tsh")]
    tsh_flip = dict(beta_base=bs_["beta"], se_base=bs_["se"],
                    beta_over=os_["beta"], se_over=os_["se"], Fw=os_["Fw"])
    # 精确反解：plim = (β·F + π·β_obs)/(F+1)，π=1 时 β_obs = β_over·(F+1) − β·F
    tsh_flip["shift_observed"] = os_["beta"] - bs_["beta"]
    tsh_flip["beta_obs_needed_pi1"] = os_["beta"] * (os_["Fw"] + 1.0) - bs_["beta"] * os_["Fw"]
    # 位移公式：plim − β = (β_obs − β)/(F+1)。取 β_obs = 0.10（OR=1.11，TSH→CAD 观察
    # 性 log-OR 的现实上限量级）与 π=1（完全重叠）→ 位移上界
    tsh_flip["beta_obs_assumed"] = 0.10
    tsh_flip["max_shift_pi1_beta_obs010"] = (tsh_flip["beta_obs_assumed"] - bs_["beta"]) / (os_["Fw"] + 1.0)
    tsh_flip["observed_over_max_shift"] = (tsh_flip["shift_observed"]
                                           / tsh_flip["max_shift_pi1_beta_obs010"])

    # ---------- 5. 写 TSV
    def wtsv(path, rows, fields):
        with io.open(path, "w", encoding="utf-8", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            wr.writeheader()
            for r in rows:
                wr.writerow({k: (("%.6g" % v) if isinstance(v, float) else v)
                             for k, v in r.items()})
        print("OK  %s  %d rows  %d B" % (path, len(rows), os.path.getsize(path)))

    wtsv(os.path.join(RES, "overlap_sim_curve.tsv"), curve_rows,
         ["cfg", "label", "trait", "Fw", "r", "pi", "R_analytic", "R_mc",
          "R_mc_ratio", "mc_lo", "mc_hi", "mc_rho_max", "note"])
    wtsv(os.path.join(RES, "overlap_sim_anchor.tsv"), inv_rows,
         ["disease", "pair", "R_obs", "Fw", "r", "pi_req", "F_req_at_pi1",
          "R_max_at_pi1", "r_req_at_pi1"])
    wtsv(os.path.join(RES, "overlap_sim_mc_validation.tsv"), mc_checks,
         ["cfg", "r", "pi", "Fw", "R_analytic", "R_mc", "R_mc_ratio", "rho_max"])

    # ---------- 6. 汇总
    ok = [x for x in mc_checks if x["R_mc_ratio"] is not None]
    maxdev = max(abs(x["R_analytic"] - x["R_mc"]) for x in ok)
    maxdev_eq = max(abs(x["R_mc"] - x["R_mc_ratio"]) for x in ok)   # 两种权重口径的等价性
    worst = max(ok, key=lambda x: abs(x["R_analytic"] - x["R_mc"]))
    Fh = idx[("cad_ukb", "hypo")]["Fw"]
    Ft = idx[("cad_ukb", "tsh")]["Fw"]
    # 头条反解值必须与 observed_pairs 一致：使用**逐对** F_w（而非单一 Fh）
    _cad = [a for a in anchor_rows if a["disease"] == "CAD"][0]
    _pad = sorted([a for a in anchor_rows if a["disease"] == "PAD"], key=lambda a: a["R_obs"])
    _pad_sakaue, _pad_r13 = _pad[0], _pad[1]
    summary = dict(
        seed=RNG_SEED, nrep=NREP, pi_grid=PIS, mc_pi_grid=MC_PIS, r_grid=RS,
        model="plim beta_hat = beta*(F_w + pi*r)/(F_w + 1);  R(pi) = 1 + pi*r/F_w",
        configs=[dict(key=c["key"], label=c["label"], trait=c["trait"], nsnp=c["nsnp"],
                      Fw=round(c["Fw"], 2), F_mean=round(c["F_mean"], 1),
                      F_min=round(c["F_min"], 1), beta=round(c["beta"], 6),
                      se=round(c["se"], 6), overlap=c["overlap"]) for c in cfgs],
        observed_pairs=[dict(disease=a["disease"],
                             pair="%s -> %s" % (a["base_label"], a["over_label"]),
                             beta_base=round(a["beta_base"], 6),
                             beta_over=round(a["beta_over"], 6),
                             R_obs=round(a["R_obs"], 4), Fw=round(a["Fw_over"], 2),
                             r_req_at_pi1=round(r_required(a["R_obs"], a["Fw_over"], 1.0), 1),
                             pi_req_at_r5=round(pi_required(a["R_obs"], 5.0, a["Fw_over"]), 2),
                             R_max_at_pi1_r5=round(R_max_at_full_overlap(5.0, a["Fw_over"]), 4),
                             R_max_at_pi1_r20=round(R_max_at_full_overlap(20.0, a["Fw_over"]), 4),
                             note=a["note"]) for a in anchor_rows],
        mc_validation=dict(max_abs_dev=round(maxdev, 6),
                           worst_case=dict(cfg=worst["cfg"], r=worst["r"], pi=worst["pi"]),
                           max_abs_dev_between_weightings=round(maxdev_eq, 8),
                           tol=0.01, passed=bool(maxdev < 0.01)),
        tsh_cad_flip=tsh_flip,
        headline=dict(
            Fw_hypo=round(Fh, 1), Fw_tsh=round(Ft, 1), Fw_ft4=round(idx[("cad_ukb", "ft4")]["Fw"], 1),
            R_max_pi1_r5=round(R_max_at_full_overlap(5.0, Fh), 4),
            R_max_pi1_r10=round(R_max_at_full_overlap(10.0, Fh), 4),
            R_max_pi1_r20=round(R_max_at_full_overlap(20.0, Fh), 4),
            # 反解值按各配对自己的 F_w 计（与 observed_pairs / 论文 Table 3 同口径）
            r_needed_for_R177=round(r_required(_cad["R_obs"], _cad["Fw_over"], 1.0), 1),
            r_needed_for_R190=round(r_required(_pad_sakaue["R_obs"], _pad_sakaue["Fw_over"], 1.0), 1),
            r_needed_for_R362=round(r_required(_pad_r13["R_obs"], _pad_r13["Fw_over"], 1.0), 1),
            pi_needed_for_R177_at_r5=round(pi_required(_cad["R_obs"], 5.0, _cad["Fw_over"]), 2),
            pi_needed_for_R177_at_r10=round(pi_required(_cad["R_obs"], 10.0, _cad["Fw_over"]), 2),
            Fw_used_for_r_needed=dict(CAD_Nikpay_to_UKB=round(_cad["Fw_over"], 2),
                                      PAD_MVP_to_Sakaue=round(_pad_sakaue["Fw_over"], 2),
                                      PAD_MVP_to_FinnGenR13=round(_pad_r13["Fw_over"], 2)),
            note_r_needed="r_needed_* 按逐对 F_w 计算（pair-specific），与 observed_pairs 及论文 Table 3 一致；R_max_* 用甲减 CAD-UKB 层 F_w。",
        ),
    )
    pj = os.path.join(RES, "overlap_sim_summary.json")
    with io.open(pj, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print("OK  %s  %d B" % (pj, os.path.getsize(pj)))

    make_figures(cfgs, idx, curve_rows, anchor_rows, mc_checks, tsh_flip,
                 maxdev, maxdev_eq)
    return summary, cfgs, idx, curve_rows, anchor_rows, mc_checks, tsh_flip


# ============================================================ 画图
COL = {1.0: "#4C72B0", 2.0: "#55A868", 5.0: "#DD8452", 10.0: "#C44E52", 20.0: "#8172B3"}


def make_figures(cfgs, idx, curve_rows, anchor_rows, mc_checks, tsh_flip,
                 maxdev=0.0, maxdev_eq=0.0):
    cf = idx[("cad_ukb", "hypo")]
    ct = idx[("cad_ukb", "tsh")]

    # ---------- 图 1：π ↔ 相对偏倚曲线
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.9))
    for ax, cfg, tag in ((axes[0], cf, "Hypothyroidism  (F_w = %.0f)" % cf["Fw"]),
                         (axes[1], ct, "TSH  (F_w = %.0f)" % ct["Fw"])):
        for r in RS:
            ys = [R_analytic(p, r, cfg["Fw"]) for p in PIS]
            ax.plot(PIS, ys, lw=1.9, color=COL[r], label="r = %g" % r)
        pts = [(x["pi"], x["R_mc"]) for x in mc_checks
               if x["cfg"] == "cad_ukb" and x["r"] == 5.0 and x["R_mc"] is not None]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=7, mfc="white",
                    mec="#222222", mew=1.4, zorder=5, label="Monte Carlo (r = 5)")
        ax.axhline(1.0, color="#888888", lw=1.0, ls="--")
        ax.axhspan(1.77, 3.70, color="#C44E52", alpha=0.10)
        ax.text(0.99, 3.58, "observed between-GWAS differences\n1.77× / 1.91× / 3.62×",
                fontsize=8.4, color="#8B1A1A", ha="right", va="top")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(0.94, 3.80)
        ax.set_xlabel("overlap proportion  π  (fraction of outcome sample shared)", fontsize=9.5)
        ax.set_ylabel("R(π) = IVW estimate relative to the no-overlap baseline", fontsize=9.5)
        ax.set_title(tag, fontsize=11)
        ax.legend(fontsize=8.2, loc="center left", title="r = β_obs / β",
                  title_fontsize=8.2, framealpha=0.95)
        ax.annotate("", xy=(1.0, 1.16), xytext=(1.0, 1.60),
                    arrowprops=dict(arrowstyle="->", color="#B22222", lw=1.4))
        ax.text(0.965, 1.62, "max reachable\nat π = 1", fontsize=7.6, color="#B22222",
                ha="right")
    fig.suptitle("Sample overlap cannot produce the observed between-GWAS differences: "
                 "R(π) = 1 + π·r/F_w, and F_w ≈ 109 caps the effect at a few per cent",
                 fontsize=11.5, y=1.02)
    p = os.path.join(FIG, "overlap_bias_curve.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  %s  %d B" % (p, os.path.getsize(p)))

    # ---------- 图 2：MC 验证（解析 vs 模拟，恒等线散点）
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.0))
    # (a) 解析 vs MC
    ax = axes[0]
    xs, ys, cs = [], [], []
    for x in mc_checks:
        if x["R_mc"] is None:
            continue
        xs.append(x["R_analytic"])
        ys.append(x["R_mc"])
        cs.append(COL[x["r"]])
    lo = min(xs + ys) - 0.02
    hi = max(xs + ys) + 0.02
    ax.plot([lo, hi], [lo, hi], color="#888888", lw=1.2, ls="--", label="identity")
    ax.scatter(xs, ys, s=38, c=cs, edgecolors="black", linewidths=0.4, zorder=4)
    ax.set_xlabel("analytic  R(π) = 1 + π·r/F_w", fontsize=9.5)
    ax.set_ylabel("Monte-Carlo  R  (%d reps)" % NREP, fontsize=9.5)
    ax.set_title("(a) Analytic solution validated by simulation\n(max |deviation| = %.4f)"
                 % maxdev, fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    # (b) 两种权重口径的等价性（数值证明）
    ax = axes[1]
    xs2 = [x["R_mc_ratio"] for x in mc_checks if x["R_mc_ratio"] is not None]
    ys2 = [x["R_mc"] for x in mc_checks if x["R_mc"] is not None]
    ax.plot([lo, hi], [lo, hi], color="#888888", lw=1.2, ls="--")
    ax.scatter(xs2, ys2, s=34, color="#8172B3", edgecolors="black", linewidths=0.4, zorder=4)
    ax.set_xlabel("ratio-mean IVW   Σw·(β_Y/β_X)/Σw,  w = 1/se_R²", fontsize=9.5)
    ax.set_ylabel("regression IVW   Σw·β_Xβ_Y / Σw·β_X²,  w = 1/se_Y²", fontsize=9.5)
    ax.set_title("(b) The two weightings are algebraically identical\n"
                 "(max |difference| = %.1e)" % maxdev_eq, fontsize=10)
    hh = [Line2D([], [], marker="o", ls="", mfc=COL[r], mec="black", ms=7, label="r = %g" % r)
          for r in RS]
    fig.legend(hh, [h.get_label() for h in hh], fontsize=8.4, ncol=5, frameon=False,
               loc="lower center", bbox_to_anchor=(0.5, -0.055))
    fig.suptitle("Monte-Carlo validation  (5 configurations × 5 overlap proportions × "
                 "5 confounding ratios = 125 cells)", fontsize=11.5, y=1.02)
    p = os.path.join(FIG, "overlap_bias_mc.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  %s  %d B" % (p, os.path.getsize(p)))

    # ---------- 图 3：要解释观测差异，需要什么
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.7))
    obs = [(1.774, "CAD  Nikpay→UKB", "#C44E52"),
           (1.913, "PAD  MVP→Sakaue", "#DD8452"),
           (3.622, "PAD  MVP→FinnGen", "#8172B3")]
    # (a) π = 1 时所需的 r  vs  F_w
    ax = axes[0]
    Fs = np.logspace(0, 3, 300)
    for R, lab_, c in obs:
        ax.plot(Fs, r_required(R, Fs, 1.0), lw=1.9, color=c, label=lab_)
    ax.axhspan(1, 20, color="#55A868", alpha=0.13)
    ax.text(620, 4.0, "plausible range of r\n(observational RR 1.3–1.7 vs\ncausal OR ≈ 1.03 → r ≈ 5–15)",
            fontsize=7.4, color="#2F6B3F", ha="right", va="center")
    ax.axvline(cf["Fw"], color="#222222", lw=1.2, ls=":")
    ax.text(cf["Fw"] * 1.15, 1200, "F_w = %.0f\n(this project)" % cf["Fw"],
            fontsize=7.8, color="#222222", va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 1000)
    ax.set_ylim(0.02, 2000)
    ax.set_xlabel("weighted instrument strength  F_w   (log scale)", fontsize=9.5)
    ax.set_ylabel("r = β_obs / β  required  (π = 1, log scale)", fontsize=9.5)
    ax.set_title("(a) To buy the observed inflation with overlap, the observational\n"
                 "association would have to be 84×–287× the causal effect", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    # (b) 所需的 π  vs  观测 R
    ax = axes[1]
    Rs = np.linspace(1.0, 4.0, 300)
    for r in (1.0, 5.0, 10.0, 20.0):
        ax.plot(Rs, pi_required(Rs, r, cf["Fw"]), lw=1.9, color=COL[r], label="r = %g" % r)
    ax.axhline(1.0, color="#B22222", ls="--", lw=1.4)
    ax.text(1.55, 1.30, "π = 1  (complete overlap)", fontsize=8.4, color="#B22222",
            ha="left", va="bottom")
    ax.axhspan(0, 1, color="#55A868", alpha=0.12)
    ax.text(1.90, 0.42, "attainable region  (π ≤ 1)", fontsize=8.6, ha="center",
            color="#2F6B3F")
    for R, lab_, c in obs:
        ax.plot([R], [1.0], "v", ms=9, color=c, mec="black", mew=0.5, zorder=6)
    ax.text(2.45, 6.0,
            "observed  R = 1.77 / 1.91 / 3.62\n"
            "→ required π = 17 / 20 / 60   (at r = 5)\n"
            "even r = 1 would need π = 84",
            fontsize=8.2, color="#8B1A1A", va="center")
    ax.set_xlim(1.0, 4.0)
    ax.set_ylim(0, 15)
    ax.set_xlabel("observed inflation  R", fontsize=9.5)
    ax.set_ylabel("overlap proportion π  required", fontsize=9.5)
    ax.set_title("(b) Required overlap exceeds the whole outcome sample by 17×–60×\n"
                 "(▼ = observed R, drawn on the π = 1 ceiling; green = achievable)",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper left", title="r = β_obs / β", title_fontsize=8)
    p = os.path.join(FIG, "overlap_bias_required.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  %s  %d B" % (p, os.path.getsize(p)))


if __name__ == "__main__":
    main()
    print("ALL DONE")
