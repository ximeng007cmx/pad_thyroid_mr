# -*- coding: utf-8 -*-
"""
Step 5 / 6-7 / 8 / 9 结果配图（投稿用）

背景：step3/figures（12 张）与 step4/figures（94 张）已覆盖主分析与复制层，
      但 Step 5（泛血管谱总览）、Step 6-7（反向 MR / MVMR）、Step 8（共定位）、
      Step 9（统计功效）四个交付此前**没有配图**。本脚本补齐这 4 组。

原则：
  · 所有数值一律从既有结果文件（JSON / TSV）读取，**不重跑 MR、不手抄、不改数**
  · 英文标注（与 s4_figures.py 同风格，规避中文字体依赖）
  · 重叠状态在图上显式区分（空心 + 灰边 = 含样本重叠）

输入
  step3/results/mr_main_<t>.json            PAD 主分析（Sakaue 2021）
  step4/results/mr_<tag>_<t>.json           tag ∈ finngen, mvp_meta, cad_nouk, cad_ukb, mi_ukb, is_laa, is_sv
  step6/results/rev_<exp>_<t>.json          exp ∈ cad_nikpay, pad_sakaue, pad_finngen
  step7/results/mvmr_<key>.json
  step8/results/coloc_results.tsv
  step9/results/power_forward.tsv

输出
  step5/figures/forest_panvascular_overview.png
  step6/figures/rev_mr_forest.png
  step7/figures/mvmr_single_vs_multi.png
  step8/figures/coloc_prob_heatmap.png
  step9/figures/power_mde80_or.png
  step9/figures/power_curve.png
"""
import io
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import csv

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
Z = 1.959963984540054
E = math.e

plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.3, "savefig.bbox": "tight"})

TRAITS = [("hypo", "Hypothyroidism"), ("tsh", "TSH"), ("ft4", "FT4")]
TCOLOR = {"hypo": "#C44E52", "tsh": "#4C72B0", "ft4": "#55A868"}


def J(p):
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def figdir(step):
    d = os.path.join(ROOT, step, "figures")
    os.makedirs(d, exist_ok=True)
    return d


def phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def power(beta_true, se, alpha=0.05):
    lam = abs(beta_true) / se
    return phi(lam - Z) + phi(-lam - Z)


# ============================================================ 图 1 · Step 5 泛血管谱总览森林图
# (显示名, 数据来源, tag, 是否含样本重叠)
OUTCOMES = [
    ("PAD (Sakaue 2021, main)", "main", None, True),
    ("PAD (FinnGen R13)", "tag", "finngen", True),
    ("PAD (MVP trans-ancestry)", "tag", "mvp_meta", False),
    ("CAD (Nikpay 2015)", "tag", "cad_nouk", False),
    ("CAD (UKB Mbatchou 2021)", "tag", "cad_ukb", True),
    ("MI (UKB Donertas 2021) ‡", "tag", "mi_ukb", True),
    ("IS-LAA (MEGASTROKE)", "tag", "is_laa", False),
    ("IS-SV (MEGASTROKE)", "tag", "is_sv", False),
]


def load_mr(kind, tag, t):
    if kind == "main":
        p = os.path.join(ROOT, "step3", "results", "mr_main_%s.json" % t)
    else:
        p = os.path.join(ROOT, "step4", "results", "mr_%s_%s.json" % (tag, t))
    if not os.path.exists(p):
        return None
    return J(p)


def _overlap(kind, tag, t):
    """暴露×结局 级的重叠判定（比结局级更准确）。
    Sakaue 主分析：暴露含 UKB，结局 EUR 亦含 UKB → 重叠；
    FinnGen：仅甲减重叠（暴露含 FinnGen DF10），TSH/FT4 无重叠；
    MVP / cad_nouk / is_laa / is_sv：无重叠；cad_ukb / mi_ukb：重叠。"""
    if kind == "main":
        return True
    if tag in ("cad_ukb", "mi_ukb"):
        return True
    if tag == "finngen":
        return (t == "hypo")
    return False


def fig1_panvascular():
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    n = len(OUTCOMES)
    # 每个结局占一个 y 槽，3 个暴露在槽内错开
    off = {"hypo": 0.24, "tsh": 0.0, "ft4": -0.24}
    yt, yl = [], []
    for i, (name, kind, tag, ov) in enumerate(OUTCOMES):
        y0 = n - 1 - i
        yt.append(y0)
        yl.append(name)
        for t, cn in TRAITS:
            d = load_mr(kind, tag, t)
            if not d:
                continue
            ov_t = _overlap(kind, tag, t)
            b = d["ivw_fixed"]["beta"]
            se = d["ivw_fixed"]["se"]
            orv = E ** b
            lo, hi = E ** (b - Z * se), E ** (b + Z * se)
            y = y0 + off[t]
            mfc = TCOLOR[t] if not ov_t else "white"
            ax.errorbar(orv, y, xerr=[[orv - lo], [hi - orv]], fmt="o",
                        ms=6.2, color=TCOLOR[t], ecolor=TCOLOR[t],
                        elinewidth=1.3, capsize=2.6,
                        markerfacecolor=mfc, markeredgecolor=TCOLOR[t], markeredgewidth=1.4,
                        zorder=3)
    ax.axvline(1.0, color="#888888", lw=1.0, ls="--", zorder=1)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=9.5)
    ax.set_xscale("log")
    ax.set_xlabel("OR per doubling of genetic liability (IVW fixed effect, 95% CI)",
                  fontsize=10)
    ax.set_title("Pan-vascular spectrum: genetically predicted thyroid function vs\n"
                 "atherosclerotic vascular outcomes (3 exposures × 8 outcomes)",
                 fontsize=11.5, pad=12)
    handles = [Line2D([], [], marker="o", ls="", color=TCOLOR[t], ms=7, label=cn)
               for t, cn in TRAITS]
    handles += [
        Line2D([], [], marker="o", ls="", mfc="white", mec="#555555", ms=7,
               label="open marker = sample overlap (reference only)"),
        Line2D([], [], color="none", label="‡ MI file is on a linear-probability scale"),
    ]
    ax.legend(handles=handles, fontsize=8.4, loc="lower right", framealpha=0.95)
    ax.set_ylim(-0.75, n - 0.25)
    p = os.path.join(figdir("step5"), "forest_panvascular_overview.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  " + p, os.path.getsize(p), "B")


# ============================================================ 图 2 · Step 6 反向 MR 森林图
REV_EXP = [
    ("cad_nikpay", "CAD (Nikpay 2015)", False),
    ("pad_sakaue", "PAD (Sakaue 2021)", True),
    ("pad_finngen", "PAD (FinnGen R13)", True),
]


def fig2_reverse():
    rows = []
    for exp, name, ov in REV_EXP:
        for t, cn in TRAITS:
            p = os.path.join(ROOT, "step6", "results", "rev_%s_%s.json" % (exp, t))
            if not os.path.exists(p):
                continue
            d = J(p)
            rows.append((name, cn, t, d["ivw_fixed"]["beta"], d["ivw_fixed"]["se"],
                         d["ivw_fixed"]["p"], ov))
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ys = []
    ylabels = []
    cur = 0
    for exp, name, ov in REV_EXP:
        for t, cn in TRAITS:
            m = [r for r in rows if r[0] == name and r[1] == cn]
            if not m:
                continue
            _, _, tt, b, se, pv, ovv = m[0]
            y = cur
            ys.append(y)
            ylabels.append("  %s → %s" % (name, cn))
            lo, hi = b - Z * se, b + Z * se
            mfc = TCOLOR[tt] if not ovv else "white"
            ax.errorbar(b, y, xerr=[[b - lo], [hi - b]], fmt="o", ms=6.2,
                        color=TCOLOR[tt], ecolor=TCOLOR[tt], elinewidth=1.3, capsize=2.6,
                        markerfacecolor=mfc, markeredgecolor=TCOLOR[tt], markeredgewidth=1.4,
                        zorder=3)
            lab = "p=%.2g" % pv
            ax.text(hi + 0.012, y, lab, va="center", fontsize=7.6, color="#333333")
            cur -= 1
        cur -= 0.7
    ax.axvline(0.0, color="#888888", lw=1.0, ls="--", zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_xlabel("β (per doubling of vascular-disease odds) — IVW fixed, 95% CI", fontsize=10)
    ax.set_title("Reverse MR: vascular disease → thyroid function (Step 6)\n"
                 "open marker = overlapping cohorts (reference only)", fontsize=11.5, pad=10)
    handles = [Line2D([], [], marker="o", ls="", color=TCOLOR[t], ms=7, label=cn)
               for t, cn in TRAITS]
    ax.legend(handles=handles, fontsize=8.4, loc="lower right", framealpha=0.95)
    p = os.path.join(figdir("step6"), "rev_mr_forest.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  " + p, os.path.getsize(p), "B")


# ============================================================ 图 3 · Step 7 MVMR：单变量 vs 多变量
MVM_ORDER = [("cad_nouk", "CAD (Nikpay 2015)"),
             ("mvp_meta", "PAD (MVP meta)"),
             ("is_sv", "IS-SV (MEGASTROKE)"),
             ("is_laa", "IS-LAA (MEGASTROKE)")]


def fig3_mvmr():
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    ys, yl = [], []
    cur = 0.0
    for key, name in MVM_ORDER:
        p = os.path.join(ROOT, "step7", "results", "mvmr_%s.json" % key)
        if not os.path.exists(p):
            continue
        d = J(p)
        for t, cn in TRAITS:
            y = cur
            ys.append(y)
            yl.append("  %s | %s" % (name, cn))
            if t in d.get("mvmr", {}):
                b = d["mvmr"][t]["beta"]
                se = d["mvmr"][t]["se"]
                ax.errorbar(b, y + 0.17, xerr=[[Z * se], [Z * se]], fmt="D", ms=5.6,
                            color=TCOLOR[t], ecolor=TCOLOR[t], elinewidth=1.2, capsize=2.4,
                            markeredgecolor="black", markeredgewidth=0.6, zorder=3)
            uv = d.get("univariable", {}).get(t, {})
            if uv and uv.get("beta") is not None:
                b2 = uv["beta"]
                ax.plot([b2], [y - 0.17], marker="o", ms=5.6, color="white",
                        markeredgecolor=TCOLOR[t], markeredgewidth=1.4, ls="", zorder=3)
            cur -= 1.0
        cur -= 0.8
    ax.axvline(0.0, color="#888888", lw=1.0, ls="--", zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(yl, fontsize=8.8)
    ax.set_xlabel("β (per doubling of genetic liability) — MV-IVW vs univariable IVW", fontsize=10)
    ax.set_title("MVMR (Step 7): hypothyroidism effect independent of TSH/FT4\n"
                 "◆ filled = multivariable (MV-IVW)    ○ open = univariable (for comparison)",
                 fontsize=11.5, pad=10)
    handles = [Line2D([], [], marker="s", ls="", color=TCOLOR[t], ms=7, label=cn)
               for t, cn in TRAITS]
    ax.legend(handles=handles, fontsize=8.4, loc="lower right", framealpha=0.95)
    p = os.path.join(figdir("step7"), "mvmr_single_vs_multi.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  " + p, os.path.getsize(p), "B")


# ============================================================ 图 4 · Step 8 共定位 H4/H3 热图
def fig4_coloc():
    p = os.path.join(ROOT, "step8", "results", "coloc_results.tsv")
    rows = list(csv.DictReader(io.open(p, encoding="utf-8"), delimiter="\t"))
    loci = []
    outs = []
    for r in rows:
        if r["region"] not in loci:
            loci.append(r["region"])
        if r["outcome"] not in outs:
            outs.append(r["outcome"])
    OURLAB = {"cad_nouk": "CAD\n(Nikpay)", "mvp_meta": "PAD\n(MVP)",
              "is_sv": "IS-SV\n(MEGASTROKE)", "is_laa": "IS-LAA\n(MEGASTROKE)"}
    LOCLAB = {}
    for r in rows:
        LOCLAB.setdefault(r["region"], r["rs"])
    import numpy as np
    H4 = np.full((len(loci), len(outs)), float("nan"))
    H3 = np.full((len(loci), len(outs)), float("nan"))
    for r in rows:
        i = loci.index(r["region"])
        j = outs.index(r["outcome"])
        H4[i, j] = float(r["H4"])
        H3[i, j] = float(r["H3"])
    fig, ax = plt.subplots(figsize=(6.6, 7.2))
    im = ax.imshow(H4, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(outs)))
    ax.set_xticklabels([OURLAB.get(o, o) for o in outs], fontsize=8.8)
    ax.set_yticks(range(len(loci)))
    ax.set_yticklabels(["%s (%s)" % (l, LOCLAB.get(l, "")) for l in loci], fontsize=8.8)
    for i in range(len(loci)):
        for j in range(len(outs)):
            v = H4[i, j]
            h3 = H3[i, j]
            if math.isnan(v):
                continue
            txt = "%.2f" % v
            if h3 > v:
                txt = "%.2f\n(H3 %.2f)" % (v, h3)
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.4,
                    color="black" if v < 0.6 else "white")
    ax.grid(False)
    ax.set_title("Colocalisation (Step 8, coloc.abf): PP.H4 across loci and outcomes",
                 fontsize=11.5, pad=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("PP.H4 (shared causal variant)", fontsize=9)
    p = os.path.join(figdir("step8"), "coloc_prob_heatmap.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  " + p, os.path.getsize(p), "B")


# ============================================================ 图 5 · Step 9 MDE80 OR
def _load_forward():
    p = os.path.join(ROOT, "step9", "results", "power_forward.tsv")
    return list(csv.DictReader(io.open(p, encoding="utf-8"), delimiter="\t"))


def fig5_power_mde():
    rows = _load_forward()
    order = ["mvp_meta", "finngen", "sakaue", "cad_nouk", "is_sv", "is_laa", "cad_ukb", "mi_ukb"]
    LAB = {"mvp_meta": "PAD(MVP)", "finngen": "PAD(FinnGen)", "sakaue": "PAD(Sakaue)",
           "cad_nouk": "CAD(Nikpay)", "is_sv": "IS-SV", "is_laa": "IS-LAA",
           "cad_ukb": "CAD(UKB)", "mi_ukb": "MI(UKB) ‡"}
    import numpy as np
    keys = [k for k in order if any(r["out_key"] == k for r in rows)]
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    x = np.arange(len(keys))
    w = 0.26
    for i, (t, cn) in enumerate(TRAITS):
        vals = []
        for k in keys:
            m = [r for r in rows if r["out_key"] == k and r["exposure_key"] == t]
            if not m or not m[0].get("mde80_or"):
                vals.append(float("nan"))
            else:
                vals.append(float(m[0]["mde80_or"]))
        hs = ["//" if k == "mi_ukb" else None for k in keys]  # MI 行加斜线提示尺度异常
        ax.bar(x + (i - 1) * w, vals, w, label=cn, color=TCOLOR[t],
               edgecolor="black", linewidth=0.5, hatch=hs)
        for xi, v in zip(x + (i - 1) * w, vals):
            if not math.isnan(v):
                ax.text(xi, v + 0.004, "%.3f" % v, ha="center", fontsize=6.6, rotation=90)
    ax.axhline(1.05, color="#666666", lw=1.1, ls="--")
    ax.text(len(keys) - 0.45, 1.051, "OR 1.05", fontsize=8, color="#666666", va="bottom", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels([LAB.get(k, k) for k in keys], fontsize=9)
    ax.set_ylabel("MDE at 80% power (OR per unit exposure)", fontsize=10)
    ax.set_title("Minimum detectable effect (80% power, α=0.05 two-sided) — forward MR\n"
                 "† MI(UKB) excluded from interpretation (linear-probability scale)",
                 fontsize=11.5, pad=10)
    ax.legend(fontsize=8.6, ncol=3)
    ax.set_ylim(0.99, 1.34)
    p = os.path.join(figdir("step9"), "power_mde80_or.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  " + p, os.path.getsize(p), "B")


# ============================================================ 图 6 · Step 9 功效曲线
def fig6_power_curve():
    rows = _load_forward()
    # 只用无样本重叠的结局（与 Step 9 结论口径一致）
    keys = [("mvp_meta", "PAD (MVP meta)"), ("cad_nouk", "CAD (Nikpay 2015)"),
            ("is_laa", "IS-LAA (MEGASTROKE)"), ("is_sv", "IS-SV (MEGASTROKE)"),
            ("finngen", "PAD (FinnGen R13)")]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4), sharey=True)
    xs = [1.00 + i * 0.005 for i in range(0, 61)]  # OR 1.00 → 1.30
    cmap = plt.get_cmap("viridis")
    for ax, (t, cn) in zip(axes, TRAITS):
        for i, (k, lab) in enumerate(keys):
            m = [r for r in rows if r["out_key"] == k and r["exposure_key"] == t]
            if not m:
                continue
            se = float(m[0]["se_fixed"])
            ys = [power(math.log(o), se) for o in xs]
            ax.plot(xs, ys, lw=1.8, color=cmap(i / max(1, len(keys) - 1)), label=lab)
        ax.axhline(0.8, color="#666666", ls="--", lw=1.1)
        ax.axvline(1.05, color="#bbbbbb", ls=":", lw=1.0)
        ax.set_title(cn, fontsize=11)
        ax.set_xlabel("True OR (per unit exposure)", fontsize=9.5)
        ax.set_ylim(0, 1.02)
        ax.set_xlim(1.0, 1.30)
    axes[0].set_ylabel("Power (α=0.05, two-sided)", fontsize=9.5)
    axes[0].text(1.004, 0.82, "80%", fontsize=8, color="#666666")
    axes[0].text(1.052, 0.03, "OR 1.05", fontsize=7.6, color="#999999", rotation=90)
    h, l = axes[-1].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8.8, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("Statistical power vs true effect size (non-overlapping outcomes only)",
                 fontsize=12.5, y=1.02)
    p = os.path.join(figdir("step9"), "power_curve.png")
    fig.savefig(p)
    plt.close(fig)
    print("OK  " + p, os.path.getsize(p), "B")


if __name__ == "__main__":
    fig1_panvascular()
    fig2_reverse()
    fig3_mvmr()
    fig4_coloc()
    fig5_power_mde()
    fig6_power_curve()
    print("ALL DONE")
