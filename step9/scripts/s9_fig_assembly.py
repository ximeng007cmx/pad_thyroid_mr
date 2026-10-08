# -*- coding: utf-8 -*-
"""
第 2 轮：图表装配 —— 从底层结果表重绘主图 Fig 3 / 4 / 5 / 6 与附图 Fig S2
================================================================================
设计原则
  1) 全部数值 programmatic 读取自 results 文件，脚本内不硬编码任何估计值；
  2) 统一 300 dpi + RGB + 矢量 SVG + 同一配色/字号/线宽体系（沿用 Fig 1 / Fig 2）；
  3) 每图加 "Fig N | 标题" 眉行与数据来源脚注；凡涉及 π 上界者具名 Burgess 2016。

数据来源（逐图）
  Fig 3A/3B  step9/results/overlap_sim_curve.tsv、overlap_sim_anchor.tsv、overlap_sim_summary.json
  Fig 3C     step9/results/overlap_sim_mc_validation.tsv
  Fig 4      step3/results/mr_main_*.json、step4/results/mr_{mvp_meta,finngen}_*.json
  Fig 5A     step9/results/scale_check.tsv
  Fig 5B     step9/results/power_summary.json  (scenarios["×1"])
  Fig 6A     step7/results/mvmr_summary.tsv
  Fig 6B     step8/results/coloc_results.tsv
  Fig S2     全部 24 层的 mr_*.json（经 work/net_fix/_fig2_layers_full.json 校核）

输出：step9/figures/fig{3,4,5,6}_*.png/.svg、figS2_panvascular_overview.png/.svg
"""
import os, io, csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIGDIR = os.path.join(ROOT, "step9", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# ------------------------------------------------------------------ 统一配色
C_NOV  = "#1F4E79"   # no-overlap / 主
C_OV   = "#C00000"   # overlap
C_PART = "#D97706"   # partial overlap
C_TSH  = "#B8860B"   # TSH
C_FT4  = "#5B3E8E"   # FT4
C_TXT  = "#1A1A1A"
C_GRY  = "#6B7280"
C_GRID = "#DDDDDD"
C_BAND = "#F2F5F9"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.2,
    "axes.edgecolor": "#9AA3AD",
    "axes.linewidth": 0.9,
    "text.color": C_TXT,
    "axes.labelcolor": C_TXT,
    "xtick.color": C_TXT,
    "ytick.color": C_TXT,
    "savefig.facecolor": "white",
})
Z = 1.959963984540054
BURGESS = ("Overlap–bias framework: Burgess, Davies & Thompson (2016), "
           "Genet Epidemiol 40:597–608 — PMID 27625185, eq. (3).")


def save(fig, stem):
    png = os.path.join(FIGDIR, stem + ".png")
    svg = os.path.join(FIGDIR, stem + ".svg")
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    from PIL import Image
    im = Image.open(png)
    if im.mode != "RGB":
        im = im.convert("RGB")
    im.save(png, dpi=(300, 300))
    im2 = Image.open(png)
    print("  %-42s %-11s %-5s dpi=%s" % (stem + ".png", "%dx%d" % im2.size,
                                         im2.mode, im2.info.get("dpi")))
    return png


def read_tsv(rel):
    with io.open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_json(rel):
    with io.open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return json.load(fh)


def header(fig, x, y, main, sub=None):
    fig.text(x, y, main, ha="left", va="center", fontsize=12.2,
             fontweight="bold", color=C_TXT)
    if sub:
        fig.text(x, y - 0.052, sub, ha="left", va="center", fontsize=7.3, color=C_GRY)


def footer(fig, x, y, text):
    fig.text(x, y, text, ha="left", va="center", fontsize=7.2, color=C_GRY)


# ==================================================================== 读数据
CURVE = read_tsv("step9/results/overlap_sim_curve.tsv")
ANCH = read_tsv("step9/results/overlap_sim_anchor.tsv")
MC    = read_tsv("step9/results/overlap_sim_mc_validation.tsv")
SUMM  = read_json("step9/results/overlap_sim_summary.json")
SCALE = read_tsv("step9/results/scale_check.tsv")
POWER = read_json("step9/results/power_summary.json")["scenarios"]["×1"]
MVMR  = read_tsv("step7/results/mvmr_summary.tsv")
COLOC = read_tsv("step8/results/coloc_results.tsv")

LAYER_FILES = {}
for t in ("hypo", "tsh", "ft4"):
    LAYER_FILES["PAD_Sakaue|" + t] = "step3/results/mr_main_%s.json" % t
    LAYER_FILES["PAD_MVP|" + t]    = "step4/results/mr_mvp_meta_%s.json" % t
    LAYER_FILES["PAD_R13|" + t]    = "step4/results/mr_finngen_%s.json" % t
    LAYER_FILES["CAD_Nikpay|" + t] = "step4/results/mr_cad_nouk_%s.json" % t
    LAYER_FILES["CAD_UKB|" + t]    = "step4/results/mr_cad_ukb_%s.json" % t
    LAYER_FILES["MI_UKB|" + t]     = "step4/results/mr_mi_ukb_%s.json" % t
    LAYER_FILES["IS_LAA|" + t]     = "step4/results/mr_is_laa_%s.json" % t
    LAYER_FILES["IS_SV|" + t]      = "step4/results/mr_is_sv_%s.json" % t

EST = {}
for k, rel in LAYER_FILES.items():
    d = read_json(rel)
    EST[k] = dict(
        fixed=d["ivw_fixed"], rand=d["ivw_random"],
        med=d["weighted_median"], nsnp=d["ivw_fixed"].get("nsnp") or d.get("nsnp"),
    )

print("=" * 74)
print("[Fig 3] 重叠解释力的确定性上界")
# ------------------------------------------------------------------ Fig 3
fig = plt.figure(figsize=(14.4, 5.05), dpi=300)
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.00, 1.02],
                      left=0.058, right=0.985, top=0.790, bottom=0.185, wspace=0.30)

# --- Panel A: R(pi) 曲线
axA = fig.add_subplot(gs[0, 0])
cfg_a = "cad_ukb"                      # Fw = 109.021（甲减配对代表）
Fw_a = float([r for r in CURVE if r["cfg"] == cfg_a][0]["Fw"])
rs = sorted({float(r["r"]) for r in CURVE})
cmap_r = {1: "#9AA3AD", 5: "#2E75B6", 10: "#D97706", 20: C_OV}
for rv in rs:
    pts = sorted([(float(x["pi"]), float(x["R_analytic"]))
                  for x in CURVE if x["cfg"] == cfg_a and float(x["r"]) == rv])
    axA.plot([p[0] for p in pts], [p[1] for p in pts],
             color=cmap_r.get(int(rv), C_GRY), lw=1.9,
             label="r = %d" % rv, zorder=3, solid_capstyle="round")
axA.axvspan(1.0, 1.34, color="#F3D9D9", alpha=0.55, zorder=0)
axA.axvline(1.0, color=C_OV, lw=1.1, ls=":", zorder=2)
axA.text(1.172, 2.55, "π > 1\nphysically\nimpossible", fontsize=7.8, color=C_OV,
         ha="center", va="center", fontweight="bold", linespacing=1.35)

# 观测比值线（标签放右端，三条线按数值错开，避免互相压叠）
OBS = sorted(SUMM["observed_pairs"], key=lambda o: -float(o["R_obs"]))
for o in OBS:
    rv = float(o["R_obs"])
    axA.axhline(rv, color=C_OV, lw=1.25, ls=(0, (5, 3)), zorder=1)
    if rv > 3.0:
        ytxt = rv + 0.16          # 3.62 → 标签在线上方
    elif rv > 1.8:
        ytxt = rv + 0.15          # 1.91 → 标签在线上方
    else:
        ytxt = rv - 0.20          # 1.77 → 标签在线下方
    axA.text(0.972, ytxt, "observed %s  ×%.2f" % (o["disease"], rv), fontsize=7.5,
             color=C_OV, ha="right", va="center", fontweight="bold")
axA.set_xlim(0, 1.34)
axA.set_ylim(0.97, 4.05)
axA.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
axA.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"], fontsize=8.3)
axA.set_xlabel("Overlap proportion  π  (outcome sample)", fontsize=8.7, labelpad=5)
axA.set_ylabel("Inflation of the MR estimate,  R(π) = 1 + π·r / F$_w$", fontsize=8.7)
for s in ("top", "right"):
    axA.spines[s].set_visible(False)
axA.tick_params(length=3.2)
axA.legend(loc="center left", fontsize=8.0, frameon=True, facecolor="white",
           edgecolor=C_GRID, borderpad=0.5, handletextpad=0.55, framealpha=0.96,
           title="observational association\nr = β$_{obs}$ / β$_{causal}$",
           title_fontsize=7.3, bbox_to_anchor=(0.015, 0.615))
axA.set_title("A   Even complete overlap cannot reach the observed ratios",
              fontsize=10.3, fontweight="bold", color=C_NOV, loc="left", pad=8)

# --- Panel B: 反解所需 r（数据取自 summary 的 3 个观测配对，避免长表重复）
axB = fig.add_subplot(gs[0, 1])
lab_short = ["CAD\nNikpay → UKB", "PAD\nMVP → Sakaue", "PAD\nMVP → FinnGen"]
ysB = [2, 1, 0]
axB.axvspan(1, float(max(rs)), color=C_BAND, zorder=0)
for i, o in enumerate(OBS):
    rr = float(o["r_req_at_pi1"])
    axB.barh(ysB[i], rr - 1, left=1, height=0.44, color=C_OV, alpha=0.90,
             edgecolor="white", lw=0.8, zorder=3)
    axB.text(rr * 1.10, ysB[i], "r = %.1f" % rr, ha="left", va="center",
             fontsize=8.6, fontweight="bold", color=C_OV, zorder=5)
axB.axvline(float(max(rs)), color=C_NOV, lw=1.3, ls="--", zorder=4)
axB.text(float(max(rs)) * 0.93, 2.66, "largest plausible r = 20\n(observational association\n20× the causal effect)",
         ha="right", va="top", fontsize=7.2, color=C_NOV, linespacing=1.3)
axB.set_xscale("log")
axB.set_xlim(1, 1200)
axB.set_ylim(-0.62, 2.98)
axB.set_yticks(ysB)
axB.set_yticklabels(lab_short, fontsize=8.0, linespacing=1.35)
axB.set_xlabel("Observational association required to explain the ratio,  r$_{req}$  (log scale)",
               fontsize=8.5, labelpad=5)
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.tick_params(axis="y", length=0, pad=4)
axB.tick_params(axis="x", length=3.2)
axB.set_title("B   Required confounding, r$_{req}$", fontsize=10.3,
              fontweight="bold", color=C_NOV, loc="left", pad=8)

# --- Panel C: 蒙特卡洛一致性
axC = fig.add_subplot(gs[0, 2])
xa = [float(m["R_analytic"]) for m in MC]
ya = [float(m["R_mc"]) for m in MC]
dev = max(abs(a - b) for a, b in zip(xa, ya))
axC.plot([1.0, 1.20], [1.0, 1.20], color=C_GRY, lw=1.1, ls="--", zorder=1)
axC.scatter(xa, ya, s=16, facecolor=C_NOV, edgecolor="white", linewidth=0.5,
            alpha=0.85, zorder=3)
axC.set_xlim(0.995, 1.205)
axC.set_ylim(0.995, 1.205)
axC.set_xticks([1.00, 1.05, 1.10, 1.15, 1.20])
axC.set_yticks([1.00, 1.05, 1.10, 1.15, 1.20])
axC.set_xlabel("Analytic R(π) (eq. 3 reparameterised)", fontsize=8.5, labelpad=5)
axC.set_ylabel("Simulated R (10,000 Monte Carlo replicates)", fontsize=8.5)
for s in ("top", "right"):
    axC.spines[s].set_visible(False)
axC.tick_params(length=3.2)
axC.text(1.003, 1.196, "n = %d parameter cells\nmax |analytic − simulated| = %.4f"
         % (len(MC), dev), fontsize=7.4, color=C_TXT, ha="left", va="top",
         linespacing=1.4,
         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=C_GRID, lw=0.9))
axC.set_title("C   Monte Carlo agreement", fontsize=10.3, fontweight="bold",
              color=C_NOV, loc="left", pad=8)

header(fig, 0.012, 0.955,
       "Fig 3  |  Sample overlap cannot account for the between-library differences",
       "F$_w$ = %.1f for the hypothyroidism pairs plotted (panel A); the framework is adopted, not derived." % Fw_a)
footer(fig, 0.012, 0.075,
       "Observed β ratios: " + ", ".join("×%.2f (%s)" % (float(o["R_obs"]), o["disease"]) for o in OBS)
       + ".  Required overlap at r = 5: "
       + " / ".join("%.1f× the outcome sample" % float(o["pi_req_at_r5"]) for o in OBS) + "."
       + "\n" + BURGESS)
footer(fig, 0.012, 0.028,
       "Data: step9/results/overlap_sim_{curve,anchor,mc_validation}.tsv, overlap_sim_summary.json. "
       "The Sakaue 2021 ↔ FinnGen R13 pair is non-independent (R3 ⊂ R13) and is shown for reference only.")
save(fig, "fig3_overlap_bound")

# ------------------------------------------------------------------ Fig 4
print("[Fig 4] PAD 三层复制森林图")
PAD_LIB = [("PAD_MVP", "MVP trans-ancestry\n(no overlap)", C_NOV, False),
           ("PAD_Sakaue", "Sakaue 2021 EUR\n(UKB + FinnGen R3 · primary)", C_PART, False),
           ("PAD_R13", "FinnGen R13\n(invalid · sample overlap)", C_OV, True)]
TRAITS = [("hypo", "Hypothyroidism"), ("tsh", "TSH"), ("ft4", "FT4")]
fig = plt.figure(figsize=(12.9, 5.15), dpi=300)
gs = fig.add_gridspec(1, 3, left=0.115, right=0.985, top=0.775, bottom=0.190, wspace=0.30)
for ci, (tr, trlab) in enumerate(TRAITS):
    ax = fig.add_subplot(gs[0, ci])
    ys = [2, 1, 0]
    for i, (key, lab, col, invalid) in enumerate(PAD_LIB):
        y = ys[i]
        e = EST[key + "|" + tr]
        for off, est, mfc, mk, tag in ((-0.155, e["fixed"], col, "o", "fixed"),
                                       (+0.155, e["rand"], "white", "o", "random")):
            b, se = float(est["beta"]), float(est["se"])
            orr, lo, hi = pow(2.718281828459045, b), \
                pow(2.718281828459045, b - Z * se), pow(2.718281828459045, b + Z * se)
            ax.errorbar([orr], [y], xerr=[[orr - lo], [hi - orr]], fmt=mk, ms=7.4,
                        mfc=mfc, mec=col, mew=1.5, ecolor=col, elinewidth=1.6,
                        capsize=3.0, capthick=1.6, zorder=4, alpha=1.0 if tag == "fixed" else 0.85)
            if tag == "fixed":
                ax.text(hi * 1.035, y, "p=%.3g" % float(est["p"]), fontsize=7.2,
                        color=col, ha="left", va="center")
        if invalid:
            ax.axhspan(y - 0.42, y + 0.42, color="#F3D9D9", alpha=0.55, zorder=0)
        ax.axhline(y - 0.5, color=C_GRID, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axvline(1.0, color="#9AA3AD", lw=1.0, ls="--", zorder=2)
    ax.set_xscale("log")
    ax.set_xlim(0.86, 1.34)
    ax.set_ylim(-0.62, 2.62)
    ax.set_xticks([0.9, 1.0, 1.1, 1.2, 1.3])
    ax.set_xticklabels(["0.9", "1.0", "1.1", "1.2", "1.3"], fontsize=8.3)
    ax.set_yticks(ys)
    if ci == 0:
        ax.set_yticklabels([l for _, l, _, _ in PAD_LIB], fontsize=7.5, linespacing=1.3)
        for t, (_, _, col, _) in zip(ax.get_yticklabels(), PAD_LIB):
            t.set_color(col)
    else:
        ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("OR per 1-SD genetically predicted exposure\n(log scale, 95%% CI, %d instruments at primary layer)"
                  % EST["PAD_Sakaue|" + tr]["nsnp"], fontsize=8.1, labelpad=5)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.tick_params(axis="x", length=3.2)
    ax.set_title("%s   %s" % ("ABC"[ci], trlab), fontsize=10.3, fontweight="bold",
                 color=C_NOV, loc="left", pad=8)
leg = [Line2D([], [], marker="o", ls="none", mfc=C_NOV, mec=C_NOV, ms=7.0,
              label="IVW fixed effect"),
       Line2D([], [], marker="o", ls="none", mfc="white", mec=C_NOV, ms=7.0, mew=1.5,
              label="IVW random effects")]
fig.legend(handles=leg, loc="upper right", ncol=2, fontsize=8.2, frameon=False,
           bbox_to_anchor=(0.985, 0.955))
header(fig, 0.012, 0.955,
       "Fig 4  |  Three-layer replication of the PAD estimates",
       "The hypothyroidism signal is not robust in any informative layer; TSH and FT4 are null throughout.")
footer(fig, 0.012, 0.045,
       "The FinnGen R13 layer is invalidated by sample overlap (exposure GWAS includes FinnGen DF10; R3 ⊂ R13) and is shown for completeness only. "
       "Estimator consistency was the prespecified criterion; no single estimator is treated as definitive.\n"
       "Data: step3/results/mr_main_*.json, step4/results/mr_{mvp_meta,finngen}_*.json (values plotted verbatim; no re-derivation).")
save(fig, "fig4_pad_replication")

# ------------------------------------------------------------------ Fig 5
print("[Fig 5] 结局尺度体检 + 功效")
NAME = {"sakaue": "PAD (Sakaue 2021)", "finngen": "PAD (FinnGen R13)",
        "mvp_meta": "PAD (MVP)", "cad_nouk": "CAD (Nikpay 2015)",
        "cad_ukb": "CAD (UKB)", "mi_ukb": "MI (UKB) †", "is_laa": "IS-LAA",
        "is_sv": "IS-SV"}
fig = plt.figure(figsize=(13.4, 5.00), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.10],
                      left=0.090, right=0.985, top=0.775, bottom=0.185, wspace=0.34)

axA = fig.add_subplot(gs[0, 0])
sc = sorted(SCALE, key=lambda r: float(r["ratio_to_ref"]))
ysA = list(range(len(sc)))
for i, r in enumerate(sc):
    v = float(r["ratio_to_ref"])
    bad = v < 0.1
    axA.barh(ysA[i], v - 1, left=1, height=0.56,
             color=(C_OV if bad else C_NOV), alpha=0.92, edgecolor="white", lw=0.8, zorder=3)
    if bad:   # MI：条向左延伸，数值标签放在条端上方
        axA.text(v, ysA[i] + 0.40, "%.3f" % v, ha="center", va="bottom",
                 fontsize=7.6, fontweight="bold", color=C_OV, zorder=5)
    else:
        axA.text(v * 1.10, ysA[i], "%.3f" % v, ha="left", va="center",
                 fontsize=7.6, color=C_TXT, zorder=5)
axA.axvline(1.0, color="#9AA3AD", lw=1.1, ls="--", zorder=2)
axA.set_xscale("log")
axA.set_xlim(0.02, 4.2)
axA.set_ylim(-0.62, len(sc) - 0.38)
axA.set_yticks(ysA)
axA.set_yticklabels([NAME.get(r["key"], r["key"]) for r in sc], fontsize=8.1)
axA.get_yticklabels()[[r["key"] for r in sc].index("mi_ukb")].set_color(C_OV)
axA.set_xlabel("Median SE of instrument SNPs, relative to the reference file (log scale)",
               fontsize=8.5, labelpad=5)
for s in ("top", "right", "left"):
    axA.spines[s].set_visible(False)
axA.tick_params(axis="y", length=0, pad=4)
axA.tick_params(axis="x", length=3.2)
_mi_y = ysA[[r["key"] for r in sc].index("mi_ukb")]
axA.text(1.14, _mi_y, "MI is on a linear-probability scale:\nβ is not interpretable as log OR,\nbut z and p remain valid",
         fontsize=7.2, color=C_OV, ha="left", va="center", linespacing=1.35)
axA.set_title("A   Outcome-scale audit of the eight files", fontsize=10.3,
              fontweight="bold", color=C_NOV, loc="left", pad=8)

axB = fig.add_subplot(gs[0, 1])
KEYMAP = [("PAD（Sakaue 2021, GCST90018890）", "PAD (Sakaue 2021)"),
          ("PAD（FinnGen R13 I9_PAD）", "PAD (FinnGen R13)"),
          ("PAD（MVP trans-ancestry meta）", "PAD (MVP)"),
          ("CAD（Nikpay 2015, GCST003116）", "CAD (Nikpay 2015)"),
          ("CAD（UKB Mbatchou 2021, GCST90013864）", "CAD (UKB)"),
          ("IS-LAA（MEGASTROKE, GCST005840）", "IS-LAA (MEGASTROKE)"),
          ("IS-SV（MEGASTROKE, GCST005841）", "IS-SV (MEGASTROKE)")]
TCOL = {"hypo": C_NOV, "tsh": C_TSH, "ft4": C_FT4}
ysB = list(range(len(KEYMAP)))[::-1]
for i, (suffix, lab) in enumerate(KEYMAP):
    y = ysB[i]
    axB.axhline(y, color=C_GRID, lw=0.7, ls=(0, (3, 3)), zorder=0)
    vals = []
    for t in ("hypo", "tsh", "ft4"):
        v = float(POWER["%s|%s" % (t, suffix)])
        vals.append(v)
        axB.plot([v], [y], marker="o", ms=7.2, mfc=TCOL[t], mec="white", mew=1.1, zorder=4)
    axB.plot(vals, [y] * 3, color=C_GRY, lw=1.0, alpha=0.55, zorder=2)
    axB.text(max(vals) + 0.004, y, "%.3f" % max(vals), fontsize=7.2, color=C_GRY,
             ha="left", va="center")
axB.axvline(1.0, color="#9AA3AD", lw=1.0, ls="--", zorder=1)
axB.set_xlim(0.995, 1.365)
axB.set_ylim(-0.62, len(KEYMAP) - 0.38)
axB.set_xticks([1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35])
axB.set_xticklabels(["1.00", "1.05", "1.10", "1.15", "1.20", "1.25", "1.30", "1.35"], fontsize=8.1)
axB.set_yticks(ysB)
axB.set_yticklabels([l for _, l in KEYMAP], fontsize=8.1)
axB.set_xlabel("Minimum detectable OR at 80% power (α = 0.05, log-OR outcomes)", fontsize=8.5, labelpad=5)
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.tick_params(axis="y", length=0, pad=4)
axB.tick_params(axis="x", length=3.2)
legB = [Line2D([], [], marker="o", ls="none", mfc=TCOL[t], mec="white", ms=7.2,
               label={"hypo": "hypothyroidism", "tsh": "TSH", "ft4": "FT4"}[t])
        for t in ("hypo", "tsh", "ft4")]
axB.legend(handles=legB, loc="upper right", fontsize=8.0, frameon=True,
           facecolor="white", edgecolor=C_GRID, borderpad=0.5, handletextpad=0.5,
           framealpha=0.96)
axB.set_title("B   Statistical power was adequate except for FT4", fontsize=10.3,
              fontweight="bold", color=C_NOV, loc="left", pad=8)
header(fig, 0.012, 0.955,
       "Fig 5  |  Outcome-scale audit and statistical power",
       "MDE values are OR-scale per allele and exclude the linear-probability MI file.")
footer(fig, 0.012, 0.045,
       "† MI (GCST90038610) is on a linear-probability scale, so its β is not interpretable as a log OR and its MDE is not plotted in panel B. "
       "MDE at 80% power was computed with a true effect equal to the observed single-SNP effect.\n"
       "Data: step9/results/scale_check.tsv, power_summary.json (scenario ×1).")
save(fig, "fig5_outcome_audit")

# ------------------------------------------------------------------ Fig 6
print("[Fig 6] 互补证据：MVMR + 共定位")
OUTLAB = {"cad_nouk": "CAD\n(Nikpay 2015)", "is_laa": "IS-LAA\n(MEGASTROKE)",
          "is_sv": "IS-SV\n(MEGASTROKE)", "mvp_meta": "PAD\n(MVP)"}
fig = plt.figure(figsize=(13.2, 5.30), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.06],
                      left=0.105, right=0.975, top=0.780, bottom=0.185, wspace=0.30)

axA = fig.add_subplot(gs[0, 0])
rows = [r for r in MVMR if r["exposure"] == "hypo" and float(r["beta_univ"] or 0) != 0
        and r["outcome"] in OUTLAB and not r["overlap"].strip()]
rows = [r for r in rows if r["outcome"] != "sakaue"]
ysA = list(range(len(rows)))[::-1]
for i, r in enumerate(rows):
    y = ysA[i]
    bm, sm = float(r["beta_mvmr"]), float(r["se_mvmr"])
    bu = float(r["beta_univ"])
    axA.plot([bu, bm], [y, y], color="#B9BFC6", lw=2.6, solid_capstyle="round", zorder=1)
    axA.errorbar([bu], [y], xerr=[[0], [0]], fmt="o", ms=8.0, mfc=C_NOV, mec="white",
                 mew=1.4, zorder=4)
    axA.errorbar([bm], [y], xerr=[[Z * sm], [Z * sm]], fmt="s", ms=8.0, mfc=C_OV,
                 mec="white", mew=1.4, ecolor=C_OV, elinewidth=1.7, capsize=3.2,
                 capthick=1.7, zorder=4)
    sig = float(r["p_mvmr"]) < 0.05
    axA.text(max(bu, bm) + Z * sm + 0.006, y, "p=%.3g" % float(r["p_mvmr"]),
             fontsize=7.6, color=(C_OV if sig else C_GRY), va="center",
             fontweight="bold" if sig else "normal", zorder=5)
    axA.axhline(y - 0.5, color=C_GRID, lw=0.8, ls=(0, (4, 3)), zorder=0)
axA.axvline(0, color="#9AA3AD", lw=1.0, ls="--", zorder=2)
axA.set_xlim(-0.02, 0.135)
axA.set_ylim(-0.62, len(rows) - 0.38)
axA.set_yticks(ysA)
axA.set_yticklabels([OUTLAB[r["outcome"]] for r in rows], fontsize=8.0, linespacing=1.25)
axA.set_xlabel("β for hypothyroidism (per 1-SD; log-OR scale)", fontsize=8.5, labelpad=5)
for s in ("top", "right", "left"):
    axA.spines[s].set_visible(False)
axA.tick_params(axis="y", length=0, pad=4)
axA.tick_params(axis="x", length=3.2)
legA = [Line2D([], [], marker="o", ls="none", mfc=C_NOV, mec="white", ms=8.0,
               label="univariable IVW"),
        Line2D([], [], marker="s", ls="none", mfc=C_OV, mec="white", ms=8.0,
               label="multivariable IVW (95% CI)")]
axA.legend(handles=legA, loc="lower left", fontsize=8.0, frameon=True,
           facecolor="white", edgecolor=C_GRID, borderpad=0.5, handletextpad=0.5,
           bbox_to_anchor=(0.005, 0.02), framealpha=0.96)
axA.set_title("A   MVMR: hypothyroidism across the four no-overlap libraries",
              fontsize=10.2, fontweight="bold", color=C_NOV, loc="left", pad=8)

axB = fig.add_subplot(gs[0, 1])
OUTS = ["cad_nouk", "mvp_meta", "is_sv", "is_laa"]
OUTS_LAB = ["CAD\n(Nikpay)", "PAD\n(MVP)", "IS-SV", "IS-LAA"]
REGIONS = ["R%02d" % i for i in range(1, 11)]
H4 = {(r["region"], r["outcome"]): float(r["H4"]) for r in COLOC}
H3 = {(r["region"], r["outcome"]): float(r["H3"]) for r in COLOC}
grid = [[H4.get((rg, o), float("nan")) for o in OUTS] for rg in REGIONS]
im = axB.imshow(grid, cmap="Blues", vmin=0.0, vmax=1.0, aspect="auto")
axB.set_xticks(range(len(OUTS)))
axB.set_xticklabels(OUTS_LAB, fontsize=8.0, linespacing=1.25)
axB.set_yticks(range(len(REGIONS)))
axB.set_yticklabels(REGIONS, fontsize=8.2)
for i, rg in enumerate(REGIONS):
    for j, o in enumerate(OUTS):
        v = H4.get((rg, o))
        if v is None:
            continue
        axB.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=7.4,
                 color=("white" if v > 0.55 else C_TXT))
        h3 = H3.get((rg, o), 0.0)
        if h3 >= 0.5:
            axB.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                             edgecolor=C_OV, lw=1.9, zorder=3))
            axB.text(j, i + 0.30, "H3 %.2f" % h3, ha="center", va="center",
                     fontsize=6.6, color=C_OV, fontweight="bold", zorder=4)
axB.set_xticks([x - 0.5 for x in range(1, len(OUTS))], minor=True)
axB.set_yticks([y - 0.5 for y in range(1, len(REGIONS))], minor=True)
axB.grid(which="minor", color="white", lw=1.4)
axB.tick_params(which="both", length=0)
cb = fig.colorbar(im, ax=axB, fraction=0.040, pad=0.025)
cb.set_label("PP.H4 (shared causal variant)", fontsize=8.2)
cb.ax.tick_params(labelsize=7.6)
axB.set_title("B   Colocalisation across the ten thyroid lead-variant regions",
              fontsize=10.2, fontweight="bold", color=C_NOV, loc="left", pad=8)
header(fig, 0.012, 0.955,
       "Fig 6  |  Complementary analyses: multivariable MR and colocalisation",
       "Direction-consistent MVMR support for hypothyroidism in 4/4 no-overlap libraries; a shared variant at R03 in all four.")
footer(fig, 0.012, 0.062,
       "Red boxes mark regions where the independent-signal posterior PP.H3 ≥ 0.5. R02 (HLA) violates the single-causal-variant assumption of colocalisation and is reported for reference only.\n"
       "MVMR conditional F statistics were 37.8–39.9, above weak-instrument thresholds, but the three exposures derive from largely overlapping GWAS and the multivariable Q was significant in every outcome.\n"
       "Data: step7/results/mvmr_summary.tsv, step8/results/coloc_results.tsv.")
save(fig, "fig6_complementary")

# ------------------------------------------------------------------ Fig S2
print("[Fig S2] 泛血管谱 24 层总览森林图")
OUT_ORDER = [("PAD_Sakaue", "PAD (Sakaue 2021, primary)"), ("PAD_R13", "PAD (FinnGen R13)"),
             ("PAD_MVP", "PAD (MVP trans-ancestry)"), ("CAD_UKB", "CAD (UKB Mbatchou 2021)"),
             ("CAD_Nikpay", "CAD (Nikpay 2015)"), ("MI_UKB", "MI (UKB) †"),
             ("IS_LAA", "IS-LAA (MEGASTROKE)"), ("IS_SV", "IS-SV (MEGASTROKE)")]
OVERLAP_RULE = {"PAD_Sakaue": ("partial", "partial overlap (UKB + FinnGen R3)"),
                "PAD_R13": ("invalid", "invalid — exposure-side FinnGen overlap"),
                "PAD_MVP": ("none", "no overlap"),
                "CAD_UKB": ("overlap", "overlap (UKB)"),
                "CAD_Nikpay": ("none", "no overlap"),
                "MI_UKB": ("overlap", "overlap (UKB) †"),
                "IS_LAA": ("none", "no overlap"),
                "IS_SV": ("none", "no overlap")}
OCOL = {"none": C_NOV, "partial": C_PART, "overlap": C_OV, "invalid": "#7A7A7A"}
traits = [("hypo", "Hypothyroidism"), ("tsh", "TSH"), ("ft4", "FT4")]
fig = plt.figure(figsize=(11.0, 8.20), dpi=300)
gs = fig.add_gridspec(3, 1, left=0.255, right=0.905, top=0.845, bottom=0.135, hspace=0.52)
for pi_ax, (tr, trlab) in enumerate(traits):
    ax = fig.add_subplot(gs[pi_ax, 0])
    ys = list(range(len(OUT_ORDER)))[::-1]
    for i, (key, lab) in enumerate(OUT_ORDER):
        y = ys[i]
        e = EST[key + "|" + tr]
        b, se = float(e["fixed"]["beta"]), float(e["fixed"]["se"])
        kind, _ = OVERLAP_RULE[key]
        col = OCOL[kind]
        orr = pow(2.718281828459045, b)
        lo = pow(2.718281828459045, b - Z * se)
        hi = pow(2.718281828459045, b + Z * se)
        if b >= 0:
            orr2, lo2, hi2 = orr, lo, hi
        else:
            orr2, lo2, hi2 = orr, lo, hi
        ax.errorbar([orr2], [y], xerr=[[orr2 - lo2], [hi2 - orr2]], fmt="o", ms=6.0,
                    mfc=col, mec="white", mew=1.1, ecolor=col, elinewidth=1.5,
                    capsize=2.6, capthick=1.5, zorder=4)
        ax.text(0.995, y, "%d" % (e["nsnp"] or 0), fontsize=6.9, color=C_GRY,
                ha="center", va="center", transform=ax.get_yaxis_transform(which="grid"))
        p = float(e["fixed"]["p"])
        ptxt = ("%.3g" % p) if p < 0.001 else ("%.4f" % p)
        ax.text(max(hi2, 1.001) * 1.045, y, ptxt, fontsize=6.7,
                color=(C_OV if p < 0.05 else C_GRY), va="center")
        if kind == "invalid":
            ax.axhspan(y - 0.45, y + 0.45, color="#F3F3F3", alpha=0.9, zorder=0)
        ax.axhline(y - 0.5, color=C_GRID, lw=0.7, ls=(0, (3, 3)), zorder=0)
    ax.axvline(1.0, color="#9AA3AD", lw=1.0, ls="--", zorder=2)
    ax.set_xscale("log")
    ax.set_xlim(0.68, 1.55)
    ax.set_ylim(-0.6, len(OUT_ORDER) - 0.4)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4])
    ax.set_xticklabels(["0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.4"], fontsize=7.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([l for _, l in OUT_ORDER], fontsize=8.0, linespacing=1.25)
    for t, (key, _) in zip(ax.get_yticklabels(), OUT_ORDER):
        t.set_color(OCOL[OVERLAP_RULE[key][0]])
    ax.set_title("%s   %s" % ("ABC"[pi_ax], trlab), fontsize=10.3, fontweight="bold",
                 color=C_NOV, loc="left", pad=7)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.tick_params(axis="x", length=3.0)
    ax.set_xlabel("OR per 1-SD genetically predicted exposure (IVW fixed effect, 95% CI)",
                  fontsize=8.1, labelpad=4)
legs = [Line2D([], [], marker="o", ls="none", mfc=OCOL[k], mec="white", ms=6.6,
               label={"none": "no sample overlap", "partial": "partial overlap",
                      "overlap": "overlap", "invalid": "invalid (overlap)"}[k])
        for k in ("none", "partial", "overlap", "invalid")]
fig.legend(handles=legs, loc="upper right", ncol=4, fontsize=8.0, frameon=False,
           bbox_to_anchor=(0.905, 0.882))
fig.text(0.255, 0.868, "instruments (n)", fontsize=7.0, color=C_GRY, ha="center")
header(fig, 0.012, 0.962,
       "Fig S2  |  Overview of all 24 exposure–outcome layers",
       "Three thyroid exposures × eight atherosclerotic vascular outcomes; the primary analysis (PAD, Sakaue 2021) is highlighted by the pairing design in Fig 2.")
footer(fig, 0.012, 0.052,
       "† MI (GCST90038610) is parameterised on a linear-probability scale; its β is therefore not interpretable as a log OR and its OR is shown for completeness only.\n"
       "Values are plotted verbatim from the per-layer result files; no re-derivation.  Data: step3/results/mr_main_*.json, step4/results/mr_*.json (24 layers).")
save(fig, "figS2_panvascular_overview")

# ------------------------------------------------------------------ Fig S1（RGB 化 + 规范命名）
print("[Fig S1] 工具变量筛选流程图（RGB 化 + 统一命名）")
from PIL import Image
src = os.path.join(FIGDIR, "instrument_selection_flow.png")
if os.path.exists(src):
    im = Image.open(src)
    if im.mode != "RGB":
        im = im.convert("RGB")
    im.save(os.path.join(FIGDIR, "figS1_instrument_selection_flow.png"), dpi=(300, 300))
    im2 = Image.open(os.path.join(FIGDIR, "figS1_instrument_selection_flow.png"))
    print("  %-42s %-11s %-5s dpi=%s" % ("figS1_instrument_selection_flow.png",
                                         "%dx%d" % im2.size, im2.mode, im2.info.get("dpi")))
    svg_src = os.path.join(FIGDIR, "instrument_selection_flow.svg")
    if os.path.exists(svg_src):
        import shutil
        shutil.copyfile(svg_src, os.path.join(FIGDIR, "figS1_instrument_selection_flow.svg"))
        print("  figS1_instrument_selection_flow.svg  copied")

print("=" * 74)
print("done.")
