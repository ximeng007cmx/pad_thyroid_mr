# -*- coding: utf-8 -*-
"""
Fig 2 —— 配对差异图 (paired-difference figure)   ★核心图之二
================================================================
Panel A  哑铃图：三组「无重叠 ↔ 有重叠」结局库配对，线上标注 β 比
Panel B  TSH -> CAD 符号翻转（-0.0135 不显著 -> +0.0394 显著）

数据来源（全部已实检，未做任何修约）：
  step3/results/mr_main_hypo.json                 Sakaue  PAD 甲减 = +0.045240 (se 0.018609)
  step4/results/mr_mvp_meta_hypo.json             MVP     PAD 甲减 = +0.023652 (se 0.011499)
  step4/results/mr_finngen_hypo.json              FinnGen PAD 甲减 = +0.085674 (se 0.012438)
  step4/results/mr_cad_nouk_hypo.json             Nikpay  CAD 甲减 = +0.030760 (se 0.012110)
  step4/results/mr_cad_ukb_hypo.json              UKB     CAD 甲减 = +0.054572 (se 0.011149)
  step4/results/mr_cad_nouk_tsh.json              Nikpay  CAD TSH  = -0.013459 (se 0.019931, p 0.500)
  step4/results/mr_cad_ukb_tsh.json               UKB     CAD TSH  = +0.039364 (se 0.018672, p 0.035)

输出：step9/figures/fig2_paired_differences.png (300 dpi, RGB) + .svg
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIGDIR = os.path.join(ROOT, "step9", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# ---------------------------------------------------------------- 配色（沿用 Fig 1）
C_NOV   = "#1F4E79"   # 无重叠结局库
C_OV    = "#C00000"   # 有重叠结局库（主分析 / 新库）
C_TXT   = "#1A1A1A"
C_GRY   = "#6B7280"
C_GRID  = "#DDDDDD"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.2,
    "axes.edgecolor": "#9AA3AD",
    "axes.linewidth": 0.9,
    "text.color": C_TXT,
    "axes.labelcolor": C_TXT,
    "xtick.color": C_TXT,
    "ytick.color": C_TXT,
})

Z = 1.959963984540054   # 95% CI


def beta_ci(beta, se):
    return beta - Z * se, beta + Z * se


# ---------------------------------------------------------------- 数据（甲减暴露）
PAIRS = [
    dict(
        group="CAD", sub="coronary artery disease",
        left_lab="Nikpay 2015\n(no UKB · no overlap)",
        left_beta=0.030760126863534933, left_se=0.012109714593707284, left_p=0.0111,
        right_lab="Mbatchou 2021 (UKB)\n(UKB overlap)",
        right_beta=0.05457231671728457, right_se=0.011148850593942402, right_p=9.84e-07,
    ),
    dict(
        group="PAD", sub="peripheral artery disease",
        left_lab="MVP trans-ancestry\n(no overlap)",
        left_beta=0.02365246843085242, left_se=0.011499008561808424, left_p=0.0397,
        right_lab="Sakaue 2021 (EUR)\n(UKB + FinnGen R3)",
        right_beta=0.04524017335839431, right_se=0.018609046915227814, right_p=0.0151,
    ),
    dict(
        group="PAD", sub="peripheral artery disease",
        left_lab="MVP trans-ancestry\n(no overlap)",
        left_beta=0.02365246843085242, left_se=0.011499008561808424, left_p=0.0397,
        right_lab="FinnGen R13\n(FinnGen overlap †)",
        right_beta=0.08567351327972754, right_se=0.01243758285993042, right_p=5.65e-12,
    ),
]

FLIP = [
    ("Nikpay 2015\n(no UKB · no overlap)", -0.013458882807350223, 0.019931027272598278, 0.500),
    ("Mbatchou 2021 (UKB)\n(UKB overlap)",  0.039364473361787285, 0.018672234803016222, 0.035),
]

# ---------------------------------------------------------------- 画布
fig = plt.figure(figsize=(11.2, 5.30), dpi=300)
gs = fig.add_gridspec(
    1, 2, width_ratios=[1.55, 1.0],
    left=0.135, right=0.980, top=0.835, bottom=0.190,
    wspace=0.36,
)

# ================================================================ Panel A
axA = fig.add_subplot(gs[0, 0])
axA.set_xlim(-0.004, 0.116)
axA.set_ylim(-0.78, 2.78)

ys = [2.0, 1.0, 0.0]
DY = 0.165

ytick_pos, ytick_lab, ytick_col = [], [], []
for i, p in enumerate(PAIRS):
    y = ys[i]
    lb, rb = p["left_beta"], p["right_beta"]
    lci = beta_ci(lb, p["left_se"])
    rci = beta_ci(rb, p["right_se"])

    # 配对连线（哑铃主体）
    axA.plot([lb, rb], [y + DY, y - DY], color="#B9BFC6", lw=2.8,
             solid_capstyle="round", zorder=1)

    # 左右误差棒 + 点
    axA.errorbar([lb], [y + DY], xerr=[[lb - lci[0]], [lci[1] - lb]],
                 fmt="o", ms=8.4, mfc=C_NOV, mec="white", mew=1.5,
                 ecolor=C_NOV, elinewidth=1.8, capsize=3.4, capthick=1.8, zorder=4)
    axA.errorbar([rb], [y - DY], xerr=[[rb - rci[0]], [rci[1] - rb]],
                 fmt="s", ms=8.4, mfc=C_OV, mec="white", mew=1.5,
                 ecolor=C_OV, elinewidth=1.8, capsize=3.4, capthick=1.8, zorder=4)

    # β 比标签（连线上方，带圆角框）
    ratio = rb / lb
    xm = (lb + rb) / 2.0
    axA.annotate(
        "", xy=(rb, y - DY - 0.075), xytext=(lb, y + DY + 0.075),
        arrowprops=dict(arrowstyle="-|>", color=C_GRY, lw=1.0,
                        shrinkA=11, shrinkB=11, connectionstyle="arc3,rad=-0.32"),
        zorder=2,
    )
    axA.text(xm + 0.004, y + 0.345, "×%.2f" % ratio, ha="center", va="center",
             fontsize=11.6, fontweight="bold", color=C_OV,
             bbox=dict(boxstyle="round,pad=0.24", fc="white", ec=C_OV,
                       lw=1.1, alpha=0.96), zorder=6)

    # β 数值标签
    axA.text(lb, y + DY + 0.095, "%.4f" % lb, ha="center", va="bottom",
             fontsize=7.4, color=C_NOV, zorder=6)
    axA.text(rb, y - DY - 0.095, "%.4f" % rb, ha="center", va="top",
             fontsize=7.4, color=C_OV, zorder=6)

    # y 轴刻度标签 = 库名（无重叠 / 有重叠 分行）
    ytick_pos += [y + DY, y - DY]
    ytick_lab  += [p["left_lab"], p["right_lab"]]
    ytick_col  += [C_NOV, C_OV]

    # 组名（右侧；第 3 组与第 2 组同为 PAD，副标题不重复）
    axA.text(0.1145, y + 0.020, p["group"], ha="right", va="center",
             fontsize=12.0, fontweight="bold", color=C_TXT, zorder=6)
    if i != 2:
        axA.text(0.1145, y - 0.225, p["sub"], ha="right", va="center",
                 fontsize=7.3, color=C_GRY, style="italic", zorder=6)

    # 组间分隔线
    if i < len(PAIRS) - 1:
        axA.axhline(y - 0.5, color=C_GRID, lw=0.9, ls=(0, (4, 3)), zorder=0)

# 零线 + 轴
axA.axvline(0, color="#9AA3AD", lw=1.0, ls="--", zorder=0)
axA.set_yticks(ytick_pos)
axA.set_yticklabels(ytick_lab, fontsize=7.5, linespacing=1.30)
for t, c in zip(axA.get_yticklabels(), ytick_col):
    t.set_color(c)
axA.set_xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.10])
axA.set_xticklabels(["0.00", "0.02", "0.04", "0.06", "0.08", "0.10"], fontsize=8.4)
axA.set_xlabel("MR effect estimate  β  (per 1-SD liability to hypothyroidism; IVW fixed effect)",
               fontsize=8.6, labelpad=6)
for s in ("top", "right", "left"):
    axA.spines[s].set_visible(False)
axA.tick_params(axis="x", length=3.2)
axA.tick_params(axis="y", length=0, pad=4)

axA.set_title("A   Paired β comparison: no-overlap vs overlap outcome library",
              fontsize=10.4, fontweight="bold", color=C_NOV, loc="left", pad=9)

# 图例（放左上空白）
leg = [
    Line2D([], [], marker="o", ls="none", mfc=C_NOV, mec="white", mew=1.3,
           ms=7.6, label="no-overlap library"),
    Line2D([], [], marker="s", ls="none", mfc=C_OV, mec="white", mew=1.3,
           ms=7.6, label="overlap library"),
]
axA.legend(handles=leg, loc="upper right", fontsize=8.0, frameon=True,
           facecolor="white", edgecolor=C_GRID, borderpad=0.55,
           handletextpad=0.5, bbox_to_anchor=(0.995, 0.995))

# ================================================================ Panel B
axB = fig.add_subplot(gs[0, 1])
axB.set_xlim(-0.062, 0.100)
axB.set_ylim(-0.66, 1.72)

yb = [1.0, 0.0]
for i, (lab, b, se, p) in enumerate(FLIP):
    y = yb[i]
    ci = beta_ci(b, se)
    col = C_NOV if i == 0 else C_OV
    mk = "o" if i == 0 else "s"
    axB.errorbar([b], [y], xerr=[[b - ci[0]], [ci[1] - b]],
                 fmt=mk, ms=9.4, mfc=col, mec="white", mew=1.5,
                 ecolor=col, elinewidth=2.0, capsize=3.6, capthick=2.0, zorder=4)
    # β / p 标签（放点的上方，避免与底部脚注打架）
    axB.text(b, y + 0.28, "β=%.4f\np=%.3f" % (b, p), ha="center", va="bottom",
             fontsize=7.6, color=col, linespacing=1.35, zorder=6)

# 零线
axB.axvline(0, color="#9AA3AD", lw=1.1, ls="--", zorder=1)
axB.text(0.0015, 1.60, "β = 0", fontsize=7.4, color=C_GRY, ha="left", va="center")

# 翻转箭头
axB.annotate(
    "", xy=(0.0240, 0.150), xytext=(-0.0060, 0.78),
    arrowprops=dict(arrowstyle="-|>", color=C_OV, lw=1.6,
                    shrinkA=7, shrinkB=7, connectionstyle="arc3,rad=-0.34"),
    zorder=5,
)
axB.text(0.0080, 0.620, "sign flip", ha="center", va="bottom",
         fontsize=10.2, fontweight="bold", color=C_OV, rotation=-19)
axB.text(-0.0585, 0.330, "TSH → CAD", ha="left", va="center",
         fontsize=8.8, fontweight="bold", color=C_OV)

# 库名用 y 轴刻度标签
axB.set_yticks(yb)
axB.set_yticklabels([FLIP[0][0], FLIP[1][0]], fontsize=7.5, linespacing=1.30)
axB.get_yticklabels()[0].set_color(C_NOV)
axB.get_yticklabels()[1].set_color(C_OV)
axB.set_xticks([-0.04, -0.02, 0.00, 0.02, 0.04, 0.06, 0.08])
axB.set_xticklabels(["-0.04", "-0.02", "0.00", "0.02", "0.04", "0.06", "0.08"],
                    fontsize=8.4)
axB.set_xlabel("MR effect estimate  β  (per 1-SD genetically proxied TSH; IVW fixed effect)",
               fontsize=8.6, labelpad=6)
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.tick_params(axis="x", length=3.2)
axB.tick_params(axis="y", length=0, pad=4)

axB.set_title("B   Exposure-specific sign instability",
              fontsize=10.4, fontweight="bold", color=C_NOV, loc="left", pad=9)

# ---------------------------------------------------------------- 总标题 + 框架署名
fig.text(0.022, 0.960,
         "Fig 2  |  Pairing the same exposure and trait across two outcome GWAS libraries",
         ha="left", va="center", fontsize=12.2, fontweight="bold", color=C_TXT)
fig.text(0.022, 0.905,
         "Effect estimates are IVW fixed-effect β with 95% CI. Pair ratios compare the overlap library with the no-overlap anchor. "
         "Overlap–bias framework: Burgess, Davies & Thompson (2016), Genet Epidemiol 40:597–608 — PMID 27625185, eq. (3).",
         ha="left", va="center", fontsize=7.3, color=C_GRY)

# 底部脚注（图外）
fig.text(0.022, 0.045,
         "† The FinnGen R13 overlap layer is non-independent (exposure GWAS includes FinnGen DF10; R3 ⊂ R13) and serves as a secondary reference only. "
         "The headline PAD ratio is the zero-overlap anchor MVP → Sakaue (×1.91).",
         ha="left", va="center", fontsize=7.2, color=C_GRY)

# ---------------------------------------------------------------- 保存
png = os.path.join(FIGDIR, "fig2_paired_differences.png")
svg = os.path.join(FIGDIR, "fig2_paired_differences.svg")
fig.savefig(png, dpi=300, facecolor="white")
fig.savefig(svg, facecolor="white")
plt.close(fig)

# PNG 规整：RGBA -> RGB，dpi 标签固定为 300
from PIL import Image
im = Image.open(png)
if im.mode != "RGB":
    im = im.convert("RGB")
im.save(png, dpi=(300, 300))

im2 = Image.open(png)
print("PNG:", png)
print("  size    =", im2.size, " mode =", im2.mode)
print("  dpi     =", im2.info.get("dpi"))
print("SVG:", svg)
