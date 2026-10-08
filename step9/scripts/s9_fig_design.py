# -*- coding: utf-8 -*-
"""Fig 1 —— 研究设计图（三面板）
(A) 3 暴露 × 8 结局分析矩阵（每格 = harmonise 后 SNP 数；底色 = 样本重叠状态；‡ = 线性概率尺度）
(B) 三组「跨结局库」配对比较（无重叠锚 -> 有重叠）
(C) π -> 偏倚示意 R(pi) = 1 + pi*r/F_w（显式标注框架来源 Burgess 2016）
输出：step9/figures/fig1_study_design.png (300 dpi) + .svg
数值全部取自 step9/results/i2_heterogeneity_all.tsv、overlap_sim_summary.json、
Step5_泛血管谱结果.md、release/.../DATA_SOURCES.md，未作人工修改。
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
OUTD = os.path.join(ROOT, "step9", "figures")
os.makedirs(OUTD, exist_ok=True)

C_MAIN = "#1F4E79"; C_FILL = "#EAF1F8"
C_ACC = "#C00000";  C_ACC_FILL = "#FCE8E6"
C_OK = "#1E6B3A";   C_OK_FILL = "#E8F3EC"
C_PART = "#B26A00"; C_PART_FILL = "#FBF3E2"
C_TXT = "#1A1A1A"

# ================================================================ 数据（实测）
EXPOSURES = [
    dict(name="Hypothyroidism", short="Hypothyroidism", gwas="GCST90572791",
         n="113,393 / 1,065,268", nsnp="222", r2="5.91%", fmin="29.4",
         cohort="UKB + FinnGen DF10 (+23andMe, Danish)", finn=True),
    dict(name="TSH", short="TSH", gwas="GCST90572789", n="482,873 (EUR)",
         nsnp="258", r2="10.75%", fmin="29.3",
         cohort="UKB (+23andMe, Danish); no Finnish data", finn=False),
    dict(name="Free T4", short="FT4", gwas="GCST90572790", n="191,449 (EUR)",
         nsnp="65", r2="2.37%", fmin="30.1",
         cohort="UKB (+23andMe, Danish); no Finnish data", finn=False),
]
OUTCOMES = [
    dict(short="PAD\n(Sakaue 2021)", acc="GCST90018890", n="7,114 / 475,964",
         group="PAD — peripheral artery disease", ov_hypo="partial", ov_other="partial"),
    dict(short="PAD\n(FinnGen R13)", acc="I9_PAD, R13", n="22,244 / 460,490",
         group="", ov_hypo="full", ov_other="none"),
    dict(short="PAD\n(MVP)", acc="pha004826", n="31,307 / 211,753",
         group="", ov_hypo="none", ov_other="none"),
    dict(short="CAD\n(Nikpay 2015)", acc="GCST003116", n="—",
         group="CAD", ov_hypo="none", ov_other="none"),
    dict(short="CAD\n(UKB, Mbatchou)", acc="GCST90013864", n="29,339 / 322,724",
         group="", ov_hypo="full", ov_other="full"),
    dict(short="MI\n(UKB) ‡", acc="GCST90038610", n="11,081 / 473,517",
         group="MI", ov_hypo="full", ov_other="full"),
    dict(short="IS-LAA\n(MEGASTROKE)", acc="GCST005840", n="4,373 / 406,111",
         group="IS — ischaemic stroke", ov_hypo="none", ov_other="none"),
    dict(short="IS-SV\n(MEGASTROKE)", acc="GCST005841", n="5,386 / 192,662",
         group="", ov_hypo="none", ov_other="none"),
]
NCLUMP = {"Hypothyroidism": [221, 220, 213, 215, 222, 221, 216, 195],
          "TSH":            [256, 248, 251, 252, 258, 256, 251, 220],
          "FT4":            [65,  60,  63,  65,  65,  65,  65,  60]}
PAIRS = [
    dict(d="CAD", base="Nikpay 2015\n(no overlap)", over="UKB · Mbatchou\n(UKB overlap)",
         b0=0.0308, b1=0.0546, R="1.77"),
    dict(d="PAD", base="MVP\n(no overlap)", over="Sakaue 2021\n(partial overlap)",
         b0=0.0237, b1=0.0452, R="1.91"),
    dict(d="PAD", base="MVP\n(no overlap)", over="FinnGen R13\n(overlap)",
         b0=0.0237, b1=0.0857, R="3.62"),
]
F_W = 109.0
R_OBS = [1.77, 1.91, 3.62]
R_MAX_TOP = 1.1835
FRAMEWORK = ("Overlap–bias framework: Burgess, Davies & Thompson (2016), "
             "Genet Epidemiol 40:597–608 — PMID 27625185, eq. (3)")

# ================================================================ 版面
fig = plt.figure(figsize=(13.6, 11.0), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.42, 1.0],
                      left=0.038, right=0.978, top=0.905, bottom=0.130,
                      wspace=0.14, hspace=0.15)
axA = fig.add_subplot(gs[0, :]); axA.set_xlim(0, 100); axA.set_ylim(0, 50); axA.axis("off")
axB = fig.add_subplot(gs[1, 0]); axB.set_xlim(0, 100); axB.set_ylim(0, 34); axB.axis("off")
axC = fig.add_subplot(gs[1, 1])

fig.text(0.5, 0.972,
         "Figure 1  |  Study design: three thyroid-function exposures × eight vascular outcome GWAS",
         ha="center", va="center", fontsize=11.4, fontweight="bold", color=C_MAIN)
fig.text(0.5, 0.946,
         "Two-sample Mendelian randomisation · 24 exposure–outcome layers · "
         "paired zero-overlap ↔ overlap design for outcome-library sensitivity",
         ha="center", va="center", fontsize=8.0, color="#4A4A4A")
fig.text(0.5, 0.925,
         "All datasets are publicly available summary statistics; no individual-level data were analysed.",
         ha="center", va="center", fontsize=7.0, color="#6A6A6A")


def txt(ax, x, y, s, fs=7.4, c=C_TXT, ha="center", va="center", bold=False,
        rot=0, ls=1.35, z=5):
    ax.text(x, y, s, fontsize=fs, color=c, ha=ha, va=va, rotation=rot,
            fontweight=("bold" if bold else "normal"), linespacing=ls, zorder=z)


def rbox(ax, x, y, w, h, fc, ec, lw=1.0, r=0.8, z=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.0,rounding_size=%s" % r,
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z))


# ================================================================ (A) 矩阵
GX0, GX1 = 27.0, 99.0
CW = (GX1 - GX0) / 8.0
ROWH = 6.0
MTOP = 32.0                     # 矩阵顶
txt(axA, 0.5, 47.6, "(A)  Analysis matrix — 3 exposures × 8 outcome GWAS",
    fs=9.6, c=C_MAIN, ha="left", bold=True)

# 分组条
GRP = [("PAD — peripheral artery disease", 0, 3), ("CAD", 3, 2),
       ("MI", 5, 1), ("IS — ischaemic stroke", 6, 2)]
for label, a, b in GRP:
    if not label:
        continue
    x0 = GX0 + a * CW; w = (b - a) * CW
    axA.add_patch(Rectangle((x0, 43.6), w, 2.4, facecolor=C_FILL,
                            edgecolor=C_MAIN, linewidth=0.8, zorder=2))
    txt(axA, x0 + w / 2, 44.8, label, fs=6.5, c=C_MAIN, bold=True)

# 列表头：竖排短名 + 横排 accession / 样本量
for j, o in enumerate(OUTCOMES):
    cx = GX0 + (j + 0.5) * CW
    txt(axA, cx - 1.6, 33.6, o["short"], fs=6.4, rot=90, va="bottom",
        ha="center", c=C_MAIN, bold=True, ls=1.2)
    txt(axA, cx + 2.6, 40.6, o["acc"], fs=5.7, c="#4A4A4A")
    txt(axA, cx + 2.6, 38.4, o["n"], fs=5.7, c="#4A4A4A")

# 行头（暴露）
for i, e in enumerate(EXPOSURES):
    y = MTOP - (i + 1) * ROWH
    txt(axA, 0.5, y + ROWH * 0.72, e["name"], fs=8.6, ha="left", bold=True, c=C_MAIN)
    txt(axA, 0.5, y + ROWH * 0.42, "%s  ·  %s" % (e["gwas"], e["n"]),
        fs=6.4, ha="left", c="#3A3A3A")
    txt(axA, 0.5, y + ROWH * 0.20, "%s SNPs · R² = %s · min F = %s"
        % (e["nsnp"], e["r2"], e["fmin"]), fs=6.4, ha="left", c="#3A3A3A")

# 格子
OV_STYLE = {"none": (C_OK_FILL, C_OK, "no overlap"),
            "partial": (C_PART_FILL, C_PART, "partial overlap"),
            "full": (C_ACC_FILL, C_ACC, "overlap")}
for i, e in enumerate(EXPOSURES):
    y = MTOP - (i + 1) * ROWH
    for j, o in enumerate(OUTCOMES):
        lvl = o["ov_hypo"] if e["finn"] else o["ov_other"]
        fc, ec, _ = OV_STYLE[lvl]
        x = GX0 + j * CW
        axA.add_patch(Rectangle((x + 0.3, y + 0.35), CW - 0.6, ROWH - 0.7,
                                facecolor=fc, edgecolor=ec, linewidth=0.9, zorder=3))
        txt(axA, x + CW / 2, y + ROWH / 2, str(NCLUMP[e["short"]][j]),
            fs=8.4, bold=True, c=ec, z=5)

# 矩阵下方：配对连接线
def bracket(j0, j1, y, label, c):
    x0 = GX0 + j0 * CW + 0.6; x1 = GX0 + (j1 + 1) * CW - 0.6
    axA.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-",
                                  linewidth=1.0, color=c, zorder=6))
    axA.add_patch(FancyArrowPatch((x0, y), (x0, y + 0.9), arrowstyle="-",
                                  linewidth=1.0, color=c, zorder=6))
    axA.add_patch(FancyArrowPatch((x1, y), (x1, y + 0.9), arrowstyle="-",
                                  linewidth=1.0, color=c, zorder=6))
    txt(axA, (x0 + x1) / 2, y - 1.25, label, fs=6.5, c=c, bold=True)

bracket(3, 4, 13.0, "×1.77   (vs Nikpay, no overlap)", "#7A3E9D")
bracket(0, 1, 9.4, "×1.91  /  ×3.62   (vs MVP, no overlap)", "#7A3E9D")

# 图例
lx, ly = GX0, 6.0
for lvl, lab in (("none", "no participant overlap"), ("partial", "partial overlap"),
                 ("full", "overlap")):
    fc, ec, _ = OV_STYLE[lvl]
    axA.add_patch(Rectangle((lx, ly - 0.85), 1.6, 1.7, facecolor=fc,
                            edgecolor=ec, linewidth=0.9, zorder=3))
    txt(axA, lx + 2.3, ly, lab, fs=6.6, ha="left", c="#3A3A3A")
    lx += 16.5
txt(axA, GX0, 3.4, "Cell values = number of instruments retained after harmonisation "
                   "(Table 1 / Table S2).   ‡ = outcome on a linear-probability scale (β/SE not a log-OR).",
    fs=6.4, ha="left", c="#6A6A6A")
txt(axA, GX0, 1.5, "Overlap: all three exposures share participants with UK-Biobank-based outcomes; "
                   "the hypothyroidism GWAS additionally contains FinnGen DF10 (⊃ FinnGen R3).",
    fs=6.4, ha="left", c="#6A6A6A")

# ================================================================ (B) 配对
txt(axB, 0.5, 32.4, "(B)  Paired cross-outcome-GWAS comparisons", fs=9.4,
    c=C_MAIN, ha="left", bold=True)
txt(axB, 0.5, 29.9, "same exposure · same phenotype · two outcome GWAS from different sources",
    fs=6.7, ha="left", c="#5A5A5A")

rowh = 6.6
for k, p in enumerate(PAIRS):
    yc = 23.4 - k * rowh
    txt(axB, 0.5, yc, p["d"], fs=8.6, ha="left", bold=True, c=C_MAIN)
    rbox(axB, 6.0, yc - 2.3, 15.0, 4.6, C_OK_FILL, C_OK)
    txt(axB, 13.5, yc, p["base"], fs=6.4, c=C_OK, ls=1.2)
    axB.add_patch(FancyArrowPatch((21.6, yc), (28.4, yc), arrowstyle="-|>",
                                  mutation_scale=8, linewidth=1.1, color=C_MAIN, zorder=4))
    rbox(axB, 29.0, yc - 2.3, 15.0, 4.6, C_ACC_FILL, C_ACC)
    txt(axB, 36.5, yc, p["over"], fs=6.4, c=C_ACC, ls=1.2)
    txt(axB, 50.5, yc, "×" + p["R"], fs=10.6, bold=True, c=C_ACC)

txt(axB, 0.5, 4.6,
    "β (hypothyroidism):  0.0308 → 0.0546   |   0.0237 → 0.0452   |   0.0237 → 0.0857\n"
    "TSH → CAD reverses sign across libraries:  −0.0135 (p = 0.500)  →  +0.0394 (p = 0.035)",
    fs=6.5, ha="left", c="#3A3A3A", va="bottom", ls=1.5)
txt(axB, 0.5, 1.4,
    "Secondary (non-independent) comparison: Sakaue 2021 ↔ FinnGen R13 = 1.89 "
    "(FinnGen R3 ⊂ R13) — reference only.",
    fs=6.2, ha="left", c="#8A6A00", va="bottom")

# ================================================================ (C) π -> 偏倚
pi = np.linspace(0, 1, 400)
cols = {5: "#4C72B0", 10: "#DD8452", 20: "#8172B3"}
for r, c in cols.items():
    axC.plot(pi, 1.0 + pi * r / F_W, lw=1.5, color=c, label="r = %d" % r)
axC.set_yscale("log")
axC.set_ylim(1.0, 4.4); axC.set_xlim(0, 1)
axC.axhspan(1.0, R_MAX_TOP, color=C_OK, alpha=0.16, zorder=0)
axC.axhspan(1.77, 3.70, color=C_ACC, alpha=0.09, zorder=0)
for R, c in zip(R_OBS, ["#C44E52", "#C44E52", "#8B1A1A"]):
    axC.axhline(R, color=c, ls=":", lw=0.9, alpha=0.9)
axC.text(0.02, 2.42, "observed between-GWAS ratios\n1.77 / 1.91 / 3.62",
         fontsize=6.6, color="#8B1A1A", ha="left", va="center", linespacing=1.35)
axC.text(0.98, 1.215, "attainable at complete overlap (π = 1):\n+4.6% / +9.2% / +18.4%",
         fontsize=6.3, color=C_OK, ha="right", va="bottom", linespacing=1.3)
axC.set_xlabel("overlap proportion  π   (fraction of the outcome sample shared)", fontsize=7.2)
axC.set_ylabel("R(π) = paired estimate relative to\nthe no-overlap baseline", fontsize=7.2)
axC.tick_params(labelsize=6.4, length=2.4, width=0.7)
for s in ("top", "right"):
    axC.spines[s].set_visible(False)
for s in ("left", "bottom"):
    axC.spines[s].set_linewidth(0.8)
axC.legend(fontsize=6.4, loc="upper left", frameon=False,
           title="r = β_obs / β", title_fontsize=6.4)
axC.set_title("(C)  Effect of sample overlap on the paired estimate\n"
              "R(π) = 1 + π·r / F_w      (F_w = %.0f, hypothyroidism)" % F_W,
              fontsize=8.0, color=C_MAIN, linespacing=1.45)

# ★ 框架来源与公式同处出现（红线要求）
fig.text(0.528, 0.062, FRAMEWORK, fontsize=6.6, color=C_ACC, ha="left", va="bottom")
fig.text(0.528, 0.038,
         "To explain the observed ×1.77 by overlap alone, the observational association would "
         "have to be 84× the causal effect (×3.62 → 287×).",
         fontsize=6.3, color="#3A3A3A", ha="left", va="bottom")
fig.text(0.038, 0.038,
         "Exposures (n = 3, one publication): Rand et al. Nat Genet 2025 — PMID 41238958.",
         fontsize=6.4, color="#6A6A6A", ha="left", va="bottom")

png = os.path.join(OUTD, "fig1_study_design.png")
svg = os.path.join(OUTD, "fig1_study_design.svg")
fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(svg, bbox_inches="tight", facecolor="white")
plt.close(fig)

# 规整为 RGB + 精确 300 dpi（期刊要求；matplotlib 输出为 RGBA 且 dpi 带浮点尾差）
from PIL import Image
im = Image.open(png)
if im.mode != "RGB":
    bg = Image.new("RGB", im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3] if im.mode == "RGBA" else None)
    im = bg
im.save(png, dpi=(300, 300))
print("OK  %s  %d B  size=%s dpi=%s mode=%s"
      % (png, os.path.getsize(png), im.size, im.info.get("dpi"), im.mode))
print("OK  %s  %d B" % (svg, os.path.getsize(svg)))
