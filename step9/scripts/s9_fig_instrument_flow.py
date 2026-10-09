# -*- coding: utf-8 -*-
"""P2-10a: 绘制工具变量筛选流程图（Fig S / 补充图），300 dpi PNG + SVG。
数据全部来自 step9/results/instrument_selection_flow.tsv（实测复核）。
图内文字用英文（入稿用）。
"""
import os, sys, io, textwrap
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
OUTD = os.path.join(ROOT, "step9", "figures")
os.makedirs(OUTD, exist_ok=True)

# ---------------- 数据（实测） ----------------
COLS = ["Hypothyroidism", "TSH", "FT4"]
STAGES = [
    ("Variants in exposure GWAS\nsummary statistics",      ["7,344,354", "6,151,997", "5,901,316"], None),
    ("p < 5e-8",                                            ["16,247",    "17,246",    "1,690"],     "Not genome-wide significant"),
    ("Autosomal biallelic SNPs",                            ["12,782",    "13,731",    "1,474"],     "Indel / multiallelic / non-autosomal"),
    ("Mapped to 1000G EUR\nreference panel",                ["11,070",    "12,007",    "1,445"],     "No panel match / ambiguous position / duplicate"),
    ("Independent instruments after\nLD clumping (r2<0.001, 10 Mb)", ["222", "258", "65"],             "LD with a stronger index variant"),
    ("Retained after\nharmonisation",                       ["221",       "256",       "65"],        "EAF discordant (>0.2) or absent in outcome"),
]
EXCL = [
    ("Not genome-wide significant",                          ["7,328,107", "6,134,751", "5,899,626"]),
    ("Indel / multiallelic / non-autosomal",                 ["3,465",     "3,515",     "216"]),
    ("No panel match / ambiguous position / duplicate",      ["1,712",     "1,724",     "29"]),
    ("LD with a stronger index variant",                     ["10,848",    "11,749",    "1,380"]),
    ("EAF discordant (>0.2) or absent in outcome",           ["1",         "2",         "0"]),
]

# ---------------- 版面 ----------------
FIG_W, FIG_H = 12.4, 8.6
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

C_MAIN = "#1F4E79"; C_FILL = "#EAF1F8"; C_ACC = "#C00000"; C_ACC_FILL = "#FCE8E6"
C_OK = "#1E6B3A"; C_OK_FILL = "#E8F3EC"
FS = 9.2

def box(x, y, w, h, text, fc, ec, fs=FS, bold=False, tc="#1A1A1A"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.1",
                                linewidth=1.15, edgecolor=ec, facecolor=fc, zorder=2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, zorder=3, fontweight=("bold" if bold else "normal"), linespacing=1.35)

def arrow(x1, y1, x2, y2, color=C_MAIN, lw=1.25, style="-|>", ms=9):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                                 linewidth=lw, color=color, zorder=1, shrinkA=1, shrinkB=1))

# 表头
ax.text(50, 97.4, "Instrument selection flow",
        ha="center", va="center", fontsize=14, fontweight="bold", color=C_MAIN)

LANE_X = [13.0, 41.5, 70.0]     # 三列主流程左边界
LANE_W = 17.0
EXC_X, EXC_W = 87.4, 12.5

top, bot = 90.0, 12.5
n = len(STAGES)
step = (top - bot) / (n - 1)
BH = 5.6

# 列标题
for j, name in enumerate(COLS):
    ax.text(LANE_X[j] + LANE_W/2, top + 4.0, name, ha="center", va="center",
            fontsize=11.5, fontweight="bold", color=C_MAIN)

ys = []
for i, (label, counts, reason) in enumerate(STAGES):
    y = top - i*step
    ys.append(y)
    for j in range(3):
        last = (i == n - 1)
        fc = C_OK_FILL if last else C_FILL
        ec = C_OK if last else C_MAIN
        box(LANE_X[j], y - BH/2, LANE_W, BH, counts[j], fc, ec, fs=11.5, bold=True)
    # 阶段说明（左侧）
    ax.text(LANE_X[0] - 1.4, y, label, ha="right", va="center", fontsize=FS,
            color="#333333", linespacing=1.3)
    # 箭头 + 排除框
    if i < n - 1:
        for j in range(3):
            arrow(LANE_X[j] + LANE_W/2, y - BH/2 - 0.25, LANE_X[j] + LANE_W/2, y - step + BH/2 + 0.25)
        if reason:
            ym = y - step/2
            # 原因自动换行（防止文字溢出框体、贴图右缘）+ 三列计数
            reason_w = "\n".join(textwrap.wrap(reason, width=22))
            txt = reason_w + "\n\u2212 " + EXCL[i][1][0] + " / " + EXCL[i][1][1] + " / " + EXCL[i][1][2]
            box(EXC_X, ym - 3.9, EXC_W, 7.8, txt, C_ACC_FILL, C_ACC, fs=6.5, tc=C_ACC)
            # 虚线连到第一列
            ax.add_patch(FancyArrowPatch((LANE_X[0] + LANE_W, ym), (EXC_X - 0.3, ym),
                                         arrowstyle="-", linestyle=(0, (3, 2)),
                                         linewidth=0.9, color=C_ACC, zorder=1))

# 图注
ax.text(50, 7.4,
        "Screening threshold p < 5e-8; LD clumping r2 < 0.001 within 10 Mb (PLINK v1.90b7); "
        "1000 Genomes Phase 3 EUR reference panel (n = 294).",
        ha="center", va="center", fontsize=8.0, color="#555555", style="italic")
ax.text(50, 4.0,
        "Exclusion counts are per trait (hypothyroidism / TSH / FT4). Harmonisation removed variants "
        "with discordant effect-allele frequency or absent from the outcome dataset.",
        ha="center", va="center", fontsize=8.0, color="#555555", style="italic")

png = os.path.join(OUTD, "instrument_selection_flow.png")
svg = os.path.join(OUTD, "instrument_selection_flow.svg")
fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(svg, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("saved:\n ", png, "\n ", svg)
print("dpi=300 | size:", os.path.getsize(png), "bytes")
