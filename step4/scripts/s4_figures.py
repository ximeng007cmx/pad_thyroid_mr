# -*- coding: utf-8 -*-
"""
Step 4 图形：FinnGen R13 I9_PAD 复制分析配图

输入  step4/results/harmonised_finngen_<t>.tsv
      step4/results/mr_finngen_<t>.json
      step4/results/loo_finngen_<t>.tsv
输出  step4/figures/<类型>_finngen_<t>.png        （每性状 4 张，共 12 张）
      step4/figures/summary_forest_finngen.png    （3 性状 × 5 方法 汇总森林图）

与 step3/s3_figures.py 同风格，另加：
  · 重叠状态警示（hypo 的 FinnGen 复制为 OVERLAP_INVALID，图上明确标注）
  · 方法学汇总森林图（OR 对数轴，跨性状对照）
"""
import argparse
import csv
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
S4 = os.path.join(ROOT, "step4")
RES = os.path.join(S4, "results")
FIG = os.path.join(S4, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.3})

TRAITS = ["hypo", "tsh", "ft4"]
CN = {"hypo": "Hypothyroidism", "tsh": "TSH", "ft4": "FT4"}

# 结局标识（--tag 切换）：决定读取/写出的文件名、图题、结局名与样本量描述
# 每项 = (文件名后缀, 图题名称, 样本量串, 结局简称 OU, 底注类型 NOTEK)
# finngen 为默认值，保证既有 13 张图重跑后几乎逐字/逐位不变（OU="PAD" 与原硬编码一致）
LABELS = {
    "finngen":  ("finngen",  "FinnGen R13 I9_PAD",        "22,244 cases / 460,490 controls", "PAD",    "ov_finngen"),
    "mvp_meta": ("mvp_meta", "MVP PAD (trans-ancestry meta)", "31,307 cases / 211,753 controls", "PAD", "clean_meta"),
    # ---- Step 5 泛血管谱 ----
    "cad_ukb":  ("cad_ukb",  "CAD (Mbatchou2021, UKB; GCST90013864)",
                 "29,339 cases / 322,724 controls", "CAD", "ov_ukb"),
    "cad_nouk": ("cad_nouk", "CAD (Nikpay2015 CARDIoGRAMplusC4D; GCST003116)",
                 "60,801 cases / 123,504 controls", "CAD", "clean"),
    "mi_ukb":   ("mi_ukb",   "MI (Donertas2021, UKB; GCST90038610)",
                 "11,081 cases / 473,517 controls", "MI", "ov_ukb"),
    "is_laa":   ("is_laa",   "IS-LAA (MEGASTROKE Malik2018; GCST005840)",
                 "4,373 cases / 406,111 controls", "IS-LAA", "clean_small"),
    "is_sv":    ("is_sv",    "IS-SV (MEGASTROKE Malik2018; GCST005841)",
                 "5,386 cases / 192,662 controls", "IS-SV", "clean_small"),
}
TAG, LABEL, NSTR, OU, NOTEK = LABELS["finngen"]

# 汇总森林图底注：按 NOTEK 分支（英文，DejaVu Sans 无 CJK 字形）
NOTES = {
    "ov_finngen": ("NOTE: Hypothyroidism -> FinnGen is flagged sample-overlap invalid "
                   "(exposure GWAS includes FinnGen); interpret hypo estimates as sensitivity only.\n"
                   "TSH / FT4 columns are overlap-valid."),
    "clean_meta": ("NOTE: the exposure GWAS (Rand 2025) does NOT include MVP, so there is no sample "
                   "overlap; all three traits are overlap-valid.\n"
                   "MVP PAD is a trans-ancestry meta (EUR/AFR/HIS), while the exposure is EUR-based — "
                   "cross-ancestry extrapolation applies."),
    "ov_ukb": ("NOTE: this outcome is UK-Biobank-derived and the exposure GWAS (Rand 2025) includes UKB, "
               "so the estimates are affected by sample overlap.\n"
               "Treat as a sensitivity analysis only; a companion UKB-free outcome of the same disease "
               "(CAD: GCST003116) is provided for comparison."),
    "clean": ("NOTE: the outcome cohort includes neither UKB nor FinnGen, so there is no sample overlap "
              "with the exposure GWAS; all three traits are overlap-valid."),
    "clean_small": ("NOTE: no sample overlap; however this stroke-subtype GWAS has a small case count, so "
                    "power is limited and the usable SNP set may be reduced.\n"
                    "Interpret null estimates with caution."),
}

METHODS = [("ivw_fixed", "IVW (fixed)"),
           ("ivw_random", "IVW (random, DL)"),
           ("mr_egger", "MR-Egger"),
           ("weighted_median", "Weighted median"),
           ("weighted_mode", "Weighted mode")]


def load(t):
    h = os.path.join(RES, "harmonised_%s_%s.tsv" % (TAG, t))
    j = os.path.join(RES, "mr_%s_%s.json" % (TAG, t))
    l = os.path.join(RES, "loo_%s_%s.tsv" % (TAG, t))
    for p in (h, j, l):
        if not os.path.exists(p):
            raise SystemExit("缺少输入文件：%s" % p)
    rows = list(csv.DictReader(open(h, encoding="utf-8"), delimiter="\t"))
    for r in rows:
        for k in ("beta_X", "se_X", "beta_Y", "se_Y", "beta_XY", "se_XY", "eaf_X", "F_X"):
            r[k] = float(r[k])
    main = json.load(open(j, encoding="utf-8"))
    loo = list(csv.DictReader(open(l, encoding="utf-8"), delimiter="\t"))
    for r in loo:
        r["beta"] = float(r["beta"]); r["se"] = float(r["se"])
    return rows, main, loo


def est_line(ax, slope, intercept, xmax, color, ls, label):
    xs = [0, xmax]
    ys = [intercept, intercept + slope * xmax]
    ax.plot(xs, ys, color=color, ls=ls, lw=1.5, label=label, zorder=3)


# 重叠状态 → 图上警示框文案（英文）
# 注意：字段名有两种历史来源 —— 旧 FinnGen 脚本写 overlap_status，
#       通用 pipeline (s4_pipeline.py) 写 overlap_note，故两者都要读。
OVERLAP_DRAW = {
    "OVERLAP_INVALID": ("SAMPLE OVERLAP INVALID",
                        "exposure GWAS includes FinnGen (endpoint DF10)"),
    "OVERLAP_UKB": ("SAMPLE OVERLAP (UKB)",
                    "exposure GWAS includes UK Biobank"),
}


def warn(ax, main):
    """重叠警示框。matplotlib 的 DejaVu Sans 不含 CJK 字形，因此一律绘制英文状态串
    （只影响图形，不改动任何结果文件）。"""
    st = main.get("overlap_status") or main.get("overlap_note") or ""
    for key, (title, detail) in OVERLAP_DRAW.items():
        if key in st:
            ax.text(0.02, 0.98, "%s\n%s" % (title, detail),
                    transform=ax.transAxes, va="top", ha="left", fontsize=7,
                    color="#cc3311", weight="bold",
                    bbox=dict(boxstyle="round,pad=0.3", fc="#fff3f0", ec="#cc3311", lw=0.6))
            return


def main():
    global TAG, LABEL, NSTR, OU, NOTEK
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="finngen", choices=sorted(LABELS),
                    help="结局标识：决定输入文件名、图题与输出文件名后缀")
    args = ap.parse_args()
    TAG, LABEL, NSTR, OU, NOTEK = LABELS[args.tag]
    cache = {}
    for t in TRAITS:
        rows, main_j, loo = load(t)
        cache[t] = (rows, main_j, loo)

        bX = [r["beta_X"] for r in rows]; sX = [r["se_X"] for r in rows]
        bY = [r["beta_Y"] for r in rows]; sY = [r["se_Y"] for r in rows]
        bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
        ivw_f = main_j["ivw_fixed"]; ivw_r = main_j["ivw_random"]; eg = main_j["mr_egger"]

        # ---- 1) 散点图 ----
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.errorbar(bX, bY, yerr=sY, xerr=sX, fmt="o", ms=2.5, lw=0.6,
                    color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5, zorder=2)
        xmax = max(abs(x) for x in bX) * 1.1
        est_line(ax, ivw_f["beta"], 0, xmax, "#cc3311", "-",
                 "IVW fixed: %.4f (p=%.2g)" % (ivw_f["beta"], ivw_f["p"]))
        est_line(ax, eg["beta"], eg["intercept"], xmax, "#0077bb", "--",
                 "MR-Egger (int p=%.2f)" % eg["intercept_pval"])
        ax.axhline(0, color="grey", lw=0.6); ax.axvline(0, color="grey", lw=0.6)
        ax.set_xlabel("SNP effect on %s (beta ± SE)" % CN[t])
        ax.set_ylabel("SNP effect on %s (log-OR ± SE)" % OU)
        ax.set_title("MR scatter — %s: %s -> %s (n=%d SNPs)" % (LABEL, CN[t], OU, len(rows)))
        ax.legend(loc="best", fontsize=8)
        warn(ax, main_j)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "scatter_%s_%s.png" % (TAG, t))); plt.close(fig)

        # ---- 2) 森林图（逐 SNP Wald 比） ----
        order = sorted(range(len(bR)), key=lambda i: bR[i])
        h = max(6.0, 0.11 * len(bR) + 2.0)
        fig, ax = plt.subplots(figsize=(8, h))
        ys = list(range(len(bR)))
        ax.errorbar([bR[i] for i in order], ys, xerr=[sR[i] for i in order],
                    fmt="o", ms=2, lw=0.6, color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5)
        ax.axvline(ivw_f["beta"], color="#cc3311", lw=1.4,
                   label="IVW fixed %.4f [%.4f, %.4f]"
                         % (ivw_f["beta"], ivw_f["beta"] - 1.96 * ivw_f["se"],
                            ivw_f["beta"] + 1.96 * ivw_f["se"]))
        ax.axvline(ivw_r["beta"], color="#ee7733", lw=1.2, ls="--",
                   label="IVW random %.4f (Q_p=%.2g)" % (ivw_r["beta"], ivw_r["Q_pval"]))
        ax.axvline(eg["beta"], color="#0077bb", lw=1.2, ls=":", label="MR-Egger %.4f" % eg["beta"])
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_yticks([])
        ax.set_xlabel("Wald ratio (log-OR %s per unit %s)" % (OU, CN[t]))
        ax.set_title("Forest — %s: %s -> %s, %d SNPs" % (LABEL, CN[t], OU, len(rows)))
        ax.legend(loc="lower right", fontsize=8)
        warn(ax, main_j)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "forest_%s_%s.png" % (TAG, t))); plt.close(fig)

        # ---- 3) 漏斗图 ----
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.errorbar(bR, [abs(s) for s in sR], xerr=sR, fmt="o", ms=2.5, lw=0.6,
                    color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5)
        ax.axvline(ivw_f["beta"], color="#cc3311", lw=1.4, label="IVW fixed")
        ax.axvline(eg["beta"], color="#0077bb", lw=1.2, ls="--", label="MR-Egger")
        ax.invert_yaxis()
        ax.set_xlabel("Wald ratio"); ax.set_ylabel("SE of ratio (inverted)")
        ax.set_title("Funnel — %s: %s -> %s" % (LABEL, CN[t], OU))
        ax.legend(loc="best", fontsize=8)
        warn(ax, main_j)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "funnel_%s_%s.png" % (TAG, t))); plt.close(fig)

        # ---- 4) 留一法 ----
        fig, ax = plt.subplots(figsize=(6, max(4.5, 0.12 * len(loo) + 1.5)))
        lbl = [d["SNP"] for d in loo]
        bs = [d["beta"] for d in loo]; ses = [d["se"] for d in loo]
        ys2 = list(range(len(loo)))
        ax.errorbar(bs, ys2, xerr=[1.96 * s for s in ses], fmt="o", ms=2, lw=0.6,
                    color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5)
        ax.axvline(ivw_f["beta"], color="#cc3311", lw=1.4, label="IVW fixed (all SNPs)")
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_yticks(ys2)
        ax.set_yticklabels(lbl, fontsize=4 if len(loo) > 80 else 6)
        ax.set_xlabel("IVW estimate after leaving out one SNP (95% CI)")
        ax.set_title("Leave-one-out — %s: %s -> %s" % (LABEL, CN[t], OU))
        ax.legend(loc="lower right", fontsize=8)
        warn(ax, main_j)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "loo_%s_%s.png" % (TAG, t))); plt.close(fig)
        print("%s: 4 张图完成（n=%d）" % (t, len(rows)))

    # ---- 5) 汇总森林图：3 性状 × 5 方法 ----
    palette = {"ivw_fixed": "#cc3311", "ivw_random": "#ee7733", "mr_egger": "#0077bb",
               "weighted_median": "#228833", "weighted_mode": "#aa3377"}
    ylabels, yvals, rows_plot, seps = [], [], [], []
    y = 0.0
    for t in TRAITS:
        rows, main_j, loo = cache[t]
        for key, lab in METHODS:
            d = main_j.get(key)
            if not d:
                continue
            rows_plot.append((d["beta"], d["se"], d["p"], palette[key], key))
            ylabels.append("%s — %s" % (CN[t], lab))
            yvals.append(y)
            y += 1.0
        seps.append(y - 0.5)
        y += 1.4

    fig, ax = plt.subplots(figsize=(9.5, 0.42 * len(yvals) + 2.2))
    for (b, se, p, c, key), yy in zip(rows_plot, yvals):
        lo, hi = b - 1.96 * se, b + 1.96 * se
        ax.errorbar([b], [yy], xerr=[[b - lo], [hi - b]], fmt="o", ms=5,
                    color=c, ecolor=c, elinewidth=1.2, capsize=2.5)
        ax.text(hi + 0.035, yy, "OR=%.3f [%.3f, %.3f]  p=%.2g"
                % (math.exp(b), math.exp(lo), math.exp(hi), p),
                va="center", fontsize=6.5, color="#333333")
    for s in seps[:-1]:
        ax.axhline(s + 0.7, color="#bbbbbb", lw=0.8, ls="-")
    ax.axvline(0, color="grey", lw=0.9)
    ax.set_yticks(yvals); ax.set_yticklabels(ylabels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(min(r[0] - 1.96 * r[1] for r in rows_plot) - 0.06,
                max(r[0] + 1.96 * r[1] for r in rows_plot) + 0.55)
    ax.set_xlabel("MR estimate (log-OR scale; right-hand text gives OR = exp(beta) with 95% CI)")
    ax.set_title("Replication in %s (%s)\nhypo / TSH / FT4 -> %s, five MR methods" % (LABEL, NSTR, OU))
    note = NOTES[NOTEK]
    fig.text(0.012, 0.012, note, fontsize=6.8, color="#cc3311", va="bottom")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=palette[k], label=l)
                       for k, l in METHODS if any(r[4] == k for r in rows_plot)],
              loc="lower right", fontsize=7.5, title="method", title_fontsize=7.5)
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.savefig(os.path.join(FIG, "summary_forest_%s.png" % TAG)); plt.close(fig)
    print("汇总森林图完成")

    print("all figures ->", FIG)


if __name__ == "__main__":
    main()
