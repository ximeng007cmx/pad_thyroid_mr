# -*- coding: utf-8 -*-
"""
Step 3.9  MR 图形：散点图 / 森林图 / 漏斗图 / 留一法 × {hypo, tsh, ft4}
输入  step3/results/harmonised_<t>.tsv, mr_main_<t>.json, loo_<t>.tsv
输出  step3/figures/<类型>_<t>.png  (12 张)
"""
import os, csv, json, math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
S3   = os.path.join(ROOT, "step3")
RES  = os.path.join(S3, "results")
FIG  = os.path.join(S3, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.3})
CN = {"hypo": "Hypothyroidism", "tsh": "TSH", "ft4": "FT4"}
TRAITS = ["hypo", "tsh", "ft4"]

def load(t):
    rows = list(csv.DictReader(open(os.path.join(RES, "harmonised_%s.tsv" % t),
                                    encoding="utf-8"), delimiter="\t"))
    for r in rows:
        r["beta_X"] = float(r["beta_X"]); r["se_X"] = float(r["se_X"])
        r["beta_Y"] = float(r["beta_Y"]); r["se_Y"] = float(r["se_Y"])
        r["beta_XY"] = float(r["beta_XY"]); r["se_XY"] = float(r["se_XY"])
    main = json.load(open(os.path.join(RES, "mr_main_%s.json" % t), encoding="utf-8"))
    loo = list(csv.DictReader(open(os.path.join(RES, "loo_%s.tsv" % t),
                                   encoding="utf-8"), delimiter="\t"))
    for r in loo:
        r["beta"] = float(r["beta"]); r["se"] = float(r["se"])
    return rows, main, loo

def est_line(ax, slope, intercept, xmax, color, ls, label):
    xs = [0, xmax]
    ys = [intercept, intercept + slope * xmax]
    ax.plot(xs, ys, color=color, ls=ls, lw=1.5, label=label, zorder=3)

for t in TRAITS:
    rows, main, loo = load(t)
    bX = [r["beta_X"] for r in rows]; sX = [r["se_X"] for r in rows]
    bY = [r["beta_Y"] for r in rows]; sY = [r["se_Y"] for r in rows]
    bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
    ivw_f = main["ivw_fixed"]; ivw_r = main["ivw_random"]
    eg = main["mr_egger"]

    # ---- 1) 散点图 b_X vs b_Y ----
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.errorbar(bX, bY, yerr=sY, xerr=sX, fmt="o", ms=2.5, lw=0.6,
                color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5, zorder=2)
    xmax = max(abs(x) for x in bX) * 1.1
    est_line(ax, ivw_f["beta"], 0, xmax, "#cc3311", "-", "IVW (fixed)")
    est_line(ax, eg["beta"], eg["intercept"], xmax, "#0077bb", "--", "MR-Egger")
    ax.axhline(0, color="grey", lw=0.6)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("SNP effect on %s (beta ± SE)" % CN[t])
    ax.set_ylabel("SNP effect on PAD (log-OR ± SE)")
    ax.set_title("MR scatter: %s -> PAD (n=%d SNPs)" % (CN[t], len(rows)))
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "scatter_%s.png" % t)); plt.close(fig)

    # ---- 2) 森林图（前 60 个 + 汇总，按 beta 排序取代表：全部画但图高自适应） ----
    order = sorted(range(len(bR)), key=lambda i: bR[i])
    h = max(6.0, 0.11 * len(bR) + 2.0)
    fig, ax = plt.subplots(figsize=(8, h))
    ys = list(range(len(bR)))
    ax.errorbar([bR[i] for i in order], ys, xerr=[sR[i] for i in order],
                fmt="o", ms=2, lw=0.6, color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5)
    ax.axvline(ivw_f["beta"], color="#cc3311", lw=1.4, label="IVW fixed: %.3f [%.3f, %.3f]"
               % (ivw_f["beta"], ivw_f["beta"]-1.96*ivw_f["se"], ivw_f["beta"]+1.96*ivw_f["se"]))
    ax.axvline(ivw_r["beta"], color="#ee7733", lw=1.2, ls="--", label="IVW random: %.3f" % ivw_r["beta"])
    ax.axvline(eg["beta"], color="#0077bb", lw=1.2, ls=":", label="Egger: %.3f" % eg["beta"])
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Wald ratio (log-OR PAD per unit %s)" % CN[t])
    ax.set_title("Forest: %s -> PAD, %d SNPs" % (CN[t], len(rows)))
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "forest_%s.png" % t)); plt.close(fig)

    # ---- 3) 漏斗图 ----
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.errorbar(bR, [abs(s) for s in sR], xerr=sR, fmt="o", ms=2.5, lw=0.6,
                color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5)
    ax.axvline(ivw_f["beta"], color="#cc3311", lw=1.4, label="IVW fixed")
    ax.axvline(eg["beta"], color="#0077bb", lw=1.2, ls="--", label="Egger")
    ax.invert_yaxis()
    ax.set_xlabel("Wald ratio"); ax.set_ylabel("SE of ratio (inverted)")
    ax.set_title("Funnel: %s -> PAD" % CN[t])
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "funnel_%s.png" % t)); plt.close(fig)

    # ---- 4) 留一法 ----
    fig, ax = plt.subplots(figsize=(6, max(4.5, 0.12*len(loo)+1.5)))
    lbl = [d["SNP"] for d in loo]
    bs = [d["beta"] for d in loo]; ses = [d["se"] for d in loo]
    ys = list(range(len(loo)))
    ax.errorbar(bs, ys, xerr=[1.96*s for s in ses], fmt="o", ms=2, lw=0.6,
                color="#4477aa", ecolor="#a0b4d6", elinewidth=0.5)
    ax.axvline(ivw_f["beta"], color="#cc3311", lw=1.4, label="IVW fixed (all SNPs)")
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_yticks(ys); ax.set_yticklabels(lbl, fontsize=4 if len(loo) > 80 else 6)
    ax.set_xlabel("IVW estimate after leaving out one SNP (95% CI)")
    ax.set_title("Leave-one-out: %s -> PAD" % CN[t])
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "loo_%s.png" % t)); plt.close(fig)
    print("%s: 4 figures done" % t)

print("all figures ->", FIG)
