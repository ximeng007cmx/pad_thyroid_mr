# -*- coding: utf-8 -*-
"""
s4_figure_crossdataset.py — 甲减 / TSH / FT4 -> PAD：主分析 + 两套复制的跨数据集对比森林图

输入
  主分析（Sakaue 2021 GCST90018890，EUR）： step3/results/mr_main_<t>.json
  FinnGen R13 I9_PAD（⚠ 样本重叠无效）：  step4/results/mr_finngen_<t>.json
  MVP trans-ancestry（dbGaP phs001672）：  step4/results/mr_mvp_meta_<t>.json
输出
  step4/figures/summary_forest_<t>_crossdataset.png   × 3
"""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
FIG = os.path.join(ROOT, "step4", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.3})

TRAITS = ["hypo", "tsh", "ft4"]
CN = {"hypo": "Hypothyroidism", "tsh": "TSH", "ft4": "FT4"}
UNIT = {"hypo": "1-SD increase in liability to hypothyroidism",
        "tsh": "1-SD increase in TSH",
        "ft4": "1-SD increase in FT4"}

# (显示名, 目录, 文件名模板, 样本量说明, 是否样本重叠无效)
DATASETS = [
    ("Sakaue 2021 GCST90018890 (EUR) - MAIN", "step3", "mr_main_%s.json",
     "7,114 / 475,964 cases/controls", False),
    ("FinnGen R13 I9_PAD (EUR) - replication", "step4", "mr_finngen_%s.json",
     "22,244 / 460,490", True),
    ("MVP trans-ancestry (phs001672) - replication", "step4", "mr_mvp_meta_%s.json",
     "31,307 / 211,753", False),
]
METHODS = [("ivw_fixed", "IVW (fixed)"), ("ivw_random", "IVW (random, DL)"),
           ("mr_egger", "MR-Egger"), ("weighted_median", "Weighted median"),
           ("weighted_mode", "Weighted mode")]
PAL = {"ivw_fixed": "#cc3311", "ivw_random": "#ee7733", "mr_egger": "#0077bb",
       "weighted_median": "#228833", "weighted_mode": "#aa3377"}


def one_trait(t):
    labels, ys, pts, seps, colors, used = [], [], [], [], [], []
    y = 0.0
    for name, sub, tpl, ns, overlap in DATASETS:
        path = os.path.join(ROOT, sub, "results", tpl % t)
        if not os.path.exists(path):
            print("  [跳过] 缺文件 %s" % path)
            continue
        d = json.load(open(path, encoding="utf-8"))
        first = True
        for key, lab in METHODS:
            m = d.get(key)
            if not m or "beta" not in m:
                continue
            if first:
                labels.append("%s\nn = %s\n%s" % (name, ns, lab))
                first = False
            else:
                labels.append(lab)
            ys.append(y)
            pts.append((m["beta"], m["se"], m["p"], PAL[key], key, overlap))
            colors.append("#cc3311" if overlap else "#222222")
            y += 1.0
        if not first:
            seps.append(y - 0.5)
            used.append((name, ns, overlap))
            y += 1.6
    if not pts:
        print("  %s：无可用结果，跳过" % t)
        return None

    fig, ax = plt.subplots(figsize=(11.5, 0.44 * len(ys) + 2.6))
    for (b, se, p, c, key, overlap), yy in zip(pts, ys):
        lo, hi = b - 1.96 * se, b + 1.96 * se
        ax.errorbar([b], [yy], xerr=[[b - lo], [hi - b]], fmt="o", ms=5,
                    color=c, ecolor=c, elinewidth=1.2, capsize=2.5,
                    alpha=0.55 if overlap else 1.0)
        ax.text(hi + 0.004, yy, "beta=%.4f  OR=%.3f [%.3f, %.3f]  p=%.3g"
                % (b, math.exp(b), math.exp(lo), math.exp(hi), p),
                va="center", fontsize=7, color="#555555" if overlap else "#222222",
                style="italic" if overlap else "normal")

    for s in seps[:-1]:
        ax.axhline(s + 0.8, color="#bbbbbb", lw=0.8)
    ax.axvline(0, color="grey", lw=1.0)
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=7.6)
    for tick, c in zip(ax.get_yticklabels(), colors):
        tick.set_color(c)
    ax.invert_yaxis()

    xmin = min(r[0] - 1.96 * r[1] for r in pts)
    xmax = max(r[0] + 1.96 * r[1] for r in pts)
    span = max(xmax - xmin, 1e-6)
    ax.set_xlim(xmin - 0.06 * span, xmax + 0.54 * span)
    ax.set_xlabel("MR estimate (log-OR of PAD per %s)" % UNIT[t])
    ax.set_title("%s -> Peripheral artery disease:\n"
                 "main analysis and replication datasets, five MR methods" % CN[t])

    has_overlap = any(r[5] for r in pts)
    note = ("FinnGen panel is SAMPLE-OVERLAP INVALID (the exposure GWAS includes FinnGen) "
            "- drawn faded, sensitivity analysis only.  Sakaue and MVP panels have no sample overlap."
            if has_overlap else
            "No sample overlap between exposure and the outcome datasets shown.")
    fig.text(0.012, 0.012, note, fontsize=7, color="#cc3311" if has_overlap else "#333333",
             va="bottom")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=PAL[k], label=l)
                       for k, l in METHODS],
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=5,
              frameon=False, fontsize=8.5, title="method", title_fontsize=8.5)
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    out = os.path.join(FIG, "summary_forest_%s_crossdataset.png" % t)
    fig.savefig(out); plt.close(fig)
    print("  已输出:", os.path.basename(out))
    return out


def main():
    for t in TRAITS:
        print("%s (%s):" % (t, CN[t]))
        one_trait(t)


if __name__ == "__main__":
    main()
