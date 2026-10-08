# -*- coding: utf-8 -*-
"""P1-12b: 从现有 Cochran Q 与 df 补算 I² 统计量。
I² = max(0, (Q - df) / Q) * 100%   (Higgins & Thompson 2002)
只做统计量派生，不改动任何原始数值。
"""
import json, os, glob, sys, io, csv
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"

# 结局显示名（中文）与英文标签
OUTCOME_LABEL = {
    "mr_main":        ("主分析 PAD (Sakaue)", "PAD (Sakaue 2021)"),
    "mr_sakaue":      ("主分析 PAD (Sakaue)", "PAD (Sakaue 2021)"),
    "mr_finngen":     ("FinnGen R13 I9_PAD", "PAD (FinnGen R13)"),
    "mr_mvp_meta":    ("MVP te.PAD (无重叠)", "PAD (MVP)"),
    "mr_cad_ukb":     ("CAD (UKB, Aragam)", "CAD (UKB)"),
    "mr_cad_nouk":    ("CAD (Nikpay, 无UKB)", "CAD (Nikpay)"),
    "mr_mi_ukb":      ("MI (UKB)", "MI (UKB)"),
    "mr_is_laa":      ("IS-LAA (MEGASTROKE)", "IS-LAA"),
    "mr_is_sv":       ("IS-SV (MEGASTROKE)", "IS-SV"),
}
TRAIT_CN = {"hypo":"甲减","tsh":"TSH","ft4":"FT4"}

def i2_from_q(Q, df):
    if Q is None or df is None or Q <= 0:
        return None
    v = (Q - df) / Q * 100.0
    return max(0.0, v)

# 归一化 prefix：mr_sakaue→mr_main（同一分析）；mr_finngen_r13→mr_finngen
PREFIX_ALIAS = {"mr_sakaue": "mr_main", "mr_finngen_r13": "mr_finngen"}

rows = []
seen = set()
search_dirs = ["step3/results", "step4/results", "step5/results", "step6/results", "step7/results", "step8/results"]
for d in search_dirs:
    p = os.path.join(ROOT, d)
    if not os.path.isdir(p):
        continue
    for f in sorted(glob.glob(os.path.join(p, "mr_*.json"))):
        base = os.path.basename(f)[:-5]  # 去掉 .json
        # 解析 前缀 + 性状
        trait = None
        for t in TRAIT_CN:
            if base.endswith("_" + t):
                trait = t; prefix = base[:-len(t)-1]; break
        if trait is None:
            continue
        prefix = PREFIX_ALIAS.get(prefix, prefix)
        key = (prefix, trait)
        if key in seen:      # 去重：同一分析同一性状只保留一次
            continue
        seen.add(key)
        try:
            j = json.load(open(f, encoding="utf-8"))
        except Exception as e:
            print("SKIP", f, e); continue
        ivw = j.get("ivw_random", {}) or {}
        Q = ivw.get("Q"); df = ivw.get("Q_df")
        tau2 = ivw.get("tau2")
        nsnp = j.get("nsnp") or ivw.get("nsnp")
        i2 = i2_from_q(Q, df)
        # 加权平均 F（F_w）：优先取 JSON 的 F_X；缺失则空
        Fw = (j.get("directionality", {}) or {}).get("F_X")
        olab_cn, olab_en = OUTCOME_LABEL.get(prefix, (prefix, prefix))
        rows.append({
            "outcome_cn": olab_cn, "outcome_en": olab_en, "prefix": prefix,
            "trait": trait, "trait_cn": TRAIT_CN[trait],
            "nsnp": nsnp, "Q": Q, "Q_df": df, "tau2": tau2,
            "I2_pct": i2, "Fw": Fw,
            "Q_pval": ivw.get("Q_pval"),
            "file": os.path.relpath(f, ROOT).replace("\\","/"),
        })

# 排序：outcome 顺序按 OUTCOME_LABEL 定义序，再 trait 顺序 hypo,tsh,ft4
order = ["mr_main","mr_sakaue","mr_finngen","mr_mvp_meta","mr_cad_ukb","mr_cad_nouk","mr_mi_ukb","mr_is_laa","mr_is_sv"]
torder = {"hypo":0,"tsh":1,"ft4":2}
rows.sort(key=lambda r: (order.index(r["prefix"]) if r["prefix"] in order else 99, torder.get(r["trait"],9)))

# 输出 TSV
out = os.path.join(ROOT, "step9/results", "i2_heterogeneity_all.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
cols = ["outcome_cn","outcome_en","trait_cn","nsnp","Q","Q_df","Q_pval","tau2","I2_pct","Fw","file"]
with open(out, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("" if r.get(k) is None else (f"{r[k]:.4f}" if isinstance(r[k],float) else r[k])) for k in cols})

print(f"saved: {out} | rows: {len(rows)}\n")
hdr = f"{'outcome':<22}{'trait':<6}{'nsnp':>5}{'Q':>10}{'df':>5}{'Q_p':>12}{'tau2':>9}{'I2%':>8}{'Fw':>8}"
print(hdr); print("-"*len(hdr))
for r in rows:
    qp = "%.3g" % r["Q_pval"] if r["Q_pval"] is not None else "-"
    fw = "%8.1f" % r["Fw"] if r.get("Fw") is not None else "%8s" % "-"
    print(f"{r['outcome_cn']:<22}{r['trait_cn']:<6}{r['nsnp'] or 0:>5}{r['Q']:>10.2f}{r['Q_df']:>5}{qp:>12}{r['tau2']:>9.4f}{(r['I2_pct'] if r['I2_pct'] is not None else float('nan')):>8.1f}{fw}")
