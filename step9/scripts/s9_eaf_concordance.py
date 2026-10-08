# -*- coding: utf-8 -*-
"""P1-10d-i: 计算各结局 harmonised 文件的 EAF 一致性（直接实证变异-暴露相似性）。
|EAF_X - EAF_Y| 的分布：中介数、>0.10、>0.15 计数。只做统计，不改原始值。
"""
import os, glob, sys, io, csv, statistics
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
dirs = ["step3/results", "step4/results"]
TRAIT_CN = {"hypo":"甲减","tsh":"TSH","ft4":"FT4"}

def label(fn):
    b = os.path.basename(fn)
    m = {
      "harmonised_hypo":("主分析 PAD (Sakaue)","hypo"), "harmonised_tsh":("主分析 PAD (Sakaue)","tsh"), "harmonised_ft4":("主分析 PAD (Sakaue)","ft4"),
      "harmonised_finngen_hypo":("FinnGen R13 PAD","hypo"), "harmonised_finngen_tsh":("FinnGen R13 PAD","tsh"), "harmonised_finngen_ft4":("FinnGen R13 PAD","ft4"),
      "harmonised_mvp_meta_hypo":("MVP PAD","hypo"), "harmonised_mvp_meta_tsh":("MVP PAD","tsh"), "harmonised_mvp_meta_ft4":("MVP PAD","ft4"),
      "harmonised_cad_ukb_hypo":("CAD (UKB)","hypo"), "harmonised_cad_ukb_tsh":("CAD (UKB)","tsh"), "harmonised_cad_ukb_ft4":("CAD (UKB)","ft4"),
      "harmonised_cad_nouk_hypo":("CAD (Nikpay)","hypo"), "harmonised_cad_nouk_tsh":("CAD (Nikpay)","tsh"), "harmonised_cad_nouk_ft4":("CAD (Nikpay)","ft4"),
      "harmonised_mi_ukb_hypo":("MI (UKB)","hypo"), "harmonised_mi_ukb_tsh":("MI (UKB)","tsh"), "harmonised_mi_ukb_ft4":("MI (UKB)","ft4"),
      "harmonised_is_laa_hypo":("IS-LAA","hypo"), "harmonised_is_laa_tsh":("IS-LAA","tsh"), "harmonised_is_laa_ft4":("IS-LAA","ft4"),
      "harmonised_is_sv_hypo":("IS-SV","hypo"), "harmonised_is_sv_tsh":("IS-SV","tsh"), "harmonised_is_sv_ft4":("IS-SV","ft4"),
    }
    stem = b[:-4]
    return m.get(stem, (stem, None))

rows = []
for d in dirs:
    for fn in sorted(glob.glob(os.path.join(ROOT, d, "harmonised_*.tsv"))):
        with open(fn, encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            idx = {c:i for i,c in enumerate(hdr)}
            if "eaf_X" not in idx or "eaf_Y" not in idx:
                continue
            gaps = []; n = 0; n_na = 0
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) <= max(idx["eaf_X"], idx["eaf_Y"]):
                    continue
                n += 1
                ex, ey = c[idx["eaf_X"]], c[idx["eaf_Y"]]
                try:
                    gx, gy = float(ex), float(ey)
                    gaps.append(abs(gx-gy))
                except (ValueError, TypeError):
                    n_na += 1
            olab, tr = label(fn)
            if tr is None: continue
            rows.append({
                "outcome": olab, "trait": tr, "trait_cn": TRAIT_CN[tr],
                "n_snp": n, "n_eaf_available": len(gaps), "n_eaf_NA": n_na,
                "median_abs_gap": (statistics.median(gaps) if gaps else None),
                "max_abs_gap": (max(gaps) if gaps else None),
                "n_gap_gt_010": sum(1 for g in gaps if g > 0.10),
                "n_gap_gt_015": sum(1 for g in gaps if g > 0.15),
                "file": os.path.relpath(fn, ROOT).replace("\\","/"),
            })

out = os.path.join(ROOT, "step9/results", "eaf_concordance_all.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
cols = ["outcome","trait_cn","n_snp","n_eaf_available","n_eaf_NA","median_abs_gap","max_abs_gap","n_gap_gt_010","n_gap_gt_015","file"]
with open(out, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
    for r in rows:
        w.writerow({k:( "" if r.get(k) is None else (f"{r[k]:.4f}" if isinstance(r[k],float) else r[k])) for k in cols})

print(f"saved: {out} | rows: {len(rows)}\n")
hdr = f"{'outcome':<18}{'trait':<6}{'nSNP':>6}{'EAFok':>7}{'NA':>5}{'med|Δ|':>9}{'max|Δ|':>9}{'>0.10':>7}{'>0.15':>7}"
print(hdr); print("-"*len(hdr))
for r in rows:
    med = "%.4f"%r["median_abs_gap"] if r["median_abs_gap"] is not None else "-"
    mx = "%.4f"%r["max_abs_gap"] if r["max_abs_gap"] is not None else "-"
    print(f"{r['outcome']:<18}{r['trait_cn']:<6}{r['n_snp']:>6}{r['n_eaf_available']:>7}{r['n_eaf_NA']:>5}{med:>9}{mx:>9}{r['n_gap_gt_010']:>7}{r['n_gap_gt_015']:>7}")
