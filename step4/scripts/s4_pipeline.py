# -*- coding: utf-8 -*-
"""
Step 4/5  通用结局复制流水线（MVP PAD 就绪版）

设计目标：MVP PAD（dbGaP phs001672）数据一到，改配置即可跑；也支持任何新结局。

用法
  --probe <path>            打印任意结局文件表头 + 前 3 行（用于确定列映射）
  --outcome <key>           运行指定结局（见 OUTCOMES 配置）
  --trait hypo|tsh|ft4|all  指定暴露（默认 all；hypo 若与结局队列重叠会给出警告）

已在 OUTCOMES 中配置（全部已完成，详见 step4/Step4_README.md）
  sakaue      : GCST90018890（Step 3 主分析）
  finngen_r13 : FinnGen R13 I9_PAD（Step 4 复制）
  mvp_meta    : MVP PAD trans-ancestry meta（Step 4 复制；列名已按实测填好，非占位）
  cad_nouk / cad_ukb / mi_ukb / is_laa / is_sv : 泛血管谱 5 结局（Step 5）

分族裔层（MVP EUR / AFR / HIS）：**不纳入本研究范围**（2026-10-07 决策）
  理由：① 分层病例数过小（AFR 5,373 例 / HIS 1,925 例 / EUR 24,009 例）；
        ② 三套暴露工具均基于 EUR 构建，用于 AFR/HIS 属跨族裔外推，偏倚不可控；
        ③ trans-ancestry meta（mvp_meta，31,307 例 / 211,753 对照，n_eff≈109,098）已综合三族裔信息。
  故 OUTCOMES 中 mvp_eur / mvp_afr / mvp_his 三块整块注释保留，仅供将来重启时参考，勿误读为「待跑」。

关键口径（与 Step 3/4 一致）
  · 对齐一律 rsID 优先、位置兜底（位置一律用面板坐标=GRCh38）
  · 等位基因：集合匹配 → 链互补翻转（非回文无歧义）→ 回文 |ΔEAF|>0.1 剔除 / 非回文 >0.2 剔除
  · 统计：mrlib（IVW 固定/随机 DL、Egger 规范形式、加权中位数/众数、Q、PRESSO 等价、方向性）
"""
import os, sys, csv, gzip, json, math, argparse, datetime

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
sys.path.insert(0, os.path.join(ROOT, "step3", "scripts"))
from mrlib import (p_norm, p_chi2, ivw, ivw_random_dl, egger, weighted_median,
                   weighted_mode, ivw_bootstrap_se, presso_equiv, p_binom_two)

WORK = os.path.join(ROOT, "work")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
TRAITS = ["hypo", "tsh", "ft4"]
TRAIT_CN = {"hypo": "甲减", "tsh": "TSH", "ft4": "FT4"}

# harmonised 中间文件列（与 s4_finngen_replication.py 的输出对齐，便于 s4_figures.py 复用）
HCOL = ["SNP", "rs_id", "chr", "A1", "A2", "beta_X", "se_X", "eaf_X", "F_X", "R2_X",
        "n_eff_X", "beta_Y", "se_Y", "eaf_Y", "p_Y", "beta_XY", "se_XY", "gene_Y", "action"]
SCOL = ["trait", "trait_cn", "overlap", "method", "nsnp", "beta", "se", "or", "p", "note"]

# ---------------------------------------------------------------- 结局配置
# cols 使用统一语义键；若目标文件列名不同，只改这里
#
# GWAS-SSF harmonised 标准列（Step 5 五个泛血管谱结局文件实测完全一致）：
#   hm_rsid / hm_chrom / hm_pos / hm_effect_allele / hm_other_allele / hm_beta
#   standard_error（注意：无 hm_ 前缀）/ hm_effect_allele_frequency / p_value
# 要点：hm_beta 已由 GWAS Catalog 对齐到 hm_effect_allele（翻转动作记在 hm_code），
#       必须与 hm_effect_allele 配套使用才自洽，不要改取原始 beta 列。
SSF_COLS = dict(rsid="hm_rsid", chr="hm_chrom", pos="hm_pos",
                ea="hm_effect_allele", oa="hm_other_allele", beta="hm_beta",
                se="standard_error", eaf="hm_effect_allele_frequency", p="p_value")

OUTCOMES = {
    "sakaue": dict(
        path=os.path.join(ROOT, "data", "outcome", "pad_GCST90018890.h.tsv.gz"),
        sep="\t", gz=True, header=True,
        cols=dict(rsid="rsid", chr="chromosome", pos="base_pair_location",
                  ea="effect_allele", oa="other_allele", beta="beta",
                  se="standard_error", eaf="effect_allele_frequency", p="p_value"),
        n_case=7114, n_control=475964, label="Sakaue2021 PAD (EUR)", overlap=None),
    "finngen_r13": dict(
        path=os.path.join(ROOT, "data", "replication", "finngen_R13_I9_PAD.gz"),
        sep="\t", gz=True, header=True,
        cols=dict(rsid="rsids", chr="#chrom", pos="pos", ea="alt", oa="ref",
                  beta="beta", se="sebeta", eaf="af_alt", p="pval"),
        n_case=22244, n_control=460490, label="FinnGen R13 I9_PAD",
        overlap={"hypo": "OVERLAP_INVALID(暴露含FinnGen DF10)"}),
    # ↓↓↓ MVP：te 已按实测填好；EUR/AFR/HIS 待文件到手后填 path（列名以 --probe 输出为准）
    # mvp_meta 列映射依据 raw/verify_mvp_te_cols.py 实测（214,813 行）：
    #   表头 SNP_ID/Chromosome/Position/Allele1/Allele2/EA/EAF/SampleSize/Effect/SE/PValue/Direction
    #   实测 EA 恒等于 Allele1（214813/214813 = 100.0000%）→ oa 取 Allele2 安全
    #   SampleSize 恒为 243,060 = 31,307 例 + 211,753 对照 ✓
    "mvp_meta": dict(
        path=os.path.join(ROOT, "raw", "mvp_pad", "MVP.te.PAD.dbGAP.txt.gz"),
        sep="\t", gz=True, header=True,
        cols=dict(rsid="SNP_ID", chr="Chromosome", pos="Position",
                  ea="EA", oa="Allele2", beta="Effect", se="SE",
                  eaf="EAF", p="PValue"),
        n_case=31307, n_control=211753, label="MVP PAD (trans-ancestry meta)",
        overlap=None),
    # ============ 分族裔层（EUR / AFR / HIS）—— 不纳入本研究范围，整块注释 ============
    # 决策（2026-10-07）：不做分族裔分层验证。理由见文件头 docstring。
    # 备案（若将来重启）：
    #   dbGaP 表型号：pha004959(EUR, 24,009 例/150,983 对照) / pha004958(AFR, 5,373/42,485)
    #                 / pha004960(HIS, 1,925/18,285)；Name 字段被误标为 "CVD"，须看 Description。
    #   其公开截断版为 14 列（SNP ID/Marker Type/P-value/Chr ID/Chr Position/Rank/
    #   Allele1/Allele2/Sample size/Effect Allele/Odds ratio/CI low/CI high/Bin ID）
    #   —— 若全量同构，则只有 OR + 95% CI、无 beta/SE，需换算：
    #   beta=ln(OR), SE=(ln(CI_high)-ln(CI_low))/(2*1.96)，并给管线加转换钩子。
    #
    # "mvp_eur": dict(
    #     path=r"<MVP_PAD_EUR_FILE>", sep="\t", gz=True, header=True,
    #     cols=dict(rsid="RSID?", chr="CHR?", pos="POS?", ea="EA?", oa="NEA?",
    #               beta="BETA?", se="SE?", eaf="EAF?", p="P?"),
    #     n_case=24009, n_control=150983, label="MVP PAD (EUR)", overlap=None),
    # "mvp_afr": dict(
    #     path=r"<MVP_PAD_AFR_FILE>", sep="\t", gz=True, header=True,
    #     cols=dict(rsid="RSID?", chr="CHR?", pos="POS?", ea="EA?", oa="NEA?",
    #               beta="BETA?", se="SE?", eaf="EAF?", p="P?"),
    #     n_case=5373, n_control=42485, label="MVP PAD (AFR)", overlap=None),
    # "mvp_his": dict(
    #     path=r"<MVP_PAD_HIS_FILE>", sep="\t", gz=True, header=True,
    #     cols=dict(rsid="RSID?", chr="CHR?", pos="POS?", ea="EA?", oa="NEA?",
    #               beta="BETA?", se="SE?", eaf="EAF?", p="P?"),
    #     n_case=1925, n_control=18285, label="MVP PAD (HIS)", overlap=None),

    # ================ Step 5 泛血管谱（EBI GWAS Catalog harmonised 版） ================
    # 重叠判定依据：暴露 GWAS（Rand 2025, PMID 41238958）含 UKB + FinnGen Freeze 10。
    #   · cad_ukb / mi_ukb 均源自 UKB → 与暴露共享 UKB 样本 → 标 OVERLAP_UKB，仅作扩展分析
    #   · cad_nouk（Nikpay 2015, PMID 26343387, CARDIoGRAMplusC4D，2015 年，UKB 未参与）
    #   · is_laa / is_sv（MEGASTROKE, Malik 2018, PMID 29531354，卒中联盟，无 UKB）
    #     → 三者均为无重叠的干净分析
    # 注意：cad_ukb 的 hm_effect_allele_frequency 整列为 NA（实测 100% 缺失），
    #       harmonise 已支持 EAF 可选；本项目工具集回文位点为 0，无判别歧义。
    "cad_ukb": dict(
        path=os.path.join(ROOT, "data", "panvascular", "cad_GCST90013864.h.tsv.gz"),
        sep="\t", gz=True, header=True, cols=SSF_COLS,
        n_case=29339, n_control=322724,
        label="CAD (Mbatchou2021, UKB; GCST90013864)",
        overlap={t: "OVERLAP_UKB(暴露GWAS含UKB)" for t in TRAITS}),
    "cad_nouk": dict(
        path=os.path.join(ROOT, "data", "panvascular", "cad_GCST003116.h.tsv.gz"),
        sep="\t", gz=True, header=True, cols=SSF_COLS,
        n_case=60801, n_control=123504,
        label="CAD (Nikpay2015 CARDIoGRAMplusC4D; GCST003116)",
        overlap=None),
    "mi_ukb": dict(
        path=os.path.join(ROOT, "data", "panvascular", "mi_GCST90038610.h.tsv.gz"),
        sep="\t", gz=True, header=True, cols=SSF_COLS,
        n_case=11081, n_control=473517,
        label="MI (Donertas2021, UKB; GCST90038610)",
        overlap={t: "OVERLAP_UKB(暴露GWAS含UKB)" for t in TRAITS}),
    "is_laa": dict(
        path=os.path.join(ROOT, "data", "panvascular", "is_laa_GCST005840.h.tsv.gz"),
        sep="\t", gz=True, header=True, cols=SSF_COLS,
        n_case=4373, n_control=406111,
        label="IS-LAA (MEGASTROKE Malik2018; GCST005840)",
        overlap=None),
    "is_sv": dict(
        path=os.path.join(ROOT, "data", "panvascular", "is_sv_GCST005841.h.tsv.gz"),
        sep="\t", gz=True, header=True, cols=SSF_COLS,
        n_case=5386, n_control=192662,
        label="IS-SV (MEGASTROKE Malik2018; GCST005841)",
        overlap=None),
}

def _open(path, gz):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if gz else \
           open(path, "r", encoding="utf-8", errors="replace")

def probe(path, gz=True):
    with _open(path, gz) as f:
        for i, line in enumerate(f):
            print(line.rstrip())
            if i >= 2: break

def load_instruments(t):
    rows = list(csv.DictReader(open(os.path.join(WORK, "instruments_%s.tsv" % t),
                                    encoding="utf-8"), delimiter="\t"))
    for r in rows:
        r["beta"] = float(r["beta"]); r["se"] = float(r["se"]); r["eaf"] = float(r["EAF"])
        r["F"] = float(r["F"]); r["R2"] = float(r["R2"]); r["n_eff"] = float(r["n_eff"])
    return rows

def extract_subset(spec, inst_all, out_path):
    """流式扫描结局文件，提取命中工具（rsID 优先，位置兜底），落盘缓存"""
    c = spec["cols"]
    rs_set, pos_set = set(), set()
    for t in TRAITS:
        for r in inst_all[t]:
            for k in (r.get("rs_id_gwas"), r.get("rsid_gwas")):
                if k: rs_set.add(k)
            pos_set.add(r["SNP"])
    hdr_needed = {k: c[k] for k in ("rsid","chr","pos","ea","oa","beta","se","eaf","p")}
    n_hit = 0
    with _open(spec["path"], spec["gz"]) as f, \
         open(out_path + ".tmp", "w", encoding="utf-8", newline="") as o:
        rd = csv.reader(f, delimiter=spec["sep"])
        hdr = next(rd)
        ix = {}
        for k, name in hdr_needed.items():
            if name in hdr: ix[k] = hdr.index(name)
            else: raise KeyError("结局文件缺少列 %r（可用列：%s）" % (name, hdr))
        w = csv.writer(o, delimiter="\t")
        w.writerow(["rsid","chr","pos","ea","oa","beta","se","eaf","p"])
        for row in rd:
            rsid = row[ix["rsid"]]
            key = None
            for x in (rsid or "").split(","):
                if x in rs_set: key = x; break
            if key is None:
                if "%s:%s" % (row[ix["chr"]], row[ix["pos"]]) in pos_set:
                    key = "POS:%s:%s" % (row[ix["chr"]], row[ix["pos"]])
            if key is None: continue
            w.writerow([key, row[ix["chr"]], row[ix["pos"]], row[ix["ea"]].upper(),
                        row[ix["oa"]].upper(), row[ix["beta"]], row[ix["se"]],
                        row[ix["eaf"]], row[ix["p"]]])
            n_hit += 1
    os.replace(out_path + ".tmp", out_path)
    return n_hit

def harmonise(t, inst, sub_rows):
    ok, st = [], dict(total=0, kept=0, no_outcome=0, allele_mismatch=0, complement=0,
                      flipped=0, palindromic_ambiguous=0, eaf_discordant=0,
                      eaf_missing=0, palindromic_no_eaf=0)
    for r in inst:
        st["total"] += 1
        y = None
        for k in (r.get("rs_id_gwas"), r.get("rsid_gwas")):
            if k and k in sub_rows and not k.startswith("POS:"): y = sub_rows[k]; break
        if y is None: y = sub_rows.get("POS:" + r["SNP"])
        if y is None: st["no_outcome"] += 1; continue
        try:
            bY, seY = float(y["beta"]), float(y["se"])
        except (TypeError, ValueError):
            st["no_outcome"] += 1; continue
        if not (seY > 0) or not math.isfinite(bY): st["no_outcome"] += 1; continue
        # EAF 可选：部分 harmonised 文件（如 GCST90013864）EAF 整列为 NA
        try:
            eafY = float(y["eaf"]); have_eaf = True
        except (TypeError, ValueError):
            eafY = None; have_eaf = False
            st["eaf_missing"] += 1
        e1, o1 = r["EA"], r["OA"]; e2, o2 = y["ea"], y["oa"]
        acted = "straight"
        if {e1, o1} != {e2, o2}:
            if e1 in COMP and o1 in COMP and {COMP[e1], COMP[o1]} == {e2, o2}:
                e2, o2 = COMP[e2], COMP[o2]; st["complement"] += 1; acted = "complement"
            else:
                st["allele_mismatch"] += 1; continue
        if e2 != e1:
            bY = -bY
            if have_eaf:
                eafY = 1.0 - eafY
            st["flipped"] += 1
            acted = "flip" if acted == "straight" else acted + "_flip"
        pal = (e1 + o1) in ("AT", "TA", "CG", "GC")
        if have_eaf:
            diff = abs(r["eaf"] - eafY)
            if pal and diff > 0.1: st["palindromic_ambiguous"] += 1; continue
            if not pal and diff > 0.2: st["eaf_discordant"] += 1; continue
        elif pal:
            # 无 EAF 时无法判别回文位点链方向 → 保守剔除（本项目工具集回文位点为 0）
            st["palindromic_no_eaf"] += 1; continue
        ok.append(dict(SNP=r["SNP"], rs_id=r.get("rs_id_gwas") or r.get("rsid_gwas") or "",
                       chr=str(y["chr"]), A1=e1, A2=o1,
                       beta_X=r["beta"], se_X=r["se"], eaf_X=r["eaf"], F_X=r["F"],
                       R2_X=r["R2"], n_eff_X=r["n_eff"],
                       beta_Y=bY, se_Y=seY, eaf_Y=(eafY if have_eaf else y["eaf"]),
                       p_Y=y["p"], beta_XY=bY / r["beta"], se_XY=seY / abs(r["beta"]),
                       gene_Y="", action=acted))
        st["kept"] += 1
    return ok, st

def mr_report(rows, n_case, n_eff_out, label, trait, overlap_note=None):
    bX = [r["beta_X"] for r in rows]; sX = [r["se_X"] for r in rows]
    bY = [r["beta_Y"] for r in rows]; sY = [r["se_Y"] for r in rows]
    bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
    k = len(rows)
    bF, seF = ivw(bR, sR); rnd = ivw_random_dl(bR, sR); eg = egger(bX, bY, sY)
    bm = weighted_median(bR, sR); sem = ivw_bootstrap_se(bR, sR, "median")
    bmo = weighted_mode(bR, sR); semo = ivw_bootstrap_se(bR, sR, "mode", nb=500)
    npos = sum(1 for x in bR if x > 0)
    r2x = sum(r["R2_X"] for r in rows)
    r2y = sum((r["beta_Y"]/r["se_Y"])**2 / ((r["beta_Y"]/r["se_Y"])**2 + n_eff_out) for r in rows)
    return {"label": label, "trait": trait, "nsnp": k, "n_case": n_case, "n_eff_out": n_eff_out,
            "overlap_note": overlap_note,
            "ivw_fixed": {"beta": bF, "se": seF, "p": p_norm(bF/seF)},
            "ivw_random": rnd, "mr_egger": eg,
            "weighted_median": {"beta": bm, "se": sem, "p": p_norm(bm/sem)},
            "weighted_mode": {"beta": bmo, "se": semo, "p": p_norm(bmo/semo)},
            "presso_equiv": presso_equiv(bR, sR),
            "sign_test": {"n_positive": npos, "n_negative": k-npos,
                          "p": p_binom_two(min(npos, k-npos), k)},
            "directionality": {"R2_X": r2x, "R2_Y": r2y, "nR2_X": rows[0]["n_eff_X"]*r2x,
                               "nR2_Y": n_eff_out*r2y,
                               "direction": "forward (X->Y)" if rows[0]["n_eff_X"]*r2x > n_eff_out*r2y else "reverse?"}}

def write_harmonised(outcome_key, t, rows):
    """落盘 harmonise 后的 MR 输入（列结构与 FinnGen 版一致，供 s4_figures.py 复用）"""
    p = os.path.join(ROOT, "step4", "results", "harmonised_%s_%s.tsv" % (outcome_key, t))
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HCOL, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    return p


def write_loo(outcome_key, t, rows):
    """留一法：逐 SNP 剔除后重估随机效应 IVW（口径与 step3/s4_finngen 完全一致）"""
    bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
    k = len(rows); loo = []
    for i in range(k):
        keep = [j for j in range(k) if j != i]
        r2 = ivw_random_dl([bR[j] for j in keep], [sR[j] for j in keep])
        loo.append({"SNP": rows[i]["SNP"], "nsnp": k - 1,
                    "beta": r2["beta"], "se": r2["se"], "p": r2["p"]})
    p = os.path.join(ROOT, "step4", "results", "loo_%s_%s.tsv" % (outcome_key, t))
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["SNP", "nsnp", "beta", "se", "p"], delimiter="\t")
        w.writeheader(); w.writerows(loo)
    return p


def summary_rows_for(t, rep, overlap_label):
    """把一份 rep 展开成 FinnGen 版同构的汇总行"""
    out = []
    for method, key, note in [
            ("IVW(固定效应)", "ivw_fixed", ""),
            ("IVW(随机效应 DL)", "ivw_random", "Q=%.2f df=%d Q_p=%.3g tau2=%.4g" % (
                rep["ivw_random"]["Q"], rep["ivw_random"]["Q_df"],
                rep["ivw_random"]["Q_pval"], rep["ivw_random"]["tau2"])),
            ("MR-Egger", "mr_egger", "intercept=%.4g p=%.3g" % (
                rep["mr_egger"]["intercept"], rep["mr_egger"]["intercept_pval"])),
            ("加权中位数", "weighted_median", "bootstrap se 2000"),
            ("加权众数", "weighted_mode", "bootstrap se 500")]:
        d = rep.get(key)
        if not d:
            continue
        out.append({"trait": t, "trait_cn": TRAIT_CN.get(t, t), "overlap": overlap_label,
                    "method": method, "nsnp": rep["nsnp"],
                    "beta": round(d["beta"], 6), "se": round(d["se"], 6),
                    "or": round(math.exp(d["beta"]), 6), "p": "%.4g" % d["p"], "note": note})
    pe = rep.get("presso_equiv") or {}
    out.append({"trait": t, "trait_cn": TRAIT_CN.get(t, t), "overlap": overlap_label,
                "method": "MR-PRESSO等价", "nsnp": rep["nsnp"], "beta": "", "se": "", "or": "",
                "p": "%.4g" % pe["global_Q_pval"] if pe.get("global_Q_pval") else "",
                "note": ("离群 %d 个, 校正beta=%s p=%s, distortion p=%s" % (
                    pe["distortion"]["n_outliers"], round(pe["distortion"]["beta_corrected"], 6),
                    "%.4g" % pe["distortion"]["p_corrected"],
                    "%.4g" % pe["distortion"]["distortion_p"]))
                if pe.get("distortion") else "未检出离群位点"})
    sg = rep.get("sign_test") or {}
    out.append({"trait": t, "trait_cn": TRAIT_CN.get(t, t), "overlap": overlap_label,
                "method": "符号一致性", "nsnp": rep["nsnp"], "beta": "", "se": "", "or": "",
                "p": "%.4g" % (sg.get("p") if sg.get("p") is not None
                               else sg.get("p_binomial_two_sided") or float("nan")),
                "note": "正 %s / 负 %s" % (sg.get("n_positive"), sg.get("n_negative"))})
    di = rep.get("directionality") or {}
    out.append({"trait": t, "trait_cn": TRAIT_CN.get(t, t), "overlap": overlap_label,
                "method": "方向性(R²对比)", "nsnp": rep["nsnp"], "beta": "", "se": "", "or": "",
                "p": "", "note": "R2X=%.4g R2Y=%.4g nR2X=%.0f nR2Y=%.0f → %s" % (
                    di["R2_X"], di["R2_Y"], di["nR2_X"], di["nR2_Y"], di["direction"])})
    return out


def meta_ivw(estimates):
    """对多族裔的因果估计做 IVW 随机效应 meta（DerSimonian-Laird）"""
    b = [e["beta"] for e in estimates]; s = [e["se"] for e in estimates]
    return ivw_random_dl(b, s) if len(b) > 1 else {"beta": b[0], "se": s[0], "p": p_norm(b[0]/s[0])}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe"); ap.add_argument("--outcome"); ap.add_argument("--trait", default="all")
    a = ap.parse_args()
    if a.probe:
        probe(a.probe, gz=a.probe.endswith(".gz")); return
    if not a.outcome: ap.error("需要 --outcome 或 --probe")
    spec = OUTCOMES[a.outcome]
    if "<" in spec["path"]:
        print("⚠ 结局文件路径尚未配置：%s（请填写 OUTCOMES[%r] 后重试）" % (spec["path"], a.outcome)); return
    traits = TRAITS if a.trait == "all" else [a.trait]
    t0 = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    logp = os.path.join(ROOT, "step4", "logs", "pipeline_%s_%s.log" % (a.outcome, t0))
    def log(m):
        print(m)
        with open(logp, "a", encoding="utf-8") as f: f.write(m + "\n")
    inst_all = {t: load_instruments(t) for t in TRAITS}
    sub_path = os.path.join(WORK, "outcome_sub_%s.tsv" % a.outcome)
    if not os.path.exists(sub_path):
        n = extract_subset(spec, inst_all, sub_path)
        log("[%s] 扫描完成，命中 %d 行 → %s" % (spec["label"], n, sub_path))
    sub_rows = {}
    for d in csv.DictReader(open(sub_path, encoding="utf-8"), delimiter="\t"):
        if not d["rsid"].startswith("POS:"): sub_rows.setdefault(d["rsid"], d)
        sub_rows.setdefault("POS:" + d["chr"] + ":" + d["pos"], d)
    n_eff_out = 4.0*spec["n_case"]*spec["n_control"]/(spec["n_case"]+spec["n_control"])
    summary_all = []
    for t in traits:
        note = (spec.get("overlap") or {}).get(t)
        rows, st = harmonise(t, inst_all[t], sub_rows)
        if not rows:
            log("%s × %s：无可用工具，跳过" % (t, spec["label"])); continue
        rep = mr_report(rows, spec["n_case"], n_eff_out, spec["label"], t, note)
        rep["harmonise"] = st
        out = os.path.join(ROOT, "step4", "results", "mr_%s_%s.json" % (a.outcome, t))
        json.dump(rep, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
        hp = write_harmonised(a.outcome, t, rows)
        lp = write_loo(a.outcome, t, rows)
        summary_all += summary_rows_for(t, rep, note or "VALID")
        log("%s × %s [%s]: %d SNP | IVW固定 β=%.4f (p=%.3g) | IVW随机 p=%.3g | Egger 截距 p=%.3g | Q p=%.3g%s"
            % (t, spec["label"], note or "OK", rep["nsnp"], rep["ivw_fixed"]["beta"],
               rep["ivw_fixed"]["p"], rep["ivw_random"]["p"],
               rep["mr_egger"]["intercept_pval"] if rep["mr_egger"] else float("nan"),
               rep["ivw_random"]["Q_pval"], ("  ⚠ " + note) if note else ""))
        log("    harmonised → %s" % os.path.basename(hp))
        log("    LOO(%d 行) → %s" % (rep["nsnp"], os.path.basename(lp)))
    if summary_all:
        sp = os.path.join(ROOT, "step4", "results", "mr_%s_summary.tsv" % a.outcome)
        with open(sp, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=SCOL, delimiter="\t")
            w.writeheader(); w.writerows(summary_all)
        log("汇总 → %s" % os.path.basename(sp))

if __name__ == "__main__":
    main()
