# -*- coding: utf-8 -*-
"""
Step 9 统计功效计算：回答「阴性结果是不是因为功效不足」。

方法（全部基于**已完成的 MR 分析的实际 SE**，不重新估计）
  IVW 估计的 Wald 统计量 = |β̂| / SE，非中心参数 λ = |β_true| / SE。
    Power(β_true) = Φ(λ − z_{1−α/2}) + Φ(−λ − z_{1−α/2})
    最小可检出效应 MDE(power) = (z_{1−α/2} + z_power) × SE
  α = 0.05 双侧 → z = 1.959964；power 0.80 → 0.841621，0.90 → 1.281552。

报告单位（关键：尺度按**暴露**那一侧定）
  · 正向 MR / MVMR —— 暴露＝甲状腺性状：
      甲减（二分）→ 每「甲减 odds 翻倍」（ln2 = 0.6931 log-OR 单位）；
      TSH / FT4（连续、GWAS 已标准化）→ 每 1 SD；
      结局（PAD/CAD/MI/卒中）均为二分 → 换算后统一报 OR。
  · 反向 MR —— 暴露＝血管疾病（二分）→ 单位「每血管疾病 odds 翻倍」：
      结局为甲减（二分）→ 报 OR；
      结局为 TSH / FT4（连续）→ 报 **SD 变化**，不可报 OR（已分开出列 mde*_sd）。
  换算：效应_报告单位 = 效应_原生 × scale（二分暴露 scale = ln2，连续暴露 scale = 1）。

数据来源（只读，不改数）
  · 正向 MR：step3/results/mr_main_*.json、step4/results/mr_*.json
  · 反向 MR：step6/results/rev_*.json
  · MVMR   ：step7/results/mvmr_*.json
输出：step9/results/power_*.tsv / power_summary.json
"""
import os
import io
import csv
import json
import math

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
RES = os.path.join(ROOT, "step9", "results")
os.makedirs(RES, exist_ok=True)

Z975 = 1.959963984540054
ZB = {0.80: 0.8416212335729143, 0.90: 1.2815515655446004, 0.95: 1.6448536269514722}
LN2 = 0.6931471805599453
TRAITS = ["hypo", "tsh", "ft4"]
TRAIT_CN = {"hypo": "甲减", "tsh": "TSH", "ft4": "FT4"}
# 报告单位换算（★ 暴露是哪一侧，尺度就按哪一侧定）
#   正向 MR / MVMR：暴露＝甲状腺性状
#       甲减（二分）→ 每「甲减 odds 翻倍」（scale = ln2）
#       TSH / FT4（连续，GWAS 已标准化）→ 每 1 SD（scale = 1）
#       结局（PAD/CAD/MI/卒中）均为二分 → 换算后报 OR。
#   反向 MR：暴露＝血管疾病（二分）→ scale = ln2、单位「每血管疾病 odds 翻倍」；
#       结局为甲减（二分）时报 OR；结局为 TSH/FT4（连续）时报 SD 变化 —— **不可报 OR**。
SCALE = {"hypo": LN2, "tsh": 1.0, "ft4": 1.0}
UNIT = {"hypo": "每 odds 翻倍", "tsh": "每 1 SD", "ft4": "每 1 SD"}
VASC_SCALE = LN2
VASC_UNIT = "每血管疾病 odds 翻倍"


def Phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def power_of(beta, se):
    if se is None or se <= 0 or beta is None:
        return None
    lam = abs(beta) / se
    return Phi(lam - Z975) + Phi(-lam - Z975)


def mde(se, pw=0.80):
    if se is None or se <= 0:
        return None
    return (Z975 + ZB[pw]) * se


def J(p):
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def rows_from(path, trait, outcome_label, extras, direction="forward"):
    if not os.path.exists(path):
        return None
    d = J(path)
    fx, rd = d.get("ivw_fixed") or {}, d.get("ivw_random") or {}
    if direction == "reverse":
        # 暴露＝血管疾病（二分）；结局＝甲状腺性状
        exp_scale, exp_unit = VASC_SCALE, VASC_UNIT
        outcome_continuous = (trait != "hypo")
    else:
        # 暴露＝甲状腺性状；结局＝血管疾病（二分）
        exp_scale, exp_unit = SCALE[trait], UNIT[trait]
        outcome_continuous = False
    se_f, se_r = fx.get("se"), rd.get("se")
    out = dict(exposure=TRAIT_CN.get(trait, trait), exposure_key=trait,
               outcome=outcome_label, nsnp=d.get("nsnp"), direction=direction,
               beta=fx.get("beta"), se_fixed=se_f, se_random=se_r,
               p_fixed=fx.get("p"), p_random=rd.get("p"),
               unit=exp_unit, scale=exp_scale,
               outcome_type="连续" if outcome_continuous else "二分")
    out.update(extras)

    def _report(se, pw):
        """MDE 换算到报告单位：二分结局→OR；连续结局→SD 变化。"""
        mv = mde(se, pw)
        if mv is None:
            return None, None, None
        rep = mv * exp_scale
        if outcome_continuous:
            return mv, None, rep
        return mv, math.exp(rep), None

    for pw, tag in ((0.80, "mde80"), (0.90, "mde90")):
        mv, orv, sdv = _report(se_f, pw)
        out[tag + "_native"], out[tag + "_or"], out[tag + "_sd"] = mv, orv, sdv
        _, orvr, sdvr = _report(se_r, pw)
        out[tag + "_or_random"], out[tag + "_sd_random"] = orvr, sdvr
    # 检出指定效应量的功效（一律换算到报告单位后再比）
    if outcome_continuous:
        for sdv in (0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
            b = sdv / exp_scale
            out["pow_sd%.2f" % sdv] = power_of(b, se_f)
            out["pow_sd%.2f_rand" % sdv] = power_of(b, se_r)
    else:
        for orv in (1.05, 1.10, 1.15, 1.20, 1.25):
            b = math.log(orv) / exp_scale
            out["pow_or%.2f" % orv] = power_of(b, se_f)
            out["pow_or%.2f_rand" % orv] = power_of(b, se_r)
    return out


# ---------------------------------------------------------------- 正向 MR
FWD = [
    ("PAD（Sakaue 2021, GCST90018890）", "step3/results/mr_main_%s.json", "主分析·含UKB", "sakaue"),
    ("PAD（FinnGen R13 I9_PAD）", "step4/results/mr_finngen_%s.json", "重叠·甲减失效", "finngen"),
    ("PAD（MVP trans-ancestry meta）", "step4/results/mr_mvp_meta_%s.json", "无重叠", "mvp_meta"),
    ("CAD（Nikpay 2015, GCST003116）", "step4/results/mr_cad_nouk_%s.json", "无重叠", "cad_nouk"),
    ("CAD（UKB Mbatchou 2021, GCST90013864）", "step4/results/mr_cad_ukb_%s.json", "重叠", "cad_ukb"),
    ("MI（UKB Donertas 2021, GCST90038610）", "step4/results/mr_mi_ukb_%s.json", "重叠", "mi_ukb"),
    ("IS-LAA（MEGASTROKE, GCST005840）", "step4/results/mr_is_laa_%s.json", "无重叠", "is_laa"),
    ("IS-SV（MEGASTROKE, GCST005841）", "step4/results/mr_is_sv_%s.json", "无重叠", "is_sv"),
]

# ================================================================ 结局尺度体检
# 同一分析里各结局文件应处于同一尺度（二分→log-OR）。若某结局的仪器 SNP
# **中位 SE** 比其他结局小 ≥5 倍，则该文件不是 log-OR 尺度（典型：线性概率模型），
# 其 β / SE 不能与其他结局并列比 OR —— 需单独标注、不得当作 OR。
# 注意：常数尺度重标定不改变 z 与 p，故显著性与方向结论不受影响，只有效应量受影响。
HARM = {
    "sakaue": "step3/results/harmonised_%s.tsv",
    "finngen": "step4/results/harmonised_finngen_%s.tsv",
    "mvp_meta": "step4/results/harmonised_mvp_meta_%s.tsv",
    "cad_nouk": "step4/results/harmonised_cad_nouk_%s.tsv",
    "cad_ukb": "step4/results/harmonised_cad_ukb_%s.tsv",
    "mi_ukb": "step4/results/harmonised_mi_ukb_%s.tsv",
    "is_laa": "step4/results/harmonised_is_laa_%s.tsv",
    "is_sv": "step4/results/harmonised_is_sv_%s.tsv",
}


def _median(v):
    v = sorted(v)
    return v[len(v) // 2] if v else None


scale_rows = []
for _k, _pat in HARM.items():
    _ses = []
    for _t in TRAITS:
        _p = os.path.join(ROOT, _pat % _t)
        if not os.path.exists(_p):
            continue
        with io.open(_p, encoding="utf-8") as _f:
            for _r in csv.DictReader(_f, delimiter="\t"):
                try:
                    _ses.append(float(_r["se_Y"]))
                except (TypeError, ValueError, KeyError):
                    pass
    scale_rows.append(dict(key=_k, n_se=len(_ses), med_se=_median(_ses)))

_ref = _median([r["med_se"] for r in scale_rows if r["med_se"]])
for _r in scale_rows:
    _r["ref_med_se"] = _ref
    _r["ratio_to_ref"] = (_r["med_se"] / _ref) if _r["med_se"] else None
    _r["flag"] = ("尺度异常·疑似线性概率（β/SE 不可当 log-OR 用）"
                  if (_r["ratio_to_ref"] and _r["ratio_to_ref"] < 0.2) else "")
SCALE_FLAG = {r["key"]: r["flag"] for r in scale_rows}
with open(os.path.join(RES, "scale_check.tsv"), "w", encoding="utf-8", newline="") as _f:
    _w = csv.DictWriter(_f, fieldnames=["key", "n_se", "med_se", "ref_med_se",
                                        "ratio_to_ref", "flag"], delimiter="\t")
    _w.writeheader()
    _w.writerows(scale_rows)

fwd_rows = []
for label, pat, note, k in FWD:
    for t in TRAITS:
        r = rows_from(os.path.join(ROOT, pat % t), t, label, dict(source_note=note))
        if r:
            r["out_key"] = k
            r["scale_flag"] = SCALE_FLAG.get(k, "")
            fwd_rows.append(r)

# ---------------------------------------------------------------- 反向 MR
REV = [
    ("CAD（Nikpay 2015, GCST003116）", "step6/results/rev_cad_nikpay_%s.json", "无重叠"),
    ("PAD（Sakaue 2021, GCST90018890）", "step6/results/rev_pad_sakaue_%s.json", "重叠"),
    ("PAD（FinnGen R13 I9_PAD）", "step6/results/rev_pad_finngen_%s.json", "重叠"),
]
rev_rows = []
for label, pat, note in REV:
    for t in TRAITS:
        r = rows_from(os.path.join(ROOT, pat % t), t, label, dict(source_note=note),
                      direction="reverse")
        if r:
            rev_rows.append(r)

# ---------------------------------------------------------------- MVMR
MVM = [
    ("CAD（Nikpay 2015）", "cad_nouk", "无重叠"),
    ("PAD（MVP trans-ancestry）", "mvp_meta", "无重叠"),
    ("IS-SV（MEGASTROKE）", "is_sv", "无重叠"),
    ("IS-LAA（MEGASTROKE）", "is_laa", "无重叠"),
    ("PAD（Sakaue 2021）", "sakaue", "重叠"),
    ("CAD（UKB Mbatchou 2021）", "cad_ukb", "重叠"),
    ("MI（UKB Donertas 2021）", "mi_ukb", "重叠"),
]
mvm_rows = []
for label, key, note in MVM:
    p = os.path.join(ROOT, "step7", "results", "mvmr_%s.json" % key)
    if not os.path.exists(p):
        continue
    d = J(p)
    nsnp = d.get("nsnp")
    condF = {x["exposure"]: x["F"] for x in d.get("cond_f", [])}
    for t in TRAITS:
        m = d["mvmr"][t]
        sc = SCALE[t]
        r = dict(exposure=TRAIT_CN[t], exposure_key=t, outcome=label, nsnp=nsnp,
                 direction="forward", outcome_type="二分",
                 beta=m["beta"], se_fixed=m["se"], se_random=None,
                 p_fixed=m["p"], p_random=None, unit=UNIT[t], scale=sc,
                 condF=condF.get(t), source_note=note, out_key=key,
                 scale_flag=SCALE_FLAG.get(key, ""))
        mv = mde(m["se"], 0.80)
        r["mde80_native"] = mv
        r["mde80_or"] = math.exp(mv * sc) if mv is not None else None
        for orv in (1.05, 1.10, 1.15, 1.20, 1.25):
            r["pow_or%.2f" % orv] = power_of(math.log(orv) / sc, m["se"])
        mvm_rows.append(r)

# ---------------------------------------------------------------- 工具强度（设计层）
inst_summary = {}
for t in TRAITS:
    p = os.path.join(ROOT, "work", "instruments_%s.tsv" % t)
    with io.open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    r2 = sum(float(x["R2"]) for x in rows)
    Fs = sorted(float(x["F"]) for x in rows)
    inst_summary[t] = dict(n=len(rows), r2=r2, f_median=Fs[len(Fs) // 2],
                           f_min=Fs[0], n_eff=float(rows[0]["n_eff"]))

# ---------------------------------------------------------------- 情景分析（结局样本量放大）
SCEN = {}
for k in (1, 2, 5, 10):
    SCEN["×%d" % k] = {}
    for r in fwd_rows:
        if r["se_fixed"]:
            SCEN["×%d" % k]["%s|%s" % (r["exposure_key"], r["outcome"])] = \
                math.exp(mde(r["se_fixed"] / math.sqrt(k), 0.80) * r["scale"])


def dump_tsv(path, rows):
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


dump_tsv(os.path.join(RES, "power_forward.tsv"), fwd_rows)
dump_tsv(os.path.join(RES, "power_reverse.tsv"), rev_rows)
dump_tsv(os.path.join(RES, "power_mvmr.tsv"), mvm_rows)
with open(os.path.join(RES, "power_summary.json"), "w", encoding="utf-8") as f:
    json.dump(dict(instruments=inst_summary, scenarios=SCEN,
                   n_forward=len(fwd_rows), n_reverse=len(rev_rows), n_mvmr=len(mvm_rows)),
              f, ensure_ascii=False, indent=2)

# ---------------------------------------------------------------- 控制台摘要
print("=" * 92)
print("结局尺度体检（仪器 SNP 的中位 SE；偏离参考 ≥5 倍即判为尺度异常）")
for r in scale_rows:
    print("  %-10s n=%-5d 中位 SE=%.5f  参考=%.5f  比=%.3f  %s" % (
        r["key"], r["n_se"], r["med_se"], r["ref_med_se"], r["ratio_to_ref"],
        ("← " + r["flag"]) if r["flag"] else "OK"))
print()
print("工具强度（设计层）")
for t in TRAITS:
    s = inst_summary[t]
    print("  %-4s 工具 %3d 个  总 R²=%.4f  中位 F=%.1f  min F=%.1f  n_eff=%.0f" % (
        TRAIT_CN[t], s["n"], s["r2"], s["f_median"], s["f_min"], s["n_eff"]))
print()
print("正向 MR：80%% 功效下的最小可检出 OR（固定效应 / 随机效应）")
print("  %-10s %-34s %5s %10s %10s %10s" % ("暴露", "结局", "nSNP", "MDE80(OR)", "MDE80随机", "p"))
for r in sorted(fwd_rows, key=lambda x: (x["exposure_key"], x["outcome"])):
    print("  %-10s %-34s %5s %10s %10s %10s" % (
        r["exposure"], r["outcome"][:34], r["nsnp"],
        "%.3f" % r["mde80_or"] if r["mde80_or"] else "-",
        "%.3f" % r["mde80_or_random"] if r["mde80_or_random"] else "-",
        "%.3g" % r["p_fixed"] if r["p_fixed"] is not None else "-"))
print()
print("反向 MR：80%% 功效下的最小可检出效应（暴露单位＝血管疾病 odds 翻倍）")
for r in sorted(rev_rows, key=lambda x: (x["outcome"], x["exposure_key"])):
    if r["outcome_type"] == "二分":
        mdetxt = "MDE80 = OR %.3f" % r["mde80_or"]
    else:
        mdetxt = "MDE80 = %.4f SD" % r["mde80_sd"]
    print("  %-10s ← %-34s nSNP=%3s  %-18s 观测 p=%-9s [%s]" % (
        r["exposure"], r["outcome"][:34], r["nsnp"], mdetxt,
        "%.3g" % r["p_fixed"], r["source_note"]))
print()
print("MVMR：甲减直接效应 80%% 功效下的最小可检出 OR")
for r in mvm_rows:
    if r["exposure_key"] == "hypo":
        print("  %-30s nSNP=%3s  condF=%5.1f  MDE80=%.3f  观测 β=%+.4f (p=%s)" % (
            r["outcome"][:30], r["nsnp"], r["condF"] or float("nan"),
            r["mde80_or"], r["beta"], "%.3g" % r["p_fixed"]))
print()
print("全部写出 → step9/results/")
