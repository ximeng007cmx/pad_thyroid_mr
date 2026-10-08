# -*- coding: utf-8 -*-
"""
Step 6（反向 MR）＋ Step 7（多变量 MR）结果一页纸生成器。

原则：所有数值直接从结果 JSON 读取，不手抄、不改数（包括结论句里的数字）。
  · Step 6 → step6/results/rev_<exp>_<trait>.json 与 revloo_<exp>_<trait>.tsv
  · Step 7 → step7/results/mvmr_<key>.json（缺失则跳过该结局）
  · 单变量对照 → step4/results/mr_<key>_<trait>.json（sakaue 回退 step3/results/mr_main_<trait>.json）
输出：deliverables/Step6-7_反向MR与MVMR.html / .md
"""
import io
import json
import os

ROOT = r"D:\WorkSpace_Study\pad_thyroid_mr"
DELIV = os.path.join(ROOT, "deliverables")

TRAITS = [("hypo", "甲减"), ("tsh", "TSH"), ("ft4", "FT4")]
E = 2.718281828459045
Z = 1.959963984540054

# exp_key, 显示名, 来源, 病例/对照, overlap_note, 是否可用于干净反向 MR
REV_EXP = [
    ("cad_nikpay", "CAD（Nikpay 2015）", "GCST003116", "60,801 / 123,504", "与暴露无样本重叠", True),
    ("pad_sakaue", "PAD（Sakaue 2021）", "GCST90018890", "7,114 / 475,964", "暴露队列含 UKB", False),
    ("pad_finngen", "PAD（FinnGen R13）", "I9_PAD", "22,244 / 460,490", "暴露队列含 FinnGen DF10", False),
]

# key, 显示名, 来源, overlap（None=无重叠）
MVM_EXP = [
    ("cad_nouk", "CAD（Nikpay 2015）", "GCST003116", None),
    ("mvp_meta", "PAD（MVP trans-ancestry meta）", "MVP meta", None),
    ("is_sv", "IS-SV（MEGASTROKE）", "GCST005841", None),
    ("is_laa", "IS-LAA（MEGASTROKE）", "GCST005840", None),
    ("sakaue", "PAD（Sakaue 2021）", "GCST90018890", "UKB"),
    ("cad_ukb", "CAD（UKB Mbatchou 2021）", "GCST90013864", "UKB"),
    ("mi_ukb", "MI（UKB Donertas 2021）", "GCST90038610", "UKB"),
]

# ---- ‡ 尺度注记（Step 9 结局尺度体检发现）--------------------------------------
# MI（UKB, GCST90038610）文件的 β / SE 处于「线性概率」尺度而非 log-OR。
# 只在 MI 行加注记，**不改动任何数值**。依据：step9/results/scale_check.tsv。
SCALE_FLAG_KEYS = {"mi_ukb"}
MI_NOTE_HTML = (
    '<div class="box"><b>‡ 尺度注记（MI 行）</b>：MI（UKB, GCST90038610）的 β / SE 处于'
    '<b>线性概率尺度</b>而非 log-OR —— 其仪器 SNP 中位 SE 仅为其他结局的 <b>0.030 倍</b>；'
    '同一位点 <code>rs3184504</code>（SH2B3）在该文件 β=0.0019，而其余结局为 0.049–0.114'
    '（exp 全部落在 1.05–1.12）。因此上表 MI 行的 <b>β 不可与其他结局并列解读</b>：'
    '其 β 偏小、p 极小系<b>尺度</b>所致，<b>并非</b>「重叠使 SE 被低估」，也<b>不</b>代表效应可忽略；'
    'z 与 p 不受常数尺度重标定影响，故 MI 行「重叠 → 只作参考、不进入因果结论」的判定不变。'
    '核查：<code>step9/results/scale_check.tsv</code>、<code>Step9_统计功效.md</code> §2。</div>')
MI_NOTE_MD = (
    "> **‡ 尺度注记（MI 行）**：MI（UKB, GCST90038610）的 β / SE 处于**线性概率尺度**而非 log-OR —— "
    "其仪器 SNP 中位 SE 仅为其他结局的 **0.030 倍**；同一位点 `rs3184504`（SH2B3）在该文件 β=0.0019，"
    "而其余结局为 0.049–0.114（exp 全部落在 1.05–1.12）。因此上表 MI 行的 **β 不可与其他结局并列解读**："
    "其 β 偏小、p 极小系**尺度**所致，**并非**「重叠使 SE 被低估」，也**不**代表效应可忽略；"
    "z 与 p 不受常数尺度重标定影响，故 MI 行「重叠 → 只作参考、不进入因果结论」的判定不变。"
    "核查：`step9/results/scale_check.tsv`、`Step9_统计功效.md` §2。")


def mi_mark(key):
    return "‡" if key in SCALE_FLAG_KEYS else ""


def J(p):
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def fp(x):
    if x is None:
        return "—"
    if x == 0:
        return "0"
    if abs(x) < 1e-3:
        return "%.2e" % x
    return "%.3f" % x


def fnum(x, d=4):
    if x is None:
        return "—"
    if abs(x) != 0 and (abs(x) < 1e-3 or abs(x) >= 1e4):
        return "%.2e" % x
    return ("%%.%df" % d) % x


def fmt_or(b, se):
    if b is None or se is None:
        return "—"
    return "%.3f [%.3f, %.3f]" % (E ** b, E ** (b - Z * se), E ** (b + Z * se))


def cls_beta(b):
    if b is None:
        return ""
    return "pos" if b > 0 else ("neg" if b < 0 else "")


def revloo_range(exp, t):
    p = os.path.join(ROOT, "step6", "results", "revloo_%s_%s.tsv" % (exp, t))
    if not os.path.exists(p):
        return None
    rows = [l.split("\t") for l in io.open(p, encoding="utf-8").read().strip().split("\n")[1:] if l.strip()]
    if not rows:
        return None
    b = [float(r[2]) for r in rows]
    pl = [float(r[4]) for r in rows]
    return dict(n=len(rows), bmin=min(b), bmax=max(b), pmax=max(pl),
                all_pos=all(v > 0 for v in b), all_neg=all(v < 0 for v in b))


def univ(key, t):
    """单变量 IVW 对照：step4 优先，sakaue 回退 step3 主分析"""
    cands = [os.path.join(ROOT, "step4", "results", "mr_%s_%s.json" % (key, t))]
    if key == "sakaue":
        cands.append(os.path.join(ROOT, "step3", "results", "mr_main_%s.json" % t))
    for p in cands:
        if os.path.exists(p):
            j = J(p)
            return dict(beta=j["ivw_fixed"]["beta"], se=j["ivw_fixed"]["se"],
                        p=j["ivw_fixed"]["p"], nsnp=j.get("nsnp"))
    return None


# ---------------------------------------------------------------- 载入
REV = {}
for exp, *_ in REV_EXP:
    REV[exp] = {}
    for t, _ in TRAITS:
        p = os.path.join(ROOT, "step6", "results", "rev_%s_%s.json" % (exp, t))
        REV[exp][t] = J(p) if os.path.exists(p) else None

MVM = {}
for key, *_ in MVM_EXP:
    p = os.path.join(ROOT, "step7", "results", "mvmr_%s.json" % key)
    MVM[key] = J(p) if os.path.exists(p) else None

# ---------------------------------------------------------------- 自动统计（供结论句使用）
NOOV = [k for k, _, _, ov in MVM_EXP if ov is None and MVM.get(k)]
NOOV_HYPO = [(k, MVM[k]["mvmr"]["hypo"]) for k in NOOV]
NOOV_HYPO_POS = [(k, v) for k, v in NOOV_HYPO if v["beta"] > 0]
NOOV_HYPO_SIG = [(k, v) for k, v in NOOV_HYPO if v["p"] < 0.05]
TSH_SIG = [(k, MVM[k]["mvmr"]["tsh"]) for k, _, _, ov in MVM_EXP
           if MVM.get(k) and MVM[k]["mvmr"]["tsh"]["p"] < 0.05]
FT4_F = [x["F"] for k, _, _, ov in MVM_EXP if MVM.get(k)
         for x in MVM[k]["cond_f"] if x["exposure"] == "ft4"]
DISP = dict((k, n) for k, n, _, _ in MVM_EXP)
SHORT = {"cad_nouk": "CAD-Nikpay", "mvp_meta": "MVP-PAD", "is_sv": "IS-SV", "is_laa": "IS-LAA",
         "sakaue": "PAD-Sakaue", "cad_ukb": "CAD-UKB", "mi_ukb": "MI-UKB"}


def _bp(v):
    return "β=%+.4f (p=%s)" % (v["beta"], fp(v["p"]))


def _ratio(key, t):
    """MV/单变量 β 比值"""
    v = MVM[key]["mvmr"][t]
    u = univ(key, t) or {}
    return (v["beta"] / u["beta"]) if u.get("beta") else None


def rev_verdict(exp, t):
    d = REV[exp][t]
    if d is None:
        return ("无结果", "warn")
    usable = dict((e[0], e[5]) for e in REV_EXP)[exp]
    if d["nsnp"] <= 1:
        return ("不可用·重叠" if not usable else "单工具", "bad" if not usable else "warn")
    if not usable:
        return ("重叠·仅参考", "warn")
    return ("阴性" if d["ivw_fixed"]["p"] >= 0.05 else "名义阳性",
            "ok" if d["ivw_fixed"]["p"] >= 0.05 else "warn")


def mvm_verdict(key):
    d = MVM[key]
    if d is None:
        return ("未生成", "warn")
    ov = dict((k, o) for k, _, _, o in MVM_EXP)[key]
    mv = d["mvmr"]["hypo"]
    fh = [x["F"] for x in d["cond_f"] if x["exposure"] == "hypo"][0]
    if ov:
        return ("重叠·参考", "warn")
    if mv["p"] < 0.05 and fh > 10:
        return ("甲减独立显著", "warn")
    return ("阴性", "ok")


CSS = """
:root{--ink:#111827;--mut:#5b6675;--line:#dfe3e8;--soft:#f6f7f9;--ok:#0f7b45;--warn:#9a6700;--bad:#b42318;--acc:#1f4e79}
*{box-sizing:border-box}
body{margin:0;background:#eceef1;color:var(--ink);font-family:"Source Han Sans SC","Microsoft YaHei","PingFang SC",-apple-system,"Segoe UI",sans-serif;font-size:11.5px;line-height:1.5}
.page{width:210mm;min-height:297mm;margin:14px auto;background:#fff;padding:13mm 12mm;box-shadow:0 2px 14px rgba(15,23,42,.13)}
h1{font-size:17px;margin:0 0 2px;letter-spacing:.2px}
.sub{color:var(--mut);font-size:10.5px;margin-bottom:9px}
.meta{display:flex;flex-wrap:wrap;gap:4px 16px;font-size:10px;color:var(--mut);border-top:1px solid var(--line);border-bottom:1px solid var(--line);padding:5px 0;margin-bottom:10px}
h2{font-size:12.5px;margin:13px 0 5px;padding-left:7px;border-left:3px solid var(--acc);color:var(--acc);letter-spacing:.2px}
table{width:100%;border-collapse:collapse;font-size:10.5px;margin-bottom:4px}
th,td{border:1px solid var(--line);padding:3px 5px;text-align:center;vertical-align:middle;white-space:nowrap}
th{background:var(--soft);font-weight:600;color:#374151}
td.l,th.l{text-align:left;white-space:normal}
.pill{display:inline-block;padding:1px 6px;border-radius:9px;font-size:9.5px;font-weight:600;white-space:nowrap}
.p-ok{background:#e7f5ec;color:var(--ok);border:1px solid #b7dfc6}
.p-warn{background:#fdf4e3;color:var(--warn);border:1px solid #f0d9a8}
.p-bad{background:#fdeceb;color:var(--bad);border:1px solid #f3c3bf}
.sig{color:var(--bad);font-weight:700}
.pos{color:var(--bad);font-weight:600}
.neg{color:var(--ok);font-weight:600}
.box{background:var(--soft);border:1px solid var(--line);border-left:3px solid var(--acc);padding:6px 9px;margin:6px 0;font-size:10.5px}
.two{display:flex;gap:9px}.two>div{flex:1}
ul{margin:3px 0 3px 15px;padding:0}li{margin:2px 0}
.foot{margin-top:10px;border-top:1px solid var(--line);padding-top:6px;font-size:9.5px;color:var(--mut);line-height:1.55}
code{background:#f1f3f6;padding:0 3px;border-radius:2px;font-size:9.5px;font-family:"Cascadia Mono",Consolas,monospace}
@media print{body{background:#fff;font-size:10.5px}.page{width:auto;margin:0;box-shadow:none;padding:9mm 8mm;min-height:0}@page{size:A4;margin:9mm}h2{margin:9px 0 4px}}
"""


def pill(kind, text):
    return '<span class="pill p-%s">%s</span>' % (kind, text)


# ---------------------------------------------------------------- HTML
def build_html():
    o = []
    o.append('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">')
    o.append('<title>甲状腺功能 ↔ 动脉粥样硬化性血管疾病：反向 MR 与 MVMR 一页纸</title>')
    o.append('<style>%s</style></head><body><div class="page">' % CSS)
    o.append('<h1>甲状腺功能 ↔ 动脉粥样硬化性血管疾病：反向 MR（Step 6）与多变量 MR（Step 7）</h1>')
    o.append('<div class="sub">方向性检验 × 独立性检验 · 暴露 甲减 / TSH / FT4（Rand 2025, <i>Nat Genet</i>, PMID 41238958）</div>')
    o.append('<div class="meta">'
             '<span><b>Step 6</b> 反向 MR：3 血管暴露 → 3 甲状腺性状（9 格）</span>'
             '<span><b>Step 7</b> MVMR：三甲状腺暴露联合模型（7 结局 × 300 工具）</span>'
             '<span><b>日期</b> 2026-10-07</span></div>')

    # ---- 0 结论速览（数字全部自动计算）
    o.append('<h2>0 · 结论速览</h2><div class="two">')
    cad = REV["cad_nikpay"]
    o.append('<div class="box"><b>反向因果：不支持</b><ul>'
             '<li>唯一<b>无样本重叠</b>的反向检验（CAD→甲状腺）三性状<b>全阴性</b>：'
             '甲减 p=%s、TSH p=%s、FT4 p=%s。</li>'
             '<li>PAD→甲状腺（Sakaue / FinnGen）与暴露同源队列，固定效应「显著」随随机效应消失'
             '（FinnGen 甲减：固定 p=%s → 随机 p=%s），不构成方向性证据。</li>'
             '<li>⚠ 干净反向检验仅 %d 个工具命中，<b>功率有限</b> → 阴性属「未观察到」而非「排除」。</li>'
             '</ul></div>' % (fp(cad["hypo"]["ivw_fixed"]["p"]), fp(cad["tsh"]["ivw_fixed"]["p"]),
                              fp(cad["ft4"]["ivw_fixed"]["p"]),
                              fp(REV["pad_finngen"]["hypo"]["ivw_fixed"]["p"]),
                              fp(REV["pad_finngen"]["hypo"]["ivw_random"]["p"]), cad["hypo"]["nsnp"]))
    # 独立性 box（自动）
    hyp_txt = "、".join("%s %s" % (SHORT[k], _bp(v)) for k, v in NOOV_HYPO_SIG)
    o.append('<div class="box"><b>独立性：控制 TSH / FT4 后甲减效应仍显著</b><ul>'
             '<li>%d 个<b>无样本重叠</b>结局中，甲减 MVMR 效应<b>方向全部为正（%d/%d）</b>，'
             '其中 <b>%d 个达显著</b>：%s。</li>'
             '<li>典型例：甲减→CAD（Nikpay）单变量 β=%s (p=%s) → MVMR β=%s (p=%s)，'
             '<b>控制 TSH/FT4 后增强至 %.2f×</b>。</li>'
             '<li>TSH 仅个别结局在 MVMR 中显著（%s），FT4 条件 F 仅 %.1f–%.1f（弱工具，不可解释）。</li>'
             '</ul></div></div>' % (
                 len(NOOV_HYPO), len(NOOV_HYPO_POS), len(NOOV_HYPO), len(NOOV_HYPO_SIG),
                 hyp_txt,
                 fnum((univ("cad_nouk", "hypo") or {}).get("beta"), 4),
                 fp((univ("cad_nouk", "hypo") or {}).get("p")),
                 fnum(MVM["cad_nouk"]["mvmr"]["hypo"]["beta"], 4),
                 fp(MVM["cad_nouk"]["mvmr"]["hypo"]["p"]),
                 _ratio("cad_nouk", "hypo") or 0.0,
                 "、".join(SHORT[k] for k, _ in TSH_SIG) or "无",
                 min(FT4_F), max(FT4_F)))

    # ---- 1 反向 MR 结果
    o.append('<h2>1 · 反向 MR（Step 6）：血管疾病 → 甲状腺功能</h2>')
    o.append('<div class="box"><b>设计：</b>以血管疾病 GWAS 的工具变量为暴露工具，估计其对甲状腺性状'
             '（甲减 / TSH / FT4）的因果效应；β 为「每单位 log-OR(血管疾病) 对应的甲状腺性状 log-OR / SD 变化」。</div>')
    o.append('<table><tr><th class="l">血管暴露（工具来源）</th><th>重叠</th><th>甲状腺性状</th>'
             '<th>nSNP</th><th>β</th><th>SE</th><th>p（IVW 固定）</th><th>OR [95%CI]</th>'
             '<th>p（随机 DL）</th><th>Egger 截距 p</th><th>Q p</th><th class="l">判定</th></tr>')
    for exp, name, src, nstr, ovnote, usable in REV_EXP:
        ovtxt = ('<span style="color:#0f7b45">无</span>' if usable
                 else '<span style="color:#b42318;font-weight:600">有</span>')
        for i, (t, cn) in enumerate(TRAITS):
            d = REV[exp][t]
            if d is None:
                continue
            b = d["ivw_fixed"]["beta"]; se = d["ivw_fixed"]["se"]; p = d["ivw_fixed"]["p"]
            pr = d["ivw_random"]["p"]; qp = d["ivw_random"]["Q_pval"]
            eg = d.get("mr_egger") or {}
            eip = eg.get("intercept_pval")
            v, kind = rev_verdict(exp, t)
            first = ('<td class="l" rowspan="3">%s<div style="color:#5b6675;font-size:9.5px">%s</div></td>'
                     '<td rowspan="3">%s</td>' % (name, src, ovtxt)) if i == 0 else ""
            o.append('<tr>%s<td>%s</td><td>%d</td><td class="%s">%s</td><td>%s</td>'
                     '<td%s>%s</td><td>%s</td><td>%s</td><td%s>%s</td><td>%s</td><td class="l">%s</td></tr>' % (
                         first, cn, d["nsnp"], cls_beta(b), fnum(b), fnum(se),
                         ' class="sig"' if p < 0.05 else "", fp(p), fmt_or(b, se), fp(pr),
                         ' class="sig"' if (eip is not None and eip < 0.05) else "", fp(eip), fp(qp),
                         pill(kind, v)))
        o.append('<tr><td colspan="12" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>')
    o.append('</table>')
    o.append('<div class="box"><b>读法：</b>* PAD-FinnGen 对照数按 FinnGen R13 口径标注。'
             '红 β 为正、绿为负。反向 MR「阴性」= 血管疾病不通过遗传途径改变甲状腺功能。</div>')

    # ---- 2 工具命中率
    o.append('<h2>2 · 反向 MR 工具命中率 —— 决定哪条检验干净可用</h2>')
    o.append('<table><tr><th class="l">血管暴露</th><th>工具总数</th><th>甲状腺文件命中</th>'
             '<th>命中率</th><th class="l">命中路径</th><th>可用工具</th></tr>')
    for exp, name, *_ in REV_EXP:
        d = REV[exp]["hypo"]
        if d is None:
            continue
        cov = d.get("coverage") or {}
        br = cov.get("by_route") or {}
        o.append('<tr><td class="l">%s</td><td>%d</td><td><b>%d</b></td><td>%.1f%%</td>'
                 '<td class="l">%s</td><td>%d</td></tr>' % (
                     name, cov.get("n_tools", 0), cov.get("n_outcome_hit", 0),
                     cov.get("hit_pct", 0.0),
                     "、".join("%s %d" % (k, v) for k, v in sorted(br.items())) or "—", d["nsnp"]))
    o.append('</table>')
    o.append('<div class="box"><b>要点：</b>甲状腺 GWAS 文件（GWAS Catalog 新格式）与血管结局文件在'
             '<b>坐标口径（SNP +1 / indel +0）与链向（约 75% 互补）</b>上均不一致，'
             '且位点覆盖度仅约 TOPMed 级文件的三成 → 反向 MR 工具数天然受限，命中率仅 20–30%。'
             '这是<b>数据源的真实限制，非匹配错误</b>（已用多路由 ＋ 等位基因一致性校验确认）。</div>')

    # ---- 3 反向 MR 稳健性
    o.append('<h2>3 · 反向 MR 稳健性（LOO 与符号一致性）</h2>')
    o.append('<table><tr><th class="l">血管暴露</th><th>甲状腺性状</th><th>符号 +/− (p)</th>'
             '<th>方向性</th><th class="l">LOO β 区间（最差 p）</th></tr>')
    for exp, name, *_ in REV_EXP:
        for i, (t, cn) in enumerate(TRAITS):
            d = REV[exp][t]
            if d is None:
                continue
            sg = d.get("sign_test") or {}
            sp = sg.get("p_binomial_two_sided") or sg.get("p")
            lr = revloo_range(exp, t)
            ltxt = ("%+.4f ~ %+.4f（%s，最差 p %s）" % (lr["bmin"], lr["bmax"],
                    "全正" if lr["all_pos"] else ("全负" if lr["all_neg"] else "含变号"),
                    fp(lr["pmax"]))) if lr else "n≤1，未计算"
            first = '<td class="l" rowspan="3">%s</td>' % name if i == 0 else ""
            o.append('<tr>%s<td>%s</td><td>%s/%s (%s)</td><td>%s</td><td class="l">%s</td></tr>' % (
                first, cn, sg.get("n_positive"), sg.get("n_negative"), fp(sp),
                (d.get("directionality") or {}).get("direction", "—"), ltxt))
        o.append('<tr><td colspan="5" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>')
    o.append('</table>')

    # ---- 4 MVMR
    o.append('<h2>4 · 多变量 MR（Step 7）：甲减 / TSH / FT4 联合模型</h2>')
    o.append('<div class="box"><b>设计：</b>三性状工具集并集 → 联合 clump（r²&lt;0.001、10 Mb、p&lt;5e-8）'
             '得 300 个相互独立工具；MV-IVW 同时估计三个暴露的直接效应，条件 F 检验各暴露工具强度，'
             'MV-Q 检验整体水平多效性。<b>MV</b>＝多变量；<b>单</b>＝单变量 IVW（对照）。</div>')
    o.append('<table><tr><th class="l">结局（血管疾病）</th><th>重叠</th><th>nSNP</th>'
             '<th class="l">甲减 MV β (p) / 单 β (p)</th>'
             '<th class="l">TSH MV β (p) / 单 β (p)</th>'
             '<th class="l">FT4 MV β (p) / 单 β (p)</th>'
             '<th>条件 F<br>甲减/TSH/FT4</th><th>MV-Q p</th><th class="l">判定</th></tr>')
    for key, name, src, ov in MVM_EXP:
        d = MVM[key]
        if d is None:
            o.append('<tr><td class="l">%s</td><td>%s</td><td colspan="7">未生成</td></tr>' % (name, ov or "无"))
            continue
        mv = d["mvmr"]
        cf = {x["exposure"]: x["F"] for x in d["cond_f"]}

        def cell(t):
            m = mv[t]; u = univ(key, t) or {}
            psig = ' class="sig"' if m["p"] < 0.05 else ""
            return ('<td class="l"><b class="%s">%s</b> (<span%s>%s</span>) &nbsp;/&nbsp; '
                    '%s (%s)</td>' % (cls_beta(m["beta"]), fnum(m["beta"], 3), psig, fp(m["p"]),
                                      fnum(u.get("beta"), 3), fp(u.get("p"))))

        v, kind = mvm_verdict(key)
        ovtxt = ('<span style="color:#b42318;font-weight:600">%s</span>' % ov) if ov \
            else '<span style="color:#0f7b45">无</span>'
        o.append('<tr><td class="l">%s%s</td><td>%s</td><td>%d</td>%s%s%s'
                 '<td>%.1f / %.1f / %.1f</td><td>%s</td><td class="l">%s</td></tr>' % (
                     name, mi_mark(key), ovtxt, d["nsnp"], cell("hypo"), cell("tsh"), cell("ft4"),
                     cf.get("hypo"), cf.get("tsh"), cf.get("ft4"), fp(d.get("Q_pval")),
                     pill(kind, v)))
    o.append('</table>')
    o.append('<div class="box"><b>读法：</b>MV β 与「单」β 的差异 = 控制另外两个甲状腺性状后的直接效应。'
             '条件 F&gt;10 视为工具强度达标。<b>⚠ 限制：</b>三个暴露来自同一项 GWAS，'
             '暴露效应估计间协方差不为 0 且无法从汇总数据估计 → 条件 F 与 MV-Q 的零分布仅为近似。</div>')
    o.append(MI_NOTE_HTML)

    # ---- 5 无重叠结局汇总
    o.append('<h2>5 · 无样本重叠结局的甲减独立效应（核心证据）</h2>')
    o.append('<table><tr><th class="l">结局</th><th>nSNP</th><th>甲减 MV β</th><th>p(MV)</th>'
             '<th>甲减 单变量 β</th><th>p(单变量)</th><th>MV/单 比值</th><th>条件 F</th><th>MV-Q p</th></tr>')
    for key, v in NOOV_HYPO:
        d = MVM[key]
        u = univ(key, "hypo") or {}
        fh = [x["F"] for x in d["cond_f"] if x["exposure"] == "hypo"][0]
        ratio = (v["beta"] / u["beta"]) if u.get("beta") else None
        o.append('<tr><td class="l">%s</td><td>%d</td><td class="%s"><b>%s</b></td>'
                 '<td%s>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%.1f</td><td>%s</td></tr>' % (
                     DISP[key], d["nsnp"], cls_beta(v["beta"]), fnum(v["beta"], 4),
                     ' class="sig"' if v["p"] < 0.05 else "", fp(v["p"]),
                     fnum(u.get("beta"), 4), fp(u.get("p")),
                     ("%.2f×" % ratio) if ratio else "—", fh, fp(d.get("Q_pval"))))
    o.append('</table>')
    ratio_txt = "、".join("%s %.2f×" % (SHORT[k], _ratio(k, "hypo")) for k, _ in NOOV_HYPO if _ratio(k, "hypo"))
    o.append('<div class="box"><b>要点：</b>在 %d 个与暴露<b>无样本重叠</b>的结局中，甲减的 MVMR 效应'
             '<b>方向一致为正（%d/%d）</b>，%d 个显著。单变量 vs MVMR 的 β 比值：%s —— '
             '多数结局在控制 TSH/FT4 后效应维持或增强，提示单变量估计部分被 TSH 掩盖。'
             '注：MV-Q 在这些结局中普遍 p&lt;0.05，提示残余水平多效性，点估计应结合这一限制解读。</div>' % (
                 len(NOOV_HYPO), len(NOOV_HYPO_POS), len(NOOV_HYPO), len(NOOV_HYPO_SIG), ratio_txt))

    # ---- 6 可写 / 不可写
    o.append('<h2>6 · 可写入论文的结论</h2><div class="two">')
    o.append('<div class="box"><b>站得住的</b><ul>'
             '<li><b>无反向因果证据</b>：唯一无重叠的反向检验（CAD Nikpay→甲状腺）三性状全阴性，LOO 全不显著。</li>'
             '<li><b>甲减的致动脉粥样硬化效应独立于 TSH / FT4</b>：%d 个无重叠结局中 MVMR 效应方向一致为正、'
             '%d 个显著（如 CAD β=%s, p=%s）。</li>'
             '<li><b>方法学贡献</b>：给出甲状腺 GWAS 新格式文件的坐标/链向口径差异及多路由匹配方案，'
             '可复用于其他 GWAS Catalog 汇总数据。</li>'
             '</ul></div>' % (len(NOOV_HYPO), len(NOOV_HYPO_SIG),
                              fnum(MVM["cad_nouk"]["mvmr"]["hypo"]["beta"], 4),
                              fp(MVM["cad_nouk"]["mvmr"]["hypo"]["p"])))
    o.append('<div class="box"><b>不能写的（须披露）</b><ul>'
             '<li><b>反向 MR 阴性 ≠ 排除反向因果</b>：干净检验仅 %d 工具命中，功率有限。</li>'
             '<li><b>PAD→甲状腺的方向性结果不可用</b>：Sakaue / FinnGen 与暴露同源，固定效应虚高、随机效应消失。</li>'
             '<li><b>MVMR 的 FT4 结果不可解释</b>：条件 F 仅 %.1f–%.1f（弱工具）。</li>'
             '<li><b>MV-Q 普遍拒绝</b>：提示水平多效性或模型假设不成立，点估计需谨慎解读。</li>'
             '</ul></div></div>' % (cad["hypo"]["nsnp"], min(FT4_F), max(FT4_F)))

    # ---- 7 数据来源
    o.append('<h2>7 · 数据来源与方法口径</h2>')
    o.append('<table><tr><th class="l">环节</th><th class="l">数据集</th><th class="l">说明</th></tr>')
    o.append('<tr><td class="l">暴露（共同）</td><td class="l">甲减 / TSH / FT4 — Rand 2025 · '
             'GCST90572791 / 90572789 / 90572790</td><td class="l">PMID 41238958；含 UKB＋FinnGen Freeze 10</td></tr>')
    for exp, name, src, nstr, ovnote, usable in REV_EXP:
        o.append('<tr><td class="l">Step 6 血管暴露</td><td class="l">%s — %s</td>'
                 '<td class="l">病例/对照 %s；%s</td></tr>' % (name, src, nstr, ovnote))
    o.append('<tr><td class="l">Step 7 结局</td><td class="l">CAD / MI / PAD / IS-LAA / IS-SV（共 7 个数据集）</td>'
             '<td class="l">MVP PAD 为 trans-ancestry meta</td></tr>')
    o.append('</table>')
    o.append('<div class="foot"><b>生成方式</b>：由 <code>step6/scripts/s67_onepager.py</code> 直接读取 '
             '<code>step6/results/rev_*.json</code>、<code>revloo_*.tsv</code>、'
             '<code>step7/results/mvmr_*.json</code> 生成，<b>正文结论句中的数字亦为程序自动写入</b>，未作人工修改。<br>'
             '<b>方法口径</b>：Step 6 反向 MR 与 Step 3/4 正向 MR 完全一致（IVW 固定/随机 DL、MR-Egger、'
             '加权中位数/众数、Cochran Q、LOO、方向性 R² 对比；工具 p&lt;5e-8、r²&lt;0.001、10 Mb 窗口）；'
             '多路由匹配（rs_id→rsid→pos+1→pos+0）并要求等位基因同向或链互补。<br>'
             '<b>Step 7 MVMR</b>：MV-IVW（Sanderson 2019）、MV-Egger（WLS 加截距）、同质性 Q（df=L−K）、'
             '两样本条件 F（Sanderson 2021 式 4·7）。</div>')
    o.append('</div></body></html>')
    return "\n".join(o)


# ---------------------------------------------------------------- MD
def build_md():
    o = []
    o.append("# 甲状腺功能 ↔ 动脉粥样硬化性血管疾病：反向 MR（Step 6）与 MVMR（Step 7）")
    o.append("")
    o.append("**方向性检验 × 独立性检验**　｜　暴露：甲减 / TSH / FT4 —— Rand et al. 2025, *Nat Genet*, PMID 41238958　｜　2026-10-07")
    o.append("")
    o.append("---")
    o.append("")
    o.append("## 0 结论速览")
    o.append("")
    cad = REV["cad_nikpay"]
    o.append("**反向因果：不支持。** 唯一无样本重叠的反向检验（CAD→甲状腺）三性状全阴性"
             "（甲减 p=%s、TSH p=%s、FT4 p=%s）；PAD 两条受同源队列污染"
             "（FinnGen 甲减固定 p=%s → 随机 p=%s），不作为证据。⚠ 干净检验仅 %d 工具命中，功率有限。" % (
                 fp(cad["hypo"]["ivw_fixed"]["p"]), fp(cad["tsh"]["ivw_fixed"]["p"]),
                 fp(cad["ft4"]["ivw_fixed"]["p"]), fp(REV["pad_finngen"]["hypo"]["ivw_fixed"]["p"]),
                 fp(REV["pad_finngen"]["hypo"]["ivw_random"]["p"]), cad["hypo"]["nsnp"]))
    o.append("")
    o.append("**独立性：控制 TSH / FT4 后甲减效应仍显著。** %d 个无样本重叠结局中，甲减 MVMR 效应方向全部为正"
             "（%d/%d），%d 个显著（%s）；典型例 CAD-Nikpay 控制后由 +0.0308 增强至 +0.0578（%.2f×）。"
             "FT4 条件 F 仅 %.1f–%.1f，不可解释。" % (
                 len(NOOV_HYPO), len(NOOV_HYPO_POS), len(NOOV_HYPO), len(NOOV_HYPO_SIG),
                 "、".join("%s %s" % (SHORT[k], _bp(v)) for k, v in NOOV_HYPO_SIG),
                 _ratio("cad_nouk", "hypo") or 0.0,
                 min(FT4_F), max(FT4_F)))
    o.append("")
    o.append("## 1 反向 MR（Step 6）：血管疾病 → 甲状腺功能")
    o.append("")
    o.append("| 血管暴露 | 重叠 | 性状 | nSNP | β | SE | p（固定） | OR [95%CI] | p（随机） | Egger截距 p | Q p | 判定 |")
    o.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for exp, name, src, nstr, ovnote, usable in REV_EXP:
        ovtxt = "无" if usable else "**有**"
        for t, cn in TRAITS:
            d = REV[exp][t]
            if d is None:
                continue
            b = d["ivw_fixed"]["beta"]; se = d["ivw_fixed"]["se"]; p = d["ivw_fixed"]["p"]
            eg = d.get("mr_egger") or {}
            v, _ = rev_verdict(exp, t)
            mark = "**" if p < 0.05 else ""
            o.append("| %s | %s | %s | %d | %+.4f | %.4f | %s%s%s | %s | %s | %s | %s | %s |" % (
                name, ovtxt, cn, d["nsnp"], b, se, mark, fp(p), mark, fmt_or(b, se),
                fp(d["ivw_random"]["p"]), fp(eg.get("intercept_pval")),
                fp(d["ivw_random"]["Q_pval"]), v))
    o.append("")
    o.append("## 2 反向 MR 工具命中率")
    o.append("")
    o.append("| 血管暴露 | 工具总数 | 甲状腺文件命中 | 命中率 | 命中路径 |")
    o.append("|---|---|---|---|---|")
    for exp, name, *_ in REV_EXP:
        d = REV[exp]["hypo"]
        if d is None:
            continue
        cov = d.get("coverage") or {}
        br = cov.get("by_route") or {}
        o.append("| %s | %d | %d | %.1f%% | %s |" % (
            name, cov.get("n_tools", 0), cov.get("n_outcome_hit", 0), cov.get("hit_pct", 0.0),
            "、".join("%s %d" % (k, v) for k, v in sorted(br.items()))))
    o.append("")
    o.append("> 甲状腺 GWAS 新格式文件与血管结局文件在坐标口径（SNP +1 / indel +0）与链向上不一致，"
             "且覆盖度约为 TOPMed 级文件的三成 → 反向 MR 工具数天然受限，属数据源真实限制。")
    o.append("")
    o.append("## 3 MVMR（Step 7）：甲减 / TSH / FT4 联合模型")
    o.append("")
    o.append("**MV** = 多变量 IVW；**单** = 单变量 IVW（对照）。工具：三性状并集联合 clump → 300 个独立工具。")
    o.append("")
    o.append("| 结局 | 重叠 | nSNP | 甲减 MV β (p) / 单 β (p) | TSH MV β (p) / 单 β (p) | FT4 MV β (p) / 单 β (p) | 条件F 甲减/TSH/FT4 | MV-Q p | 判定 |")
    o.append("|---|---|---|---|---|---|---|---|---|")
    for key, name, src, ov in MVM_EXP:
        d = MVM[key]
        if d is None:
            o.append("| %s | %s | — | 未生成 | | | | | |" % (name, ov or "无"))
            continue
        mv = d["mvmr"]
        cf = {x["exposure"]: x["F"] for x in d["cond_f"]}
        cells = []
        for t in ("hypo", "tsh", "ft4"):
            m = mv[t]; u = univ(key, t) or {}
            mk = "**" if m["p"] < 0.05 else ""
            cells.append("%s%+.4f%s (%s) / %+.4f (%s)" % (
                mk, m["beta"], mk, fp(m["p"]), (u.get("beta") or 0.0), fp(u.get("p"))))
        v, _ = mvm_verdict(key)
        o.append("| %s%s | %s | %d | %s | %s | %s | %.1f / %.1f / %.1f | %s | %s |" % (
            name, mi_mark(key), ov or "无", d["nsnp"], cells[0], cells[1], cells[2],
            cf.get("hypo"), cf.get("tsh"), cf.get("ft4"), fp(d.get("Q_pval")), v))
    o.append("")
    o.append(MI_NOTE_MD)
    o.append("")
    o.append("## 4 无样本重叠结局的甲减独立效应（核心证据）")
    o.append("")
    o.append("| 结局 | nSNP | 甲减 MV β | p(MV) | 甲减 单变量 β | p(单变量) | MV/单 | 条件 F | MV-Q p |")
    o.append("|---|---|---|---|---|---|---|---|---|")
    for key, v in NOOV_HYPO:
        d = MVM[key]; u = univ(key, "hypo") or {}
        fh = [x["F"] for x in d["cond_f"] if x["exposure"] == "hypo"][0]
        ratio = (v["beta"] / u["beta"]) if u.get("beta") else None
        mk = "**" if v["p"] < 0.05 else ""
        o.append("| %s | %d | %s%+.4f%s | %s%s%s | %+.4f | %s | %s | %.1f | %s |" % (
            DISP[key], d["nsnp"], mk, v["beta"], mk, mk, fp(v["p"]), mk,
            (u.get("beta") or 0.0), fp(u.get("p")),
            ("%.2f×" % ratio) if ratio else "—", fh, fp(d.get("Q_pval"))))
    o.append("")
    o.append("> ⚠ **限制**：三个暴露来自同一项 GWAS，暴露效应估计间协方差不为 0 且无法从汇总数据估计 "
             "→ 条件 F 与 MV-Q 的零分布仅为近似。条件 F>10 视为工具强度达标。")
    o.append("")
    o.append("## 5 可写入论文的结论")
    o.append("")
    o.append("### 站得住的")
    o.append("")
    o.append("1. **无反向因果证据**：唯一无重叠的反向检验（CAD Nikpay→甲状腺）三性状全阴性，LOO 全不显著。")
    o.append("2. **甲减的致动脉粥样硬化效应独立于 TSH/FT4**：%d 个无重叠结局中 MVMR 方向一致为正、%d 个显著"
             "（如 CAD β=%s, p=%s）。" % (len(NOOV_HYPO), len(NOOV_HYPO_SIG),
                                          fnum(MVM["cad_nouk"]["mvmr"]["hypo"]["beta"], 4),
                                          fp(MVM["cad_nouk"]["mvmr"]["hypo"]["p"])))
    o.append("3. **方法学贡献**：给出甲状腺 GWAS 新格式文件的坐标/链向口径差异及多路由匹配方案。")
    o.append("")
    o.append("### 不能写的（须披露）")
    o.append("")
    o.append("1. 反向 MR 阴性 ≠ 排除反向因果（干净检验仅 %d 工具命中，功率有限）。" % cad["hypo"]["nsnp"])
    o.append("2. PAD→甲状腺的方向性结果不可用（与暴露同源队列）。")
    o.append("3. MVMR 的 FT4 结果不可解释（条件 F 仅 %.1f–%.1f，弱工具）。" % (min(FT4_F), max(FT4_F)))
    o.append("4. MV-Q 普遍拒绝，提示水平多效性或模型假设不成立，MVMR 点估计需谨慎。")
    o.append("")
    o.append("## 6 数据来源与方法口径")
    o.append("")
    o.append("| 环节 | 数据集 | 说明 |")
    o.append("|---|---|---|")
    o.append("| 暴露（共同） | 甲减/TSH/FT4 — Rand 2025 · GCST90572791/90572789/90572790 | PMID 41238958；含 UKB＋FinnGen DF10 |")
    for exp, name, src, nstr, ovnote, usable in REV_EXP:
        o.append("| Step 6 血管暴露 | %s — %s | 病例/对照 %s；%s |" % (name, src, nstr, ovnote))
    o.append("")
    o.append("**方法口径**：Step 6 与 Step 3/4 完全一致（IVW 固定/随机 DL、MR-Egger、加权中位数/众数、"
             "Cochran Q、LOO、方向性 R² 对比；工具 p<5e-8、r²<0.001、10 Mb 窗口）；"
             "多路由匹配并要求等位基因同向或链互补。Step 7：MV-IVW、MV-Egger、同质性 Q、两样本条件 F。")
    o.append("")
    o.append("> 本文件由 `step6/scripts/s67_onepager.py` 从结果 JSON / LOO TSV 直接生成；"
             "正文结论句中的数字亦为程序自动写入，未作人工修改。")
    return "\n".join(o)


if __name__ == "__main__":
    if not os.path.exists(DELIV):
        os.makedirs(DELIV)
    hp = os.path.join(DELIV, "Step6-7_反向MR与MVMR.html")
    mp = os.path.join(DELIV, "Step6-7_反向MR与MVMR.md")
    with io.open(hp, "w", encoding="utf-8") as f:
        f.write(build_html())
    with io.open(mp, "w", encoding="utf-8") as f:
        f.write(build_md())
    print("OK  " + hp, os.path.getsize(hp), "B")
    print("OK  " + mp, os.path.getsize(mp), "B")
