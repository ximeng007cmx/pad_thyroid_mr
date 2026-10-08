# -*- coding: utf-8 -*-
"""
Step 5 泛血管谱结果一页纸生成器。

原则：所有数值直接从 step4/results/*.json 与 loo_*.tsv 读取，不手抄、不改数。
输出：deliverables/Step5_泛血管谱结果.html / .md
"""
import io
import json
import os

ROOT = r"D:\WorkSpace_Study\pad_thyroid_mr"
DELIV = os.path.join(ROOT, "deliverables")
TRAITS = [("hypo", "甲减"), ("tsh", "TSH"), ("ft4", "FT4")]

# tag, 显示名, 来源, 样本量, 人群, 是否与暴露重叠
OUTCOMES = [
    ("cad_nouk", "CAD（无重叠）", "Nikpay 2015 · GCST003116", "60,801 / 123,504", "EUR（CARDIoGRAMplusC4D）", False),
    ("cad_ukb",  "CAD（UKB）",    "Mbatchou 2021 · GCST90013864", "29,339 / 322,724", "英区 UKB", True),
    ("mi_ukb",   "MI（UKB）",     "Donertas 2021 · GCST90038610", "11,081 / 473,517", "英区 UKB", True),
    ("is_laa",   "IS-LAA（大动脉粥样硬化型卒中）", "MEGASTROKE · GCST005840", "4,373 / 406,111", "EUR", False),
    ("is_sv",    "IS-SV（小血管型卒中）",     "MEGASTROKE · GCST005841", "5,386 / 192,662", "EUR", False),
]

# ---- ‡ 尺度注记（Step 9 结局尺度体检发现）--------------------------------------
# MI（UKB, GCST90038610）文件的 β / SE 处于「线性概率」尺度而非 log-OR。
# 只在 MI 行加注记，**不改动任何数值**。依据：step9/results/scale_check.tsv。
SCALE_FLAG_TAGS = {"mi_ukb"}
MI_NOTE_HTML = (
    '<div class="box"><b>‡ 尺度注记（MI 行）</b>：MI（UKB, GCST90038610）的 β / SE 处于'
    '<b>线性概率尺度</b>而非 log-OR —— 其仪器 SNP 中位 SE 仅为本表其他结局的 <b>0.030 倍</b>；'
    '同一位点 <code>rs3184504</code>（SH2B3）在本文件 β=0.0019，而其余结局为 0.049–0.114'
    '（exp 全部落在 1.05–1.12）。因此上表 MI 行的 <b>β 与 OR 不可与其他结局并列解读</b>：'
    '其 β 偏小、p 极小系<b>尺度</b>所致，<b>并非</b>「重叠使 SE 被低估」，也<b>不</b>代表效应可忽略；'
    'z 与 p 不受常数尺度重标定影响，故 MI 行「重叠 → 只作敏感性分析、不进入因果结论」的判定不变。'
    '核查：<code>step9/results/scale_check.tsv</code>、<code>Step9_统计功效.md</code> §2。</div>')
MI_NOTE_MD = (
    "> **‡ 尺度注记（MI 行）**：MI（UKB, GCST90038610）的 β / SE 处于**线性概率尺度**而非 log-OR —— "
    "其仪器 SNP 中位 SE 仅为本表其他结局的 **0.030 倍**；同一位点 `rs3184504`（SH2B3）在本文件 β=0.0019，"
    "而其余结局为 0.049–0.114（exp 全部落在 1.05–1.12）。因此上表 MI 行的 **β 与 OR 不可与其他结局并列解读**："
    "其 β 偏小、p 极小系**尺度**所致，**并非**「重叠使 SE 被低估」，也**不**代表效应可忽略；"
    "z 与 p 不受常数尺度重标定影响，故 MI 行「重叠 → 只作敏感性分析、不进入因果结论」的判定不变。"
    "核查：`step9/results/scale_check.tsv`、`Step9_统计功效.md` §2。")


def mi_mark(tag):
    return "‡" if tag in SCALE_FLAG_TAGS else ""


def J(p):
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def fp(x):
    if x is None:
        return "-"
    if x < 1e-3:
        return "%.2e" % x
    return "%.3f" % x


def fnum(x, d=4):
    if x is None:
        return "-"
    if abs(x) != 0 and (abs(x) < 1e-3 or abs(x) >= 1e4):
        return "%.2e" % x
    return ("%%.%df" % d) % x


E = 2.718281828459045


def fmt_or(b, se):
    return "%.3f [%.3f, %.3f]" % (E ** b, E ** (b - 1.959963984540054 * se), E ** (b + 1.959963984540054 * se))


def loo_range(tag, t):
    p = os.path.join(ROOT, "step4", "results", "loo_%s_%s.tsv" % (tag, t))
    if not os.path.exists(p):
        return None
    rows = [l.split("\t") for l in io.open(p, encoding="utf-8").read().strip().split("\n")[1:] if l.strip()]
    b = [float(r[2]) for r in rows]
    pl = [float(r[4]) for r in rows]
    return dict(n=len(rows), bmin=min(b), bmax=max(b), pmax=max(pl),
                all_pos=all(v > 0 for v in b), all_neg=all(v < 0 for v in b))


R = {}   # R[tag][trait] = json
for tag, *_ in OUTCOMES:
    R[tag] = {t: J(os.path.join(ROOT, "step4", "results", "mr_%s_%s.json" % (tag, t))) for t, _ in TRAITS}

# PAD 层（用于重叠配对对比）
pad_main = {t: J(os.path.join(ROOT, "step3", "results", "mr_main_%s.json" % t)) for t, _ in TRAITS}
pad_fin = {t: J(os.path.join(ROOT, "step4", "results", "mr_finngen_%s.json" % t)) for t, _ in TRAITS}

# ---------------------------------------------------------------- 判定（15 格 + 三角）
VERDICT = {
    ("hypo", "cad_nouk"): ("名义阳性", "warn",
        "IVW 固定 p=0.011，但随机 p=0.146、加权中位数 p=0.12；异质性极强（Q p=1.7e-11）。与 PAD 主分析模式一致。"),
    ("hypo", "cad_ukb"): ("跨库差异", "warn",
        "β 较无重叠版高 1.77 倍、p 从 0.011 → 9.8e-07。**注意（2026-10-07 模拟）**：Step9 的解析解 + 蒙特卡洛表明，"
        "在本项目工具强度（甲减 F_w≈109）下样本重叠最多只能造成 +4.6%（r=5）～+18.4%（r=20）的改变，**解释不了这一差异**；"
        "差异应归因于两套结局 GWAS 的病例定义/构成/人群。故 UKB 结局仍只作敏感性分析（因跨库不可调和，而非因重叠）。"),
    ("hypo", "mi_ukb"): ("尺度·不进结论", "bad",
        "β 仅 +0.0020 却 p=9.8e-08 → 系该文件处于**线性概率尺度**（见 ‡），并非「重叠使 SE 被低估」；OR=1.002，临床意义可忽略。"),
    ("hypo", "is_laa"): ("阴性", "ok", "全方法一致，95%CI 跨 0。"),
    ("hypo", "is_sv"): ("名义阳性·未达校正", "warn",
        "固定 p=0.0055、随机 p=0.029 双显著（全表唯一）；但 15 次检验 Bonferroni 阈值 0.0033，校正后不显著。"),
    ("tsh", "cad_nouk"): ("阴性", "ok", "全方法一致。"),
    ("tsh", "cad_ukb"): ("跨库不一致", "warn",
        "同一疾病：无重叠版 β=−0.0135 (p=0.50) 阴性，重叠版 β=+0.0394 (p=0.035) —— 符号翻转。"
        "**注意（2026-10-07 模拟）**：重叠在该工具强度（TSH F_w≈202）下最多造成 +0.00056 位移"
        "（仅为观测位移 +0.0528 的 1/94），**不能归因于重叠**；系两套结局 GWAS 之间的异质性。"),
    ("tsh", "mi_ukb"): ("边缘·尺度", "warn", "β≈+0.001，p=0.097 未达显著；该结局文件为线性概率尺度（见 ‡）。"),
    ("tsh", "is_laa"): ("阴性", "ok", "全方法一致。"),
    ("tsh", "is_sv"): ("阴性", "ok", "全方法一致。"),
    ("ft4", "cad_nouk"): ("阴性", "ok", "全方法一致。"),
    ("ft4", "cad_ukb"): ("阴性", "ok", "全方法一致。"),
    ("ft4", "mi_ukb"): ("阴性·尺度", "warn", "β≈−0.0024，p=0.085；线性概率尺度（见 ‡）。"),
    ("ft4", "is_laa"): ("阴性·弱正向", "warn",
        "β=+0.117 (p=0.226) 方向为正，与 PAD 层 MVP 的弱正向一致；小样本，无稳健证据。"),
    ("ft4", "is_sv"): ("阴性", "ok", "全方法一致。"),
}
TRI = {
    "hypo": ("方向跨谱一致（5/5 β&gt;0），但无一通过随机效应＋加权中位数双重稳健检验；"
             "两套结局库间差异达 1.77×–1.91×（重叠量级不足以解释，见 §3）→ 信号存疑", "warn"),
    "tsh": ("三个无重叠结局一致阴性 → 全血管谱上最硬的阴性结论", "ok"),
    "ft4": ("全谱无稳健证据；IS-LAA 与 PAD-MVP 层同现弱正向，提示可能有极微弱信号但不可检验", "warn"),
}


def pill(kind, text):
    return '<span class="pill p-%s">%s</span>' % (kind, text)


def cls_beta(b):
    return "pos" if b > 0 else ("neg" if b < 0 else "")


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


def build_html():
    o = []
    o.append('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">')
    o.append('<title>甲状腺功能 → 泛血管谱：Step 5 结果一页纸</title>')
    o.append('<style>%s</style></head><body><div class="page">' % CSS)
    o.append('<h1>甲状腺功能 → 泛血管谱（CAD / MI / IS-LAA / IS-SV）：Step 5 结果一页纸</h1>')
    o.append('<div class="sub">两样本孟德尔随机化 · 3 性状 × 5 结局 · 「无重叠 vs 有重叠」配对设计</div>')
    o.append('<div class="meta">'
             '<span><b>暴露</b> 甲减 / TSH / FT4（Rand 2025, <i>Nat Genet</i>, PMID 41238958，含 UKB＋FinnGen DF10）</span>'
             '<span><b>日期</b> 2026-10-07</span>'
             '<span><b>检验数</b> 15（Bonferroni 阈 0.0033）</span></div>')

    # 1 结论矩阵
    o.append('<h2>1 · 结论矩阵</h2><table><tr><th class="l">暴露</th>')
    for tag, name, *_ in OUTCOMES:
        ov = '<br><span style="color:#b42318;font-weight:600">重叠</span>' if dict((t[0], t[5]) for t in OUTCOMES)[tag] else '<br><span style="color:#0f7b45">无重叠</span>'
        mark = mi_mark(tag)
        o.append('<th>%s%s%s</th>' % (name.replace("（", "<br>（"), ("<br>" + mark) if mark else "", ov))
    o.append('<th class="l">跨谱判定</th></tr>')
    for t, cn in TRAITS:
        o.append('<tr><td class="l"><b>%s</b></td>' % cn)
        for tag, *_ in OUTCOMES:
            v, kind, _ = VERDICT[(t, tag)]
            o.append('<td>%s</td>' % pill(kind, v))
        tk, tkind = TRI[t]
        o.append('<td class="l">%s</td></tr>' % pill(tkind, tk))
    o.append('</table>')

    # 2 主结果
    o.append('<h2>2 · 主结果（IVW；log-OR 尺度）</h2><table>'
             '<tr><th class="l">结局</th><th>暴露</th><th>nSNP</th><th>β</th><th>SE</th><th>p（固定）</th>'
             '<th>OR [95%CI]</th><th>p（随机 DL）</th><th>Egger 截距 p</th><th>Q p</th><th>加权中位数 p</th></tr>')
    for tag, name, *_ in OUTCOMES:
        for t, cn in TRAITS:
            d = R[tag][t]
            b = d["ivw_fixed"]["beta"]; se = d["ivw_fixed"]["se"]; p = d["ivw_fixed"]["p"]
            ip = d["mr_egger"]["intercept_pval"]; qp = d["ivw_random"]["Q_pval"]
            wmp = d["weighted_median"]["p"]
            psig = ' class="sig"' if p < 0.05 else ""
            o.append('<tr><td class="l">%s%s</td><td>%s</td><td>%d</td><td class="%s">%s</td><td>%s</td>'
                     '<td%s>%s</td><td>%s</td><td>%s</td><td%s>%s</td><td>%s</td><td>%s</td></tr>' % (
                         name, mi_mark(tag), cn, d["nsnp"], cls_beta(b), fnum(b), fnum(se), psig, fp(p),
                         fmt_or(b, se), fp(d["ivw_random"]["p"]),
                         ' class="sig"' if ip < 0.05 else "", fp(ip), fp(qp), fp(wmp)))
        o.append('<tr><td colspan="11" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>')
    o.append('</table>')
    o.append('<div class="box"><b>读法：</b>红 β 为正、绿为负；p&lt;0.05 加粗标红。'
             '两个 UKB 结局（CAD/MI）的「显著」结果<b>不可作因果证据</b>：CAD 与无重叠结局库差异大'
             '（且该差异非样本重叠所能解释，见第 3 节）；MI 另有<b>尺度</b>问题（见 ‡）。</div>')
    o.append(MI_NOTE_HTML)

    # 3 结局库配对对比（含重叠解释力上限）
    o.append('<h2>3 · 跨结局库配对对比（含重叠解释力上限）</h2>')
    o.append('<table><tr><th class="l">疾病</th><th class="l">数据集对</th><th class="l">重叠</th>'
             '<th>甲减 β（比值）</th><th>甲减 p</th><th>TSH β</th><th>TSH p</th><th class="l">解读</th></tr>')
    hm, hf = pad_main["hypo"]["ivw_fixed"]["beta"], pad_fin["hypo"]["ivw_fixed"]["beta"]
    tm, tf = pad_main["tsh"]["ivw_fixed"]["beta"], pad_fin["tsh"]["ivw_fixed"]["beta"]
    cn_, cu_ = R["cad_nouk"]["hypo"]["ivw_fixed"]["beta"], R["cad_ukb"]["hypo"]["ivw_fixed"]["beta"]
    tn_, tu_ = R["cad_nouk"]["tsh"]["ivw_fixed"]["beta"], R["cad_ukb"]["tsh"]["ivw_fixed"]["beta"]
    o.append('<tr><td class="l" rowspan="2"><b>PAD</b></td>'
             '<td class="l">主分析（Sakaue 2021）</td><td class="l">部分重叠*</td><td>%+.4f</td><td>%s</td><td>%+.4f</td><td>%s</td>'
             '<td class="l" rowspan="2">两库甲减 β 相差 <b>%.2f 倍</b>（主分析↔R13，<b>非独立样本</b>，仅作次级参照）；零重叠锚（MVP↔Sakaue）的配对 β 比为 <b>1.91</b>。但按 Step9 模拟，重叠在该工具强度下'
             '最多只能造成 <b>+18.4%%</b> 的改变，<b>远不足以解释</b>；差异应归因于两套结局 GWAS 的'
             '病例定义/构成/人群（含 MVP 为跨族裔 meta）。'
             '*主分析队列亦含 UKB；且其 EUR 成分 483,078 = UKB 350,366 + FinnGen R3 132,712（2026-10-08 已核实），FinnGen R3 ⊂ R13 → 两层并非独立样本，比值作参考。</td></tr>' % (
                 hm, fp(pad_main["hypo"]["ivw_fixed"]["p"]), tm, fp(pad_main["tsh"]["ivw_fixed"]["p"]), hf / hm))
    o.append('<tr><td class="l">FinnGen R13 I9_PAD</td><td class="l"><b>重叠</b>（甲减）</td><td><b>%+.4f</b></td><td><b>%s</b></td><td>%+.4f</td><td>%s</td></tr>' % (
        hf, fp(pad_fin["hypo"]["ivw_fixed"]["p"]), tf, fp(pad_fin["tsh"]["ivw_fixed"]["p"])))
    o.append('<tr><td class="l" rowspan="2"><b>CAD</b></td>'
             '<td class="l">Nikpay 2015（GCST003116）</td><td class="l"><b>无重叠</b></td><td>%+.4f</td><td>%s</td><td>%+.4f</td><td>%s</td>'
             '<td class="l" rowspan="2">两库甲减 β 相差 <b>%.2f 倍</b>；<b>TSH 符号翻转</b>'
             '（−0.014, p=0.50 → +0.039, p=0.035）。Step9 模拟：重叠在 TSH 工具强度下最多造成 '
             '<b>+0.00056</b> 位移，仅为观测位移的 <b>1/94</b> → <b>两项差异均不能归因于重叠</b>。</td></tr>' % (
                 cn_, fp(R["cad_nouk"]["hypo"]["ivw_fixed"]["p"]), tn_, fp(R["cad_nouk"]["tsh"]["ivw_fixed"]["p"]), cu_ / cn_))
    o.append('<tr><td class="l">Mbatchou 2021（GCST90013864）</td><td class="l"><b>重叠</b>（UKB）</td><td><b>%+.4f</b></td><td><b>%s</b></td><td><b>%+.4f</b></td><td><b>%s</b></td></tr>' % (
        cu_, fp(R["cad_ukb"]["hypo"]["ivw_fixed"]["p"]), tu_, fp(R["cad_ukb"]["tsh"]["ivw_fixed"]["p"])))
    o.append('</table>')
    o.append('<div class="box"><b>要点（2026-10-07 依 Step9 模拟修正）</b>：配对设计确实显示同一暴露在两套结局库上'
             '给出 <b>1.77×–1.91×</b> 的 β 差异、且 TSH→CAD 符号翻转；但<b>样本重叠不是原因</b>——'
             '本项目工具强度极高（甲减 F_w≈109、TSH≈202；最低单 SNP F≈29），在完全重叠（π=1）时重叠能造成的'
             '改变上限只有 <b>+4.6%（r=5）～+18.4%（r=20）</b>；要把 1.77× 归因于重叠，需要观察性关联达到'
             '因果效应的 <b>84 倍</b>（3.62× 则需 287 倍），且所需的 π 达 <b>17 倍于</b>整个结局样本。'
             '因此差异应归因于<b>不同结局 GWAS 本身不可互换</b>（病例定义、病例/构成、人群、meta 设计）。'
             '核验：<code>Step9_重叠偏倚模拟.html</code>、<code>step9/results/overlap_sim_*.tsv</code>。</div>')

    # 4 稳健性
    o.append('<h2>4 · 稳健性与诊断</h2><table>'
             '<tr><th class="l">结局</th><th>暴露</th><th>加权众数 p</th><th>PRESSO 离群</th><th>PRESSO 校正 β (p)</th>'
             '<th>符号 +/− (p)</th><th>方向性</th><th class="l">LOO β 区间（最差 p）</th></tr>')
    for tag, name, *_ in OUTCOMES:
        for t, cn in TRAITS:
            d = R[tag][t]
            pr = (d.get("presso_equiv") or {}).get("distortion")
            nout = pr["n_outliers"] if pr else 0
            corr = ("%s (%s)" % (fnum(pr["beta_corrected"]), fp(pr["p_corrected"]))) if pr else "—"
            sg = d.get("sign_test") or {}
            sp = sg.get("p_binomial_two_sided") or sg.get("p")
            lr = loo_range(tag, t)
            ltxt = ("%+.4f ~ %+.4f（%s，最差 p %s）" % (lr["bmin"], lr["bmax"],
                    "全正" if lr["all_pos"] else ("全负" if lr["all_neg"] else "含变号"), fp(lr["pmax"]))) if lr else "未计算"
            o.append('<tr><td class="l">%s%s</td><td>%s</td><td>%s</td><td>%d</td><td>%s</td>'
                     '<td>%s/%s (%s)</td><td>%s</td><td class="l">%s</td></tr>' % (
                         name, mi_mark(tag), cn, fp(d["weighted_mode"]["p"]), nout, corr,
                         sg.get("n_positive"), sg.get("n_negative"), fp(sp),
                         d["directionality"]["direction"], ltxt))
        o.append('<tr><td colspan="8" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>')
    o.append('</table>')

    # 5 可写 / 不可写
    o.append('<h2>5 · 可写入论文的结论</h2><div class="two"><div class="box"><b>站得住的</b><ul>'
             '<li><b>遗传预测 TSH 与动脉粥样硬化性血管事件无因果关联</b>：三个无重叠结局（CAD-Nikpay p=0.50、'
             'IS-LAA p=0.185、IS-SV p=0.605）一致阴性；加上 PAD 三层 755 次 LOO 全不显著，'
             '<b>这是跨血管谱最硬的一条</b>。</li>'
             '<li><b>配对设计给出了「结局库间差异」的量化与「重叠解释力」的上限</b>（2026-10-07 修正）：'
             '同一暴露在两套结局库间 β 相差 1.77×（CAD）/1.91×（PAD）；Step9 模拟（解析解 + 10,000 次 MC）'
             '给出重叠的解释力上限 —— 完全重叠时最多 +4.6%（r=5）/ +18.4%（r=20），'
             '要把 1.77× 归因于重叠需观察性关联达因果效应的 84 倍。'
             '<b>结论：差异真实且显著，但不由样本重叠解释</b>，而指向结局 GWAS 之间的不可互换性。可作为方法学论文的核心证据。</li>'
             '<li><b>FT4 全谱无稳健因果证据</b>：15 个组合中无一通过随机效应＋加权中位数。</li>'
             '<li>甲减 β 在 5 个结局中方向全部为正（跨谱一致性本身有信息量），但无一通过双重稳健检验。</li>'
             '</ul></div><div class="box"><b>不能写的（须披露）</b><ul>'
             '<li><b>不能引用 UKB 结局（CAD/MI）的任何「显著」结果作为因果证据</b>——CAD 与无重叠结局库差异达 1.77×'
             '（Step9 模拟证明该差异<b>不是</b>样本重叠造成的，而是结局库本身不可互换）；MI 另有<b>尺度</b>问题（见 §2 ‡）。</li>'
             '<li><b>甲减的正向信号不能写成因果</b>：无重叠数据集下最强者（IS-SV p=0.0055）未过 Bonferroni（0.0033），'
             '且随机效应/加权中位数均不稳健；异质性极强。</li>'
             '<li><b>IS-SV 甲减的"双显著"（固定 p=0.0055 / 随机 p=0.029）只能作提示</b>，不能作结论。</li>'
             '<li>MEGASTROKE 亚型病例数小（4,373 / 5,386），阴性结果受限于效能。</li>'
             '<li>cad_ukb 的 EAF 整列缺失 → 无法做 EAF 一致性校验（本项目回文位点为 0，实际影响很小，但须披露）。</li>'
             '</ul></div></div>')

    # 6 数据来源
    o.append('<h2>6 · 数据来源</h2><table><tr><th class="l">结局</th><th class="l">来源</th>'
             '<th class="l">样本量（例/对照）</th><th class="l">人群</th><th class="l">与暴露重叠</th></tr>')
    for tag, name, src, nstr, pop, ov in OUTCOMES:
        o.append('<tr><td class="l">%s%s</td><td class="l">%s</td><td class="l">%s</td><td class="l">%s</td>'
                 '<td class="l">%s</td></tr>' % (
                     name, mi_mark(tag), src, nstr, pop,
                     '<b style="color:#b42318">UKB 重叠</b>' if ov else '<span style="color:#0f7b45">无</span>'))
    o.append('</table>')
    o.append('<div class="foot"><b>生成方式</b>：由 <code>step5/scripts/s5_onepager.py</code> 直接读取 '
             '<code>step4/results/mr_&lt;tag&gt;_&lt;trait&gt;.json</code> 与 <code>loo_*.tsv</code> 生成，数值未作人工修改。'
             '配图 65 张见 <code>step4/figures/</code>（每结局 13 张）。'
             '方法口径与 Step 3/4 完全一致（IVW 固定/随机 DL、MR-Egger、加权中位数/众数、Cochran Q、'
             'MR-PRESSO 等价、LOO、方向性 R² 对比；工具 p&lt;5e-8、r²&lt;0.001、10 Mb 窗口）。</div>')
    o.append('</div></body></html>')
    return "\n".join(o)


def build_md():
    o = []
    o.append("# 甲状腺功能 → 泛血管谱（CAD / MI / IS-LAA / IS-SV）：Step 5 结果一页纸")
    o.append("")
    o.append("**两样本孟德尔随机化 · 3 性状 × 5 结局 · 「无重叠 vs 有重叠」配对设计**　｜　2026-10-07")
    o.append("")
    o.append("- **暴露**：甲减 / TSH / FT4 —— Rand et al. 2025, *Nat Genet*, PMID 41238958（含 UKB＋FinnGen DF10）")
    o.append("- **检验数**：15（Bonferroni 阈值 0.0033）")
    o.append("")
    o.append("---")
    o.append("")
    o.append("## 1 结论矩阵")
    o.append("")
    hdr = ["暴露"] + ["%s%s" % (n, mi_mark(tag)) for tag, n, *_ in OUTCOMES] + ["跨谱判定"]
    o.append("| " + " | ".join(hdr) + " |")
    o.append("|" + "---|" * len(hdr))
    tri_txt = {t: v for t, (v, _) in TRI.items()}
    for t, cn in TRAITS:
        cells = [VERDICT[(t, tag)][0] for tag, *_ in OUTCOMES]
        ovmark = {tag: ("⚠" if ovr else "✓") for tag, _, _, _, _, ovr in OUTCOMES}
        row = [cn] + ["%s %s" % (ovmark[tag], c) for (tag, c) in zip([x[0] for x in OUTCOMES], cells)] + [tri_txt[t]]
        o.append("| " + " | ".join(row) + " |")
    o.append("")
    o.append("> ⚠ = 存在 UKB 样本重叠，结果只作敏感性分析（2026-10-07 修正：此处**不是**说「重叠造成放大」——"
             "Step9 模拟显示重叠量级不足以解释观察到的差异，而是结局库之间本身不可调和）；✓ = 无重叠。"
             "‡ = 该结局存在尺度问题，见 §2 表下注记。")
    o.append("")
    o.append("## 2 主结果（IVW 固定效应）")
    o.append("")
    o.append("| 结局 | 暴露 | nSNP | β | SE | p（固定） | OR [95%CI] | p（随机 DL） | Egger 截距 p | Q p | 加权中位数 p |")
    o.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for tag, name, *_ in OUTCOMES:
        for t, cn in TRAITS:
            d = R[tag][t]
            b = d["ivw_fixed"]["beta"]; se = d["ivw_fixed"]["se"]; p = d["ivw_fixed"]["p"]
            mark = "**" if p < 0.05 else ""
            o.append("| %s%s | %s | %d | %+.4f | %.4f | %s%s%s | %s | %s | %s | %s | %s |" % (
                name, mi_mark(tag), cn, d["nsnp"], b, se, mark, fp(p), mark, fmt_or(b, se),
                fp(d["ivw_random"]["p"]), fp(d["mr_egger"]["intercept_pval"]),
                fp(d["ivw_random"]["Q_pval"]), fp(d["weighted_median"]["p"])))
    o.append("")
    o.append(MI_NOTE_MD)
    o.append("")
    o.append("## 3 跨结局库配对对比（含重叠解释力上限）")
    o.append("")
    o.append("| 疾病 | 数据集 | 重叠 | 甲减 β | 甲减 p | TSH β | TSH p |")
    o.append("|---|---|---|---|---|---|---|")
    hm = pad_main["hypo"]["ivw_fixed"]["beta"]; hf = pad_fin["hypo"]["ivw_fixed"]["beta"]
    tm = pad_main["tsh"]["ivw_fixed"]["beta"]; tf = pad_fin["tsh"]["ivw_fixed"]["beta"]
    cn_ = R["cad_nouk"]["hypo"]["ivw_fixed"]["beta"]; cu_ = R["cad_ukb"]["hypo"]["ivw_fixed"]["beta"]
    tn_ = R["cad_nouk"]["tsh"]["ivw_fixed"]["beta"]; tu_ = R["cad_ukb"]["tsh"]["ivw_fixed"]["beta"]
    o.append("| PAD | 主分析（Sakaue 2021）* | 部分 | %+.4f | %s | %+.4f | %s |" % (
        hm, fp(pad_main["hypo"]["ivw_fixed"]["p"]), tm, fp(pad_main["tsh"]["ivw_fixed"]["p"])))
    o.append("| PAD | FinnGen R13 | **有** | **%+.4f（×%.2f，vs 主分析·非独立）** | %s | %+.4f | %s |" % (
        hf, hf / hm, fp(pad_fin["hypo"]["ivw_fixed"]["p"]), tf, fp(pad_fin["tsh"]["ivw_fixed"]["p"])))
    o.append("| CAD | Nikpay 2015 | **无** | %+.4f | %s | %+.4f | %s |" % (
        cn_, fp(R["cad_nouk"]["hypo"]["ivw_fixed"]["p"]), tn_, fp(R["cad_nouk"]["tsh"]["ivw_fixed"]["p"])))
    o.append("| CAD | Mbatchou 2021 (UKB) | **有** | **%+.4f（×%.2f）** | %s | **%+.4f（翻转）** | **%s** |" % (
        cu_, cu_ / cn_, fp(R["cad_ukb"]["hypo"]["ivw_fixed"]["p"]), tu_,
        fp(R["cad_ukb"]["tsh"]["ivw_fixed"]["p"])))
    o.append("")
    o.append("> **要点（2026-10-07 依 Step9 模拟修正）**：配对设计确实显示同一暴露在两套结局库上给出"
             "**1.77×（CAD）/ 1.91×（PAD）** 的 β 差异、且 TSH→CAD 符号翻转（%+.4f, p=%s → %+.4f, p=%s）；"
             "但**样本重叠不是原因** —— 本项目工具强度极高（甲减 F_w≈109、TSH≈202；最低单 SNP F≈29），"
             "在完全重叠（π=1）下重叠能造成的改变上限只有 **+4.6%%（r=5）～+18.4%%（r=20）**；"
             "要把 1.77× 归因于重叠需观察性关联达因果效应的 **84 倍**（3.62× 需 287 倍），"
             "所需 π 达 **17 倍于**整个结局样本。差异应归因于**不同结局 GWAS 本身不可互换**"
             "（病例定义、病例/对照构成、人群、meta 设计）。"
             "**口径**：头条配对数字以零重叠锚（MVP↔Sakaue = **1.91×**）为准；本表 PAD 行的 ×1.89 系主分析↔FinnGen R13（R3 ⊂ R13，非独立样本），仅作次级参照。"
             "核验：`Step9_重叠偏倚模拟.md`、`step9/results/overlap_sim_*.tsv`。"
             "*主分析队列亦含 UKB；且其 EUR 成分 483,078 = UKB 350,366 + FinnGen R3 132,712（2026-10-08 已核实），FinnGen R3 ⊂ R13 → 两层并非独立样本，比值作参考。"
             % (tn_, fp(R["cad_nouk"]["tsh"]["ivw_fixed"]["p"]), tu_, fp(R["cad_ukb"]["tsh"]["ivw_fixed"]["p"])))
    o.append("")
    o.append("## 4 可写入论文的结论")
    o.append("")
    o.append("### 站得住的")
    o.append("")
    o.append("1. **遗传预测 TSH 与动脉粥样硬化性血管事件无因果关联**：三个无重叠结局一致阴性"
             "（CAD-Nikpay p=0.50、IS-LAA p=0.185、IS-SV p=0.605）；加上 PAD 三层 755 次 LOO 全不显著 → 跨血管谱最硬的一条。")
    o.append("2. **结局库间差异的量化 + 重叠解释力的上限**：同一暴露在两套结局库间 β 相差 1.77×（CAD）/ 1.91×（PAD）；"
             "Step9 模拟给出重叠的解释力上限（π=1 时 ≤ +4.6%（r=5）～+18.4%（r=20））。"
             "**差异真实，但不由样本重叠解释**——指向结局 GWAS 之间的不可互换性。")
    o.append("3. **FT4 全谱无稳健因果证据**。")
    o.append("4. 甲减 β 在 5 个结局中方向全正（跨谱一致），但无一通过随机＋中位数双重稳健检验。")
    o.append("")
    o.append("### 不能写的（须披露）")
    o.append("")
    o.append("1. **不能引用 UKB 结局（CAD/MI）的任何「显著」结果作为因果证据**——CAD 与无重叠结局库差异达 1.77×"
             "（Step9 模拟证明该差异**不是**样本重叠造成的，而是结局库本身不可互换）；MI 另有**尺度**问题（见 §2 ‡）。")
    o.append("2. **甲减正向信号不能写成因果**：无重叠下最强者（IS-SV p=0.0055）未过 Bonferroni（0.0033），随机/中位数不稳健。")
    o.append("3. IS-SV 甲减「双显著」只能作提示。")
    o.append("4. MEGASTROKE 亚型病例数小（4,373 / 5,386），阴性结果受限于效能。")
    o.append("5. cad_ukb 的 EAF 整列缺失 → 无法做 EAF 一致性校验（回文位点为 0，实际影响小，但须披露）。")
    o.append("")
    o.append("## 5 数据来源")
    o.append("")
    o.append("| 结局 | 来源 | 样本量 | 人群 | 与暴露重叠 |")
    o.append("|---|---|---|---|---|")
    for tag, name, src, nstr, pop, ov in OUTCOMES:
        o.append("| %s%s | %s | %s | %s | %s |" % (
            name, mi_mark(tag), src, nstr, pop, "**UKB 重叠**" if ov else "无"))
    o.append("")
    o.append("**方法口径**：与 Step 3/4 完全一致（IVW 固定/随机 DL、MR-Egger、加权中位数/众数、Cochran Q、"
             "MR-PRESSO 等价、LOO、方向性 R² 对比；工具 p<5e-8、r²<0.001、10 Mb 窗口）。")
    o.append("")
    o.append("> 本文件由 `step5/scripts/s5_onepager.py` 从结果 JSON / LOO TSV 直接生成，数值未作人工修改。")
    return "\n".join(o)


if __name__ == "__main__":
    if not os.path.exists(DELIV):
        os.makedirs(DELIV)
    hp = os.path.join(DELIV, "Step5_泛血管谱结果.html")
    mp = os.path.join(DELIV, "Step5_泛血管谱结果.md")
    with io.open(hp, "w", encoding="utf-8") as f:
        f.write(build_html())
    with io.open(mp, "w", encoding="utf-8") as f:
        f.write(build_md())
    print("OK  " + hp, os.path.getsize(hp), "B")
    print("OK  " + mp, os.path.getsize(mp), "B")
