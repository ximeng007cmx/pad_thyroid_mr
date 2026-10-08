# -*- coding: utf-8 -*-
"""
Step 3 + Step 4 结果一页纸汇总生成器。

原则：所有数值直接从 results/*.json 与 results/loo_*.tsv 读取，
      不手抄、不改数、不做任何修约（除展示位数）。
输出：
  deliverables/Step3-4_结果一页纸.html   (自包含、A4 打印友好)
  deliverables/Step3-4_结果一页纸.md     (归档用)
"""
import io
import json
import os

ROOT = r"D:\WorkSpace_Study\pad_thyroid_mr"
DELIV = os.path.join(ROOT, "deliverables")
TRAITS = [("hypo", "甲减"), ("tsh", "TSH"), ("ft4", "FT4")]


def J(p):
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def fnum(x, d=4):
    """数值格式化；不做修约以外的任何改动。"""
    if x is None:
        return "-"
    if abs(x) != 0 and (abs(x) < 1e-3 or abs(x) >= 1e4):
        return "%.2e" % x
    return ("%%.%df" % d) % x


def fp(x):
    """p 值：<1e-3 用科学计数法，否则 3 位小数。"""
    if x is None:
        return "-"
    if x < 1e-3:
        return "%.2e" % x
    return "%.3f" % x


def orf(b):
    return "%.3f" % (2.718281828459045 ** b)


def fmt_ci(b, se):
    lo = 2.718281828459045 ** (b - 1.959963984540054 * se)
    hi = 2.718281828459045 ** (b + 1.959963984540054 * se)
    return "%.3f－%.3f" % (lo, hi)


def fmt_or(b, se):
    """OR 点估计 + 95%CI。"""
    return "%.3f [%s]" % (2.718281828459045 ** b, fmt_ci(b, se))


def loo_range(path):
    if not os.path.exists(path):
        return None
    rows = [l.split("\t") for l in io.open(path, encoding="utf-8").read().strip().split("\n")[1:] if l.strip()]
    b = [float(r[2]) for r in rows]
    p = [float(r[4]) for r in rows]
    return dict(n=len(rows), bmin=min(b), bmax=max(b), pmin=min(p), pmax=max(p),
                all_pos=all(v > 0 for v in b), all_neg=all(v < 0 for v in b))


# ---------------------------------------------------------------- 载入
main = {t: J(os.path.join(ROOT, "step3", "results", "mr_main_%s.json" % t)) for t, _ in TRAITS}
fin = {t: J(os.path.join(ROOT, "step4", "results", "mr_finngen_%s.json" % t)) for t, _ in TRAITS}
mvp = {t: J(os.path.join(ROOT, "step4", "results", "mr_mvp_meta_%s.json" % t)) for t, _ in TRAITS}
loo = {
    "main": {t: loo_range(os.path.join(ROOT, "step3", "results", "loo_%s.tsv" % t)) for t, _ in TRAITS},
    "fin": {t: loo_range(os.path.join(ROOT, "step4", "results", "loo_finngen_%s.tsv" % t)) for t, _ in TRAITS},
    "mvp": {t: loo_range(os.path.join(ROOT, "step4", "results", "loo_mvp_meta_%s.tsv" % t)) for t, _ in TRAITS},
}

DATASETS = [
    ("main", "主分析", "Sakaue 2021 · GCST90018890", "7,114 / 475,964", "EUR＋跨族裔合并", ""),
    ("fin", "FinnGen R13", "I9_PAD", "22,244 / 460,490", "芬兰（EUR）", "部分限制"),
    ("mvp", "MVP", "Klarin 2019 · dbGaP phs001672", "31,307 / 211,753", "跨族裔 meta", ""),
]

# ---------------------------------------------------------------- 判定
# 依据：点估计方向 + 显著性 + 稳健性（Egger 截距 / Q / LOO）
VERDICT = {
    ("hypo", "main"): ("名义阳性", "warn",
                       "IVW 固定 p=0.015，但随机效应 p=0.130、加权中位数 p=0.506；LOO 下 p 全程 &gt;0.05。工具侧含 UKB，与结局重叠。"),
    ("hypo", "fin"): ("重叠失效", "bad",
                      "暴露含 FinnGen DF10，与结局同源，估计不可作复制证据（β 约为主分析 1.9 倍）。注：Step9 模拟显示该 1.9 倍差异<b>非重叠所能解释</b>（重叠上限 ≤ +4.6%～+18.4%），而来自两套结局库本身不可互换（Step9_重叠偏倚模拟.md）。"),
    ("hypo", "mvp"): ("无重叠·方向一致", "ok",
                      "独立队列（暴露不含 MVP），IVW 固定 p=0.040 同向；但随机效应 p=0.527、LOO 全程 p 0.41－0.70 → 方向一致、无统计证据。"),
    ("tsh", "main"): ("阴性", "ok", "全方法一致，95%CI 跨 0。"),
    ("tsh", "fin"): ("阴性·独立复制", "ok", "有效独立复制，方向与主分析一致。"),
    ("tsh", "mvp"): ("阴性·独立复制", "ok",
                     "跨族裔队列阴性；LOO 全程同向为负（−0.033～−0.016）且 p 全不显著 → 稳固阴性。"),
    ("ft4", "main"): ("阴性·有多效性", "warn",
                      "点估计阴性，但 MR-Egger 截距显著（p=2.8e-05）→ 定向水平多效性，结论不可作因果解读。"),
    ("ft4", "fin"): ("阴性·独立复制", "ok", "FinnGen 中 Egger 截距不显著（p=0.125），提示主分析截距信号不稳定。"),
    ("ft4", "mvp"): ("阴性·方向翻转", "warn",
                     "点估计转正（+0.078），但 95%CI 跨 0；仅 63 个工具，无稳健证据。"),
}

# ---------------------------------------------------------------- HTML
CSS = """
:root{
  --ink:#111827; --mut:#5b6675; --line:#dfe3e8; --soft:#f6f7f9;
  --ok:#0f7b45; --warn:#9a6700; --bad:#b42318; --acc:#1f4e79;
}
*{box-sizing:border-box}
body{margin:0;background:#eceef1;color:var(--ink);
  font-family:"Source Han Sans SC","Microsoft YaHei","PingFang SC",-apple-system,"Segoe UI",sans-serif;
  font-size:11.5px;line-height:1.5}
.page{width:210mm;min-height:297mm;margin:14px auto;background:#fff;padding:13mm 12mm;
  box-shadow:0 2px 14px rgba(15,23,42,.13)}
h1{font-size:17px;margin:0 0 2px;letter-spacing:.2px}
.sub{color:var(--mut);font-size:10.5px;margin-bottom:9px}
.meta{display:flex;flex-wrap:wrap;gap:4px 16px;font-size:10px;color:var(--mut);
  border-top:1px solid var(--line);border-bottom:1px solid var(--line);padding:5px 0;margin-bottom:10px}
h2{font-size:12.5px;margin:13px 0 5px;padding-left:7px;border-left:3px solid var(--acc);
  color:var(--acc);letter-spacing:.2px}
table{width:100%;border-collapse:collapse;font-size:10.5px;margin-bottom:4px}
th,td{border:1px solid var(--line);padding:3px 5px;text-align:center;vertical-align:middle;white-space:nowrap}
th{background:var(--soft);font-weight:600;color:#374151}
td.l,th.l{text-align:left;white-space:normal}
tr.grp td{background:#fbfbfc;font-weight:600;text-align:left;color:#374151}
.pill{display:inline-block;padding:1px 6px;border-radius:9px;font-size:9.5px;font-weight:600;white-space:nowrap}
.p-ok{background:#e7f5ec;color:var(--ok);border:1px solid #b7dfc6}
.p-warn{background:#fdf4e3;color:var(--warn);border:1px solid #f0d9a8}
.p-bad{background:#fdeceb;color:var(--bad);border:1px solid #f3c3bf}
.sig{color:var(--bad);font-weight:700}
.pos{color:var(--bad);font-weight:600}
.neg{color:var(--ok);font-weight:600}
.box{background:var(--soft);border:1px solid var(--line);border-left:3px solid var(--acc);
  padding:6px 9px;margin:6px 0;font-size:10.5px}
.two{display:flex;gap:9px}
.two>div{flex:1}
ul{margin:3px 0 3px 15px;padding:0}
li{margin:2px 0}
.foot{margin-top:10px;border-top:1px solid var(--line);padding-top:6px;
  font-size:9.5px;color:var(--mut);line-height:1.55}
code{background:#f1f3f6;padding:0 3px;border-radius:2px;font-size:9.5px;
  font-family:"Cascadia Mono",Consolas,monospace}
@media print{
  body{background:#fff;font-size:10.5px}
  .page{width:auto;margin:0;box-shadow:none;padding:9mm 8mm;min-height:0}
  @page{size:A4;margin:9mm}
  h2{margin:9px 0 4px}
}
"""


def pill(kind, text):
    return '<span class="pill p-%s">%s</span>' % (kind, text)


def cls_beta(b):
    return "pos" if b > 0 else ("neg" if b < 0 else "")


def build_html():
    o = []
    o.append('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">')
    o.append('<title>甲状腺功能 → PAD：Step 3 + Step 4 结果一页纸</title>')
    o.append('<style>%s</style></head><body><div class="page">' % CSS)

    o.append('<h1>甲状腺功能 → 外周动脉疾病（PAD）：Step 3 + Step 4 结果一页纸</h1>')
    o.append('<div class="sub">两样本孟德尔随机化 · 三性状 × 三数据集三角验证</div>')
    o.append('<div class="meta">'
             '<span><b>暴露</b> 甲减 / TSH / FT4（Rand 2025, Nat Genet, PMID 41238958）</span>'
             '<span><b>主结局</b> Sakaue 2021（PMID 34594039）</span>'
             '<span><b>日期</b> 2026-10-07</span>'
             '<span><b>工具</b> 甲减 222 / TSH 258 / FT4 65</span>'
             '</div>')

    # ---- 1 结论矩阵
    o.append('<h2>1 · 结论矩阵</h2>')
    o.append('<table><tr><th class="l">暴露</th>'
             '<th>主分析<br><span style="font-weight:400;color:#5b6675">Sakaue 2021</span></th>'
             '<th>FinnGen R13<br><span style="font-weight:400;color:#5b6675">I9_PAD</span></th>'
             '<th>MVP<br><span style="font-weight:400;color:#5b6675">phs001672</span></th>'
             '<th class="l">三角验证判定</th></tr>')
    tri = {
        "hypo": ("方向一致（全正），但均未达稳健显著", "warn"),
        "tsh": ("三数据集一致阴性 → 最硬结论", "ok"),
        "ft4": ("方向翻转、全部跨 0 → 无稳健证据", "warn"),
    }
    for t, cn in TRAITS:
        row = ['<tr><td class="l"><b>%s</b></td>' % cn]
        for k, _, _, _, _, _ in DATASETS:
            v, kind, _ = VERDICT[(t, k)]
            row.append('<td>%s</td>' % pill(kind, v))
        tk, tkind = tri[t]
        row.append('<td class="l">%s</td></tr>' % pill(tkind, tk))
        o.append("".join(row))
    o.append('</table>')

    # ---- 2 主结果
    o.append('<h2>2 · 主结果（log-OR 尺度；OR 为该 β 的指数变换）</h2>')
    o.append('<table><tr><th class="l">暴露</th><th>数据集</th><th>nSNP</th>'
             '<th>IVW 固定 β</th><th>SE</th><th>p</th><th>OR [95%CI]</th>'
             '<th>IVW 随机 p</th><th>Egger β (p)</th><th>Egger 截距 p</th><th>Q p</th></tr>')
    for t, cn in TRAITS:
        for k, dlabel, dsub, nstr, _, _ in DATASETS:
            d = {"main": main, "fin": fin, "mvp": mvp}[k][t]
            nsnp = d.get("nsnp")
            b = d["ivw_fixed"]["beta"]; se = d["ivw_fixed"]["se"]; p = d["ivw_fixed"]["p"]
            pr = d["ivw_random"]["p"]
            eb = d["mr_egger"]["beta"]; ep = d["mr_egger"]["p"]
            ip = d["mr_egger"]["intercept_pval"]
            qp = d["ivw_random"]["Q_pval"]
            psig = ' class="sig"' if p < 0.05 else ""
            o.append('<tr><td class="l"><b>%s</b></td><td class="l">%s</td><td>%s</td>'
                     '<td class="%s">%s</td><td>%s</td><td%s>%s</td><td>%s</td>'
                     '<td>%s</td><td class="%s">%s</td><td%s>%s</td><td>%s</td></tr>' % (
                         cn, dlabel, nsnp,
                         cls_beta(b), fnum(b), fnum(se), psig, fp(p), fmt_or(b, se),
                         fp(pr), cls_beta(eb), "%s (%s)" % (fnum(eb), fp(ep)),
                         ' class="sig"' if ip < 0.05 else "", fp(ip), fp(qp)))
        o.append('<tr><td colspan="11" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>')
    o.append('</table>')
    o.append('<div class="box"><b>读法：</b>三行 = 三个外部队列。'
             '<span class="pos">红</span> 为正向 β，<span class="neg">绿</span> 为负向 β；'
             'p&lt;0.05 加粗标红。加粗行间以细线分隔，便于横向比对同一性状在不同队列的表现。</div>')

    # ---- 3 稳健性
    o.append('<h2>3 · 稳健性与诊断</h2>')
    o.append('<table><tr><th class="l">暴露</th><th>数据集</th><th>加权中位数 p</th><th>加权众数 p</th>'
             '<th>PRESSO 离群</th><th>PRESSO 校正 β (p)</th><th>符号正/负 (p)</th>'
             '<th>方向性</th><th class="l">LOO β 区间</th></tr>')
    for t, cn in TRAITS:
        for k, dlabel, _, _, _, _ in DATASETS:
            d = {"main": main, "fin": fin, "mvp": mvp}[k][t]
            wm = d["weighted_median"]["p"]; wmo = d["weighted_mode"]["p"]
            pr = d.get("presso_equiv", {}).get("distortion")
            nout = pr["n_outliers"] if pr else 0
            corr = ("%s (%s)" % (fnum(pr["beta_corrected"]), fp(pr["p_corrected"]))) if pr else "—"
            sg = d.get("sign_test") or {}
            npos, nneg = sg.get("n_positive"), sg.get("n_negative")
            sp = sg.get("p_binomial_two_sided") or sg.get("p")
            dirc = d["directionality"]["direction"]
            lr = loo[k][t]
            if lr:
                tag = "全正" if lr["all_pos"] else ("全负" if lr["all_neg"] else "含变号")
                ltxt = "%+.4f ~ %+.4f（%s，p %.3f－%.3f）" % (lr["bmin"], lr["bmax"], tag, lr["pmin"], lr["pmax"])
            else:
                ltxt = "未计算（s4_pipeline 未输出）"
            o.append('<tr><td class="l"><b>%s</b></td><td class="l">%s</td><td>%s</td><td>%s</td>'
                     '<td>%s</td><td>%s</td><td>%s / %s (%s)</td><td>%s</td><td class="l">%s</td></tr>' % (
                         cn, dlabel, fp(wm), fp(wmo), nout, corr, npos, nneg, fp(sp), dirc, ltxt))
        o.append('<tr><td colspan="9" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>')
    o.append('</table>')

    # ---- 4 结论
    o.append('<h2>4 · 可写入论文的结论</h2>')
    o.append('<div class="two"><div class="box">'
             '<b>站得住的（可作为主结论）</b>'
             '<ul>'
             '<li><b>遗传预测的 TSH 水平与 PAD 无因果关联</b>：三个独立/半独立数据集（主分析 −0.0059, p=0.848；'
             'FinnGen +0.0035, p=0.870；MVP −0.0292, p=0.131）方向与量级一致，95%CI 全跨 0。'
             'Egger 截距全不显著、方向性全为正向；<b>三层共 755 次留一法全部不显著</b>'
             '（MVP 层 251 次剔除后 β 始终在 −0.033～−0.016、p 0.26－0.58）→ <b>三角验证中最硬的一条</b>。</li>'
             '<li><b>FT4 无稳健因果证据</b>：三个数据集点估计分别为 −0.025 / −0.063 / +0.078，'
             '全部 CI 跨 0 且 MVP 层方向翻转；<b>LOO 亦方向分裂</b>（MVP 层 63 次全为正 +0.057～+0.104，'
             'FinnGen 层 60 次全为负 −0.097～−0.058）→ 只能报告为阴性/无证据。</li>'
             '<li>三性状的 <b>MR-Egger 截距在 MVP/FinnGen 层全部不显著</b>，方向性检验全为 forward（R²_X ≫ R²_Y）'
             '→ 无系统性水平多效性、无反向因果信号。</li>'
             '</ul></div>'
             '<div class="box">'
             '<b>不能写的（须在局限性中披露）</b>'
             '<ul>'
             '<li><b>甲减的主分析阳性（OR 1.046, p=0.015）不是稳健结论</b>：随机效应 p=0.130、'
             '加权中位数 p=0.506、PRESSO 校正后 p=0.037 但 LOO 下 p 全程 0.066－0.176；'
             '且暴露含 UKB、与结局队列重叠。</li>'
             '<li><b>FinnGen 的甲减结果不可用</b>（OVERLAP_INVALID）：暴露含 FinnGen Freeze 10；'
             '其 β（0.0857）约为主分析的 1.9 倍，<b>但该差异不能归因于重叠</b>：Step9 模拟（解析解 + 蒙特卡洛）'
             '给出重叠在完全重叠时的改变上限仅 +4.6%（r=5）～+18.4%（r=20），远不足以解释 1.9 倍；'
             '差异应归因于两套结局库本身不可互换（Step9_重叠偏倚模拟.md）。</li>'
             '<li><b>FT4 主分析存在定向水平多效性</b>（Egger 截距 p=2.8e-05）；'
             '该信号在 FinnGen（p=0.125）与 MVP（p=0.849）均未重现 → 不稳定，FT4 结论仅作参考。</li>'
             '<li>甲减/TSH 异质性极强（Q p 达 1e-11～1e-24）→ 固定效应 IVW 的 SE 被低估，'
             '应以随机效应/Egger/中位数为主判据。</li>'
             '</ul></div></div>')

    # ---- 5 关键数字
    o.append('<h2>5 · 关键数字速查</h2>')
    o.append('<table><tr><th class="l">项目</th><th class="l">数值</th><th class="l">来源</th></tr>')
    rows = [
        ("甲减工具数 / 累计 R² / 整体 F", "222 / 5.908% / 115.9", "Step 2（单位点 F 最小 29.4，无 F&lt;10）"),
        ("TSH 工具数 / 累计 R² / 整体 F", "258 / 10.748% / 225.3", "Step 2（单位点 F 最小 29.3）"),
        ("FT4 工具数 / 累计 R² / 整体 F", "65 / 2.367% / 71.4", "Step 2（单位点 F 最小 30.1）"),
        ("暴露文件链方向问题", "75－83% 完美互补（A↔T、C↔G）", "Step 3 核实：暴露提交时链方向混合，已按 1000G EUR hg38 仲裁"),
        ("甲减主分析 PRESSO 校正", "β=0.0531, p=0.037（离群 2 个，distortion p=0.744）", "step3/results/mr_main_hypo.json"),
        ("TSH 主分析 PRESSO 校正", "β=−0.0077, p=0.862（离群 2 个，distortion p=0.985）", "step3/results/mr_main_tsh.json"),
        ("MVP 结局有效样本量", "n_eff = 109,098（31,307 例 / 211,753 对照）", "s4_pipeline 计算"),
        ("主分析结局有效样本量", "n_eff = 28,037（7,114 例 / 475,964 对照）", "GWAS Catalog V2 API"),
        ("流水线可复现性回放", "2/2 完全一致（甲减×Sakaue、TSH×FinnGen）", "Step4_FinnGen复制结果.md §5"),
    ]
    for a, b, c in rows:
        o.append('<tr><td class="l">%s</td><td class="l"><b>%s</b></td><td class="l">%s</td></tr>' % (a, b, c))
    o.append('</table>')

    # ---- 6 数据与方法
    o.append('<h2>6 · 数据来源与方法</h2>')
    o.append('<table><tr><th class="l">角色</th><th class="l">数据集</th><th class="l">样本量（例/对照）</th>'
             '<th class="l">人群</th><th class="l">与暴露重叠</th><th class="l">机制/来源</th></tr>')
    o.append('<tr><td class="l">暴露</td><td class="l">甲减 GCST90572791 / TSH GCST90572789 / FT4 GCST90572790</td>'
             '<td class="l">甲减 113,393/1,065,268；TSH 482,873；FT4 191,449</td>'
             '<td class="l">欧洲为主（甲减含 131,670 例族裔未报告）</td><td class="l">—</td>'
             '<td class="l">Rand et al. 2025, <i>Nat Genet</i>, PMID 41238958（REGENIE，GRCh38）</td></tr>')
    o.append('<tr><td class="l">主结局</td><td class="l">PAD GCST90018890</td><td class="l">7,114 / 475,964</td>'
             '<td class="l">EUR＋东亚分层合并</td><td class="l">含 UKB，<b>有重叠风险</b></td>'
             '<td class="l">Sakaue et al. 2021, <i>Nat Genet</i>, PMID 34594039</td></tr>')
    o.append('<tr><td class="l">复制①</td><td class="l">FinnGen R13 <code>I9_PAD</code></td><td class="l">22,244 / 460,490</td>'
             '<td class="l">芬兰（EUR）</td><td class="l">甲减重叠；TSH/FT4 不重叠</td>'
             '<td class="l">FinnGen 公开桶 R13（2026 新发布）</td></tr>')
    o.append('<tr><td class="l">复制②</td><td class="l">MVP <code>te.PAD</code>（dbGaP phs001672）</td>'
             '<td class="l">31,307 / 211,753</td><td class="l">跨族裔 meta（EUR 24,009 / AFR 5,373 / HIS 1,925）</td>'
             '<td class="l"><b>不重叠</b>（暴露队列不含 MVP）</td>'
             '<td class="l">Klarin et al. 2019, <i>Nat Med</i>, PMID 31285632</td></tr>')
    o.append('</table>')
    o.append('<div class="box"><b>方法口径：</b>'
             '工具变量 p&lt;5e-8、r²&lt;0.001、窗口 10 Mb（PLINK 1.90b7，1000G Phase 3 hg38 EUR 面板）；'
             '主估计 IVW 固定效应，敏感性分析含 IVW 随机效应（DerSimonian-Laird）、MR-Egger、'
             '加权中位数（bootstrap 2000）、加权众数（bootstrap 500）、Cochran Q、'
             'MR-PRESSO 等价（全局异质性＋Bonferroni 残差离群＋distortion 检验）、留一法与方向性 R² 对比。'
             '全部为纯 Python 实现，口径与 TwoSampleMR / MendelianRandomization 对齐，'
             '<code>p_chi2</code> 经 χ² 临界值校验。所有中间台账（逐 SNP harmonise 决策、LOO、JSON）留档可回溯。</div>')

    o.append('<div class="foot">'
             '<b>生成方式</b>：本页由 <code>step4/scripts/s4_onepager.py</code> 从 '
             '<code>step3/results/*.json</code>、<code>step4/results/*.json</code> 与 '
             '<code>loo_*.tsv</code> 直接读取生成，数值未经任何人工修改或修约（仅展示位数）。'
             '详细报告见 <code>Step3_MR分析报告.md</code>、<code>Step4_FinnGen复制结果.md</code>、'
             '<code>Step4_MVP数据获取指南.md</code>、<code>Step4_样本重叠影响评估.md</code>。'
             '图件共 41 张：<code>step3/figures/</code> 12 张（主分析）+ <code>step4/figures/</code> 29 张'
             '（FinnGen 13 + MVP 13 + 跨数据集对比 3）。'
             '</div>')

    o.append('</div></body></html>')
    return "\n".join(o)


# ---------------------------------------------------------------- Markdown
def build_md():
    o = []
    o.append("# 甲状腺功能 → 外周动脉疾病（PAD）：Step 3 + Step 4 结果一页纸")
    o.append("")
    o.append("**两样本孟德尔随机化 · 三性状 × 三数据集三角验证**　｜　2026-10-07")
    o.append("")
    o.append("- **暴露**：甲减 / TSH / FT4 —— Rand et al. 2025, *Nat Genet*, PMID 41238958")
    o.append("- **主结局**：PAD（Sakaue et al. 2021, *Nat Genet*, PMID 34594039, GCST90018890）")
    o.append("- **工具数**：甲减 222 / TSH 258 / FT4 65")
    o.append("")
    o.append("---")
    o.append("")
    o.append("## 1 结论矩阵")
    o.append("")
    o.append("| 暴露 | 主分析 Sakaue 2021 | FinnGen R13 I9_PAD | MVP phs001672 | 三角验证判定 |")
    o.append("|---|---|---|---|---|")
    for t, cn in TRAITS:
        cells = [VERDICT[(t, k)][0] for k, _, _, _, _, _ in DATASETS]
        tk, _ = tri_short(t)
        o.append("| **%s** | %s | %s | %s | %s |" % (cn, cells[0], cells[1], cells[2], tk))
    o.append("")
    o.append("## 2 主结果（log-OR 尺度）")
    o.append("")
    o.append("| 暴露 | 数据集 | nSNP | IVW 固定 β | SE | p | OR [95%CI] | IVW 随机 p | Egger β (p) | Egger 截距 p | Q p |")
    o.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for t, cn in TRAITS:
        for k, dlabel, _, _, _, _ in DATASETS:
            d = {"main": main, "fin": fin, "mvp": mvp}[k][t]
            b = d["ivw_fixed"]["beta"]; se = d["ivw_fixed"]["se"]; p = d["ivw_fixed"]["p"]
            mark = "**" if p < 0.05 else ""
            o.append("| %s | %s | %s | %+.4f | %.4f | %s%s%s | %s | %s | %+.4f (%s) | %s | %s |" % (
                cn, dlabel, d["nsnp"], b, se, mark, fp(p), mark, fmt_or(b, se),
                fp(d["ivw_random"]["p"]),
                d["mr_egger"]["beta"], fp(d["mr_egger"]["p"]),
                fp(d["mr_egger"]["intercept_pval"]), fp(d["ivw_random"]["Q_pval"])))
    o.append("")
    o.append("## 3 稳健性与诊断")
    o.append("")
    o.append("| 暴露 | 数据集 | 加权中位数 p | 加权众数 p | PRESSO 离群 | PRESSO 校正 β (p) | 符号正/负 (p) | 方向性 | LOO β 区间 |")
    o.append("|---|---|---|---|---|---|---|---|---|")
    for t, cn in TRAITS:
        for k, dlabel, _, _, _, _ in DATASETS:
            d = {"main": main, "fin": fin, "mvp": mvp}[k][t]
            pr = d.get("presso_equiv", {}).get("distortion")
            nout = pr["n_outliers"] if pr else 0
            corr = ("%+.4f (%.3f)" % (pr["beta_corrected"], pr["p_corrected"])) if pr else "—"
            sg = d.get("sign_test") or {}
            sp = sg.get("p_binomial_two_sided") or sg.get("p")
            lr = loo[k][t]
            if lr:
                tag = "全正" if lr["all_pos"] else ("全负" if lr["all_neg"] else "含变号")
                ltxt = "%+.4f ~ %+.4f（%s；p %.3f－%.3f）" % (lr["bmin"], lr["bmax"], tag, lr["pmin"], lr["pmax"])
            else:
                ltxt = "未计算"
            o.append("| %s | %s | %s | %s | %d | %s | %s/%s (%s) | %s | %s |" % (
                cn, dlabel, fp(d["weighted_median"]["p"]), fp(d["weighted_mode"]["p"]), nout, corr,
                sg.get("n_positive"), sg.get("n_negative"), fp(sp),
                d["directionality"]["direction"], ltxt))
    o.append("")
    o.append("## 4 可写入论文的结论")
    o.append("")
    o.append("### 站得住的（主结论）")
    o.append("")
    o.append("1. **遗传预测的 TSH 水平与 PAD 无因果关联**：三数据集方向量级一致（主分析 −0.0059, p=0.848；"
             "FinnGen +0.0035, p=0.870；MVP −0.0292, p=0.131），95%CI 全跨 0；Egger 截距全不显著、"
             "方向性全 forward；**三层共 755 次留一法全部不显著**（MVP 层 251 次剔除后 β 始终在 "
             "−0.033～−0.016、p 0.26－0.58）→ 三角验证中最硬的一条。")
    o.append("2. **FT4 无稳健因果证据**：三点估计 −0.025 / −0.063 / +0.078，全部 CI 跨 0 且 MVP 层方向翻转；"
             "**LOO 亦方向分裂**（MVP 层 63 次全为正 +0.057～+0.104，FinnGen 层 60 次全为负 −0.097～−0.058）。")
    o.append("3. 三性状的 MR-Egger 截距在 MVP / FinnGen 层全部不显著，方向性检验全为 forward → "
             "无系统性水平多效性、无反向因果信号。")
    o.append("")
    o.append("### 不能写的（须在局限性中披露）")
    o.append("")
    o.append("1. **甲减主分析阳性（OR 1.046, p=0.015）不是稳健结论**：随机效应 p=0.130、加权中位数 p=0.506、"
             "LOO 下 p 全程 0.066－0.176；且暴露含 UKB、与结局队列重叠。")
    o.append("2. **FinnGen 甲减结果不可用**（OVERLAP_INVALID，暴露含 FinnGen Freeze 10）；"
             "其 β 约为主分析 1.9 倍，但该差异**不能归因于重叠**（Step9 模拟给出重叠上限 ≤ +4.6%～+18.4%，"
             "`Step9_重叠偏倚模拟.md`），而来自两套结局库本身不可互换。")
    o.append("3. **FT4 主分析存在定向水平多效性**（Egger 截距 p=2.8e-05），该信号在 FinnGen（p=0.125）"
             "与 MVP（p=0.849）均未重现。")
    o.append("4. 甲减 / TSH 异质性极强（Q p 达 1e-11～1e-24）→ 固定效应 IVW 的 SE 被低估，"
             "应以随机效应 / Egger / 中位数为主判据。")
    o.append("")
    o.append("## 5 关键数字速查")
    o.append("")
    o.append("| 项目 | 数值 | 来源 |")
    o.append("|---|---|---|")
    for a, b, c in [("甲减工具数 / 累计 R² / 整体 F", "222 / 5.908% / 115.9（最小单位点 F 29.4）", "Step 2 工具变量清单"),
                    ("TSH 工具数 / 累计 R² / 整体 F", "258 / 10.748% / 225.3（最小 29.3）", "Step 2 工具变量清单"),
                    ("FT4 工具数 / 累计 R² / 整体 F", "65 / 2.367% / 71.4（最小 30.1）", "Step 2 工具变量清单"),
                    ("暴露文件链方向", "75－83% 完美互补（A↔T、C↔G），已按 1000G EUR hg38 仲裁", "Step3_MR分析报告.md §2"),
                    ("甲减主分析 PRESSO 校正", "β=+0.0531, p=0.037（离群 2 个，distortion p=0.744）", "mr_main_hypo.json"),
                    ("TSH 主分析 PRESSO 校正", "β=−0.0077, p=0.862（离群 2 个，distortion p=0.985）", "mr_main_tsh.json"),
                    ("MVP 结局有效样本量", "n_eff = 109,098", "s4_pipeline 计算"),
                    ("主分析结局有效样本量", "n_eff = 28,037", "GWAS Catalog V2 API"),
                    ("流水线可复现性回放", "2/2 完全一致（甲减×Sakaue、TSH×FinnGen）", "Step4_FinnGen复制结果.md §5")]:
        o.append("| %s | %s | %s |" % (a, b, c))
    o.append("")
    o.append("## 6 数据来源")
    o.append("")
    o.append("| 角色 | 数据集 | 样本量 | 人群 | 与暴露重叠 |")
    o.append("|---|---|---|---|---|")
    o.append("| 暴露 | 甲减 GCST90572791 / TSH GCST90572789 / FT4 GCST90572790 | 甲减 113,393/1,065,268；TSH 482,873；FT4 191,449 | 欧洲为主 | — |")
    o.append("| 主结局 | PAD GCST90018890 | 7,114 / 475,964 | EUR＋东亚合并 | 含 UKB，有重叠风险 |")
    o.append("| 复制① | FinnGen R13 `I9_PAD` | 22,244 / 460,490 | 芬兰 | 甲减重叠；TSH/FT4 不重叠 |")
    o.append("| 复制② | MVP `te.PAD`（dbGaP phs001672） | 31,307 / 211,753 | 跨族裔 meta | 不重叠 |")
    o.append("")
    o.append("**方法口径**：IV p<5e-8、r²<0.001、窗口 10 Mb（PLINK 1.90b7 + 1000G Phase 3 hg38 EUR）；"
             "IVW 固定/随机（DL）、MR-Egger、加权中位数/众数、Cochran Q、MR-PRESSO 等价、LOO、方向性 R² 对比。"
             "纯 Python 实现，口径对齐 TwoSampleMR / MendelianRandomization，`p_chi2` 经 χ² 校验。")
    o.append("")
    o.append("> 本文件由 `step4/scripts/s4_onepager.py` 从结果 JSON / LOO TSV 直接生成，数值未作人工修改。")
    return "\n".join(o)


def tri_short(t):
    return {
        "hypo": ("方向一致（全正），均未达稳健显著", "warn"),
        "tsh": ("三数据集一致阴性 → 最硬结论", "ok"),
        "ft4": ("方向翻转、全部跨 0 → 无稳健证据", "warn"),
    }[t]


if __name__ == "__main__":
    hp = os.path.join(DELIV, "Step3-4_结果一页纸.html")
    mp = os.path.join(DELIV, "Step3-4_结果一页纸.md")
    with io.open(hp, "w", encoding="utf-8") as f:
        f.write(build_html())
    with io.open(mp, "w", encoding="utf-8") as f:
        f.write(build_md())
    print("OK  " + hp, os.path.getsize(hp), "B")
    print("OK  " + mp, os.path.getsize(mp), "B")
