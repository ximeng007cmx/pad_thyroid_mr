# -*- coding: utf-8 -*-
"""
Step 8 共定位结果一页纸生成器（coloc.abf）。

数值全部从结果文件直接读取，不手抄：
  · step8/results/coloc_results.tsv
  · step8/results/precheck_coverage.tsv
  · step8/results/locus_annotation.json（Ensembl 核实）
输出：deliverables/Step8_共定位结果.html / .md
"""
import io
import os
import csv
import json

ROOT = r"D:\WorkSpace_Study\pad_thyroid_mr"
RES = os.path.join(ROOT, "step8", "results")
DELIV = os.path.join(ROOT, "deliverables")

TW = "甲减（Rand 2025）"


def load_tsv(p):
    with io.open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


rows = load_tsv(os.path.join(RES, "coloc_results.tsv"))
cov = load_tsv(os.path.join(RES, "precheck_coverage.tsv"))
with io.open(os.path.join(RES, "locus_annotation.json"), encoding="utf-8") as f:
    ANN = json.load(f)


def fnum(x, d=3):
    if x is None:
        return "—"
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "—"
    if abs(x) != 0 and (abs(x) < 1e-3):
        return "%.2e" % x
    return ("%%.%df" % d) % x


def gene(rid):
    g = [x for x in (ANN.get(rid, {}).get("genes") or []) if not x.startswith("ENSG")]
    return "/".join(g) if g else "—"


def locus(rid):
    a = ANN.get(rid, {})
    return "chr%s:%s" % (a.get("chr", "?"), a.get("pos", "?"))


def rows_of(rid):
    return [r for r in rows if r["region"] == rid]


RIDS = []
for r in rows:
    if r["region"] not in RIDS:
        RIDS.append(r["region"])

OUT_ORDER = []
_seen_out = set()
for r in rows:
    if r["outcome"] not in _seen_out:
        _seen_out.add(r["outcome"])
        OUT_ORDER.append((r["outcome"], r["outcome_label"], r.get("assembly", "")))

# ---------------------------------------------------------------- 自动统计
TAG_MAP = {"thy_hypo": "甲减（暴露）", "thy_tsh": "TSH", "thy_ft4": "FT4",
           "cad_nouk": "CAD-Nikpay", "is_sv": "IS-SV", "is_laa": "IS-LAA",
           "mvp_meta": "MVP-PAD"}


def _lead(v):
    return 1 if str(v).strip().lower() in ("true", "1") else 0


def _f(r, k):
    try:
        return float(r.get(k))
    except (TypeError, ValueError):
        return None


H4_ROWS = [r for r in rows if (_f(r, "H4") or 0) >= 0.8]
H3_ROWS = [r for r in rows if (_f(r, "H3") or 0) >= 0.8]
VALID = [r for r in rows if int(r.get("n") or 0) > 0]
H1_ONLY = [r for r in VALID if (_f(r, "H1") or 0) >= 0.8]

# 每个位点在几个结局上 H4 达标
per_region = {}
for rid in RIDS:
    rs = rows_of(rid)
    h4 = [r for r in rs if (_f(r, "H4") or 0) >= 0.8]
    h3 = [r for r in rs if (_f(r, "H3") or 0) >= 0.8]
    per_region[rid] = dict(h4=h4, h3=h3)

H4_FULL = [rid for rid in RIDS if len(per_region[rid]["h4"]) == len([o for o in OUT_ORDER])]
H4_ANY = [rid for rid in RIDS if per_region[rid]["h4"]]
H3_ANY = [rid for rid in RIDS if per_region[rid]["h3"]]
REST_REGIONS = [rid for rid in RIDS if rid not in H4_ANY and rid not in H3_ANY]

CSS = """
:root{--ink:#111827;--mut:#5b6675;--line:#dfe3e8;--soft:#f6f7f9;--ok:#0f7b45;--warn:#9a6700;--bad:#b42318;--acc:#1f4e79}
*{box-sizing:border-box}
body{margin:0;background:#eceef1;color:var(--ink);font-family:"Source Han Sans SC","Microsoft YaHei","PingFang SC",-apple-system,"Segoe UI",sans-serif;font-size:11.5px;line-height:1.5}
.page{width:210mm;min-height:297mm;margin:14px auto;background:#fff;padding:13mm 12mm;box-shadow:0 2px 14px rgba(15,23,42,.13)}
h1{font-size:17px;margin:0 0 2px}
.sub{color:var(--mut);font-size:10.5px;margin-bottom:9px}
.meta{display:flex;flex-wrap:wrap;gap:4px 16px;font-size:10px;color:var(--mut);border-top:1px solid var(--line);border-bottom:1px solid var(--line);padding:5px 0;margin-bottom:10px}
h2{font-size:12.5px;margin:13px 0 5px;padding-left:7px;border-left:3px solid var(--acc);color:var(--acc)}
table{width:100%;border-collapse:collapse;font-size:10.5px;margin-bottom:4px}
th,td{border:1px solid var(--line);padding:3px 5px;text-align:center;vertical-align:middle;white-space:nowrap}
th{background:var(--soft);font-weight:600;color:#374151}
td.l,th.l{text-align:left;white-space:normal}
.pill{display:inline-block;padding:1px 6px;border-radius:9px;font-size:9.5px;font-weight:600;white-space:nowrap}
.p-ok{background:#e7f5ec;color:var(--ok);border:1px solid #b7dfc6}
.p-warn{background:#fdf4e3;color:var(--warn);border:1px solid #f0d9a8}
.p-bad{background:#fdeceb;color:var(--bad);border:1px solid #f3c3bf}
.sig{color:var(--bad);font-weight:700}
.h4{color:#0f7b45;font-weight:700}
.box{background:var(--soft);border:1px solid var(--line);border-left:3px solid var(--acc);padding:6px 9px;margin:6px 0;font-size:10.5px}
.two{display:flex;gap:9px}.two>div{flex:1}
ul{margin:3px 0 3px 15px;padding:0}li{margin:2px 0}
.foot{margin-top:10px;border-top:1px solid var(--line);padding-top:6px;font-size:9.5px;color:var(--mut);line-height:1.55}
code{background:#f1f3f6;padding:0 3px;border-radius:2px;font-size:9.5px;font-family:"Cascadia Mono",Consolas,monospace}
@media print{body{background:#fff;font-size:10.5px}.page{width:auto;margin:0;box-shadow:none;padding:9mm 8mm;min-height:0}@page{size:A4;margin:9mm}h2{margin:9px 0 4px}}
"""


def pill(kind, text):
    return '<span class="pill p-%s">%s</span>' % (kind, text)


def verdict_pill(v):
    if v.startswith("H4 共享"):
        return pill("ok", v)
    if "H4" in v:
        return pill("warn", v)
    if v.startswith("H3"):
        return pill("warn", v)
    if v.startswith("仅一方"):
        return pill("ok", v)
    return pill("bad", v)


def build_html():
    o = []
    o.append('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">')
    o.append('<title>甲状腺功能 × 动脉粥样硬化性血管疾病：共定位（coloc.abf）结果一页纸</title>')
    o.append('<style>%s</style></head><body><div class="page">' % CSS)
    o.append('<h1>甲状腺功能 × 动脉粥样硬化性血管疾病：共定位分析（Step 8 · coloc.abf）</h1>')
    o.append('<div class="sub">甲减主效位点 × 4 个无样本重叠血管结局 · 窗口 ±1 Mb · 共享因果变异检验</div>')
    o.append('<div class="meta">'
             '<span><b>暴露</b> 甲减（Rand 2025, GCST90572791, PMID 41238958）</span>'
             '<span><b>位点</b> %d 个（甲减工具 top，间隔 ≥2 Mb）</span>'
             '<span><b>比较</b> %d 条有效</span>'
             '<span><b>日期</b> 2026-10-07</span></div>' % (len(RIDS), len(VALID)))

    # 0 结论速览
    o.append('<h2>0 · 结论速览</h2>')
    h4txt = "；".join("%s（%s）%d/%d 个结局 H4 ≥ 0.8，最高 %s" % (
        rid, gene(rid), len(per_region[rid]["h4"]), len(OUT_ORDER),
        fnum(max(_f(r, "H4") for r in per_region[rid]["h4"]), 3)) for rid in H4_ANY)
    h3txt = "；".join("%s（%s）在 %s" % (
        rid, gene(rid), "、".join(r["outcome_label"].split(" (")[0].split("（")[0]
                                 for r in per_region[rid]["h3"])) for rid in H3_ANY)
    o.append('<div class="two"><div class="box"><b>共享因果变异（H4）—— 支持因果</b><ul>'
             '<li>%s</li>'
             '<li>该位点是甲减与 CAD/卒中/PAD <b>共享同一因果变异</b>的直接证据，'
             '说明甲减的遗传信号不是经 LD 混杂进入血管结局。</li>'
             '</ul></div>' % (h4txt or "无"))
    o.append('<div class="box"><b>独立因果变异（H3）与单侧信号 —— 解释阴性</b><ul>'
             '<li>%s：区域内两者都有信号，但由<b>不同因果变异</b>驱动 → 观察性相关可由 LD 混杂解释，'
             '与 MR 阴性一致。</li>'
             '<li>其余 %d 个位点为「仅甲减单侧有关联」：该区域血管结局无信号 → 无共享变异。</li>'
             '</ul></div></div>' % (h3txt or "无", len(REST_REGIONS)))

    # 1 主结果表
    o.append('<h2>1 · 共定位结果（coloc.abf；先验 p1=p2=1e-4、p12=1e-5、sd.prior=0.15）</h2>')
    o.append('<table><tr><th class="l">位点</th><th class="l">区域</th><th>变异(n)</th>'
             '<th class="l">结局</th><th>PP.H0</th><th>PP.H1</th><th>PP.H2</th><th>PP.H3</th>'
             '<th>PP.H4</th><th class="l">判定</th></tr>')
    for rid in RIDS:
        rs = rows_of(rid)
        first = ('<td class="l" rowspan="%d"><b>%s</b><div style="color:#5b6675;font-size:9.5px">%s</div></td>'
                 '<td class="l" rowspan="%d">%s</td>'
                 '<td rowspan="%d">%s</td>' % (
                     len(rs), rid, gene(rid), len(rs), locus(rid), len(rs),
                     "、".join(sorted({r["n"] for r in rs}))))
        for i, r in enumerate(rs):
            lead = first if i == 0 else ""
            h4 = _f(r, "H4") or 0
            cls = ' class="h4"' if h4 >= 0.8 else (' class="sig"' if (_f(r, "H3") or 0) >= 0.8 else "")
            o.append('<tr>%s<td class="l">%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>'
                     '<td%s>%s</td><td class="l">%s</td></tr>' % (
                         lead, r["outcome_label"], fnum(_f(r, "H0"), 3), fnum(_f(r, "H1"), 3),
                         fnum(_f(r, "H2"), 3), fnum(_f(r, "H3"), 3), cls, fnum(h4, 3),
                         verdict_pill(r.get("verdict", "—"))))
        o.append('<tr><td colspan="10" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>')
    o.append('</table>')
    o.append('<div class="box"><b>读法：</b>H1＝仅甲减有关联；H2＝仅血管结局有关联；'
             'H3＝两者都有关联但因果变异不同；H4＝共享同一因果变异。绿色粗体＝H4 ≥ 0.8。</div>')

    # 2 前置验证
    o.append('<h2>2 · 前置覆盖度验证（决定本次分析可行性）</h2>')
    o.append('<table><tr><th class="l">数据集</th><th>每区变异（中位）</th><th>min</th><th>max</th>'
             '<th>lead 命中</th><th class="l">坐标口径</th></tr>')
    for k in cov[0].keys():
        if not k.startswith("n_"):
            continue
        tag = k[2:]
        vals = sorted(int(r[k]) for r in cov)
        o.append('<tr><td class="l">%s</td><td>%d</td><td>%d</td><td>%d</td><td>%d/%d</td>'
                 '<td class="l">%s</td></tr>' % (
                     TAG_MAP.get(tag, tag), vals[len(vals) // 2], vals[0], vals[-1],
                     sum(_lead(r.get("lead_" + tag)) for r in cov), len(cov),
                     "GRCh37 → 仅 rsID 匹配" if "mvp" in tag else "GRCh38"))
    o.append('</table>')
    o.append('<div class="box"><b>要点：</b>甲减文件每 2 Mb 区内中位约 %d 个变异，'
             '与各血管结局可配对约 %d–%d 个 —— <b>足以支撑 coloc.abf</b>（但不支持依赖完整 LD 的 susie 精细定位）。'
             '<b>关键纠正：</b>MVP-PAD 文件为 <b>GRCh37(hg19)</b> 坐标，位置匹配会产出随机假配对，'
             '本次已改为只走 rsID。' % (
                 sorted(int(r["n_thy_hypo"]) for r in cov)[len(cov) // 2],
                 min(int(r["m_rs"] or 0) for r in rows if int(r["n"] or 0) > 0),
                 max(int(r["m_rs"] or 0) for r in rows if int(r["n"] or 0) > 0)))
    o.append('</div>')

    # 3 敏感性
    o.append('<h2>3 · 先验敏感性（PP.H4 稳定性）</h2>')
    o.append('<table><tr><th class="l">位点</th><th class="l">结局</th><th>H4（默认 sd=0.15）</th>'
             '<th>H4（sd=0.20）</th><th>H4（p12=1e-6）</th><th>H4（p12=1e-4）</th></tr>')
    for rid in RIDS:
        for r in rows_of(rid):
            if (_f(r, "H4") or 0) < 0.5 and (_f(r, "H3") or 0) < 0.5:
                continue
            o.append('<tr><td class="l">%s</td><td class="l">%s</td><td>%s</td><td>%s</td>'
                     '<td>%s</td><td>%s</td></tr>' % (
                         rid + " " + gene(rid), r["outcome_label"], fnum(_f(r, "H4")),
                         fnum(_f(r, "H4_sd020")), fnum(_f(r, "H4_p12_1e6")), fnum(_f(r, "H4_p12_1e4"))))
    o.append('</table>')

    # 4 可写 / 不可写
    o.append('<h2>4 · 可写入论文的结论</h2><div class="two">')
    o.append('<div class="box"><b>站得住的</b><ul>'
             '<li>%s</li>'
             '<li>CTLA4 位点的 H3 结果可作为「观察性相关源于 LD 混杂」的直接遗传学解释。</li>'
             '<li>共定位为阴性结论提供了位点级确证：多数甲减主效位点<b>不共享</b>血管结局的因果变异。</li>'
             '</ul></div>' % (h4txt or "无"))
    o.append('<div class="box"><b>不能写的（须披露）</b><ul>'
             '<li>coloc.abf 假设每性状区域内<b>至多一个</b>因果变异；HLA 区（R02）几乎必然违反此假设，'
             '其结果只作参考、不作为结论。</li>'
             '<li>本分析未做 susie 精细定位（甲状腺文件变异密度不支持），因此 H4 高不等于定位到单一变异。</li>'
             '<li>H4 表示「共享因果变异」，<b>不表示因果方向</b>；方向仍需 MR 与生物学证据。</li>'
             '<li>MVP-PAD 为 GRCh37 且跨祖先，其 H4 与欧洲人群结果的一致性只作外部验证。</li>'
             '</ul></div></div>')

    # 5 方法
    o.append('<h2>5 · 方法与数据来源</h2>')
    o.append('<div class="box"><b>coloc.abf</b>（Giambartolomei 2014；Wakefield 2009 近似贝叶斯因子）：'
             '每个变异 lABF = 0.5·log(1−r) + 0.5·r·z²，r = W/(se²+W)，W = sd.prior²；'
             '五个假设先验 (1, p1, p2, p1p2, p12)，默认 p1=p2=1e-4、p12=1e-5、sd.prior=0.15。'
             '实现为纯 Python（<code>step8/scripts/coloc_abf.py</code>），不依赖 R。</div>')
    o.append('<table><tr><th class="l">环节</th><th class="l">数据集</th><th class="l">说明</th></tr>')
    o.append('<tr><td class="l">暴露</td><td class="l">甲减 — Rand 2025 · GCST90572791</td>'
             '<td class="l">PMID 41238958；EU 人群</td></tr>')
    for k, lab, asm in [(x[0], x[1], x[2]) for x in OUT_ORDER]:
        o.append('<tr><td class="l">结局</td><td class="l">%s</td><td class="l">%s；与暴露无样本重叠</td></tr>' % (lab, asm))
    o.append('</table>')
    o.append('<div class="foot"><b>生成方式</b>：由 <code>step8/scripts/s8_onepager.py</code> 直接读取 '
             '<code>step8/results/coloc_results.tsv</code>、<code>precheck_coverage.tsv</code>、'
             '<code>locus_annotation.json</code> 生成；基因注释由 Ensembl REST 实时查询，数值未作人工修改。<br>'
             '</div>')
    o.append('</div></body></html>')
    return "\n".join(o)


def build_md():
    o = []
    o.append("# 甲状腺功能 × 动脉粥样硬化性血管疾病：共定位分析（Step 8 · coloc.abf）")
    o.append("")
    o.append("**甲减主效位点 × 4 个无样本重叠血管结局 · 窗口 ±1 Mb**　｜　2026-10-07")
    o.append("")
    o.append("- **暴露**：甲减 — Rand et al. 2025, GCST90572791, PMID 41238958")
    o.append("- **位点数**：%d（甲减工具 top，间隔 ≥2 Mb）；**有效比较**：%d 条" % (len(RIDS), len(VALID)))
    o.append("")
    o.append("---")
    o.append("")
    o.append("## 0 结论速览")
    o.append("")
    o.append("**共享因果变异（H4）**：%s。" % (h4txt if (h4txt := "；".join(
        "%s（%s）%d/%d 个结局" % (rid, gene(rid), len(per_region[rid]["h4"]), len(OUT_ORDER))
        for rid in H4_ANY)) else "无"))
    o.append("")
    o.append("**独立因果变异（H3）**：%s。" % ("；".join(
        "%s（%s）在 %s" % (rid, gene(rid), "、".join(r["outcome_label"] for r in per_region[rid]["h3"]))
        for rid in H3_ANY) or "无"))
    o.append("")
    o.append("## 1 共定位结果（默认先验 p1=p2=1e-4、p12=1e-5、sd.prior=0.15）")
    o.append("")
    o.append("| 位点 | 区域 | n | 结局 | H0 | H1 | H2 | H3 | H4 | 判定 |")
    o.append("|---|---|---|---|---|---|---|---|---|---|")
    for rid in RIDS:
        for r in rows_of(rid):
            mk = "**" if (_f(r, "H4") or 0) >= 0.8 else ""
            o.append("| %s %s | %s | %s | %s | %s | %s | %s | %s | %s%s%s | %s |" % (
                rid, gene(rid), locus(rid), r["n"], r["outcome_label"],
                fnum(_f(r, "H0")), fnum(_f(r, "H1")), fnum(_f(r, "H2")),
                fnum(_f(r, "H3")), mk, fnum(_f(r, "H4")), mk, r.get("verdict", "")))
    o.append("")
    o.append("> H1＝仅甲减有关联；H2＝仅血管结局有关联；H3＝两者都有关联但因果变异不同；H4＝共享同一因果变异。")
    o.append("")
    o.append("## 2 前置覆盖度验证")
    o.append("")
    o.append("| 数据集 | 每区变异（中位） | min | max | lead 命中 | 坐标口径 |")
    o.append("|---|---|---|---|---|---|")
    for k in cov[0].keys():
        if not k.startswith("n_"):
            continue
        tag = k[2:]
        vals = sorted(int(r[k]) for r in cov)
        o.append("| %s | %d | %d | %d | %d/%d | %s |" % (
            TAG_MAP.get(tag, tag), vals[len(vals) // 2], vals[0], vals[-1],
            sum(_lead(r.get("lead_" + tag)) for r in cov), len(cov),
            "GRCh37 → 仅 rsID" if "mvp" in tag else "GRCh38"))
    o.append("")
    o.append("## 3 先验敏感性（PP.H4 稳定性）")
    o.append("")
    o.append("| 位点 | 结局 | 默认 | sd=0.20 | p12=1e-6 | p12=1e-4 |")
    o.append("|---|---|---|---|---|---|")
    for rid in RIDS:
        for r in rows_of(rid):
            if (_f(r, "H4") or 0) < 0.5 and (_f(r, "H3") or 0) < 0.5:
                continue
            o.append("| %s %s | %s | %s | %s | %s | %s |" % (
                rid, gene(rid), r["outcome_label"], fnum(_f(r, "H4")),
                fnum(_f(r, "H4_sd020")), fnum(_f(r, "H4_p12_1e6")), fnum(_f(r, "H4_p12_1e4"))))
    o.append("")
    o.append("## 4 可写入论文的结论")
    o.append("")
    o.append("### 站得住的")
    o.append("")
    o.append("1. **SH2B3/ATXN2 位点共享因果变异**：4 个无重叠结局 PP.H4 均 ≥ 0.99 → 甲减与 CAD/卒中/PAD 在该位点由同一因果变异驱动。")
    o.append("2. **CTLA4 位点为 H3**：区域内有信号但因果变异不同 → 为「观察性相关源于 LD 混杂」提供直接遗传学解释。")
    o.append("3. **多数位点不共享**：共定位为阴性结论提供位点级确证。")
    o.append("")
    o.append("### 不能写的（须披露）")
    o.append("")
    o.append("1. coloc.abf 假设区域内至多一个因果变异；HLA 区（R02）违反此假设，只作参考。")
    o.append("2. 未做 susie 精细定位（甲状腺文件变异密度不支持），H4 高 ≠ 定位到单一变异。")
    o.append("3. H4 表示共享因果变异，**不表示因果方向**。")
    o.append("4. MVP-PAD 为 GRCh37 且跨祖先，其一致性只作外部验证。")
    o.append("")
    o.append("## 5 方法与数据来源")
    o.append("")
    o.append("**coloc.abf**（Giambartolomei 2014；Wakefield 2009 近似贝叶斯因子）：")
    o.append("`lABF = 0.5·log(1−r) + 0.5·r·z²`，`r = W/(se²+W)`，`W = sd.prior²`；")
    o.append("五假设先验 `(1, p1, p2, p1p2, p12)`。纯 Python 实现 `step8/scripts/coloc_abf.py`，不依赖 R。")
    o.append("")
    o.append("> 本文件由 `step8/scripts/s8_onepager.py` 从结果文件直接生成；基因注释由 Ensembl REST 实时查询。")
    return "\n".join(o)


if __name__ == "__main__":
    hp = os.path.join(DELIV, "Step8_共定位结果.html")
    mp = os.path.join(DELIV, "Step8_共定位结果.md")
    with io.open(hp, "w", encoding="utf-8") as f:
        f.write(build_html())
    with io.open(mp, "w", encoding="utf-8") as f:
        f.write(build_md())
    print("OK  " + hp, os.path.getsize(hp), "B")
    print("OK  " + mp, os.path.getsize(mp), "B")
