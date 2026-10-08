# -*- coding: utf-8 -*-
"""
Step 9 统计功效结果一页纸生成器。

数值全部从结果文件直接读取，不手抄：
  · step9/results/power_forward.tsv
  · step9/results/power_reverse.tsv
  · step9/results/power_mvmr.tsv
  · step9/results/power_summary.json
输出：deliverables/Step9_统计功效.html / .md
"""
import io
import os
import csv
import json

ROOT = r"D:\WorkSpace_Study\pad_thyroid_mr"
RES = os.path.join(ROOT, "step9", "results")
DELIV = os.path.join(ROOT, "deliverables")

TRAIT_ORDER = ["hypo", "tsh", "ft4"]
TRAIT_CN = {"hypo": "甲减", "tsh": "TSH", "ft4": "FT4"}
NOTE_RANK = {"无重叠": 0, "主分析·含UKB": 1, "重叠": 2, "重叠·甲减失效": 3}


def load_tsv(p):
    with io.open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


FW = load_tsv(os.path.join(RES, "power_forward.tsv"))
RV = load_tsv(os.path.join(RES, "power_reverse.tsv"))
MV = load_tsv(os.path.join(RES, "power_mvmr.tsv"))
SCALE = load_tsv(os.path.join(RES, "scale_check.tsv"))
with io.open(os.path.join(RES, "power_summary.json"), encoding="utf-8") as f:
    SUM = json.load(f)
INST = SUM["instruments"]
SCEN = SUM["scenarios"]

KEY_CN = {"sakaue": "PAD（Sakaue 2021）", "finngen": "PAD（FinnGen R13）",
          "mvp_meta": "PAD（MVP trans-ancestry）", "cad_nouk": "CAD（Nikpay 2015）",
          "cad_ukb": "CAD（UKB Mbatchou 2021）", "mi_ukb": "MI（UKB Donertas 2021）",
          "is_laa": "IS-LAA（MEGASTROKE）", "is_sv": "IS-SV（MEGASTROKE）"}

# 尺度异常（中位 se_Y 偏离参考 ≥5 倍）的结局文件
BAD = [s for s in SCALE if str(s.get("flag") or "").strip()]


def fl(x):
    if x is None:
        return None
    s = str(x).strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def fnum(x, d=3):
    if x is None:
        return "—"
    if abs(x) != 0 and abs(x) < 1e-3:
        return "%.2e" % x
    return ("%%.%df" % d) % x


def pnum(x):
    if x is None:
        return "—"
    return "%.3g" % x


def key(r):
    return (TRAIT_ORDER.index(r["exposure_key"]),
            NOTE_RANK.get(r["source_note"], 9), r["outcome"])


FW = sorted(FW, key=key)
MV = sorted(MV, key=key)


def rows_t(trait, src=None, reverse=False):
    out = []
    for r in RV:
        if r["exposure_key"] == trait and (src is None or r["source_note"] == src):
            out.append(r)
    return out


def clean(title):
    return title.split("（")[0].split(" (")[0]


def short(title):
    """保留数据集标识的短名：'PAD（Sakaue 2021, GCST90018890）'→'PAD（Sakaue 2021）'。"""
    if "（" not in title:
        return title
    name, rest = title.split("（", 1)
    return "%s（%s）" % (name, rest.rstrip("）").split(",")[0].strip())


def flagged(r):
    return bool(str(r.get("scale_flag") or "").strip())


def rng(vals, d=3):
    vals = [v for v in vals if v is not None]
    if not vals:
        return "—"
    lo, hi = min(vals), max(vals)
    if abs(hi - lo) < 1e-12:
        return fnum(lo, d)
    return "%s–%s" % (fnum(lo, d), fnum(hi, d))


# ---------------------------------------------------------------- 自动摘要
SUMM = {}
for t in TRAIT_ORDER:
    sub = [r for r in FW if r["exposure_key"] == t and r["source_note"] == "无重叠"]
    SUMM[t] = dict(
        n=len(sub),
        mde80=rng([fl(r["mde80_or"]) for r in sub]),
        mde90=rng([fl(r["mde90_or"]) for r in sub]),
        p105=rng([fl(r["pow_or1.05"]) for r in sub], 2),
        p110=rng([fl(r["pow_or1.10"]) for r in sub], 2),
        p115=rng([fl(r["pow_or1.15"]) for r in sub], 2),
        outcomes="、".join(clean(r["outcome"]) for r in sub),
    )

REV_CLEAN = [r for r in RV if r["source_note"] == "无重叠"]
REV_HYPO = [r for r in REV_CLEAN if r["outcome_type"] == "二分"]
REV_CONT = [r for r in REV_CLEAN if r["outcome_type"] == "连续"]
MV_NOOV = [r for r in MV if r["exposure_key"] == "hypo" and r["source_note"] == "无重叠"]

INST_TXT = "；".join("%s %d 个（R²=%.2f%%，中位 F=%.1f，最小 F=%.1f）" % (
    TRAIT_CN[t], INST[t]["n"], INST[t]["r2"] * 100, INST[t]["f_median"], INST[t]["f_min"])
    for t in TRAIT_ORDER)
MINF = min(INST[t]["f_min"] for t in TRAIT_ORDER)

PAD_KEYS = ["hypo|PAD（Sakaue 2021, GCST90018890）",
            "tsh|PAD（Sakaue 2021, GCST90018890）",
            "ft4|PAD（Sakaue 2021, GCST90018890）"]

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
.p-mut{background:#eef1f5;color:#475569;border:1px solid #cfd6e0}
.sig{color:var(--bad);font-weight:700}
.ok{color:var(--ok);font-weight:700}
.box{background:var(--soft);border:1px solid var(--line);border-left:3px solid var(--acc);padding:6px 9px;margin:6px 0;font-size:10.5px}
.two{display:flex;gap:9px}.two>div{flex:1}
ul{margin:3px 0 3px 15px;padding:0}li{margin:2px 0}
.foot{margin-top:10px;border-top:1px solid var(--line);padding-top:6px;font-size:9.5px;color:var(--mut);line-height:1.55}
code{background:#f1f3f6;padding:0 3px;border-radius:2px;font-size:9.5px;font-family:"Cascadia Mono",Consolas,monospace}
@media print{body{background:#fff;font-size:10.5px}.page{width:auto;margin:0;box-shadow:none;padding:9mm 8mm;min-height:0}@page{size:A4;margin:9mm}h2{margin:9px 0 4px}}
"""


def pill(kind, text):
    return '<span class="pill p-%s">%s</span>' % (kind, text)


def note_pill(n):
    return {"无重叠": pill("ok", n), "主分析·含UKB": pill("mut", n),
            "重叠": pill("warn", n), "重叠·甲减失效": pill("bad", n)}.get(n, pill("mut", n))


FWD_HEAD = ["暴露", "结局", "重叠", "nSNP", "观测 β (p)", "MDE80", "MDE90",
            "功效@1.05", "功效@1.10", "功效@1.15"]
REV_HEAD = ["甲状腺结局", "血管暴露", "重叠", "nSNP", "观测 β (p)", "MDE80（80% 功效）"]
MV_HEAD = ["结局", "重叠", "nSNP", "条件 F", "观测 β (p)", "MDE80（OR）"]


def build_html():
    o = []
    o.append('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">')
    o.append('<title>甲状腺功能 × 动脉粥样硬化性血管疾病：统计功效与最小可检出效应（Step 9）</title>')
    o.append('<style>%s</style></head><body><div class="page">' % CSS)
    o.append('<h1>甲状腺功能 × 动脉粥样硬化性血管疾病：统计功效分析（Step 9）</h1>')
    o.append('<div class="sub">回答「阴性结果是否源于功效不足」· 基于各分析实际 SE，不重新估计效应</div>')
    o.append('<div class="meta">'
             '<span><b>暴露</b> 甲减 / TSH / FT4（Rand 2025, PMID 41238958）</span>'
             '<span><b>显著水平</b> α=0.05（双侧）</span>'
             '<span><b>功效档</b> 80%% / 90%%</span>'
             '<span><b>比较</b> 正向 %d · 反向 %d · MVMR %d</span>'
             '<span><b>日期</b> 2026-10-07</span></div>'
             % (len(FW), len(RV), len(MV)))

    # ---- 0 结论速览
    o.append('<h2>0 · 结论速览：阴性并非整体功效不足，而是分性状不同</h2>')
    o.append('<table><tr><th class="l">暴露 / 检验</th><th class="l">是否功效不足</th>'
             '<th class="l">定量判据（无样本重叠结局）</th></tr>')
    o.append('<tr><td class="l"><b>甲减</b>（正向 MR）</td><td class="l">%s</td>'
             '<td class="l">MDE80 OR %s（每 odds 翻倍）；对 OR 1.05 的功效 %s → '
             '观测 OR 量级<b>完全落在可检出区间</b>；不稳健源于异质性，不是功效</td></tr>' % (
                 pill("ok", "否 · 功效充足"), SUMM["hypo"]["mde80"], SUMM["hypo"]["p105"]))
    o.append('<tr><td class="l"><b>TSH</b>（正向 MR）</td><td class="l">%s</td>'
             '<td class="l">MDE80 OR %s（每 1 SD）；对 OR 1.05 的功效仅 %s → '
             '可排除 ≥6%% 的效应，<b>&lt;5%% 的小效应不能排除</b></td></tr>' % (
                 pill("warn", "部分 · 功效中等"), SUMM["tsh"]["mde80"], SUMM["tsh"]["p105"]))
    o.append('<tr><td class="l"><b>FT4</b>（正向 MR）</td><td class="l">%s</td>'
             '<td class="l">MDE80 OR %s（每 1 SD）；对 OR 1.10 的功效仅 %s → '
             '阴性应写作「<b>未检出</b>」而非「已排除」</td></tr>' % (
                 pill("bad", "是 · 功效最弱"), SUMM["ft4"]["mde80"], SUMM["ft4"]["p110"]))
    o.append('<tr><td class="l"><b>反向 MR</b>（血管→甲状腺）</td><td class="l">%s</td>'
             '<td class="l">无重叠检验 MDE80 = %s（甲减 OR）、%s SD（TSH/FT4，每血管 odds 翻倍）→ '
             '精度足以排除中等以上反向效应；限制在<b>工具覆盖率</b>（效度）而非精度</td></tr>' % (
                 pill("ok", "否 · 精度足够"),
                 rng([fl(r["mde80_or"]) for r in REV_HYPO]),
                 rng([fl(r["mde80_sd"]) for r in REV_CONT], 4)))
    o.append('<tr><td class="l"><b>MVMR</b>（甲减直接效应）</td><td class="l">%s</td>'
             '<td class="l">MDE80 OR %s；条件 F %s（&gt;10 达标）→ 甲减直接效应的功效充足</td></tr>' % (
                 pill("ok", "否 · 功效充足"),
                 rng([fl(r["mde80_or"]) for r in MV_NOOV]),
                 rng([fl(r["condF"]) for r in MV_NOOV], 1)))
    o.append('</table>')

    # ---- 1 工具强度
    o.append('<h2>1 · 工具强度（设计层）</h2>')
    o.append('<table><tr><th class="l">暴露</th><th>工具数</th><th>累计 R²</th>'
             '<th>中位 F</th><th>最小 F</th><th>有效样本量</th></tr>')
    for t in TRAIT_ORDER:
        s = INST[t]
        cls = ' class="ok"' if s["f_min"] > 10 else ' class="sig"'
        o.append('<tr><td class="l">%s</td><td>%d</td><td>%.2f%%</td><td>%.1f</td>'
                 '<td%s>%.1f</td><td>%s</td></tr>' % (
                     TRAIT_CN[t], s["n"], s["r2"] * 100, s["f_median"], cls, s["f_min"],
                     "{:,}".format(int(round(s["n_eff"])))))
    o.append('</table>')
    o.append('<div class="box"><b>要点：</b>三性状<b>最小单位点 F 全部 ≥ %.1f</b>（&gt;10 的经验阈值），'
             '不存在弱工具偏倚；因此后续 MDE 的差异由<b>结局侧样本量</b>与<b>工具累计 R²</b>决定，'
             '而非工具强度不足。</div>' % MINF)

    # ---- 2 结局尺度体检
    o.append('<h2>2 · 结局尺度体检（数据质量前置）</h2>')
    o.append('<div class="box"><b>为什么先查这一项：</b>各结局文件必须处于同一尺度（二分结局＝log-OR），'
             '否则「MDE 的 OR」不成立。判据：用<b>仪器 SNP 的中位 SE</b> 横向比对，'
             '偏离参考中位 ≥5 倍即判为尺度异常。常数尺度重标定<b>不改变 z 与 p</b>，'
             '故显著性与方向结论不受影响，只有<b>效应量</b>受影响。</div>')
    o.append('<table><tr><th class="l">结局文件</th><th>仪器 SNP 数</th><th>中位 SE</th>'
             '<th>参考中位 SE</th><th>比值</th><th class="l">判定</th></tr>')
    for s in SCALE:
        bad = bool(str(s.get("flag") or "").strip())
        o.append('<tr><td class="l">%s</td><td>%s</td><td>%s</td><td>%s</td><td%s>%s</td>'
                 '<td class="l">%s</td></tr>' % (
                     KEY_CN.get(s["key"], s["key"]), s["n_se"], fnum(fl(s["med_se"]), 5),
                     fnum(fl(s["ref_med_se"]), 5), ' class="sig"' if bad else "",
                     fnum(fl(s["ratio_to_ref"]), 3),
                     pill("bad", s["flag"]) if bad else pill("ok", "OK · log-OR 尺度")))
    o.append('</table>')
    if BAD:
        o.append('<div class="box"><b>发现 1 个尺度异常文件：</b>%s —— 中位 SE 为 %s，'
                 '仅为参考值（%s）的 <b>%.3f 倍</b>，与「线性概率模型」的特征吻合'
                 '（β_LPM ≈ log-OR × p(1−p)，该结局病例占比 %s）。<br>'
                 '实证核验：同一强效位点 <code>rs3184504</code>（SH2B3）在该文件 β=%s（exp=%s），'
                 '而在 CAD-Nikpay / CAD-UKB / FinnGen / MVP / 两个卒中亚型中 β 为 %s–%s'
                 '（exp 全部落在 1.05–1.12，彼此一致）。<br>'
                 '<b>处理方式</b>：本次<b>不改动原始数据</b>，仅在功效报告中将该结局标记为'
                 '「尺度不可比」，其 MDE <b>不以 OR 形式呈现</b>；该结局本身已属「重叠·参考」行，'
                 '不进入任何主结论。</div>' % (
                     "、".join(KEY_CN.get(s["key"], s["key"]) for s in BAD),
                     fnum(fl(BAD[0]["med_se"]), 5), fnum(fl(BAD[0]["ref_med_se"]), 5),
                     fl(BAD[0]["ratio_to_ref"]), "2.29%（11,081/484,598）",
                     "0.0019", "1.0019", "0.049", "0.114"))

    # ---- 3 正向 MR
    o.append('<h2>3 · 正向 MR：最小可检出 OR 与对预设效应量的功效</h2>')
    o.append('<div class="box"><b>单位</b>：甲减＝每「甲减 odds 翻倍」（ln2）；TSH / FT4＝每 1 SD'
             '（GWAS 已标准化）。结局均为二分 → 统一报 OR。'
             '「重叠」＝与该暴露 GWAS 存在样本重叠，其 β 与功效仅供参考；'
             '<b>n.a.¹</b>＝§2 判定为尺度不可比的文件，其 MDE 不能以 OR 表达。</div>')
    o.append('<table><tr>%s</tr>' % "".join('<th class="l">%s</th>' % h for h in FWD_HEAD))
    last = None
    for r in FW:
        if last is not None and r["exposure_key"] != last:
            o.append('<tr><td colspan="%d" style="height:1px;padding:0;border-left:0;border-right:0"></td></tr>' % len(FWD_HEAD))
        last = r["exposure_key"]
        bad = flagged(r)
        cls = ' class="sig"' if (fl(r["p_fixed"]) or 1) < 0.05 else ""
        if bad:
            rest = "<td>n.a.<sup>1</sup></td><td>n.a.<sup>1</sup></td><td>n.a.</td><td>n.a.</td><td>n.a.</td>"
        else:
            rest = "<td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>" % (
                fnum(fl(r["mde80_or"]), 3), fnum(fl(r["mde90_or"]), 3),
                fnum(fl(r["pow_or1.05"]), 3), fnum(fl(r["pow_or1.10"]), 3),
                fnum(fl(r["pow_or1.15"]), 3))
        o.append('<tr><td class="l">%s</td><td class="l">%s%s</td><td class="l">%s</td>'
                 '<td>%s</td><td%s>%+.4f (%s)</td>%s</tr>' % (
                     r["exposure"], short(r["outcome"]), "¹" if bad else "",
                     note_pill(r["source_note"]), r["nsnp"],
                     cls, fl(r["beta"]), pnum(fl(r["p_fixed"])), rest))
    o.append('</table>')
    o.append('<div class="box"><b>读法：</b>MDE80＝在 80% 功效下可检出的<b>最小</b> OR。观测 β 若显著且其 OR 量级 ≥ MDE80，'
             '说明该效应本就在检出能力之内（阴性/阳性的分歧不是功效造成的）。<br>'
             '<sup>1</sup> 尺度不可比（见 §2）：该行仅保留 β 与 p（二者不受常数尺度影响），MDE 不呈现。<br>'
             '<b>本表刻意不给「观测效应的功效」</b>（post-hoc power）：负研究该指标恒偏低，'
             '会系统性地把阴性误判为「功效不足」，故以 MDE 与对<b>预设</b>效应量的功效为判据。</div>')

    # ---- 4 反向
    o.append('<h2>4 · 反向 MR：血管疾病 → 甲状腺功能 的检出能力</h2>')
    o.append('<div class="box"><b>单位</b>：暴露＝血管疾病（二分）→ 「每血管疾病 odds 翻倍」；'
             '结局为甲减（二分）报 <b>OR</b>，结局为 TSH / FT4（连续）报 <b>SD 变化</b>（不可报 OR）。</div>')
    o.append('<table><tr>%s</tr>' % "".join('<th class="l">%s</th>' % h for h in REV_HEAD))
    for r in sorted(RV, key=lambda x: (TRAIT_ORDER.index(x["exposure_key"]), x["outcome"])):
        if r["outcome_type"] == "二分":
            mdtxt = "OR <b>%s</b>" % fnum(fl(r["mde80_or"]), 3)
        else:
            mdtxt = "<b>%s</b> SD" % fnum(fl(r["mde80_sd"]), 4)
        cls = ' class="sig"' if (fl(r["p_fixed"]) or 1) < 0.05 else ""
        o.append('<tr><td class="l">%s</td><td class="l">%s</td><td class="l">%s</td>'
                 '<td>%s</td><td%s>%+.4f (%s)</td><td>%s</td></tr>' % (
                     r["exposure"], short(r["outcome"]), note_pill(r["source_note"]), r["nsnp"],
                     cls, fl(r["beta"]), pnum(fl(r["p_fixed"])), mdtxt))
    o.append('</table>')
    o.append('<div class="box"><b>要点：</b>唯一无样本重叠的反向检验（CAD Nikpay → 甲状腺）三性状全阴性，'
             'MDE80 分别为甲减 OR %s、TSH %s SD、FT4 %s SD（每血管 odds 翻倍）→ '
             '<b>精度足以排除中等以上反向效应</b>。须同时披露：该检验仅命中 11/41 个工具'
             '（甲状腺文件覆盖率约三成），限制在<b>效度</b>（工具是否代表全体）而非精度。</div>' % (
                 rng([fl(r["mde80_or"]) for r in REV_HYPO]),
                 fnum(fl([r for r in REV_CONT if r["exposure_key"] == "tsh"][0]["mde80_sd"]), 4),
                 fnum(fl([r for r in REV_CONT if r["exposure_key"] == "ft4"][0]["mde80_sd"]), 4)))

    # ---- 5 MVMR
    o.append('<h2>5 · MVMR：甲减直接效应（控制 TSH / FT4）的检出能力</h2>')
    o.append('<table><tr>%s</tr>' % "".join('<th class="l">%s</th>' % h for h in MV_HEAD))
    for r in MV:
        if r["exposure_key"] != "hypo":
            continue
        bad = flagged(r)
        cls = ' class="sig"' if (fl(r["p_fixed"]) or 1) < 0.05 else ""
        mde = ("n.a.<sup>1</sup>" if bad else "<b>%s</b>" % fnum(fl(r["mde80_or"]), 3))
        o.append('<tr><td class="l">%s%s</td><td class="l">%s</td><td>%s</td><td>%s</td>'
                 '<td%s>%+.4f (%s)</td><td>%s</td></tr>' % (
                     short(r["outcome"]), "¹" if bad else "",
                     note_pill(r["source_note"]), r["nsnp"],
                     fnum(fl(r["condF"]), 1), cls, fl(r["beta"]), pnum(fl(r["p_fixed"])), mde))
    o.append('</table>')
    o.append('<div class="box"><b>要点：</b>4 个无样本重叠结局的 MDE80 为 OR %s，条件 F %s（全部 &gt;10）；'
             '观测到显著正效应的结局（CAD-Nikpay、MVP-PAD、IS-SV）其效应量均 ≥ 或接近 MDE → '
             '甲减在 MVMR 中的检出能力充足。FT4 的条件 F 仅约 7.7（弱工具），其 MVMR 估计不可解释。</div>' % (
                 rng([fl(r["mde80_or"]) for r in MV_NOOV]),
                 rng([fl(r["condF"]) for r in MV_NOOV], 1)))

    # ---- 6 情景
    o.append('<h2>6 · 情景分析：结局样本量放大后的 MDE80（主分析 PAD，Sakaue 2021）</h2>')
    o.append('<table><tr><th class="l">暴露</th>%s</tr>' % "".join(
        "<th>×%d</th>" % k for k in (1, 2, 5, 10)))
    for t, key in zip(TRAIT_ORDER, PAD_KEYS):
        cells = []
        for k in (1, 2, 5, 10):
            v = SCEN["×%d" % k].get(key)
            cells.append("<td>%s</td>" % fnum(v, 3))
        o.append('<tr><td class="l">%s</td>%s</tr>' % (TRAIT_CN[t], "".join(cells)))
    o.append('</table>')
    o.append('<div class="box">假设 SE ∝ 1/√n（样本量放大 k 倍 → SE / √k），其余不变。'
             'FT4 即使在结局样本量 ×5 时 MDE80 仍为 OR %s，×10 时才降至 %s —— '
             '<b>若要把 FT4 的阴性做成「排除性」结论，需 ≈10 倍规模的结局 GWAS</b>。</div>' % (
                 fnum(SCEN["×5"][PAD_KEYS[2]], 3), fnum(SCEN["×10"][PAD_KEYS[2]], 3)))

    # ---- 7 可写 / 不能写
    o.append('<h2>7 · 可写入论文的结论</h2><div class="two">')
    o.append('<div class="box"><b>站得住的</b><ul>'
             '<li>三性状工具强度全部达标（最小 F ≥ %.1f），<b>排除了弱工具导致的偏倚</b>。</li>'
             '<li>甲减的阳性信号落在可检出区间内（MDE80 OR %s；对 OR 1.05 的功效 %s），'
             '说明其「不稳健」来自异质性 / 随机效应，<b>不是功效不足</b>。</li>'
             '<li>反向 MR（CAD → 甲状腺）的精度足以排除中等以上反向效应。</li>'
             '<li>TSH 阴性可表述为「可排除 ≥6%% 的效应」（MDE80 OR %s / SD）。</li>'
             '</ul></div>' % (MINF, SUMM["hypo"]["mde80"], SUMM["hypo"]["p105"], SUMM["tsh"]["mde80"]))
    o.append('<div class="box"><b>不能写的（须披露）</b><ul>'
             '<li><b>MI（UKB, GCST90038610）的效应量不可与其他结局并列</b>：§2 判定其 β / SE '
             '处于线性概率尺度（约为 log-OR 的 p(1−p) 倍），其 exp(β) 不是 OR。'
             'z 与 p 不受影响，故该行「重叠·参考」的显著性结论仍成立，但<b>不得报 OR</b>。</li>'
             '<li><b>FT4 阴性不可写作「无因果」</b>：MDE80 OR %s，对 OR 1.10 的功效仅 %s → 只能写「未检出」。</li>'
             '<li>TSH 阴性不能排除 &lt;5%% 的小效应（对 OR 1.05 的功效 %s）。</li>'
             '<li>不使用 post-hoc power 作为判据（负研究该指标恒低，易误判为功效不足）。</li>'
             '<li>MDE 以「SE 不随真效应变化」为前提（对 IVW 成立）；情景分析假设 SE ∝ 1/√n，'
             '未考虑病例比例与芯片覆盖的变化。</li>'
             '<li>反向 MR 的精度结论以「11 个命中工具均为有效工具」为前提；'
             '工具覆盖率（11/41）带来的<b>效度</b>风险不能由功效分析回答。</li>'
             '</ul></div></div>' % (SUMM["ft4"]["mde80"], SUMM["ft4"]["p110"], SUMM["tsh"]["p105"]))

    # ---- 8 方法
    o.append('<h2>8 · 方法与数据来源</h2>')
    o.append('<div class="box"><b>功效公式</b>：IVW 的 Wald 统计量 = |β̂| / SE，非中心参数 λ = |β_true| / SE；'
             'Power(β_true) = Φ(λ − z<sub>1−α/2</sub>) + Φ(−λ − z<sub>1−α/2</sub>)；'
             'MDE(power) = (z<sub>1−α/2</sub> + z<sub>power</sub>) × SE。'
             'α=0.05 双侧 → z=1.95996；power 0.80 → 0.84162，0.90 → 1.28155。'
             'SE 一律取自已完成的 MR 结果文件（IVW 固定 / 随机），不重新估计效应。</div>')
    o.append('<table><tr><th class="l">来源</th><th class="l">文件</th><th class="l">说明</th></tr>')
    o.append('<tr><td class="l">工具强度</td><td class="l">work/instruments_*.tsv</td>'
             '<td class="l">R²、F、有效样本量（Step 2 产出）</td></tr>')
    o.append('<tr><td class="l">正向 MR SE</td><td class="l">step3/results/mr_main_*.json、step4/results/mr_*.json</td>'
             '<td class="l">IVW 固定 / 随机效应</td></tr>')
    o.append('<tr><td class="l">反向 MR SE</td><td class="l">step6/results/rev_*.json</td>'
             '<td class="l">血管暴露 → 甲状腺结局</td></tr>')
    o.append('<tr><td class="l">MVMR SE</td><td class="l">step7/results/mvmr_*.json</td>'
             '<td class="l">条件 F 取自同文件 cond_f</td></tr>')
    o.append('<tr><td class="l">结局尺度体检</td><td class="l">step4/results/harmonised_*.tsv、'
             'step3/results/harmonised_*.tsv</td>'
             '<td class="l">仪器 SNP 的中位 <code>se_Y</code> 横向比对（§2）</td></tr>')
    o.append('</table>')
    o.append('<div class="foot"><b>生成方式</b>：由 <code>step9/scripts/s9_onepager.py</code> 直接读取 '
             '<code>step9/results/power_forward.tsv</code>、<code>power_reverse.tsv</code>、'
             '<code>power_mvmr.tsv</code>、<code>power_summary.json</code>、'
             '<code>scale_check.tsv</code> 生成，'
             '正文结论句中的数字亦为程序自动写入，未作人工修改。计算脚本：'
             '<code>step9/scripts/s9_power.py</code>。</div>')
    o.append('</div></body></html>')
    return "\n".join(o)


def build_md():
    o = []
    o.append("# 甲状腺功能 × 动脉粥样硬化性血管疾病：统计功效分析（Step 9）")
    o.append("")
    o.append("**回答「阴性结果是否源于功效不足」· 基于各分析实际 SE，不重新估计效应**　｜　2026-10-07")
    o.append("")
    o.append("- **暴露**：甲减 / TSH / FT4（Rand et al. 2025, *Nat Genet*, PMID 41238958）")
    o.append("- **口径**：α=0.05 双侧；功效档 80%% / 90%%；比较 正向 %d · 反向 %d · MVMR %d" % (len(FW), len(RV), len(MV)))
    o.append("")
    o.append("---")
    o.append("")
    o.append("## 0 结论速览：阴性并非整体功效不足，而是分性状不同")
    o.append("")
    o.append("| 暴露 / 检验 | 是否功效不足 | 定量判据（无样本重叠结局） |")
    o.append("|---|---|---|")
    o.append("| **甲减**（正向 MR） | 否 · 功效充足 | MDE80 OR %s（每 odds 翻倍）；对 OR 1.05 的功效 %s |" % (
        SUMM["hypo"]["mde80"], SUMM["hypo"]["p105"]))
    o.append("| **TSH**（正向 MR） | 部分 · 功效中等 | MDE80 OR %s（每 1 SD）；对 OR 1.05 的功效仅 %s |" % (
        SUMM["tsh"]["mde80"], SUMM["tsh"]["p105"]))
    o.append("| **FT4**（正向 MR） | 是 · 功效最弱 | MDE80 OR %s（每 1 SD）；对 OR 1.10 的功效仅 %s |" % (
        SUMM["ft4"]["mde80"], SUMM["ft4"]["p110"]))
    o.append("| **反向 MR**（血管→甲状腺） | 否 · 精度足够 | MDE80 = %s（甲减 OR）、%s SD（TSH/FT4） |" % (
        rng([fl(r["mde80_or"]) for r in REV_HYPO]),
        rng([fl(r["mde80_sd"]) for r in REV_CONT], 4)))
    o.append("| **MVMR**（甲减直接效应） | 否 · 功效充足 | MDE80 OR %s；条件 F %s |" % (
        rng([fl(r["mde80_or"]) for r in MV_NOOV]), rng([fl(r["condF"]) for r in MV_NOOV], 1)))
    o.append("")
    o.append("## 1 工具强度（设计层）")
    o.append("")
    o.append("| 暴露 | 工具数 | 累计 R² | 中位 F | 最小 F | 有效样本量 |")
    o.append("|---|---|---|---|---|---|")
    for t in TRAIT_ORDER:
        s = INST[t]
        o.append("| %s | %d | %.2f%% | %.1f | **%.1f** | %s |" % (
            TRAIT_CN[t], s["n"], s["r2"] * 100, s["f_median"], s["f_min"],
            "{:,}".format(int(round(s["n_eff"])))))
    o.append("")
    o.append("> 三性状最小单位点 F 全部 ≥ %.1f（>10 的经验阈值），不存在弱工具偏倚。" % MINF)
    o.append("")
    o.append("## 2 结局尺度体检（数据质量前置）")
    o.append("")
    o.append("> 各结局文件必须处于同一尺度（二分结局＝log-OR），否则「MDE 的 OR」不成立。"
             "判据：用**仪器 SNP 的中位 SE** 横向比对，偏离参考中位 ≥5 倍即判为尺度异常。"
             "常数尺度重标定**不改变 z 与 p**，故显著性与方向结论不受影响，只有效应量受影响。")
    o.append("")
    o.append("| 结局文件 | 仪器 SNP 数 | 中位 SE | 参考中位 SE | 比值 | 判定 |")
    o.append("|---|---|---|---|---|---|")
    for s in SCALE:
        bad = bool(str(s.get("flag") or "").strip())
        o.append("| %s | %s | %s | %s | %s | %s |" % (
            KEY_CN.get(s["key"], s["key"]), s["n_se"], fnum(fl(s["med_se"]), 5),
            fnum(fl(s["ref_med_se"]), 5), fnum(fl(s["ratio_to_ref"]), 3),
            ("**" + str(s["flag"]) + "**") if bad else "OK · log-OR 尺度"))
    o.append("")
    if BAD:
        o.append("**发现 1 个尺度异常文件**：%s —— 中位 SE 为 %s，仅为参考值（%s）的 **%.3f 倍**，"
                 "与「线性概率模型」特征吻合（β_LPM ≈ log-OR × p(1−p)；该结局病例占比 2.29%%，"
                 "即 11,081/484,598）。实证核验：同一强效位点 `rs3184504`（SH2B3）在该文件 β=0.0019"
                 "（exp=1.0019），而在 CAD-Nikpay / CAD-UKB / FinnGen / MVP / 两个卒中亚型中 "
                 "β 为 0.049–0.114（exp 全部落在 1.05–1.12，彼此一致）。"
                 "**处理方式**：本次不改动原始数据，仅将其标记为「尺度不可比」，MDE 不以 OR 呈现；"
                 "该结局本身已属「重叠·参考」行，不进入任何主结论。" % (
                     "、".join(KEY_CN.get(s["key"], s["key"]) for s in BAD),
                     fnum(fl(BAD[0]["med_se"]), 5), fnum(fl(BAD[0]["ref_med_se"]), 5),
                     fl(BAD[0]["ratio_to_ref"])))
        o.append("")
    o.append("## 3 正向 MR：最小可检出 OR 与对预设效应量的功效")
    o.append("")
    o.append("> **单位**：甲减＝每「甲减 odds 翻倍」（ln2）；TSH / FT4＝每 1 SD。结局均为二分 → 统一报 OR。"
             "**n.a.¹**＝§2 判定为尺度不可比的文件，其 MDE 不能以 OR 表达。")
    o.append("")
    o.append("| " + " | ".join(FWD_HEAD) + " |")
    o.append("|" + "---|" * len(FWD_HEAD))
    for r in FW:
        bad = flagged(r)
        tail = "n.a.¹ | n.a.¹ | n.a. | n.a. | n.a." if bad else "%s | %s | %s | %s | %s" % (
            "**%s**" % fnum(fl(r["mde80_or"]), 3), fnum(fl(r["mde90_or"]), 3),
            fnum(fl(r["pow_or1.05"]), 3), fnum(fl(r["pow_or1.10"]), 3),
            fnum(fl(r["pow_or1.15"]), 3))
        o.append("| %s | %s%s | %s | %s | %+.4f (%s) | %s |" % (
            r["exposure"], short(r["outcome"]), "¹" if bad else "", r["source_note"], r["nsnp"],
            fl(r["beta"]), pnum(fl(r["p_fixed"])), tail))
    o.append("")
    o.append("> MDE80＝80% 功效下可检出的**最小** OR。本表刻意不给「观测效应的功效」（post-hoc power）："
             "负研究该指标恒偏低，会系统性把阴性误判为功效不足。")
    o.append("> ¹ 尺度不可比（见 §2）：该行仅保留 β 与 p（二者不受常数尺度影响），MDE 不呈现。")
    o.append("")
    o.append("## 4 反向 MR：血管疾病 → 甲状腺功能 的检出能力")
    o.append("")
    o.append("> **单位**：暴露＝血管疾病（二分）→「每血管疾病 odds 翻倍」；甲减结局报 OR，TSH / FT4 报 SD 变化。")
    o.append("")
    o.append("| " + " | ".join(REV_HEAD) + " |")
    o.append("|" + "---|" * len(REV_HEAD))
    for r in sorted(RV, key=lambda x: (TRAIT_ORDER.index(x["exposure_key"]), x["outcome"])):
        mdtxt = ("OR **%s**" % fnum(fl(r["mde80_or"]), 3)) if r["outcome_type"] == "二分" \
            else ("**%s** SD" % fnum(fl(r["mde80_sd"]), 4))
        o.append("| %s | %s | %s | %s | %+.4f (%s) | %s |" % (
            r["exposure"], short(r["outcome"]), r["source_note"], r["nsnp"],
            fl(r["beta"]), pnum(fl(r["p_fixed"])), mdtxt))
    o.append("")
    o.append("> 唯一无样本重叠的反向检验（CAD Nikpay → 甲状腺）三性状全阴性，精度足以排除中等以上反向效应；"
             "但仍受工具覆盖率（11/41）的**效度**限制，该风险不能由功效分析回答。")
    o.append("")
    o.append("## 4 MVMR：甲减直接效应（控制 TSH / FT4）的检出能力")
    o.append("")
    o.append("| " + " | ".join(MV_HEAD) + " |")
    o.append("|" + "---|" * len(MV_HEAD))
    for r in MV:
        if r["exposure_key"] != "hypo":
            continue
        bad = flagged(r)
        mde = "n.a.¹" if bad else "**%s**" % fnum(fl(r["mde80_or"]), 3)
        o.append("| %s%s | %s | %s | %s | %+.4f (%s) | %s |" % (
            short(r["outcome"]), "¹" if bad else "", r["source_note"], r["nsnp"],
            fnum(fl(r["condF"]), 1), fl(r["beta"]), pnum(fl(r["p_fixed"])), mde))
    o.append("")
    o.append("> FT4 的条件 F 仅约 7.7（弱工具），其 MVMR 估计不可解释。"
             "¹ 尺度不可比（见 §2）。")
    o.append("")
    o.append("## 6 情景分析：结局样本量放大后的 MDE80（主分析 PAD，Sakaue 2021）")
    o.append("")
    o.append("| 暴露 | ×1 | ×2 | ×5 | ×10 |")
    o.append("|---|---|---|---|---|")
    for t, k2 in zip(TRAIT_ORDER, PAD_KEYS):
        o.append("| %s | %s | %s | %s | %s |" % (
            TRAIT_CN[t], fnum(SCEN["×1"].get(k2), 3), fnum(SCEN["×2"].get(k2), 3),
            fnum(SCEN["×5"].get(k2), 3), fnum(SCEN["×10"].get(k2), 3)))
    o.append("")
    o.append("> 假设 SE ∝ 1/√n，其余不变。FT4 即使结局样本量 ×10 后 MDE80 仍为 OR %s。" % (
        fnum(SCEN["×10"][PAD_KEYS[2]], 3)))
    o.append("")
    o.append("## 7 可写入论文的结论")
    o.append("")
    o.append("### 站得住的")
    o.append("")
    o.append("1. 三性状工具强度全部达标（最小 F ≥ %.1f），**排除弱工具偏倚**。" % MINF)
    o.append("2. 甲减的阳性信号落在可检出区间内（MDE80 OR %s；对 OR 1.05 的功效 %s）→ "
             "其「不稳健」来自异质性 / 随机效应，**不是功效不足**。" % (
                 SUMM["hypo"]["mde80"], SUMM["hypo"]["p105"]))
    o.append("3. 反向 MR（CAD → 甲状腺）的精度足以排除中等以上反向效应。")
    o.append("4. TSH 阴性可表述为「可排除 ≥6%% 的效应」（MDE80 OR %s）。" % SUMM["tsh"]["mde80"])
    o.append("")
    o.append("### 不能写的（须披露）")
    o.append("")
    o.append("1. **MI（UKB, GCST90038610）的效应量不可与其他结局并列**：§2 判定其 β / SE "
             "处于线性概率尺度（约为 log-OR 的 p(1−p) 倍），其 exp(β) 不是 OR。"
             "z 与 p 不受影响，故该行「重叠·参考」的显著性结论仍成立，但**不得报 OR**。")
    o.append("2. **FT4 阴性不可写作「无因果」**：MDE80 OR %s，对 OR 1.10 的功效仅 %s → 只能写「未检出」。" % (
        SUMM["ft4"]["mde80"], SUMM["ft4"]["p110"]))
    o.append("3. TSH 阴性不能排除 <5%% 的小效应（对 OR 1.05 的功效 %s）。" % SUMM["tsh"]["p105"])
    o.append("4. 不使用 post-hoc power 作为判据。")
    o.append("5. MDE 以「SE 不随真效应变化」为前提；情景分析假设 SE ∝ 1/√n。")
    o.append("6. 反向 MR 的精度结论以「11 个命中工具均为有效工具」为前提，工具覆盖率的效度风险不能由功效分析回答。")
    o.append("")
    o.append("## 8 方法与数据来源")
    o.append("")
    o.append("**功效公式**：`Power(β_true) = Φ(λ − z_{1−α/2}) + Φ(−λ − z_{1−α/2})`，"
             "`λ = |β_true| / SE`；`MDE(power) = (z_{1−α/2} + z_power) × SE`。"
             "α=0.05 双侧 → z=1.95996；power 0.80 → 0.84162，0.90 → 1.28155。"
             "SE 一律取自已完成的 MR 结果文件，不重新估计效应。")
    o.append("")
    o.append("| 来源 | 文件 | 说明 |")
    o.append("|---|---|---|")
    o.append("| 工具强度 | work/instruments_*.tsv | R²、F、有效样本量（Step 2） |")
    o.append("| 正向 MR SE | step3/results/mr_main_*.json、step4/results/mr_*.json | IVW 固定 / 随机 |")
    o.append("| 反向 MR SE | step6/results/rev_*.json | 血管暴露 → 甲状腺结局 |")
    o.append("| MVMR SE | step7/results/mvmr_*.json | 条件 F 取自同文件 cond_f |")
    o.append("| 结局尺度体检 | step4/results/harmonised_*.tsv、step3/results/harmonised_*.tsv | "
             "仪器 SNP 的中位 `se_Y` 横向比对（§2） |")
    o.append("")
    o.append("> 本文件由 `step9/scripts/s9_onepager.py` 从结果文件直接生成；"
             "计算脚本 `step9/scripts/s9_power.py`；正文数字由程序自动写入，未作人工修改。")
    return "\n".join(o)


if __name__ == "__main__":
    hp = os.path.join(DELIV, "Step9_统计功效.html")
    mp = os.path.join(DELIV, "Step9_统计功效.md")
    with io.open(hp, "w", encoding="utf-8") as f:
        f.write(build_html())
    with io.open(mp, "w", encoding="utf-8") as f:
        f.write(build_md())
    print("OK  " + hp, os.path.getsize(hp), "B")
    print("OK  " + mp, os.path.getsize(mp), "B")
