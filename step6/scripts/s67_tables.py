# -*- coding: utf-8 -*-
"""
Step 6（反向 MR）＋ Step 7（MVMR）结果数据表导出。
数值全部来自 step6/results/*.json 与 step7/results/*.json，不手抄。
输出：deliverables/Step6-7_反向MR与MVMR_数据表.xlsx
"""
import io
import json
import os

from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side

ROOT = r"D:\WorkSpace_Study\pad_thyroid_mr"
DELIV = os.path.join(ROOT, "deliverables")
TRAITS = [("hypo", "甲减"), ("tsh", "TSH"), ("ft4", "FT4")]

REV_EXP = [
    ("cad_nikpay", "CAD (Nikpay 2015; GCST003116)", "无重叠"),
    ("pad_sakaue", "PAD (Sakaue 2021; GCST90018890)", "含 UKB"),
    ("pad_finngen", "PAD (FinnGen R13 I9_PAD)", "含 FinnGen DF10"),
]
MVM_EXP = [
    ("cad_nouk", "CAD (Nikpay 2015; GCST003116)", "无"),
    ("cad_ukb", "CAD (Mbatchou 2021 UKB; GCST90013864)", "UKB"),
    ("mi_ukb", "MI (Donertas 2021 UKB; GCST90038610)", "UKB"),
    ("is_laa", "IS-LAA (MEGASTROKE; GCST005840)", "无"),
    ("is_sv", "IS-SV (MEGASTROKE; GCST005841)", "无"),
    ("sakaue", "PAD (Sakaue 2021; GCST90018890)", "UKB"),
    ("mvp_meta", "PAD (MVP trans-ancestry meta)", "无"),
]

HDR_FILL = PatternFill("solid", fgColor="1F4E79")
HDR_FONT = Font(bold=True, color="FFFFFF", size=10)
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
BAD_FILL = PatternFill("solid", fgColor="FDEBEB")

# ---- ‡ 尺度注记（Step 9 结局尺度体检发现）--------------------------------------
# MI（UKB, GCST90038610）文件的 β / SE 为「线性概率」尺度而非 log-OR。
# 只在 MI 行加注记，**不改动任何数值**。依据：step9/results/scale_check.tsv。
SCALE_NOTE = ("‡ 尺度注记：MI（UKB, GCST90038610）的 β/SE 处于线性概率尺度而非 log-OR"
              "（仪器 SNP 中位 SE 仅为其他结局的 0.030 倍；同一位点 rs3184504 在该文件 β=0.0019、"
              "其余结局 0.049–0.114）。故本行 β 不可与其他结局并列解读，其 β 偏小 / p 极小系尺度所致，"
              "并非「重叠使 SE 被低估」，也不代表效应可忽略；z 与 p 不受常数尺度影响，"
              "「重叠·仅参考、不进入因果结论」的判定不变。核查：step9/results/scale_check.tsv")
NOTES = {"mi_ukb": SCALE_NOTE}


def J(p):
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def style_header(ws, ncol, row=1):
    for c in range(1, ncol + 1):
        cell = ws.cell(row=row, column=c)
        cell.fill = HDR_FILL
        cell.font = HDR_FONT
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BORDER
    ws.freeze_panes = ws.cell(row=row + 1, column=1)


def autowidth(ws):
    for col in ws.columns:
        w = 10
        for cell in col:
            v = cell.value
            if v is None:
                continue
            ln = sum(2 if ord(ch) > 127 else 1 for ch in str(v))
            w = max(w, min(ln + 2, 60))
        ws.column_dimensions[col[0].column_letter].width = w


wb = Workbook()

# ---------------- 表1：反向 MR 结果
ws = wb.active
ws.title = "反向MR结果"
cols = ["血管暴露（X）", "重叠状态", "甲状腺性状（Y）", "nSNP", "β", "SE", "p(IVW固定)",
        "p(IVW随机)", "Egger截距p", "Q p", "OR[95%CI]", "符号 +/-", "LOO β区间", "判定"]
ws.append(cols)
for exp, name, ov in REV_EXP:
    for t, cn in TRAITS:
        p = os.path.join(ROOT, "step6", "results", "rev_%s_%s.json" % (exp, t))
        if not os.path.exists(p):
            ws.append([name, ov, cn, None, None, None, None, None, None, None, None, None, None, "未生成"])
            continue
        d = J(p)
        eg = d.get("mr_egger") or {}
        sg = d.get("sign_test") or {}
        # LOO
        lp = os.path.join(ROOT, "step6", "results", "revloo_%s_%s.tsv" % (exp, t))
        ltxt = ""
        if os.path.exists(lp):
            rows = [l.split("\t") for l in io.open(lp, encoding="utf-8").read().strip().split("\n")[1:] if l.strip()]
            if rows:
                b = [float(r[2]) for r in rows]
                ltxt = "%.4f ~ %.4f" % (min(b), max(b))
        usable = (ov == "无重叠")
        if d["nsnp"] <= 1:
            verdict = "不可用(单工具/重叠)"
        elif not usable:
            verdict = "重叠·仅参考"
        elif d["ivw_fixed"]["p"] < 0.05:
            verdict = "名义阳性"
        else:
            verdict = "阴性"
        ws.append([name, ov, cn, d["nsnp"], d["ivw_fixed"]["beta"], d["ivw_fixed"]["se"],
                   d["ivw_fixed"]["p"], d["ivw_random"]["p"], eg.get("intercept_pval"),
                   d["ivw_random"]["Q_pval"],
                   "%.3f [%.3f, %.3f]" % (2.718281828459045 ** d["ivw_fixed"]["beta"],
                                          2.718281828459045 ** (d["ivw_fixed"]["beta"] - 1.959963984540054 * d["ivw_fixed"]["se"]),
                                          2.718281828459045 ** (d["ivw_fixed"]["beta"] + 1.959963984540054 * d["ivw_fixed"]["se"])),
                   "%d / %d" % (sg.get("n_positive") or 0, sg.get("n_negative") or 0),
                   ltxt, verdict])
style_header(ws, len(cols))
autowidth(ws)

# ---------------- 表2：反向 MR 工具命中率
ws2 = wb.create_sheet("反向MR工具命中率")
cols2 = ["血管暴露", "工具总数", "甲状腺文件命中", "命中率(%)", "命中路径", "可用工具(nSNP)"]
ws2.append(cols2)
for exp, name, ov in REV_EXP:
    p = os.path.join(ROOT, "step6", "results", "rev_%s_hypo.json" % exp)
    if not os.path.exists(p):
        continue
    d = J(p)
    cov = d.get("coverage") or {}
    br = cov.get("by_route") or {}
    ws2.append([name, cov.get("n_tools"), cov.get("n_outcome_hit"), cov.get("hit_pct"),
                "、".join("%s %d" % (k, v) for k, v in sorted(br.items())), d["nsnp"]])
style_header(ws2, len(cols2))
autowidth(ws2)

# ---------------- 表3：MVMR 结果
ws3 = wb.create_sheet("MVMR结果")
cols3 = ["结局（Y）", "重叠", "nSNP",
         "甲减 MV β", "甲减 MV p", "甲减 单变量 β", "甲减 单变量 p",
         "TSH MV β", "TSH MV p", "TSH 单变量 β", "TSH 单变量 p",
         "FT4 MV β", "FT4 MV p", "FT4 单变量 β", "FT4 单变量 p",
         "条件F 甲减", "条件F TSH", "条件F FT4", "MV-Egger截距 p", "MV-Q", "MV-Q p", "弱工具检查",
         "尺度注记"]
ws3.append(cols3)
for key, name, ov in MVM_EXP:
    p = os.path.join(ROOT, "step7", "results", "mvmr_%s.json" % key)
    if not os.path.exists(p):
        ws3.append([name, ov, None] + [None] * 18 + ["未生成", ""])
        continue
    d = J(p)
    mv = d["mvmr"]; un = d.get("univariable") or {}
    cf = {x["exposure"]: x["F"] for x in d["cond_f"]}
    row = [name + ("‡" if key in NOTES else ""), ov, d["nsnp"]]
    for t in ("hypo", "tsh", "ft4"):
        u = un.get(t) or {}
        row += [mv[t]["beta"], mv[t]["p"], u.get("beta"), u.get("p")]
    row += [cf.get("hypo"), cf.get("tsh"), cf.get("ft4"),
            (d.get("mv_egger") or {}).get("intercept_pval"),
            d.get("Q"), d.get("Q_pval"), d.get("weak_check"), NOTES.get(key, "")]
    ws3.append(row)
    if key in NOTES:
        for c in range(1, len(cols3) + 1):
            ws3.cell(row=ws3.max_row, column=c).fill = BAD_FILL
style_header(ws3, len(cols3))
autowidth(ws3)

if not os.path.exists(DELIV):
    os.makedirs(DELIV)
out = os.path.join(DELIV, "Step6-7_反向MR与MVMR_数据表.xlsx")
wb.save(out)
print("OK  " + out, os.path.getsize(out), "B")
