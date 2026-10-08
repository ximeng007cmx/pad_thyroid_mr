# -*- coding: utf-8 -*-
"""
Step 3.2-3.8  MR 主分析（甲减/TSH/FT4 → PAD）

输入  step3/results/harmonised_{hypo,tsh,ft4}.tsv   (Step 3.1 产出)
输出  step3/results/mr_main_<trait>.json            每性状全部统计量
      step3/results/loo_<trait>.tsv                 留一法
      step3/results/mr_main_summary.tsv             三性状 × 方法汇总表
      step3/logs/s3_mr_main_<时间戳>.log

方法（纯 Python 实现，与 R TwoSampleMR/MendelianRandomization 口径对齐）
  3.2 Wald ratio 逐 SNP + 符号一致性二项检验
  3.3 IVW（固定 + DL 随机效应）
  3.4 MR-Egger（加权回归，含截距）、加权中位数、加权众数（高斯核 KDE）
  3.5 Cochran's Q（IVW 与 Egger）、Egger 截距检验
  3.6 MR-PRESSO 等价：全局异质性（=Q）+ 基于 Bonferroni 校正的残差离群检出
      + distortion 检验（剔除离群前后 IVW 估计变化）
  3.7 留一法（随机效应 IVW，逐 SNP 剔除）
  3.8 方向性检验：比较暴露/结局方差解释量（R² 与 n·R²；二值结局用 n_eff 近似）
      —— 注：为 Steiger 方向性检验的透明近似实现（比较两侧 F/R²），
         非精确 Steiger Z 公式，报告中如实标注。

结局 GCST90018890（Sakaue 2021, PMID 34594039）样本量（GWAS Catalog V2 API 2026-10-07 取回）：
  EUR 7,114 cases / 475,964 controls（n_eff_EUR = 4·n1·n0/(n1+n0)）
"""
import os, csv, json, math, time, datetime, random

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
S3   = os.path.join(ROOT, "step3")
RES  = os.path.join(S3, "results")
LOG  = os.path.join(S3, "logs")
os.makedirs(RES, exist_ok=True)
STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOGF  = os.path.join(LOG, "s3_mr_main_%s.log" % STAMP)

TRAITS = ["hypo", "tsh", "ft4"]
TRAIT_CN = {"hypo": "甲减", "tsh": "TSH", "ft4": "FT4"}
# 结局 PAD：EUR 病例/对照（用于 Steiger 方向性近似 + 二值 n_eff）
OUT_NCASE_EUR, OUT_NCTRL_EUR = 7114, 475964
OUT_N_EFF_EUR = 4.0 * OUT_NCASE_EUR * OUT_NCTRL_EUR / (OUT_NCASE_EUR + OUT_NCTRL_EUR)

def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line)
    with open(LOGF, "a", encoding="utf-8") as f:
        f.write(line + "\n")

# ---------------- 统计基础 ----------------
def p_norm(z):                      # 双侧正态 p
    return math.erfc(abs(z) / math.sqrt(2.0))

def _gser(a, x):                    # P(a,x) 级数（x 小）
    ap, s, d = a, 1.0 / a, 1.0 / a
    for _ in range(300):
        ap += 1.0
        d *= x / ap
        s += d
        if abs(d) < abs(s) * 1e-12:
            break
    return s * math.exp(-x + a * math.log(x) - math.lgamma(a)) if x > 0 else 1.0

def _gcf(a, x):                     # Q(a,x) 连分式（x 大）
    b, c, d = x + 1.0 - a, 1e300, 1.0 / (x + 1.0 - a)
    h = d
    for i in range(1, 300):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300: d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300: c = 1e-300
        d = 1.0 / d
        delt = d * c
        h *= delt
        if abs(delt - 1.0) < 1e-12:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h

def p_chi2(x, df):                  # 上侧 p
    if x <= 0: return 1.0
    a, xx = df / 2.0, x / 2.0
    if xx < a + 1.0:
        return 1.0 - _gser(a, xx) if xx > 0 else 1.0
    return _gcf(a, xx)

def p_binom_two(k, n, p=0.5):       # 精确二项双侧 p
    from math import comb
    tail = sum(comb(n, i) * p**i * (1-p)**(n-i) for i in range(0, min(k, n-k)+1))
    return min(1.0, 2.0 * tail)

# ---------------- IVW / Egger ----------------
def ivw(b, s, weights=None):
    if weights is None: weights = [1.0 / si**2 for si in s]
    W = sum(weights)
    beta = sum(wi * bi for wi, bi in zip(weights, b)) / W
    se = math.sqrt(1.0 / W)
    return beta, se

def q_stat(b, s, beta=None):
    k = len(b)
    if beta is None:
        beta, _ = ivw(b, s)
    Q = sum((bi - beta)**2 / si**2 for bi, si in zip(b, s))
    return Q, max(k - 1, 0)

def ivw_random_dl(b, s):
    """DerSimonian-Laird 随机效应 IVW"""
    bF, seF = ivw(b, s)
    Q, df = q_stat(b, s, bF)
    w = [1.0 / si**2 for si in s]
    Sw, Sww = sum(w), sum(wi**2 for wi in w)
    tau2 = max(0.0, (Q - df) / (Sw - Sww / Sw)) if Sw > Sww / Sw else 0.0
    wr = [1.0 / (si**2 + tau2) for si in s]
    beta, se = ivw(b, s, wr)
    p = p_norm(beta / se)
    return {"beta": beta, "se": se, "p": p, "Q": Q, "Q_df": df,
            "Q_pval": p_chi2(Q, df), "tau2": tau2, "nsnp": len(b)}

def egger(bX, bY, sY):
    """MR-Egger（Bowden 2016 规范形式）：b_Y ~ b_X，权重 1/se_Y^2
       斜率 B1 = 因果估计；截距 B0 = 水平多效性。
       se 按残差方差缩放（sigma^2 = Qe/(n-2)，与 R lm(weighted) 一致，偏保守）"""
    n = len(bX)
    w = [1.0 / si**2 for si in sY]
    W = sum(w)
    Sx  = sum(wi*x for wi, x in zip(w, bX))
    Sy  = sum(wi*y for wi, y in zip(w, bY))
    Sxx = sum(wi*x*x for wi, x in zip(w, bX))
    Sxy = sum(wi*x*y for wi, x, y in zip(w, bX, bY))
    D = Sxx - Sx*Sx/W
    if abs(D) < 1e-300: return None
    B1 = (Sxy - Sx*Sy/W) / D            # 斜率（因果）
    B0 = Sy/W - B1*Sx/W                 # 截距（多效性）
    fit = [B0 + B1*x for x in bX]
    Qe = sum(wi*(y-f)**2 for wi, y, f in zip(w, bY, fit))
    dfe = n - 2
    s2 = Qe/dfe if dfe > 0 else 1.0
    xbar = Sx/W
    se1 = math.sqrt(s2 / D)
    se0 = math.sqrt(s2 * (1.0/W + xbar*xbar/D))
    return {"beta": B1, "se": se1, "p": p_norm(B1/se1),
            "intercept": B0, "intercept_se": se0, "intercept_pval": p_norm(B0/se0),
            "Q": Qe, "Q_df": dfe, "Q_pval": p_chi2(Qe, dfe), "nsnp": n}

def weighted_median(b, s):
    idx = sorted(range(len(b)), key=lambda i: b[i])
    w = [1.0 / si**2 for si in s]
    W = sum(w); cum = 0.0
    for i in idx:
        cum += w[i]
        if cum >= 0.5 * W:
            return b[i]
    return b[idx[-1]]

def weighted_mode(b, s):
    """加权众数：加权高斯核 KDE 峰值（对齐 MRLegacy density('nrd0') 思路）"""
    n = len(b)
    mean = sum(b)/n
    sd = math.sqrt(sum((x-mean)**2 for x in b)/(n-1)) if n > 1 else 1.0
    bs = sorted(b)
    q1, q3 = bs[int(0.25*(n-1))], bs[int(0.75*(n-1))]
    iqr = max(q3 - q1, 1e-12)
    bw = 0.9 * min(sd, iqr/1.34) * n**(-0.2)
    bw = max(bw, 1e-10)
    w = [1.0/si**2 for si in s]
    W = sum(w)
    lo, hi = min(b) - 2*bw, max(b) + 2*bw
    best_phi, best_d = None, -1.0
    G = 1024
    for g in range(G):
        phi = lo + (hi-lo)*g/(G-1)
        d = sum(wi*math.exp(-0.5*((bi-phi)/bw)**2) for wi, bi in zip(w, b)) / W
        if d > best_d: best_d, best_phi = d, phi
    return best_phi

def ivw_bootstrap_se(b, s, method="median", nb=2000, seed=42):
    """加权中位数/众数的 bootstrap se（重抽样 SNP 行，种子固定）"""
    rnd = random.Random(seed)
    n = len(b); vals = []
    for _ in range(nb):
        idx = [rnd.randrange(n) for _ in range(n)]
        bb = [b[i] for i in idx]; ss = [s[i] for i in idx]
        try:
            v = weighted_median(bb, ss) if method == "median" else weighted_mode(bb, ss)
        except Exception:
            continue
        vals.append(v)
    if not vals: return float("nan")
    m = sum(vals)/len(vals)
    return math.sqrt(sum((v-m)**2 for v in vals)/(len(vals)-1))

def presso_equiv(b, s):
    """MR-PRESSO 等价：全局异质性 + Bonferroni 残差离群 + distortion"""
    n = len(b)
    res = ivw_random_dl(b, s)
    tau2, bIV = res["tau2"], res["beta"]
    zs, outliers = [], []
    for i, (bi, si) in enumerate(zip(b, s)):
        si2 = si**2 + tau2
        zi = (bi - bIV) / math.sqrt(si2)
        pi = p_norm(zi)
        zs.append(zi)
        if pi < 0.05 / n:
            outliers.append(i)
    dist = None
    if outliers:
        keep = [i for i in range(n) if i not in outliers]
        res2 = ivw_random_dl([b[i] for i in keep], [s[i] for i in keep])
        zdiff = (res["beta"] - res2["beta"]) / math.sqrt(res["se"]**2 + res2["se"]**2)
        dist = {"n_outliers": len(outliers),
                "outlier_snps_idx": outliers,
                "beta_corrected": res2["beta"], "se_corrected": res2["se"],
                "p_corrected": res2["p"],
                "distortion_z": zdiff, "distortion_p": p_norm(zdiff)}
    return {"global_Q": res["Q"], "global_Q_df": res["Q_df"],
            "global_Q_pval": res["Q_pval"], "distortion": dist}

# ---------------- 方向性（Steiger 近似） ----------------
def directionality(rows):
    """比较两侧方差解释量。R2X 来自 Step 2（工具强度），R2Y 用 n_eff 近似。"""
    r2x = sum(r["R2_X"] for r in rows)
    k = len(rows)
    zy = [(r["beta_Y"] / r["se_Y"])**2 for r in rows]
    r2y = sum(z / (z + OUT_N_EFF_EUR) for z in zy)
    # 两侧各自的整体 F（同 Step 2 口径）
    nX = rows[0]["n_eff_X"]
    fX = (r2x/(1-r2x)) * ((nX - k - 1)/k) if r2x < 1 else float("inf")
    fY = (r2y/(1-r2y)) * ((OUT_N_EFF_EUR - k - 1)/k) if r2y < 1 else float("inf")
    return {"R2_X": r2x, "R2_Y": r2y, "n_eff_X": nX, "n_eff_Y_obs_scale": OUT_N_EFF_EUR,
            "F_X": fX, "F_Y": fY, "nR2_X": nX*r2x, "nR2_Y": OUT_N_EFF_EUR*r2y,
            "direction": "forward (X->Y)" if nX*r2x > OUT_N_EFF_EUR*r2y else "reverse?"}

# ---------------- 主流程 ----------------
allrows = {}
for t in TRAITS:
    path = os.path.join(RES, "harmonised_%s.tsv" % t)
    rows = []
    for d in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        rows.append({
            "SNP": d["SNP"], "rs_id": d["rs_id"],
            "beta_X": float(d["beta_X"]), "se_X": float(d["se_X"]),
            "beta_Y": float(d["beta_Y"]), "se_Y": float(d["se_Y"]),
            "F_X": float(d["F_X"]), "R2_X": float(d.get("R2_X", 0) or 0),
            "n_eff_X": float(d.get("n_eff_X", 0) or 0),
            "beta_XY": float(d["beta_XY"]), "se_XY": float(d["se_XY"])})
    allrows[t] = rows
    log("%s: 读入 harmonised 工具 %d 个" % (t, len(rows)))

# 从 instruments 文件合并 R2_X / n_eff_X（harmonised 文件未携带）
for t in TRAITS:
    inst = {}
    for d in csv.DictReader(open(os.path.join(ROOT, "work", "instruments_%s.tsv" % t),
                                 encoding="utf-8"), delimiter="\t"):
        inst[d["SNP"]] = (float(d["R2"]), float(d["n_eff"]))
    for r in allrows[t]:
        r["R2_X"], r["n_eff_X"] = inst.get(r["SNP"], (0.0, 0.0))
    miss = sum(1 for r in allrows[t] if r["n_eff_X"] == 0)
    log("%s: 合并 R2/n_eff，缺失 %d 个" % (t, miss))

summary_rows = []
for t in TRAITS:
    rows = allrows[t]
    bX = [r["beta_X"] for r in rows]; sX = [r["se_X"] for r in rows]
    bY = [r["beta_Y"] for r in rows]; sY = [r["se_Y"] for r in rows]
    bR = [r["beta_XY"] for r in rows]; sR = [r["se_XY"] for r in rows]
    k  = len(rows)
    out = {"trait": t, "trait_cn": TRAIT_CN[t], "nsnp": k}

    # 3.2 逐 SNP Wald ratio + 符号一致性
    zsign = [1 if x > 0 else -1 for x in bR]
    npos = sum(1 for z in zsign if z > 0)
    out["sign_test"] = {"n_positive": npos, "n_negative": k - npos,
                        "p_binomial_two_sided": p_binom_two(min(npos, k-npos), k)}

    # 3.3 IVW
    bF, seF = ivw(bR, sR)
    out["ivw_fixed"] = {"beta": bF, "se": seF, "p": p_norm(bF/seF)}
    res_rand = ivw_random_dl(bR, sR)
    out["ivw_random"] = res_rand

    # 3.4/3.5 Egger + 中位数 + 众数
    eg = egger(bX, bY, sY)
    out["mr_egger"] = eg
    bm = weighted_median(bR, sR)
    sem = ivw_bootstrap_se(bR, sR, "median")
    out["weighted_median"] = {"beta": bm, "se": sem, "p": p_norm(bm/sem)}
    bmo = weighted_mode(bR, sR)
    semo = ivw_bootstrap_se(bR, sR, "mode", nb=500)
    out["weighted_mode"] = {"beta": bmo, "se": semo, "p": p_norm(bmo/semo)}

    # 3.6 MR-PRESSO 等价
    out["presso_equiv"] = presso_equiv(bR, sR)

    # 3.7 留一法
    loo = []
    for i in range(k):
        keep = [j for j in range(k) if j != i]
        r2 = ivw_random_dl([bR[j] for j in keep], [sR[j] for j in keep])
        loo.append({"SNP": rows[i]["SNP"], "nsnp": k-1,
                    "beta": r2["beta"], "se": r2["se"], "p": r2["p"]})
    with open(os.path.join(RES, "loo_%s.tsv" % t), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["SNP", "nsnp", "beta", "se", "p"], delimiter="\t")
        w.writeheader(); w.writerows(loo)

    # 3.8 方向性
    out["directionality"] = directionality(rows)

    with open(os.path.join(RES, "mr_main_%s.json" % t), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    # 汇总行
    def row(method, d, note=""):
        summary_rows.append({"trait": t, "trait_cn": TRAIT_CN[t], "method": method,
                             "nsnp": d.get("nsnp", k),
                             "beta": round(d["beta"], 6), "se": round(d["se"], 6),
                             "or": round(math.exp(d["beta"]), 6),
                             "p": "%.4g" % d["p"], "note": note})
    row("IVW(固定效应)", out["ivw_fixed"])
    row("IVW(随机效应 DL)", out["ivw_random"],
        "Q=%s df=%d Q_p=%.3g tau2=%.3g" % ("%.2f" % res_rand["Q"], res_rand["Q_df"],
                                            res_rand["Q_pval"], res_rand["tau2"]))
    if eg:
        row("MR-Egger", eg, "intercept=%s se=%s p=%.3g Q_p=%.3g"
            % (round(eg["intercept"], 6), round(eg["intercept_se"], 6),
               eg["intercept_pval"], eg["Q_pval"]))
    row("加权中位数", out["weighted_median"], "bootstrap se (2000)")
    row("加权众数", out["weighted_mode"], "bootstrap se (500)")
    pe = out["presso_equiv"]
    dnote = ""
    if pe["distortion"]:
        dc = pe["distortion"]
        dnote = "离群 %d 个, 校正beta=%s p=%.3g, distortion p=%.3g" % (
            dc["n_outliers"], round(dc["beta_corrected"], 6), dc["p_corrected"], dc["distortion_p"])
    summary_rows.append({"trait": t, "trait_cn": TRAIT_CN[t], "method": "MR-PRESSO等价(全局异质性)",
                         "nsnp": k, "beta": "", "se": "", "or": "", "p": "%.4g" % pe["global_Q_pval"],
                         "note": dnote})
    dj = out["directionality"]
    summary_rows.append({"trait": t, "trait_cn": TRAIT_CN[t], "method": "方向性(R²对比, Steiger近似)",
                         "nsnp": k, "beta": "", "se": "", "or": "",
                         "p": "", "note": "R2X=%.4f R2Y=%.5f nR2X=%.0f nR2Y=%.0f → %s"
                         % (dj["R2_X"], dj["R2_Y"], dj["nR2_X"], dj["nR2_Y"], dj["direction"])})
    st = out["sign_test"]
    summary_rows.append({"trait": t, "trait_cn": TRAIT_CN[t], "method": "符号一致性",
                         "nsnp": k, "beta": "", "se": "", "or": "",
                         "p": "%.4g" % st["p_binomial_two_sided"],
                         "note": "正 %d / 负 %d" % (st["n_positive"], st["n_negative"])})

    log("%s 完成: IVW固定 β=%.4f (p=%.3g) | IVW随机 β=%.4f (p=%.3g) | Egger β=%.4f 截距p=%.3g | Q_p=%.3g"
        % (t, out["ivw_fixed"]["beta"], out["ivw_fixed"]["p"],
           res_rand["beta"], res_rand["p"], eg["beta"] if eg else float("nan"),
           eg["intercept_pval"] if eg else float("nan"), res_rand["Q_pval"]))

# 汇总表
with open(os.path.join(RES, "mr_main_summary.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["trait", "trait_cn", "method", "nsnp", "beta", "se", "or", "p", "note"],
                       delimiter="\t")
    w.writeheader(); w.writerows(summary_rows)
log("全部完成。汇总表: %s" % os.path.join(RES, "mr_main_summary.tsv"))
