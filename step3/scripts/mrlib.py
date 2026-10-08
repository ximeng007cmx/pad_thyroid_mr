# -*- coding: utf-8 -*-
"""共享 MR 统计函数库（自 s3_mr_main.py 程序化抽取，Step 3 已验证口径）

包含：p_norm / p_chi2（含 Γ(a) 归一化，已用 χ² 临界值校验）/ ivw / q_stat /
ivw_random_dl(DL) / egger（Bowden 2016 规范形式 b_Y~b_X）/ weighted_median /
weighted_mode（加权高斯核 KDE）/ ivw_bootstrap_se / presso_equiv / p_binom_two
"""
import os, csv, json, math, time, datetime, random

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

