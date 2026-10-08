# -*- coding: utf-8 -*-
"""
coloc.abf 的纯 Python 实现（Giambartolomei 2014, PLoS Genet；Wakefield 2009 近似贝叶斯因子）。

不依赖 R 的 coloc 包，便于在既有纯 Python 流程中复现。

核心公式
  每个变异 i 的近似贝叶斯因子（ABF）：
      V_i = se_i²           （似然方差）
      W   = sd.prior²       （先验方差）
      r_i = W / (V_i + W)
      lABF_i = 0.5*log(1 - r_i) + 0.5 * r_i * z_i²        z_i = beta_i / se_i
  五个假设的似然（L1=Σ ABF1, L2=Σ ABF2, L12=Σ ABF1*ABF2）：
      H0 = 1
      H1 = L1                        （仅性状1 有关联）
      H2 = L2                        （仅性状2 有关联）
      H3 = L1*L2 - L12               （两者都关联，但为不同因果变异）
      H4 = L12                       （两者共享同一因果变异）
  后验 ∝ 先验 × 似然，先验 = (H0:1, H1:p1, H2:p2, H3:p1*p2, H4:p12)

判读惯例（本文采用）
  PP.H4 ≥ 0.8 → 强证据支持共享因果变异
  PP.H3 ≥ 0.8 → 强证据支持「两个独立因果变异因 LD 相关」
  PP.H4 与 PP.H3 均 < 0.5 且 PP.H0/H1/H2 之一占优 → 无共享、仅一方有关联
  其余 → 不确定（功效不足）
"""
import math

LOG10 = math.log(10.0)


def _logsumexp(vals):
    """log10(Σ 10^vals)"""
    if not vals:
        return -math.inf
    m = max(vals)
    if m == -math.inf:
        return -math.inf
    return m + math.log10(sum(10 ** (v - m) for v in vals))


def labf(beta, se, sd_prior=0.15):
    """单个变异的 log10(ABF)；参数非法返回 None"""
    if beta is None or se is None:
        return None
    if not (se > 0) or not math.isfinite(beta) or not math.isfinite(se):
        return None
    V = se * se
    W = sd_prior * sd_prior
    r = W / (V + W)
    z2 = (beta / se) ** 2
    return (0.5 * math.log(1.0 - r) + 0.5 * r * z2) / LOG10


def coloc_abf(rows, p1=1e-4, p2=1e-4, p12=1e-5, sd_prior=0.15):
    """
    rows: 可迭代的 (beta1, se1, beta2, se2) —— 已对齐等位基因、同一批变异
    返回 dict(n, H0..H4 后验, L1/L2/L12, 判定)
    """
    l1, l2, l12 = [], [], []
    n = 0
    for b1, s1, b2, s2 in rows:
        a = labf(b1, s1, sd_prior)
        b = labf(b2, s2, sd_prior)
        if a is None or b is None:
            continue
        n += 1
        l1.append(a); l2.append(b); l12.append(a + b)
    if n == 0:
        return dict(n=0, verdict="无可用变异")
    L0 = 0.0                                    # log10(1)
    L1 = _logsumexp(l1)
    L2 = _logsumexp(l2)
    L12 = _logsumexp(l12)
    # H3 = L1*L2 - L12（数值上用 log 差）
    L3 = _logsubexp(L1 + L2, L12)
    L4 = L12
    lp = [0.0, math.log10(p1), math.log10(p2), math.log10(p1) + math.log10(p2), math.log10(p12)]
    lk = [L0, L1, L2, L3, L4]
    tot = [_logsumexp([lp[k] + lk[k]]) for k in range(5)]
    denom = _logsumexp(tot)
    pp = [10 ** (tot[k] - denom) for k in range(5)]
    return dict(n=n, H0=pp[0], H1=pp[1], H2=pp[2], H3=pp[3], H4=pp[4],
                verdict=verdict(pp))


def _logsubexp(a, b):
    """log10(10^a - 10^b)，b 可能为 -inf"""
    if b == -math.inf:
        return a
    if a <= b:
        return -math.inf
    return a + math.log10(1.0 - 10 ** (b - a))


def verdict(pp):
    _, h1, h2, h3, h4 = pp
    if h4 >= 0.8:
        return "H4 共享因果变异（强）"
    if h3 >= 0.8:
        return "H3 独立因果变异·LD 相关（强）"
    if h4 >= 0.5:
        return "H4 倾向共享（中等）"
    if h3 >= 0.5:
        return "H3 倾向独立（中等）"
    if max(h1, h2) >= 0.8:
        return "仅一方有关联"
    if pp[0] >= 0.8:
        return "双方均无关联"
    return "不确定（功效不足）"


if __name__ == "__main__":
    # 自检：构造一个明显的共享信号
    rows = []
    for i in range(200):
        b1 = 0.3 if i == 50 else 0.0
        b2 = 0.3 if i == 50 else 0.0
        rows.append((b1, 0.05, b2, 0.05))
    print(coloc_abf(rows))
    # 独立信号
    rows = [(0.3 if i == 50 else 0.0, 0.05, 0.3 if i == 120 else 0.0, 0.05) for i in range(200)]
    print(coloc_abf(rows))
