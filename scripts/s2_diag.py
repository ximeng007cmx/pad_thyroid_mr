# -*- coding: utf-8 -*-
"""诊断面板 .bed 是否损坏 / 候选位点匹配率"""
import os, gzip, random, math

ROOT  = r"D:/WorkSpace_Study/pad_thyroid_mr"
PANEL = os.path.join(ROOT, "ref", "1kg_eur_hg38", "allchr.EUR.biallelicsnps.cbb")
WORK  = os.path.join(ROOT, "work")
NSAMP = 294

# ---- 读 bim
print("读取 bim ...")
bim = []
with open(PANEL + ".bim", "r", encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        bim.append((p[0], int(p[3]), p[1], p[4], p[5]))
print("  bim 变异数:", len(bim))

pos2idx = {}
for i, (c, p, rs, a1, a2) in enumerate(bim):
    pos2idx.setdefault("%s:%d" % (c, p), i)

# ---- 读 bed
raw = open(PANEL + ".bed", "rb").read()
print("  bed 头 3 字节:", raw[:3].hex(), "(期望 6c1b01)")
bsize = (NSAMP + 3) // 4
print("  每 SNP 字节数:", bsize, " 预期总 SNP:", (len(raw) - 3) / bsize)

def get_geno(idx):
    off = 3 + idx * bsize
    blk = raw[off:off + bsize]
    out = []
    for k in range(NSAMP):
        byte = blk[k >> 2]
        sh = (k & 3) << 1
        code = (byte >> sh) & 3
        # PLINK1: 00=A1/A1, 01=missing, 10=A1/A2, 11=A2/A2  → A2 剂量 0/miss/1/2
        out.append(None if code == 1 else (0 if code == 0 else (1 if code == 2 else 2)))
    return out

def r2(g1, g2):
    n = 0; s1 = 0; s2 = 0; s11 = 0; s22 = 0; s12 = 0
    for a, b in zip(g1, g2):
        if a is None or b is None:
            continue
        n += 1; s1 += a; s2 += b; s11 += a * a; s22 += b * b; s12 += a * b
    if n < 3:
        return None
    v1 = s11 / n - (s1 / n) ** 2
    v2 = s22 / n - (s2 / n) ** 2
    if v1 <= 0 or v2 <= 0:
        return None
    cov = s12 / n - (s1 / n) * (s2 / n)
    return (cov * cov) / (v1 * v2)

def maf(g):
    c = [x for x in g if x is not None]
    if not c:
        return None, 0
    f = sum(c) / (2 * len(c))
    return min(f, 1 - f), len(c)

print("\n=== 1) 随机 12 个面板 SNP 的 MAF / 缺失率 ===")
random.seed(7)
for idx in random.sample(range(len(bim)), 12):
    g = get_geno(idx)
    m, nn = maf(g)
    nm = sum(1 for x in g if x is None)
    print("  %-16s chr%s:%-10d MAF=%-8s n=%d missing=%d" % (bim[idx][2], bim[idx][0], bim[idx][1],
          ("%.4f" % m) if m is not None else "NA", nn, nm))

print("\n=== 2) 距离与 r² 的关系（chr1 相邻 / 远距）===")
c1 = sorted([b for b in bim if b[0] == "1"], key=lambda x: x[1])
import bisect
pos1 = [b[1] for b in c1]
for gap in (1000, 50_000, 1_000_000, 5_000_000):
    pairs = []
    for i in range(0, len(c1) - 1, max(1, len(c1) // 40)):
        j = bisect.bisect_left(pos1, c1[i][1] + gap)
        if j < len(c1) and abs(pos1[j] - pos1[i] - gap) < gap * 0.05:
            pairs.append((i, j))
    if not pairs:
        print("  gap %-9d : 无样本" % gap); continue
    rs = []
    for i, j in pairs[:25]:
        v = r2(get_geno(bim.index(c1[i]) if False else [k for k, b in enumerate(bim) if b is c1[i]][0]),
               get_geno([k for k, b in enumerate(bim) if b is c1[j]][0]))
        if v is not None:
            rs.append(v)
    if rs:
        print("  gap %-9d : n=%d  r² 均值=%.4f  中位=%.4f  >0.99 占比=%.1f%%"
              % (gap, len(rs), sum(rs) / len(rs), sorted(rs)[len(rs) // 2],
                 100 * sum(1 for x in rs if x > 0.99) / len(rs)))

print("\n=== 3) 候选位点与面板的匹配率（按 MAF 分层）===")
for tag in ("hypo", "tsh", "ft4"):
    tot = hit = 0
    bins = {"<1%": [0, 0], "1-5%": [0, 0], ">=5%": [0, 0]}
    with open(os.path.join(WORK, "cand_%s.tsv" % tag), encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            eaf = float(p[ci["eaf"]]); m = min(eaf, 1 - eaf)
            tot += 1
            k = "%s:%s" % (p[ci["chr"]], p[ci["pos"]])
            ok = k in pos2idx
            hit += ok
            b = "<1%" if m < 0.01 else ("1-5%" if m < 0.05 else ">=5%")
            bins[b][0] += 1; bins[b][1] += ok
    print("  %-5s 候选 %6d  匹配 %5d (%.1f%%)" % (tag, tot, hit, 100 * hit / tot), end="  ")
    print(" | ".join("%s %d/%d(%.0f%%)" % (b, v[1], v[0], 100 * v[1] / max(v[0], 1)) for b, v in bins.items()))
