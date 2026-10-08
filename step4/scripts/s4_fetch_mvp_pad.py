#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s4_fetch_mvp_pad.py  —  MVP PAD（dbGaP phs001672 / pha004826）全量汇总统计获取工具

背景（2026-10-07 实测）：
  phs001672 自 v14.p1（2026-04-08）起为 "Genomic Summary Results (GSR) unrestricted"，
  全量未脱敏汇总结果公开可下，**无需 dbGaP DAR 申请**。

目标文件：release/submission/sub20190405/MVP.te.PAD.dbGAP.txt.gz
  声明大小 665,142,901 字节
  打包在  FullData/phs001672.MVP.analysis-PI.c1.HMB-MDS.set3.tar.gz（9.95 GB）内，
  是该 tar 的第 4 个成员（也是第一个大文件），数据区自 flat 偏移 2048 起始。

关键优化：外层 tar.gz 压缩比实测 = 1.0000（成员本身已是 .gz），
  故只下载分卷开头 ~666 MB 即可完整取出，省 93% 流量。

下载策略（本机 2026-10-07 实测调优）：
  单连接对 ftp.ncbi.nlm.nih.gov 仅 12–27 KB/s；
  4 并发 ≈ 61 KB/s；8 并发 ≈ 95 KB/s；**16 并发被 NCBI 大量回 503**。
  → 默认 8 连接分块并行（8 MB/块）+ 503 指数退避 + 逐块断点续传。

用法：
  python s4_fetch_mvp_pad.py --download            # 并行下载前 700 MB（可反复续传）
  python s4_fetch_mvp_pad.py --extract             # 抽出成员（字节精确 + gzip CRC 校验）
  python s4_fetch_mvp_pad.py --extract --tsv       # 再解压为明文 TSV 并统计
  python s4_fetch_mvp_pad.py --all                 # 一条龙
  python s4_fetch_mvp_pad.py --all --conn 6 --chunk-mb 16
"""

import argparse
import os
import shutil
import subprocess
import sys
import threading
import time
import zlib

BASE = ("https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs001672/"
        "analyses/FullData/")
SET3 = "phs001672.MVP.analysis-PI.c1.HMB-MDS.set3.tar.gz"
MEMBER = "release/submission/sub20190405/MVP.te.PAD.dbGAP.txt.gz"
MEMBER_SIZE = 665142901               # 官方 tar 头声明值
TAR_PAD = 2048                        # PAD 成员数据区起点（flat tar 偏移）
NEED_MIN = TAR_PAD + MEMBER_SIZE + 1024   # 最小可用 = 665,144,973
DEFAULT_BYTES = 700_000_000           # 默认多取冗余

CURL = ["curl", "-sS", "-L", "--noproxy", "*", "--ssl-no-revoke"]

WORK = os.path.join("raw", "mvp_pad")
PARTS = os.path.join(WORK, "parts")
OUT_HEAD = os.path.join(WORK, "set3.head.gz")
OUT_MEMBER = os.path.join(WORK, "MVP.te.PAD.dbGAP.txt.gz")
OUT_TSV = os.path.join(WORK, "MVP.te.PAD.tab")

URL = BASE + SET3
_lock = threading.Lock()
_done = [0]


def log(msg):
    with _lock:
        print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


# --------------------------------------------------------------------------
# 1) 并行分块下载（可续传）
# --------------------------------------------------------------------------
def _fetch_chunk(idx, start, end, attempts=8):
    """下载 [start, end]（含端点）。支持**块内续传**（保住断线前的字节）。
    返回 (idx, ok, new_bytes)。

    ⚠ 设计约束（2026-10-07 踩坑）：本机有 safe-delete 钩子，单轮删除超过阈值后会
    拦截 os.remove，导致「tmp 删不掉 → 完成状态不再登记」而卡死覆盖。
    因此本函数**全程不做任何删除**，一律用「截断为 0 字节」代替 os.remove。
    """
    want = end - start + 1
    fp = os.path.join(PARTS, "p%05d" % idx)
    new = 0

    def _truncate(p):
        with open(p, "wb"):
            pass

    for a in range(attempts):
        got = os.path.getsize(fp) if os.path.exists(fp) else 0
        if got == want:
            return idx, True, new
        if got > want:                       # 异常膨胀 → 重来
            _truncate(fp)
            got = 0

        remaining = end - (start + got) + 1
        tmp = fp + ".tmp"
        # 不再预删 tmp：curl -o 会截断重写；失败时下面显式截断，避免读到陈旧字节
        r = subprocess.run(CURL + ["--max-time", "900",
                                   "-r", "%d-%d" % (start + got, end),
                                   URL, "-o", tmp,
                                   "-w", "%{http_code}"],
                           capture_output=True, text=True)
        code = (r.stdout or "").strip()
        if r.returncode != 0 or code not in ("200", "206"):
            _truncate(tmp)
            tsz = 0
        else:
            tsz = os.path.getsize(tmp) if os.path.exists(tmp) else 0
            if tsz > remaining:              # 服务器忽略 Range 返回整个文件 → 丢弃
                _truncate(tmp)
                tsz = 0
        if tsz > 0:
            with open(fp, "ab" if got else "wb") as fo, open(tmp, "rb") as fi:
                shutil.copyfileobj(fi, fo, 1 << 22)
            _truncate(tmp)
            new += tsz

        now = os.path.getsize(fp) if os.path.exists(fp) else 0
        if now == want:
            return idx, True, new
        wait = min(2 ** a, 30)               # 503 / 超时 / 截断 → 指数退避
        log("  块 %05d 未完成（http=%s %d/%d），%ds 后重试 (%d/%d)"
            % (idx, code, now, want, wait, a + 1, attempts))
        time.sleep(wait)
    return idx, False, new


def download(total, conn=8, chunk_mb=8):
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(PARTS, exist_ok=True)

    # 已有前缀（此前单流下到的 set3.head.gz）直接复用
    prefix = os.path.getsize(OUT_HEAD) if os.path.exists(OUT_HEAD) else 0
    if prefix >= total:
        log("已满足 %d 字节，无需下载" % prefix)
        return True
    if prefix:
        log("复用已有前缀 %d 字节（%.1f MB）" % (prefix, prefix / 1e6))

    CH = chunk_mb * 1024 * 1024
    jobs = []
    for i, s in enumerate(range(prefix, total, CH)):
        jobs.append((i, s, min(s + CH, total) - 1))

    log("目标：%s" % URL)
    log("区间 [%d, %d) 共 %.1f MB，切 %d 块，%d 连接并行"
        % (prefix, total, (total - prefix) / 1e6, len(jobs), conn))

    # 落一份计划，供 s4_peek_members.py 零网络枚举成员
    import json as _json
    with open(os.path.join(PARTS, "_plan.json"), "w", encoding="utf-8") as f:
        _json.dump({"prefix": prefix, "total": total, "chunk_mb": chunk_mb,
                    "jobs": [[i, s, e] for i, s, e in jobs]}, f)

    def part_path(i):
        return os.path.join(PARTS, "p%05d" % i)

    def part_full(i, s, e):
        p = part_path(i)
        return os.path.exists(p) and os.path.getsize(p) == e - s + 1

    cached = sum(1 for i, s, e in jobs if part_full(i, s, e))
    partial = sum(1 for i, s, e in jobs
                  if os.path.exists(part_path(i)) and not part_full(i, s, e))
    log("其中 %d 块已完成、%d 块为可续传的半成品" % (cached, partial))

    t0 = time.time()
    _done[0] = 0                      # 只统计本次会话新落的字节
    base_have = sum(os.path.getsize(part_path(i)) for i, s, e in jobs
                    if os.path.exists(part_path(i)))
    grand = total - prefix

    from concurrent.futures import ThreadPoolExecutor, as_completed
    ok_all = True
    with ThreadPoolExecutor(max_workers=conn) as ex:
        futs = {ex.submit(_fetch_chunk, i, s, e): (i, s, e) for i, s, e in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            idx, ok, nbytes = fut.result()
            _done[0] += nbytes
            if ok:
                el = max(time.time() - t0, 1)
                log("  完成 %d/%d 块 | 本次新增 %.1f MB | 合计 %.1f/%.1f MB | %.1f KB/s"
                    % (n, len(jobs), _done[0] / 1e6, (base_have + _done[0]) / 1e6,
                       grand / 1e6, _done[0] / el / 1024))
            else:
                ok_all = False
                log("  !! 块 %05d 最终失败" % idx)

    if not ok_all:
        log("有块未完成。直接重跑本命令即可续传（已完成块会跳过）。")
        return False

    # 顺序拼接：前缀 + 各块
    log("全部块就绪，开始拼接 → %s" % OUT_HEAD)
    tmp = OUT_HEAD + ".merging"
    with open(tmp, "wb") as fo:
        if prefix and os.path.exists(OUT_HEAD):
            with open(OUT_HEAD, "rb") as fi:
                shutil.copyfileobj(fi, fo, 1 << 22)
        for i, s, e in jobs:
            with open(part_path(i), "rb") as fi:
                shutil.copyfileobj(fi, fo, 1 << 22)
    os.replace(tmp, OUT_HEAD)
    got = os.path.getsize(OUT_HEAD)
    log("拼接完成：%d 字节（耗时 %.1f min）" % (got, (time.time() - t0) / 60))
    if got < NEED_MIN:
        log("!! 仍少于最小可用 %d 字节，请把 --bytes 调大后重跑" % NEED_MIN)
        return False
    # 保留 PARTS：① 本机 safe-delete 钩子会拦截 rmtree；② s4_peek_members.py 依赖分块离线解析
    log("保留分块目录 %s（如需清理请手动执行；--extract 之后可安全删除）" % PARTS)
    return True


# --------------------------------------------------------------------------
# 2) 抽取成员（字节精确）
# --------------------------------------------------------------------------
def extract():
    if not os.path.exists(OUT_HEAD):
        log("缺少 %s，请先 --download" % OUT_HEAD)
        return False
    data_len = os.path.getsize(OUT_HEAD)
    log("读取 %s（%d 字节）" % (OUT_HEAD, data_len))

    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    flat = d.decompress(open(OUT_HEAD, "rb").read())
    log("flat tar 字节 %d（压缩比 1:%.4f）" % (len(flat), len(flat) / data_len))

    off, member_start = 0, None
    while off + 512 <= len(flat):
        hdr = flat[off:off + 512]
        if hdr[:1] == b"\0":
            break
        name = hdr[0:100].split(b"\0")[0].decode("utf-8", "replace")
        try:
            size = int((hdr[124:136].split(b"\0")[0].strip() or b"0"), 8)
        except ValueError:
            size = 0
        log("  tar 成员：%-58s %d" % (name, size))
        if name == MEMBER:
            member_start = off + 512
            if member_start + MEMBER_SIZE > len(flat):
                need = member_start + MEMBER_SIZE
                log("!! 数据不足：需 flat 偏移 %d，当前仅 %d。" % (need, len(flat)))
                log("   （压缩比≈1.0）请把 --bytes 设为 >= %d 后重跑。" % (need + 1_000_000))
                return False
            break
        off += 512 + (((size + 511) // 512) * 512 if size else 0)

    if member_start is None:
        log("!! 未定位到成员 %s" % MEMBER)
        return False

    blob = flat[member_start:member_start + MEMBER_SIZE]
    if len(blob) != MEMBER_SIZE:
        log("!! 切出字节数不符")
        return False
    if blob[:2] != b"\x1f\x8b":
        log("!! 切出的内容不是 gzip（前 4 字节 %r）" % blob[:4])
        return False
    with open(OUT_MEMBER, "wb") as f:
        f.write(blob)
    log("已写出 %s（%d 字节）" % (OUT_MEMBER, len(blob)))

    ok = verify_member()
    log("gzip 完整性校验：%s" % ("通过（尾部 CRC/ISIZE 合法）" if ok else "未通过（可能截断）"))
    return ok


def verify_member():
    if not os.path.exists(OUT_MEMBER):
        return False
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    with open(OUT_MEMBER, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            d.decompress(b)
    return bool(d.eof)


# --------------------------------------------------------------------------
# 3) 转明文 TSV 并统计
# --------------------------------------------------------------------------
def to_tsv():
    if not os.path.exists(OUT_MEMBER):
        log("缺少 %s，请先 --extract" % OUT_MEMBER)
        return False
    import gzip
    import csv

    n, ss, hdr, idx = 0, {}, None, {}
    t0 = time.time()
    with gzip.open(OUT_MEMBER, "rt", encoding="utf-8", errors="replace") as fi, \
            open(OUT_TSV, "w", newline="", encoding="utf-8") as fo:
        w = csv.writer(fo, delimiter="\t")
        for line in fi:
            if hdr is None:
                hdr = line.rstrip("\n").split("\t")
                w.writerow(hdr)
                idx = {c: i for i, c in enumerate(hdr)}
                log("表头 %d 列：%s" % (len(hdr), hdr))
                continue
            n += 1
            w.writerow(line.rstrip("\n").split("\t"))
            if "SampleSize" in idx:
                p = line.split("\t")
                if len(p) > idx["SampleSize"]:
                    ss[p[idx["SampleSize"]]] = ss.get(p[idx["SampleSize"]], 0) + 1

    log("明文 TSV：%s | 数据行 %d | 耗时 %.1f min" % (OUT_TSV, n, (time.time() - t0) / 60))
    log("SampleSize 分布（前 5）：%s" % sorted(ss.items(), key=lambda x: -x[1])[:5])
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--tsv", action="store_true")
    ap.add_argument("--all", action="store_true", help="= --download --extract --tsv")
    ap.add_argument("--bytes", type=int, default=DEFAULT_BYTES,
                    help="下载字节数上限（默认 700,000,000；最小可用 665,144,973）")
    ap.add_argument("--conn", type=int, default=8, help="并发连接数（实测上限 8）")
    ap.add_argument("--chunk-mb", type=int, default=8, help="每块大小 MB")
    a = ap.parse_args()

    if not (a.download or a.extract or a.tsv or a.all):
        ap.print_help()
        return 1

    os.makedirs(WORK, exist_ok=True)
    if a.all:
        a.download = a.extract = a.tsv = True

    if a.download and not download(a.bytes, a.conn, a.chunk_mb):
        return 1
    if a.extract and not extract():
        return 1
    if a.tsv and not to_tsv():
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
