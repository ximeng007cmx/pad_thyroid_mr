#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s4_peek_members.py  —  边下边看：零网络开销枚举 set3.tar.gz 里已覆盖的 tar 成员

背景：set3.tar.gz 是**固实 gzip**，无法跳读；要拿到后面的成员必须顺序流过前面的。
      但下载是按块并行的，已落盘的分块可以直接本地解压——于是「已下到哪、下一个成员是谁」
      可以随时离线问出来，不用等整卷下完。

依赖：`s4_fetch_mvp_pad.py --download` 运行时写的 `raw/mvp_pad/parts/_plan.json`。

用法：
  python s4_peek_members.py            # 列出连续覆盖范围内的全部 tar 成员
  python s4_peek_members.py --grep PAD # 只显示名字含 PAD 的
"""

import argparse
import json
import os
import sys
import zlib

WORK = os.path.join("raw", "mvp_pad")
PARTS = os.path.join(WORK, "parts")
PREFIX = os.path.join(WORK, "set3.head.gz")


def flat_stream(paths):
    """把若干按序拼接的文件当一条 gzip 流解压，逐段 yield 解压后的字节。"""
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    for p in paths:
        if not os.path.exists(p):
            return
        try:
            with open(p, "rb") as f:
                while True:
                    b = f.read(1 << 20)
                    if not b:
                        break
                    out = d.decompress(b)
                    if out:
                        yield out
        except OSError:
            return
        if d.eof:
            return


class Reader:
    def __init__(self, it):
        self.it = it
        self.buf = b""
        self.n = 0          # 已消费的解压字节数

    def _fill(self):
        try:
            self.buf += next(self.it)
            return True
        except StopIteration:
            return False

    def read(self, k):
        while len(self.buf) < k:
            if not self._fill():
                break
        r, self.buf = self.buf[:k], self.buf[k:]
        self.n += len(r)
        return r

    def skip(self, k):
        while k > 0:
            if self.buf:
                t = min(k, len(self.buf))
                self.buf = self.buf[t:]
                k -= t
                self.n += t
            elif not self._fill():
                return

    def copy_to(self, k, fo):
        """流式写出 k 字节（不把大成员缓进内存），返回实际写出字节数。"""
        written = 0
        while k > 0:
            if not self.buf and not self._fill():
                break
            if self.buf:
                t = min(k, len(self.buf))
                fo.write(self.buf[:t])
                self.buf = self.buf[t:]
                k -= t
                written += t
                self.n += t
        return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grep", default=None, help="只显示名字含该串的成员")
    ap.add_argument("--all", action="store_true", help="不过滤")
    ap.add_argument("--save-member", default=None,
                    help="顺带把这一个成员完整抠出来（写 --out）")
    ap.add_argument("--out", default=None, help="--save-member 的输出路径")
    a = ap.parse_args()
    save_name, out_path = a.save_member, a.out
    if save_name and not out_path:
        out_path = os.path.join(WORK, os.path.basename(save_name))

    plan_p = os.path.join(PARTS, "_plan.json")
    if not os.path.exists(plan_p):
        print("找不到 %s —— 请先让下载脚本跑起来" % plan_p)
        return 1
    plan = json.load(open(plan_p, encoding="utf-8"))
    prefix, jobs = plan["prefix"], plan["jobs"]

    # 找连续覆盖：前缀 + 从第 0 块起所有「存在且完整」的块
    files = []
    if prefix and os.path.exists(PREFIX):
        files.append(PREFIX)
    covered = prefix
    for i, s, e in jobs:
        p = os.path.join(PARTS, "p%05d" % i)
        want = e - s + 1
        if os.path.exists(p) and os.path.getsize(p) == want:
            files.append(p)
            covered = e + 1
        else:
            got = os.path.getsize(p) if os.path.exists(p) else 0
            print("（第 %d 块未完成：%d/%d 字节 → 后面的块先不参与）" % (i, got, want))
            break

    print("计划总量 %d 字节（%.1f MB）| 连续可用 %.1f MB"
          % (plan["total"], plan["total"] / 1e6, (covered - prefix) / 1e6))

    r = Reader(flat_stream(files))
    members = []
    while True:
        hdr = r.read(512)
        if len(hdr) < 512:
            break
        if hdr[:1] == b"\0":
            break
        name = hdr[0:100].split(b"\0")[0].decode("utf-8", "replace")
        try:
            size = int((hdr[124:136].split(b"\0")[0].strip() or b"0"), 8)
        except ValueError:
            size = 0
        base = r.n - 512                      # 该成员头部的 flat 偏移
        data_off = base + 512
        members.append((name, size, data_off))
        consumed = 0
        if save_name and name == save_name:
            os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
            with open(out_path, "wb") as fo:
                consumed = r.copy_to(size, fo)
            print(">>> 已保存成员：%s\n    → %s（%d / %d 字节）%s"
                  % (name, out_path, consumed, size,
                     "  ✅ 完整" if consumed == size else "  !! 截断"))
        if consumed < size:
            r.skip(size - consumed)
        r.skip((size + 511) // 512 * 512 - size)

    flat_seen = r.n
    print("flat tar 实际解出 %.1f MB（含头部与填充）" % (flat_seen / 1e6))
    print()
    print("%-62s %14s %16s  %s" % ("tar 成员", "声明大小", "数据区起点", "完整性"))
    for name, size, off in members:
        if size == 0:
            status = "dir"
        elif off + size <= flat_seen:
            status = "完整 ✅"
        else:
            status = "截断（差 %.1f MB）" % ((off + size - flat_seen) / 1e6)
        if a.grep and a.grep.lower() not in name.lower():
            continue
        print("%-62s %14d %16d  %s" % (name, size, off, status))
    return 0


if __name__ == "__main__":
    sys.exit(main())
