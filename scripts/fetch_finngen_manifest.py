#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FinnGen 公共桶端点清单（manifest）抓取器 —— 通用版，支持任意 release。

背景
----
FinnGen 各版本的样本量不在 summary-stats 文件头里，只在 bucket 内的 manifest
中。路径规律：
    gs://finngen-public-data-<rel>/summary_stats/<rel>_manifest.tsv
其中 <rel> 为小写 release 号，如 r3 / r13（注意 R3 是 `r3_manifest.tsv`）。
manifest 列：phenocode / name / category / n_cases / n_controls / path_bucket / path_https

用法
----
    python fetch_finngen_manifest.py r3
    python fetch_finngen_manifest.py r13 --grep PAD
    python fetch_finngen_manifest.py r3 --out D:/.../finngen_r3_manifest.tsv

实测（2026-10-08）：
  - R3  bucket = finngen-public-data-r3，端点数 1,801，最大 N 端点 AB1_ARTHROPOD = 135,638（= R3 总样本）
  - R3  `I9_PAD` = 4,287 cases / 128,425 controls = 132,712
  - 枚举整桶对象用 JSON API（见 list_bucket_objects），可发现非 summary_stats/ 下的其它资源

坑
--
  - r3.finngen.fi 站点在本机不可达（DNS/HTTP 000）；storage.googleapis.com 可达。
  - bucket 对象枚举用 `https://storage.googleapis.com/storage/v1/b/<bucket>/o?maxResults=1000&pageToken=...`
"""
import argparse
import csv
import io
import os
import subprocess
import sys

BUCKET_TMPL = "finngen-public-data-{rel}"


def curl(url, binary=False, timeout="90"):
    out = subprocess.run(["curl", "-s", "--noproxy", "*", "--max-time", timeout, url],
                         capture_output=True)
    if out.returncode != 0:
        raise RuntimeError("curl failed: %s" % out.stderr.decode("utf-8", "replace")[:200])
    return out.stdout if binary else out.stdout.decode("utf-8", "replace")


def fetch_manifest(rel):
    rel = rel.lower().lstrip("r").join(["r", ""]) if False else "r" + rel.lower().lstrip("r")
    bucket = BUCKET_TMPL.format(rel=rel)
    url = "https://storage.googleapis.com/%s/summary_stats/%s_manifest.tsv" % (bucket, rel)
    print("GET", url, flush=True)
    txt = curl(url)
    if txt.lstrip().startswith("<?xml"):
        raise RuntimeError("NoSuchKey — 该 release 无 manifest：\n" + txt[:300])
    return rel, bucket, txt


def list_bucket_objects(rel, prefix=""):
    bucket = BUCKET_TMPL.format(rel="r" + rel.lower().lstrip("r"))
    names, token = [], None
    while True:
        url = "https://storage.googleapis.com/storage/v1/b/%s/o?maxResults=1000" % bucket
        if prefix:
            url += "&prefix=" + prefix
        if token:
            url += "&pageToken=" + token
        import json
        d = json.loads(curl(url))
        for it in d.get("items", []):
            names.append((it["name"], it.get("size")))
        token = d.get("nextPageToken")
        if not token:
            break
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("release", help="FinnGen release，如 r3 / r13")
    ap.add_argument("--grep", default=None, help="按 phenocode/name 过滤（大小写不敏感）")
    ap.add_argument("--out", default=None, help="保存 manifest 的路径（默认不保存）")
    ap.add_argument("--list-bucket", action="store_true", help="改为枚举整桶对象")
    args = ap.parse_args()

    if args.list_bucket:
        for n, s in list_bucket_objects(args.release):
            print("%-90s %s" % (n, s))
        return

    rel, bucket, txt = fetch_manifest(args.release)
    rows = list(csv.DictReader(io.StringIO(txt), delimiter="\t"))
    print("[%s] bucket=%s endpoints=%d" % (rel.upper(), bucket, len(rows)))

    mx = max(rows, key=lambda r: int(r["n_cases"]) + int(r["n_controls"]))
    print("  最大 N 端点: %s (%s) = %s cases + %s controls = %d"
          % (mx["phenocode"], mx["name"], mx["n_cases"], mx["n_controls"],
             int(mx["n_cases"]) + int(mx["n_controls"])))

    if args.grep:
        kw = args.grep.upper()
        hit = [r for r in rows
               if kw in r["phenocode"].upper() or kw in (r.get("name") or "").upper()]
        print("  匹配 %r 的端点 %d 个：" % (args.grep, len(hit)))
        for r in hit:
            print("    %-28s %-55s %8s / %8s = %9d  [%s]"
                  % (r["phenocode"], (r.get("name") or "")[:55], r["n_cases"],
                     r["n_controls"], int(r["n_cases"]) + int(r["n_controls"]),
                     r.get("category")))

    if args.out:
        with open(args.out, "w", encoding="utf-8", newline="") as f:
            f.write(txt)
        print("  已保存 ->", os.path.abspath(args.out))


if __name__ == "__main__":
    main()
