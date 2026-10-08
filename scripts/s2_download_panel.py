# -*- coding: utf-8 -*-
"""
Step 2 前置：下载 1000G Phase 3 hg38 EUR 参考面板（Zenodo record 19068199, CC-BY-4.0）

设计要点
  1. Zenodo 对本机是「按连接限速」（实测单连接 ~12 KB/s，32 连接聚合 ~426 KB/s），
     因此采用 **多线程 + HTTP Range 分块** 下载。
  2. 每个分块单独落盘为 .parts/<文件名>/<序号>.part，**天然断点续传**：
     下次重跑只要 .part 大小正确就跳过，中断也不会丢进度。
  3. 全部下载完才拼接成正式文件并做 **官方 md5 校验**，校验不过不落正式文件。
  4. 支持把「顺序下载到一半的旧文件」自动接管为第 0 块。

用法
  python s2_download_panel.py                 # 下载 + 校验
  python s2_download_panel.py --verify-only   # 只校验
  python s2_download_panel.py --workers 48    # 调并发
"""
import os, sys, time, json, hashlib, threading, queue, argparse
import urllib.request, urllib.error

RECORD = "19068199"
URL = "https://zenodo.org/api/records/%s/files/%s/content" % (RECORD, "%s")
OUTDIR = r"D:/WorkSpace_Study/pad_thyroid_mr/ref/1kg_eur_hg38"
PARTS = os.path.join(OUTDIR, ".parts")

# (文件名, 字节数, 官方 md5) —— 全部来自 Zenodo API records/19068199
FILES = [
    ("allchr.EUR.biallelicsnps.cbb.fam",     7350, "4ddff9ee53ed3972cbfb317f8c17c5ed"),
    ("allchr.EUR.biallelicsnps.cbb.bim", 247454327, "d2c30f66bff12ea85d0bff54f2135a24"),
    ("allchr.EUR.biallelicsnps.cbb.bed", 644660327, "5e4ba32a0b1baedc50f0695003eee0e9"),
]

CHUNK = 4 * 1024 * 1024          # 每块 4 MiB
OPEN_PROXY = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def md5_of(path, bufsize=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def http_range(url, start, end, retries=8, timeout=120):
    """取 [start, end] 闭区间的字节；返回 bytes。"""
    last = None
    for a in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Range": "bytes=%d-%d" % (start, end),
                "User-Agent": "panel-downloader/1.0",
            })
            with OPEN_PROXY.open(req, timeout=timeout) as r:
                if r.status not in (200, 206):
                    raise IOError("HTTP %s" % r.status)
                return r.read()
        except Exception as e:                       # noqa
            last = e
            time.sleep(min(2 * (a + 1), 20))
    raise IOError("range %d-%d failed: %s" % (start, end, last))


def seed_legacy_part(name, size):
    """把旧的顺序下载残片接管为第 0 块。"""
    fin = os.path.join(OUTDIR, name)
    p0 = os.path.join(PARTS, name, "000000.part")
    if os.path.exists(fin) and not os.path.exists(p0):
        cur = os.path.getsize(fin)
        if 0 < cur < size:
            os.makedirs(os.path.dirname(p0), exist_ok=True)
            os.replace(fin, p0)
            print("   [seed] %s 的 %d 字节已接管为第 0 块" % (name, cur))


def download_one(name, size, expect_md5, workers):
    os.makedirs(os.path.join(PARTS, name), exist_ok=True)
    nchunk = (size + CHUNK - 1) // CHUNK
    url = URL % name

    def csize(i):
        return min(CHUNK, size - i * CHUNK)

    todo = []
    for i in range(nchunk):
        p = os.path.join(PARTS, name, "%06d.part" % i)
        if os.path.exists(p) and os.path.getsize(p) == csize(i):
            continue
        todo.append(i)

    done0 = nchunk - len(todo)
    print("   %s : %d 块，已完成 %d，待下 %d" % (name, nchunk, done0, len(todo)))
    if not todo:
        return True

    q = queue.Queue()
    for i in todo:
        q.put(i)
    lock = threading.Lock()
    state = {"done": done0, "bytes": done0 * CHUNK}

    def worker():
        while True:
            try:
                i = q.get_nowait()
            except queue.Empty:
                return
            p = os.path.join(PARTS, name, "%06d.part" % i)
            need = csize(i)
            # 块内续传
            have = os.path.getsize(p) if os.path.exists(p) else 0
            while have < need:
                end = min(have + CHUNK, need) - 1
                data = http_range(url, i * CHUNK + have, i * CHUNK + end)
                mode = "ab" if have else "wb"
                with open(p, mode) as f:
                    f.write(data)
                have += len(data)
                with lock:
                    state["bytes"] += len(data)
            with lock:
                state["done"] += 1
            q.task_done()

    ts = [threading.Thread(target=worker, daemon=True) for _ in range(workers)]
    t0 = time.time()
    for t in ts:
        t.start()
    last = 0
    while any(t.is_alive() for t in ts):
        time.sleep(10)
        el = time.time() - t0
        b = state["bytes"]
        sp = (b - last) / 10.0
        last = b
        print("     %s  %d/%d 块  %.1f/%.1f MB  %.0f KB/s  ETA %.1f min"
              % (name.split(".")[-2], state["done"], nchunk,
                 b / 1048576, size / 1048576, sp / 1024,
                 (size - b) / max(sp, 1) / 60))
    for t in ts:
        t.join()

    # 拼接
    final = os.path.join(OUTDIR, name)
    tmp = final + ".asm"
    with open(tmp, "wb") as out:
        for i in range(nchunk):
            p = os.path.join(PARTS, name, "%06d.part" % i)
            if not (os.path.exists(p) and os.path.getsize(p) == csize(i)):
                print("   [FAIL] 缺块 %d" % i)
                out.close()
                os.remove(tmp)
                return False
            with open(p, "rb") as f:
                while True:
                    b = f.read(1 << 22)
                    if not b:
                        break
                    out.write(b)
    got = os.path.getsize(tmp)
    if got != size:
        print("   [FAIL] 拼接后大小 %d != %d" % (got, size))
        os.remove(tmp)
        return False
    if expect_md5:
        real = md5_of(tmp)
        if real != expect_md5:
            print("   [FAIL] md5 %s != 官方 %s" % (real, expect_md5))
            os.remove(tmp)
            return False
        print("   [OK] md5 %s 与官方一致" % real)
    os.replace(tmp, final)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--verify-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(PARTS, exist_ok=True)

    ok = True
    for name, size, md5 in FILES:
        final = os.path.join(OUTDIR, name)
        print(">>>", name)
        if os.path.exists(final) and os.path.getsize(final) == size:
            real = md5_of(final)
            if real == md5:
                print("   [OK] 已完整且 md5 一致")
                continue
            print("   [WARN] 大小对但 md5 不符，重新下载")
            os.remove(final)
        if a.verify_only:
            print("   [MISSING] 尚未完整")
            ok = False
            continue
        seed_legacy_part(name, size)
        if not download_one(name, size, md5, a.workers):
            ok = False
    print("\n结论：", "全部就绪 ✅" if ok else "仍有未完成项 ❌")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
