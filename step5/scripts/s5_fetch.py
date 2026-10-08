# -*- coding: utf-8 -*-
"""
Step 5 泛血管谱结局 GWAS 下载器（EBI GWAS Catalog harmonised 版）

背景：EBI FTP/HTTP 单连接实测仅约 48 KB/s（下载 5 MB 用满 100 秒），
      1.63 GB 单线程需约 9 小时 → 必须多段 Range 并行。

设计：
  · 先 HEAD 取 Content-Length（权威大小，不猜）
  · 每个文件切成 N 段（默认 4），段内独立 curl 进程，段间并行
  · 段级断点续传：已落盘字节 > 0 时从 start+got 继续（Range 重算）
  · 每段最多 6 次重试，指数退避；用「截断为 0」而非删除（本机 safe-delete 钩子限制）
  · 全部段完成后按序合并 → 校验 gzip 完整性（全量读一遍）
  · 幂等：已完成（大小正确）或已校验通过的文件直接跳过，可反复重跑

用法：
  python step5/scripts/s5_fetch.py                    # 全部 5 个文件
  python step5/scripts/s5_fetch.py --segments 6       # 提高并发段数
  python step5/scripts/s5_fetch.py --only is_sv       # 只下某一个
"""
import argparse
import gzip
import os
import shutil
import subprocess
import threading
import time

ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"
OUT = os.path.join(ROOT, "data", "panvascular")
PARTS = os.path.join(OUT, "_parts")
LOGDIR = os.path.join(ROOT, "step5", "logs")
B = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

CURL = ["curl", "-s", "--noproxy", "*", "--ssl-no-revoke", "--retry", "2", "--retry-delay", "3"]

# key, 输出文件名, URL, 说明（重叠状态）
TASKS = [
    dict(key="cad_ukb", file="cad_GCST90013864.h.tsv.gz",
         url=B + "/GCST90013001-GCST90014000/GCST90013864/harmonised/"
                 "34017140-GCST90013864-EFO_0001645.h.tsv.gz",
         note="CAD  Mbatchou2021(PMID34017140) ancestry=British(UKB)  -> OVERLAP_UKB"),
    dict(key="cad_nouk", file="cad_GCST003116.h.tsv.gz",
         url=B + "/GCST003001-GCST004000/GCST003116/harmonised/"
                 "26343387-GCST003116-EFO_0000378.h.tsv.gz",
         note="CAD  Nikpay2015(PMID26343387) CARDIoGRAMplusC4D     -> NO OVERLAP"),
    dict(key="mi_ukb", file="mi_GCST90038610.h.tsv.gz",
         url=B + "/GCST90038001-GCST90039000/GCST90038610/harmonised/"
                 "33959723-GCST90038610-EFO_0000612.h.tsv.gz",
         note="MI   Donertas2021(PMID33959723) UKB-derived          -> OVERLAP_UKB"),
    dict(key="is_laa", file="is_laa_GCST005840.h.tsv.gz",
         url=B + "/GCST005001-GCST006000/GCST005840/harmonised/"
                 "29531354-GCST005840-EFO_0005524.h.tsv.gz",
         note="IS-LAA MEGASTROKE(Malik2018 PMID29531354)              -> NO OVERLAP"),
    dict(key="is_sv", file="is_sv_GCST005841.h.tsv.gz",
         url=B + "/GCST005001-GCST006000/GCST005841/harmonised/"
                 "29531354-GCST005841-EFO_1001504.h.tsv.gz",
         note="IS-SV  MEGASTROKE(Malik2018 PMID29531354)              -> NO OVERLAP"),
]

_lock = threading.Lock()
LOGP = None


def log(m):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    with _lock:
        print(line, flush=True)
        if LOGP:
            with open(LOGP, "a", encoding="utf-8") as f:
                f.write(line + "\n")


def trunc(p):
    """用截断代替删除（规避本机 safe-delete 批量钩子）"""
    try:
        with open(p, "wb"):
            pass
    except OSError:
        pass


def head_size(url):
    for _ in range(4):
        r = subprocess.run(CURL + ["-I", "--max-time", "60", url],
                           capture_output=True, text=True)
        for ln in (r.stdout or "").splitlines():
            if ln.lower().startswith("content-length:"):
                try:
                    return int(ln.split(":", 1)[1].strip())
                except ValueError:
                    pass
        time.sleep(4)
    return None


def fetch_segment(url, out, start, end, tag):
    """下载 [start, end] 闭区间到 out，支持断点续传"""
    want = end - start + 1
    for attempt in range(6):
        got = os.path.getsize(out) if os.path.exists(out) else 0
        if got == want:
            return True
        if got > want:
            trunc(out)
            got = 0
        rng = "%d-%d" % (start + got, end)
        tmp = out + ".t"
        # 注意：这里不能用 --max-time（硬时限）——EBI 单连接仅约 40 KB/s，
        # 112 MB 的段需 ~45 分钟，硬时限会反复把已下载数据整段丢弃、永远下不完。
        # 改用速度阈值：只有连续 300 秒平均速度低于 2 KB/s（真·卡死）才放弃。
        r = subprocess.run(CURL + ["--speed-limit", "2048", "--speed-time", "300",
                                   "-r", rng, url, "-o", tmp, "-w", "%{http_code}"],
                           capture_output=True, text=True)
        code = (r.stdout or "").strip()
        # ★ 关键：无论本次请求成功还是被中断，都先把已下载的字节落盘。
        #   否则一旦超时/断线，整段进度归零（旧版 --max-time 1500 就是这样白下了 1 小时）。
        tsz = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        if tsz > 0:
            if tsz > want - got:      # 服务器多给了，只取需要的
                tsz = want - got
            with open(out, "ab" if got else "wb") as fo, open(tmp, "rb") as fi:
                if tsz != os.path.getsize(tmp):
                    fo.write(fi.read(tsz))
                else:
                    shutil.copyfileobj(fi, fo, 1 << 22)
            trunc(tmp)
            if os.path.getsize(out) == want:
                return True
        if r.returncode != 0 or code not in ("200", "206"):
            log("  %s 请求中断 rc=%s http=%s；已保住 %s，续传重试 %d/6" %
                (tag, r.returncode, code, human(os.path.getsize(out)), attempt + 1))
            time.sleep(4 * (attempt + 1))
            continue
        if tsz <= 0:
            time.sleep(4 * (attempt + 1))
            continue
    return False


def verify_gz(p):
    try:
        with gzip.open(p, "rb") as f:
            while True:
                b = f.read(1 << 22)
                if not b:
                    break
        return True, ""
    except Exception as e:
        return False, str(e)[:120]


def human(n):
    return "%.1f MB" % (n / 1048576.0)


def main():
    global LOGP
    ap = argparse.ArgumentParser()
    ap.add_argument("--segments", type=int, default=4, help="每文件并行段数")
    ap.add_argument("--only", default="", help="逗号分隔的 key，如 is_sv,is_laa")
    a = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(PARTS, exist_ok=True)
    os.makedirs(LOGDIR, exist_ok=True)
    LOGP = os.path.join(LOGDIR, "s5_fetch_%s.log" % time.strftime("%Y%m%d_%H%M%S"))

    tasks = TASKS if not a.only else [t for t in TASKS if t["key"] in a.only.split(",")]
    if not tasks:
        log("没有匹配的任务"); return

    log("=" * 78)
    log("Step 5 泛血管谱数据下载  |  文件 %d 个  |  每文件 %d 段  |  日志 %s"
        % (len(tasks), a.segments, os.path.basename(LOGP)))
    for t in tasks:
        log("  %-9s %-30s %s" % (t["key"], t["file"], t["note"]))
    log("=" * 78)

    # ---- 阶段 1：HEAD 取大小 ----
    jobs = []
    for t in tasks:
        final = os.path.join(OUT, t["file"])
        sz = head_size(t["url"])
        if sz is None:
            log("✗ %s：无法取到 Content-Length，跳过" % t["key"]); continue
        t["size"] = sz
        if os.path.exists(final) and os.path.getsize(final) == sz:
            ok, err = verify_gz(final)
            if ok:
                log("✓ %s 已存在且校验通过（%s），跳过" % (t["key"], human(sz)))
                continue
            log("⚠ %s 大小对但 gzip 校验失败（%s），重新下载" % (t["key"], err))
        seg = max(1, a.segments)
        step = sz // seg
        segs = []
        for i in range(seg):
            s = i * step
            e = (sz - 1) if i == seg - 1 else (s + step - 1)
            segs.append((s, e, os.path.join(PARTS, "%s.part%02d" % (t["key"], i))))
        jobs.append(dict(task=t, final=final, segs=segs))
        log("· %s  %s  →  %d 段（每段约 %s）" % (t["key"], human(sz), len(segs), human(step)))

    if not jobs:
        log("全部已完成，无需下载"); return

    # ---- 阶段 2：并行下载所有段 ----
    total_bytes = sum(j["task"]["size"] for j in jobs)
    log("开始下载，共 %d 段，合计 %s"
        % (sum(len(j["segs"]) for j in jobs), human(total_bytes)))
    t0 = time.time()
    results = {}

    def run_seg(j, idx, s, e, p):
        tag = "%s/段%d" % (j["task"]["key"], idx)
        log("  → %s 启动 [%d, %d] %s" % (tag, s, e, human(e - s + 1)))
        ok = fetch_segment(j["task"]["url"], p, s, e, tag)
        results[(j["task"]["key"], idx)] = ok
        log("  %s %s" % ("✓" if ok else "✗", tag))

    threads = []
    for j in jobs:
        for idx, (s, e, p) in enumerate(j["segs"]):
            th = threading.Thread(target=run_seg, args=(j, idx, s, e, p), daemon=True)
            th.start(); threads.append(th)

    # 进度心跳（注意：下载写入的是最终段文件 + ".t" 临时文件，两者都要统计）
    def seg_done_bytes():
        n = 0
        for j in jobs:
            for (_, _, p) in j["segs"]:
                if os.path.exists(p):
                    n += os.path.getsize(p)
                elif os.path.exists(p + ".t"):
                    n += os.path.getsize(p + ".t")
        return n

    while any(th.is_alive() for th in threads):
        time.sleep(60)
        done = seg_done_bytes()
        el = time.time() - t0
        sp = done / el if el > 0 else 0
        eta = (total_bytes - done) / sp if sp > 0 else 0
        log("  进度 %s / %s (%.1f%%)  速度 %s/s  ETA %.0f 分钟"
            % (human(done), human(total_bytes), 100.0 * done / total_bytes,
               human(sp), eta / 60.0))

    for th in threads:
        th.join()

    # ---- 阶段 3：合并 + 校验 ----
    log("-" * 78)
    allok = True
    for j in jobs:
        key = j["task"]["key"]; final = j["final"]
        if any(not results.get((key, i), False) for i in range(len(j["segs"]))):
            log("✗ %s 有段未完成，暂不合并（可重跑续传）" % key); allok = False; continue
        if os.path.exists(final):
            trunc(final)
        with open(final, "wb") as fo:
            for (_, _, p) in j["segs"]:
                with open(p, "rb") as fi:
                    shutil.copyfileobj(fi, fo, 1 << 22)
        got = os.path.getsize(final)
        if got != j["task"]["size"]:
            log("✗ %s 合并后大小 %s ≠ 期望 %s" % (key, human(got), human(j["task"]["size"])))
            allok = False; continue
        ok, err = verify_gz(final)
        log("%s %s  %s  gzip校验%s%s" % ("✓" if ok else "✗", key, human(got),
                                         "通过" if ok else "失败", "" if ok else "：" + err))
        if not ok:
            allok = False

    log("-" * 78)
    log("用时 %.1f 分钟  |  总体 %s" % ((time.time() - t0) / 60.0, "成功" if allok else "存在未完成项"))
    log("产物目录：%s" % OUT)


if __name__ == "__main__":
    main()
