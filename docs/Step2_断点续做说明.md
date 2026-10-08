# Step 2 断点续做说明

> 记录时间：2026-10-06 22:10 ｜ 状态：**已安全暂停，可无限期中断，重跑自动续上**

---

## 一、一句话结论

Step 2 的**算法链路已全部写好并通过语法自检**，唯一卡点是「参考面板下载慢」。
下次只需一条命令：

```bash
bash D:/WorkSpace_Study/pad_thyroid_mr/run_step2.sh
```

它会自动：续传面板 → 校验官方 md5 → 抽查坐标版本 → 改写 bim ID → 跑 PLINK clumping → 输出工具变量清单。
**中断任意次都没关系，已下载的分块会被跳过。**

其它子命令：
```bash
bash run_step2.sh status     # 看断点状态（面板完成度、候选数、已有交付）
bash run_step2.sh download   # 只续传面板
bash run_step2.sh clump      # 只跑 clumping
```

---

## 二、已经做好的事（不需要重做）

### 1. 候选工具变量已提取完毕 ✅
`work/cand_{hypo,tsh,ft4}.tsv`（脚本 `scripts/s2a_extract.py`）

| 性状 | p<5e-8 总数 | 其中常染色体双等位 SNP |
|---|---|---|
| 甲减 Hypothyroidism | 16,247 | **12,782** |
| TSH | 17,246 | **13,731** |
| FT4 | 1,690 | **1,474** |

> 被排除的主要是非 SNP 变异（indel / 多碱基 / 等位基因相同），已写入 `cand_*.diag.txt`。

### 2. PLINK 已就位 ✅
`tools/plink.exe` = **PLINK v1.90b7 64-bit (16 Jan 2023)**，来自
`https://s3.amazonaws.com/plink1-assets/plink_win64_20230116.zip`（本机可直连）。

### 3. 参考面板已选定并拿到官方 md5 ✅
- 来源：**Zenodo record 19068199 —— "1000G Phase 3 hg38 EUR plink files"**（CC-BY-4.0，DOI 10.5281/zenodo.19068199）
- **hg38 版本**（与我们的 GWAS 同为 GRCh38，可直接按位置匹配，无需 liftover）
- 294 个 EUR 样本、双等位 SNP 子集
- 官方校验值已存档 `ref/1kg_eur_hg38/OFFICIAL_MD5.txt`：

| 文件 | 字节数 | 官方 md5 | 本地现状 |
|---|---|---|---|
| `allchr.EUR.biallelicsnps.cbb.fam` | 7,350 | `4ddff9ee53ed3972cbfb317f8c17c5ed` | **已完整，md5 已核对一致 ✅** |
| `allchr.EUR.biallelicsnps.cbb.bim` | 247,454,327 | `d2c30f66bff12ea85d0bff54f2135a24` | 3,411,968（1.4%） |
| `allchr.EUR.biallelicsnps.cbb.bed` | 644,660,327 | `5e4ba32a0b1baedc50f0695003eee0e9` | 3,096,576（0.5%） |

现有残片会被新下载器**接管为第 0 块并接着往下续**，不会浪费。

### 4. 全流程脚本已写好并通过语法自检 ✅
- `scripts/s2_download_panel.py` —— 32 线程分块并行下载 + 分块级断点续传 + 官方 md5 校验
- `scripts/s2_clump.py` —— 改写 bim ID → 生成 clump 输入 → 调 PLINK → 解析 → F/R² 统计 → 输出

---

## 三、为什么慢（已排查到底，下次不用再摸一遍）

本机出口带宽本身不差，问题出在 **Zenodo 按连接限速**：

| 目标 | 单连接实测 | 结论 |
|---|---|---|
| Zenodo（面板所在） | **~12 KB/s** | 6 连接 72 KB/s、16 连接 223 KB/s、**32 连接 426 KB/s**（近似线性） |
| EBI ftp（1000genomes） | **942 KB/s** | 快 13 倍，但只有 b37，无 hg38 EUR PLINK 面板 |
| Broad GCS 桶 | **757 KB/s** | 快，但桶内 `1000G.EUR.QC.*` **只有 .bim，没有 .bed/.fam**，无法算 LD |

所以：**只能走 Zenodo 多连接并行**。按 426 KB/s 估算，剩余 885 MB 约需 **35–40 分钟**（单连接则要 3.4 小时）。

### 已经排掉、不用再试的死路（省得下次重复踩）
- `gwas-api.mrcieu.ac.uk`（IEU/OpenGWAS 在线 clumping）→ **502**；`api.opengwas.io/api/ld/clump` → **502**，且 HTTPS 证书已过期
- `ctg.cncr.nl/.../g1000_eur.zip`（MAGMA 面板）→ 连接被重置（000）
- Dropbox 全部链接（PLINK2 官方 1000G 资源都挂 Dropbox）→ 不可达
- GitHub / raw.githubusercontent → 不可达
- PLINK2 的 s3 下载路径（`plink2-assets/alpha6/...`）→ 404（PLINK **1.9** 的 s3 路径可用，已下好）
- Zenodo 12725512「RAISS hg38 面板」→ 只有 .bim，无基因型
- Zenodo 19127594「PRScope 参考包」→ 15.7 GB，太大
- Broad `LDSCORE/1000G_Phase3_plinkfiles.tgz` → 404

---

## 四、本次新增的技术决定（下次直接沿用）

1. **面板 SNP ID 统一改写为 `chr:pos`**（原文件备份为 `*.bim.rsid.orig`）。
   好处：完全绕开 rsID 版本差异、重名、`.` 缺失 ID；因为面板与 GWAS 同为 hg38，按位置匹配最稳。
2. **不需要 numpy / pandas**。clumping 交给 PLINK（C 实现，快），统计与出表用标准库 + openpyxl（本机已装）。
   —— 顺带规避了本机 PyPI 走国外源极慢、uv 装包屡屡超时的问题。
3. **R² 口径**（会在交付物里明确标注）：
   - `F = z² = (beta/se)²`，逐位点精确值，MR 惯例要求 F > 10
   - `R² = z²/(z² + n − 2)`（TwoSampleMR `get_r_from_pn` 的代数等价式）
   - 甲减为二值性状，用 `n_eff = 4·n_case·n_ctrl/(n_case+n_ctrl)` 代入，属**观测尺度近似**
   - 整体 `F = (R²/(1−R²))·((n_eff − k − 1)/k)`
4. **clumping 参数**：`r² < 0.001`，窗口 `10,000 kb`，p1 = 5×10⁻⁸（MR 领域通用口径）。
5. 回文位点（A/T、C/G）与 MAF<0.01 位点**本步不剔除**，只在清单里高亮标记，留到 Step 3 协调时统一处理。

---

## 五、本机 Python 环境 —— 已闭环，不再是遗留项（2026-10-06 22:12 更新）

原先挂着的后台任务 `uv pip install pandas numpy scipy openpyxl` 已**取消**（不是暂停），原因：

- 它在 PyPI 上**挂了 27 分钟、零输出**，属于「卡死」而非「推进中」；
- **Step 2 已经刻意不依赖 numpy / pandas**（clumping 交给 PLINK，统计用标准库 + openpyxl），所以没有任何东西等它；
- 中断是安全的：`uv` 安装是原子式的，取消后 venv **无半装残包**，已逐包核验。

### 实测出来的正确装包方式（下次照这个来）

本机走国外 PyPI 极慢，**必须加国内镜像**。实测对比：

| 方式 | 结果 |
|---|---|
| 默认 PyPI（`uv pip install`） | 27 分钟仍无输出，被杀 |
| 清华镜像 | **numpy 27 秒装完** |

正确命令（`--link-mode=copy` 用于消除跨盘 hardlink 警告）：

```bash
"C:/Users/Chenmx7/.local/bin/uv.exe" pip install \
  --index-url https://pypi.tuna.tsinghua.edu.cn/simple --link-mode=copy \
  --python "D:/WorkSpace_Study/.venv/Scripts/python.exe" <包名>
```

备选镜像：`https://mirrors.aliyun.com/pypi/simple/`、`https://mirrors.cloud.tencent.com/pypi/simple/`。

### venv 当前清单（已核对，全部可用）

| 包 | 版本 |
|---|---|
| numpy | 2.5.3 |
| pandas | 3.0.6 |
| scipy | 1.18.1 |
| openpyxl | 3.1.5 |
| xlrd | 2.0.2 |
| python-docx | 1.2.0 |
| pymupdf | 1.28.2 |

---

## 六、下次恢复的具体动作

```bash
# 1）看当前断点
bash D:/WorkSpace_Study/pad_thyroid_mr/run_step2.sh status

# 2）一条命令跑完（面板续传 + clumping + 出清单）
bash D:/WorkSpace_Study/pad_thyroid_mr/run_step2.sh
```

预计产出：
- `deliverables/Step2_工具变量清单.xlsx` —— 汇总表 + 三套工具明细（回文/低 MAF 位点自动高亮）
- `deliverables/Step2_工具变量清单.md` —— 汇总 + 各性状前 15 个工具 + 口径说明
- `work/instruments_{hypo,tsh,ft4}.tsv` —— 完整工具表
- `work/clump_{hypo,tsh,ft4}.clumped` —— PLINK 原始输出（留档备查）

**自检点**：甲减工具数应落在 **300–400** 区间（原文报告 350 个独立位点），若明显偏离需回头查面板坐标版本。

---

## 七、Step 3 之前需要你拍板的一件事（不影响现在暂停）

结局主分析用哪一套：
- **GCST90018890**（Sakaue 2021，7,114 例）——与 2024 年那篇 MR 同源，可直接对照复现其 OR 1.47（95%CI 0.64–3.36），**已下载好，零成本**；
- **MVP PAD**（Klarin 2019，31,307 例，n_eff≈109,098）——精度约 ×3.9，但需走 dbGaP 受控申请；
- FinnGen R11 `I9_ARTEMBTHRLOW` 作独立复制集，随时可下（0.8 GB）。
