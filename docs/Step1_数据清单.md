# Step 1 数据获取清单 — 甲状腺功能 / 甲减 → PAD 的 MR

> 生成时间：2026-10-06 20:35 ｜ 所有 GWAS 登记号、样本量、FTP 路径均**当场实测验证**（HTTP 200 + 目录列表 + md5）
> 项目目录：`D:/WorkSpace_Study/pad_thyroid_mr/`

---

## 0. 目录结构

```
pad_thyroid_mr/
└── data/
    ├── exposure/     # 暴露：甲减 / TSH / FT4（harmonised，GRCh38，GWAS-SSF 标准列）
    ├── outcome/      # 结局：PAD 及泛血管谱
    └── meta/         # md5 校验文件、数据来源记录
```

---

## 1. 暴露（三套，均来自同一篇 Nat Genet 2025，PMID 41238958）

同一篇论文同时提供三套汇总统计 → 天然支持**三重暴露三角验证**。

| 角色 | 性状 | GWAS Catalog | 样本量 | harmonised 文件 | 大小 |
|---|---|---|---|---|---|
| 暴露①（主） | Hypothyroidism（二值） | **GCST90572791** | 98,805 + 14,588 例 / 948,186 + 117,082 对照（合 113,393 / 1,065,268） | `GCST90572791.h.tsv.gz` | 230,912,838 B (0.23 GB) |
| 暴露② | TSH（连续） | **GCST90572789** | 482,873 欧洲 | `GCST90572789.h.tsv.gz` | 204,321,482 B (0.20 GB) |
| 暴露③ | Thyroxine / FT4（连续） | **GCST90572790** | 191,449 欧洲 | `GCST90572790.h.tsv.gz` | 184,291,932 B (0.18 GB) |

**原文 Data availability 原文引用**：
> "GWAS summary statistics from the meta-analysis of hypothyroidism (excluding 23andMe), TSH and T4 are publicly available at the GWAS Catalog under accession IDs GCST90572791, GCST90572789 and GCST90572790."
> "The corresponding hypothyroidism PRS … is available at the PGS Catalog under accession ID PGS005218."
> 23andMe 部分需向 23andMe 申请；UKB 个体数据需申请；FinnGen 需注册。

**FTP 基址**：`https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90572001-GCST90573000/<GCST>/harmonised/`

**md5（harmonised，已存档 `data/meta/`）**
```
8b47a15de4e2220898809146777084df  GCST90572791.h.tsv.gz   ← 甲减
59689c1b938681871eb863a0bbffee8b  GCST90572789.h.tsv.gz   ← TSH
dda8c77bd540f223bcaa1734f782c226  GCST90572790.h.tsv.gz   ← FT4
```

> 说明：同目录还有非 harmonised 版本（`GCST90572791.tsv.gz` 704 MB 等），体积约为 harmonised 的 3 倍，已弃用，只取 harmonised。

---

## 2. 结局：PAD

| 角色 | 来源 | GWAS Catalog | 例 / 对照 | 有效样本量 n_eff* | 文件 |
|---|---|---|---|---|---|
| **主分析** | Sakaue 2021 *Nat Genet*（PMID 34594039） | **GCST90018890** | 7,114 / 475,964 | ≈ 28,037 | `GCST90018890.h.tsv.gz` (636,990,393 B) |
| **复制集** | FinnGen **R13**（公共桶直下，2026 新发布） | — | 端点 **`I9_PAD`**（22,244 例 / 460,490 对照） | — | `finngen_R13_I9_PAD.gz` (805,766,381 B) |
| **待申请（更强）** | Klarin 2019 *Nat Med* MVP（PMID 31285632） | 需查 dbGaP | 31,307 / 211,753 | ≈ 109,098 | 受控访问，需申请 |

\* n_eff = 4·n_case·n_control/(n_case+n_control)，用于功效估算

> **FinnGen 端点说明 —— 含 2026-10-07 事后更正，以「最终采用」为准**
> - 原始记录（2026-10-06 撰写时）：R11/R12 公共桶中**没有** `I9_PAD` 这一名（实测 404），当时据此选用 `I9_ARTEMBTHRLOW`。
> - **更正（2026-10-07，Step 4 阶段实测）**：`I9_ARTEMBTHRLOW` 经核实实为「下肢动脉栓塞与血栓形成」(I74.3)，R11 仅 1,170 例 / R13 1,400 例 → **表型不符且样本量过小**；备选 `I9_ATHSCLE` 官方定义为「Atherosclerosis, excluding cerebral, coronary and PAD」→ **明确排除 PAD**。
> - **最终采用**：FinnGen **R13** 公开桶的 **`I9_PAD`**（Peripheral artery disease，**22,244 例 / 460,490 对照**，manifest 实测），落地文件 `data/replication/finngen_R13_I9_PAD.gz`（805,766,381 B，gzip 完整性已校验）。纠错全过程见 `Step4_FinnGen复制结果.md` §1。

---

## 3. ⚠ 关键情报修正（直接影响研究设计，务必注意）

核查 2024 年那篇甲减→ASCVD 的 MR（PMID 39512367）**表 1** 后确认：

> **它用的 PAD 结局 GWAS 是 `ebi-a-GCST90018890`（Sakaue 2021），仅 7,114 例 / 475,964 对照**，
> **不是** MVP 的 31,307 例。结果：**OR 1.47，95%CI 0.64–3.36，P = 0.358**。

| | 2024 旧研究 | 本方案 |
|---|---|---|
| 工具 SNP 数 | 122 | **350**（×2.9） |
| PAD 例数 | 7,114 | 7,114（同源，可比）／ **31,307**（MVP，×4.4） |
| 有效样本量 | 28,037 | 28,037 ／ **109,098**（×3.9） |
| PAD 结果 | OR 1.47（**0.64–3.36**）P=0.358 — CI 极宽，明显功效不足 | 待测 |

**结论**：
- 若主分析沿用 `GCST90018890`（与旧研究同源、可直接比较），精度提升主要来自**工具数 ×2.9**；
- 若以 MVP PAD 为主分析、`GCST90018890` 为复制，则**工具数 ×2.9 × 结局样本 ×3.9 ≈ 11 倍**的精度提升，CI 将显著收窄；
- 旧的 PAD 阴性结果很可能只是**功效不足**——这正是本研究的立据。（精确功效将在 Step 2 用真实数据重算，不再用解析近似。）

**该文引用的全部结局 GWAS（表 1，均公开可下）**：

| Trait | GWAS ID | 例 / 对照 | 作者 / PMID |
|---|---|---|---|
| CAD | GCST90013864 | 29,339 / 322,724 | Mbatchou J, 34017140 |
| AP | GCST90018793 | 30,025 / 440,906 | Elsworth B, 34594039 |
| MI | GCST90038610 | 11,081 / 473,517 | Dönertaş HM, 33959723 |
| IS | GCST90018864 | 11,929 / 472,192 | Sakaue S, 34594039 |
| IS-LAA | GCST005840 | 4,373 / 406,111 | Malik R, 29531354 |
| IS-SV | GCST005841 | 5,386 / 192,662 | Malik R, 29531354 |
| **PAD** | **GCST90018890** | **7,114 / 475,964** | Sakaue S, 34594039 |

---

## 4. 泛血管谱结局（Step 5 用，同一套工具平行检验）

FTP 已验证可下（未全部下载，按需取）：

| 结局 | GCST | 文件 | 分组目录 |
|---|---|---|---|
| CAD | GCST90013864 | `GCST90013864_buildGRCh38.tsv.gz` | GCST90013001-GCST90014000 |
| MI | GCST90038610 | `GCST90038610_buildGRCh37.tsv`（未压缩） | GCST90038001-GCST90039000 |
| IS | GCST90018864 | `GCST90018864_buildGRCh37.tsv.gz` | GCST90018001-GCST90019000 |
| IS-LAA | GCST005840 | `MEGASTROKE.3.LAS.TRANS.out` + `harmonised/` | GCST005001-GCST006000 |
| IS-SV | GCST005841 | 同上（MEGASTROKE） | GCST005001-GCST006000 |

---

## 5. 中介变量 GWAS —— **已出研究范围，不再使用**

> **范围决策（2026-10-07）**：本研究**不做中介分析**。
> Step 6/7 实际执行的是**反向 MR + MVMR**（以 TSH / FT4 作为共暴露，检验甲减效应是否独立于两者），
> 而非以血脂 / BMI / 吸烟 / CRP 为中介的正式中介分析。
> 理由：① MVMR 已直接回答「甲减致动脉粥样硬化效应是否独立于 TSH/FT4」这一核心问题；
> ② 正式中介分析需另下 4 组 GWAS 并新建「两步 MR + 中介比例」模块，超出本文范围。
> 下表为 Step 1 规划期的候选清单，**保留备查**；若将来补做中介分析，可从此表取源。

| 中介 | 来源 | PMID / DOI |
|---|---|---|
| 血脂（LDL-C/TC/TG/HDL） | TOPMed/GLGC, *Nature* 2021 | 34887591 / 10.1038/s41586-021-04064-3 |
| BMI | Yengo 2018 *Hum Mol Genet*（~700,000） | 30124842 / 10.1093/hmg/ddy271 |
| 吸烟 | GSCAN, *Nat Genet* 2019（~1.2 M） | 30643251 / 10.1038/s41588-018-0307-5 |
| CRP | *Nat Commun* 2025（513,273） | 40274876 / 10.1038/s41467-025-59155-w |

---

## 6. 下载状态 —— ✅ 全部完成并校验（2026-10-06 21:30）

| 文件 | 字节数 | 目标字节数 | md5 校验 | 变异数 | p<5×10⁻⁸ SNP 数 |
|---|---|---|---|---|---|
| hypo_GCST90572791.h.tsv.gz | 230,912,838 | 230,912,838 ✅ | **8b47a15de4e2220898809146777084df ✅ 与官方一致** | **7,344,354** | **16,247** |
| tsh_GCST90572789.h.tsv.gz | 204,321,482 | 204,321,482 ✅ | **59689c1b938681871eb863a0bbffee8b ✅ 与官方一致** | **6,151,997** | **17,246** |
| ft4_GCST90572790.h.tsv.gz | 184,291,932 | 184,291,932 ✅ | **dda8c77bd540f223bcaa1734f782c226 ✅ 与官方一致** | **5,901,316** | **1,690** |
| pad_GCST90018890.h.tsv.gz | 636,990,393 | 636,990,393 ✅ | **676209242bcb8d15fbbc3be1d71382f8 ✅ 与 meta.yaml 一致** | **23,492,920** | （结局，无需） |

> PAD 的 md5 官方参考来自 `GCST90018890.h.tsv.gz-meta.yaml` 的 `data_file_md5sum` 字段（其 `md5sum.txt` 只有文件名没有哈希，不可用）。

**表头（实测，供 Step 2 直接对接）**
- 暴露三套（GWAS-SSF v1.0，列序一致）：
  `chromosome | base_pair_location | effect_allele | other_allele | beta | standard_error | effect_allele_frequency | p_value | rsid | rs_id | n | hm_coordinate_conversion | hm_code | variant_id`
  → p_value 在第 **8** 列；variant_id 形如 `1_10518_T_C`
- 结局 PAD（pre-GWAS-SSF，**列序不同**）：
  `chromosome | base_pair_location | effect_allele | other_allele | beta | standard_error | effect_allele_frequency | p_value | variant_id | hm_coordinate_conversion | hm_code | rsid`
  → **无 n 列**；rsid 在最后一列（第 12 列）；读入时需显式映射列名
- meta.yaml 确认：`analysis_software: REGENIE`、`imputation_software: IMPUTE4|Eagle|Beagle|Graphtyper`、`sex: combined`、EFO_0004705

---

## 7. Step 1 完成情况与进入 Step 2 的衔接

**已完成**
- [x] 三套暴露 + 一套结局下载完成，字节数与官方完全一致
- [x] 四个文件 md5 全部校验通过（暴露三套 vs 官方 md5sum.txt；PAD vs meta.yaml）
- [x] 表头列名确认（暴露 GWAS-SSF；PAD 为 pre-GWAS-SSF，列序不同需映射）
- [x] 变异数统计：hypo 7.34M / tsh 6.15M / ft4 5.90M / pad 23.49M
- [x] 工具变量供应量确认：hypo 16,247 个、tsh 17,246 个、ft4 1,690 个 p<5e-8 SNP
      - **事后更正（2026-10-07）**：按 r²<0.001、10 Mb 窗口、EUR 1000G 标准 clump 后，**实际独立工具数为 甲减 222 / TSH 258 / FT4 65**（见 `deliverables/Step2_工具变量清单.md`、`deliverables/Step9_统计功效.md`）。原「≈300–400」为撰写时的预估，偏高；主因是未计入 10 Mb 长程 LD 下的合并强度。

**进入 Step 2 的衔接**
- 结局待决：主分析用 `GCST90018890`（与旧研究同源，可直接比较，已下好）；
  若要 MVP PAD（31,307 例，更强）需走 dbGaP 申请；FinnGen `I9_ARTEMBTHRLOW`（0.80 GB）作复制集，按需下载。
- Step 2 动作：从暴露三套分别取 p<5e-8、LD clump（r²<0.001, 10Mb，EUR 1000G）、F>10 → 生成三套工具清单。
- 注意：PAD 文件**无 n 列**（case-control），TwoSampleMR 读入时需显式指定列映射。
