# Step 4 目录说明与复制分析手册

> 本文件于 **2026-10-07 20:20** 按实际运行状态重写（原版停留在「MVP 待到位」的中间态，已过时）。

## 目录约定

```
step4/
├── scripts/
│   ├── s4_finngen_replication.py   FinnGen R13 I9_PAD 复制（已完成）
│   ├── s4_fetch_mvp_pad.py         MVP PAD（dbGaP phs001672）抓取（分卷定位 + 分段续传）
│   ├── s4_pipeline.py              通用结局流水线（配置驱动；覆盖 8 个结局）
│   ├── s4_peek_members.py          分卷成员嗅探
│   ├── s4_figures.py               森林/漏斗/LOO/散点图
│   └── s4_figure_crossdataset.py   跨数据集汇总图
├── results/   72 个文件：mr_<tag>_<trait>.json ｜ harmonised_<tag>_<trait>.tsv ｜ loo_<tag>_<trait>.tsv
├── checks/    3 个文件：harmonise 逐 SNP 决策台账
├── figures/   94 张图：各结局 forest_/funnel_/loo_/scatter_ + summary_forest_*
└── logs/      13 个日志：s4_finngen_*.log ｜ pipeline_*.log
```

数据位置（实际）：

| 数据 | 路径 |
|---|---|
| FinnGen R13 I9_PAD | `data/replication/finngen_R13_I9_PAD.gz`（806 MB） |
| MVP trans-ancestry meta | `raw/mvp_pad/MVP.te.PAD.dbGAP.txt.gz`（665 MB） |
| 泛血管谱 5 个结局 | `data/panvascular/*.h.tsv.gz`（EBI GWAS Catalog harmonised） |

## 已完成的复制分析（截至 2026-10-07）

| 结局 | tag | 状态 | 报告 |
|---|---|---|---|
| FinnGen R13 I9_PAD | `finngen` | ✅ 完成 | `deliverables/Step4_FinnGen复制结果.md` |
| MVP PAD（trans-ancestry） | `mvp_meta` | ✅ 完成 | `deliverables/Step4_MVP复制结果.md` |
| CAD（Nikpay 2015，无重叠） | `cad_nouk` | ✅ 完成 | `deliverables/Step5_泛血管谱结果.{html,md}` |
| CAD（UKB Mbatchou 2021） | `cad_ukb` | ✅ 完成（重叠·参考） | 同上 |
| MI（UKB Donertas 2021） | `mi_ukb` | ✅ 完成（重叠·参考 + ‡ 尺度注记） | 同上 |
| IS-LAA（MEGASTROKE） | `is_laa` | ✅ 完成 | 同上 |
| IS-SV（MEGASTROKE） | `is_sv` | ✅ 完成 | 同上 |

复跑方式（幂等，成品存在即跳过）：

```
python step4/scripts/s4_pipeline.py --outcome <tag> --trait all
```

## MVP 数据现状

- **trans-ancestry meta 已就位并跑完**（`mvp_meta` 配置在 `s4_pipeline.py` 中已按实测列名填好：
  `EA ≡ Allele1`（100.0000% 实测）→ `oa` 取 `Allele2`；`SampleSize` 恒为 243,060 = 31,307 例 + 211,753 对照）。
- **分族裔层（EUR / AFR / HIS）不纳入本研究范围**——见下方「范围决策」。`s4_pipeline.py` 中三个 `mvp_eur / mvp_afr / mvp_his` 配置块已整块注释，仅保留说明，避免误读为"待跑"。

## 范围决策（2026-10-07 确立）

| 事项 | 决策 | 理由 |
|---|---|---|
| MVP 分族裔（EUR/AFR/HIS）分层验证 | **不做** | ① 分层病例数过小（AFR 5,373、HIS 1,925、EUR 24,009）；② 三套暴露工具均为 EUR 构建，对 AFR/HIS 属跨族裔外推，会引入不可控偏倚；③ trans-ancestry meta（31,307 例 / 211,753 对照，n_eff≈109,098）已综合三族裔信息，为更优的主分析 |
| 中介分析（血脂/BMI/吸烟/CRP） | **不做** | 研究设计以 **MVMR 独立性检验**（Step 7）回答「甲减效应是否独立于 TSH/FT4」；中介分析另需 4 组 GWAS 与独立分析模块，超出本文范围 |

## 关键口径（勿改）

- 对齐：rsID 优先、位置兜底（位置用面板坐标 = GRCh38）
- 等位基因：集合匹配 → **链互补翻转**（暴露文件混合链，非回文无歧义）→ 回文 |ΔEAF|>0.1 剔除、非回文 >0.2 剔除
- 统计口径：`step3/scripts/mrlib.py`（IVW 固定/随机 DL、Egger 规范形式 b_Y~b_X、加权中位数/众数 bootstrap、Q、MR-PRESSO 等价、方向性 R² 对比）
- 已在多套结局上完成回放验证，改配置不会改变口径

## 重叠提醒

| 暴露 | FinnGen R13 复制 | MVP 复制 | 泛血管谱（CAD/MI/IS） |
|---|---|---|---|
| 甲减 | ❌ 暴露含 FinnGen DF10（重叠） | ✅ 可用（暴露不含 MVP） | cad_nouk / is_laa / is_sv 可用；cad_ukb / mi_ukb 重叠·仅参考 |
| TSH | ✅ 已完成（阴性一致） | ✅ 可用 | 同上 |
| FT4 | ✅ 已完成（阴性一致） | ✅ 可用 | 同上 |
