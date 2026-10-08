# Analysis pipeline — run order, inputs and outputs

Every script is pure Python (numpy / matplotlib / openpyxl only) except `s2_clump.py`, which calls PLINK 1.9.

---

## Stage 0 — Reference panel and instrument selection (`scripts/`)

| # | Script | Purpose | Key inputs | Key outputs |
|---|---|---|---|---|
| 0.1 | `s2_download_panel.py` | Download 1000G Phase 3 GRCh38 EUR reference panel | Zenodo record 19068199 | `ref/1kg_eur_hg38/*.bed/.bim/.fam` |
| 0.2 | `s2a_extract.py` | Extract candidate instruments (p < 5×10⁻⁸) from the three exposure files | `data/exposure/*.h.tsv.gz` | candidate SNP lists |
| 0.3 | `s2_clump.py` | PLINK LD clumping (r² < 0.001, 10 Mb) → instrument list; compute z, F, R² per variant and the overall F | reference panel + candidate lists | `Step2` instrument lists (222 / 258 / 65 instruments) |
| 0.4 | `s2_diag.py` | Diagnostic: reference panel integrity and candidate match rate | panel `.bed` | console report |
| — | `fetch_finngen_manifest.py` | Retrieve any FinnGen release endpoint manifest from the public GCS bucket | release tag | `*_manifest.tsv` |
| — | `s_fig_results.py` | Publication figures for Steps 5–9 | `step*/results/*` | `step*/figures/*.png` |

---

## Stage 1 — Harmonisation and primary MR (`step3/`)

| # | Script | Purpose | Key inputs | Key outputs |
|---|---|---|---|---|
| 1.1 | `s3_harmonise.py` | Extract outcome SNPs and harmonise exposure × outcome (allele + strand complement handling) | exposure instruments, `data/outcome/pad_*` | `step3/results/harmonised_{hypo,tsh,ft4}.tsv`, `harmonise_summary.json` |
| 1.2 | `s3_mr_main.py` | **Primary MR**: Wald ratio + sign-consistency binomial; IVW (fixed + DerSimonian–Laird random); MR-Egger (weighted, with intercept); weighted median; weighted mode (Gaussian KDE); Cochran's Q; MR-PRESSO-equivalent global test + Bonferroni outlier detection + distortion test; leave-one-out; directionality (R² comparison) | `step3/results/harmonised_*.tsv` | `mr_main_<trait>.json`, `loo_<trait>.tsv`, `mr_main_summary.tsv` |
| 1.3 | `mrlib.py` | Shared MR library (normal/χ² p-values validated against tabulated critical values; IVW; Q; Egger; weighted median/mode) | — | (module) |
| 1.4 | `s3_figures.py` | Scatter, forest, funnel, leave-one-out figures | Step 3 results | `step3/figures/*.png` |

---

## Stage 2 — Independent replication and pan-vascular extension (`step4/`, `step5/`)

| # | Script | Purpose | Key inputs | Key outputs |
|---|---|---|---|---|
| 2.1 | `s4_fetch_mvp_pad.py` | Download MVP PAD summary statistics from dbGaP (multi-connection, chunked, resumable) | dbGaP phs001672 | `raw/mvp_pad/MVP.te.PAD.dbGAP.txt.gz` |
| 2.2 | `s4_peek_members.py` | Enumerate members of a solid `tar.gz` without downloading the whole archive; extract a single member | local `tar.gz` chunks | member listing / extracted file |
| 2.3 | `s4_pipeline.py` | **Generic outcome replication pipeline**: harmonise + MR for any configured outcome (FinnGen R13, MVP meta, CAD no-UKB, CAD UKB, MI, IS-LAA, IS-SV) × any trait | `data/replication/*`, `data/panvascular/*` | `step4/results/harmonised_<outcome>_<trait>.tsv`, MR JSONs |
| 2.4 | `s4_finngen_replication.py` | FinnGen R13 `I9_PAD` replication (traits hypo/TSH/FT4) | `data/replication/finngen_R13_I9_PAD.gz` | `step4/results/harmonised_finngen_*.tsv` |
| 2.5 | `s4_figures.py` | FinnGen replication figures | Step 4 results | `step4/figures/*.png` |
| 2.6 | `s4_figure_crossdataset.py` | Cross-dataset forest plot: primary + two replications | Step 3/4 results | `step4/figures/summary_forest_*.png` |
| 2.7 | `s4_onepager.py` | Steps 3+4 one-page summary (numbers read directly from result files) | `step3|4/results/*` | `Step3-4_结果一页纸.md/.html` |
| 2.8 | `s5_fetch.py` | Download pan-vascular outcome GWAS (GWAS Catalog harmonised) | GWAS Catalog FTP | `data/panvascular/*` |
| 2.9 | `s5_onepager.py` | Pan-vascular spectrum one-page summary | `step4/results/*` | `Step5_泛血管谱结果.md/.html` |

---

## Stage 3 — Reverse MR and MVMR (`step6/`, `step7/`)

| # | Script | Purpose | Key outputs |
|---|---|---|---|
| 3.1 | `s6_rev_pipeline.py` | **Reverse MR**: atherosclerotic vascular outcomes as exposure → thyroid function as outcome (multi-route SNP matching; accepts only same-strand or strand-complement alleles) | `step6/results/*.json` |
| 3.2 | `s7_mvmr.py` | **MVMR**: all three thyroid traits in one model (MV-IVW, MV-Egger, conditional F, heterogeneity Q), to test whether the hypothyroidism effect is independent of TSH and FT4 | `step7/results/*.json` |
| 3.3 | `s67_tables.py` | Export Step 6/7 result tables to XLSX | `Step6-7_..._数据表.xlsx` |
| 3.4 | `s67_onepager.py` | Steps 6+7 one-page summary | `Step6-7_反向MR与MVMR.md/.html` |

---

## Stage 4 — Colocalisation (`step8/`)

| # | Script | Purpose | Key outputs |
|---|---|---|---|
| 4.1 | `s8_precheck.py` | Variant-coverage check in each candidate locus (determines whether `coloc.abf` is feasible) | `precheck_coverage.tsv/.json`, `precheck_raw.json` |
| 4.2 | `coloc_abf.py` | Pure-Python implementation of `coloc.abf` (Giambartolomei 2014; Wakefield 2009 approximate Bayes factors). ⚠ Internal log-sum-exp must be **log10** to match the lABF base | (module) |
| 4.3 | `s8_coloc.py` | Colocalisation across 10 hypothyroidism loci against four non-overlapping outcomes | `coloc_results.tsv`, `locus_annotation.json` |
| 4.4 | `s8_tables.py`, `s8_onepager.py` | XLSX tables and one-page summary | `Step8_共定位结果*` |

---

## Stage 5 — Power analysis and sample-overlap framework (`step9/`)

| # | Script | Purpose | Key outputs |
|---|---|---|---|
| 5.1 | `s9_power.py` | Power analysis based on the **actual SEs** of the completed MR analyses (forward, reverse, MVMR); MDE at 80% power; **outcome-scale screening** (median instrument se_Y compared across outcome files) | `power_*.tsv`, `power_summary.json`, `scale_check.tsv` |
| 5.2 | **`s9_overlap_sim.py`** | **Sample-overlap framework.** Adopts the analytic relationship of **Burgess et al. (2016)** — reparameterised here as `R(π) = 1 + π·r/F_w` — and **applies it in reverse** (observed paired ratio → required π or r). Includes a 10,000-replicate × 125-cell Monte Carlo consistency check. **No methodological novelty is claimed.** | `overlap_sim_curve.tsv`, `overlap_sim_anchor.tsv`, `overlap_sim_mc_validation.tsv`, `overlap_sim_summary.json` |
| 5.3 | `s9_tables.py`, `s9_onepager.py`, `s9_overlap_onepager.py` | XLSX tables and one-page summaries | `Step9_统计功效*`, `Step9_重叠偏倚模拟*` |

---

## Notes on reproducibility

1. **Root path.** Scripts contain a hard-coded project root. Run `python relocate.py` first (see README §5).
2. **Byte-level provenance.** Scripts are shipped **unmodified** so that the archived code is identical to the code that produced the published results.
3. **Numbers are never hand-copied.** All summary documents (`*_onepager.py`) read their numbers directly from the result JSON/TSV files.
4. **Raw data are not included** — see `DATA_SOURCES.md`.
