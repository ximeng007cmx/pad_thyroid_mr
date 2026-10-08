# Data sources

All datasets used in this study are **publicly available summary statistics**. No individual-level data were analysed. Raw files are not redistributed in this repository (≈5.5 GB); every accession and download location is listed below.

**All PMIDs and DOIs in this file were verified against Europe PMC on 2026-10-08.**

---

## 1. Exposures — thyroid function (3 traits)

Single publication providing all three trait GWAS, so the three exposures share the same cohort structure (**triple-exposure triangulation by design**).

| Trait | GWAS Catalog | n (cases / controls, or n) | Local file |
|---|---|---|---|
| Hypothyroidism (binary) | **GCST90572791** | 113,393 cases / 1,065,268 controls | `data/exposure/hypo_GCST90572791.h.tsv.gz` |
| TSH (continuous) | **GCST90572789** | 482,873 (EUR) | `data/exposure/tsh_GCST90572789.h.tsv.gz` |
| Free T4 (continuous) | **GCST90572790** | 191,449 (EUR) | `data/exposure/ft4_GCST90572790.h.tsv.gz` |

- **Reference**: Rand et al. *Nat Genet* 2025. **PMID 41238958**; DOI 10.1038/s41588-025-02410-z
- **Download**: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90572001-GCST90573000/<GCST>/harmonised/`
- **Format**: GWAS-SSF v1.0, GRCh38
- **md5 (verified byte-identical to the official manifest)**:
  - `8b47a15de4e2220898809146777084df` — GCST90572791
  - `59689c1b938681871eb863a0bbffee8b` — GCST90572789
  - `dda8c77bd540f223bcaa1734f782c226` — GCST90572790
- **Cohort composition** (from the source publication): UK Biobank + FinnGen (**DF10**) + 23andMe + Danish cohorts. The 23andMe contribution is not in the public release; TSH and FT4 do **not** include Finnish data.
- **Protein-score companion**: PGS Catalog **PGS005218** (hypothyroidism PRS)

---

## 2. Primary outcome — peripheral artery disease (PAD)

| Role | Source | GWAS Catalog | Cases / controls | Local file |
|---|---|---|---|---|
| **Primary** | Sakaue et al. 2021 | **GCST90018890** | 7,114 / 475,964 | `data/outcome/pad_GCST90018890.h.tsv.gz` |
| **Replication** | FinnGen **release 13** | endpoint **`I9_PAD`** | 22,244 / 460,490 | `data/replication/finngen_R13_I9_PAD.gz` |
| **Replication** | MVP (dbGaP) | **pha004826** (phs001672.v14.p1) | 31,307 / 211,753 | `raw/mvp_pad/MVP.te.PAD.dbGAP.txt.gz` |

- **Sakaue 2021**: *Nat Genet* 2021. **PMID 34594039**; DOI 10.1038/s41588-021-00931-x
  - ⚠ **The European component is an intra-European meta-analysis of UK Biobank and FinnGen release 3**: EUR 483,078 = **UKB 350,366 + FinnGen R3 132,712** (verified against the source publication and the FinnGen R3 manifest, 2026-10-08). This is central to the sample-overlap argument.
- **FinnGen R13**: `https://storage.googleapis.com/finngen-public-data-r13/summary_stats/` (public bucket; the R13 endpoint manifest is `data/replication/finngen_r3_manifest.tsv` for the R3 release used in the cross-check)
- **MVP PAD**: Klarin et al. *Nat Med* 2019. **PMID 31285632**; DOI 10.1038/s41591-019-0492-5
  - Accession: dbGaP **phs001672**, analysis **pha004826**; unrestricted summary statistics from release **v14.p1** (2026-04-08)
  - ⚠ The `Position` column in this file is **GRCh37 (hg19)**, not GRCh38 → matching must be by rsID only.

---

## 3. Additional outcomes — pan-vascular spectrum

| Outcome | GWAS Catalog | Cases / controls | Reference | PMID | DOI |
|---|---|---|---|---|---|
| CAD (no UKB) | **GCST003116** | — | Nikpay et al. *Nat Genet* 2015 | **26343387** | 10.1038/ng.3396 |
| CAD (UKB) | **GCST90013864** | 29,339 / 322,724 | Mbatchou et al. *Nat Genet* 2021 | **34017140** | 10.1038/s41588-021-00870-7 |
| Myocardial infarction | **GCST90038610** | 11,081 / 473,517 | Dönertaş et al. *Nat Aging* 2021 | **33959723** | 10.1038/s43587-021-00051-5 |
| Ischaemic stroke — large artery | **GCST005840** | 4,373 / 406,111 | MEGASTROKE (Malik et al. *Nat Genet* 2018) | **29531354** | 10.1038/s41588-018-0058-3 |
| Ischaemic stroke — small vessel | **GCST005841** | 5,386 / 192,662 | MEGASTROKE (Malik et al. *Nat Genet* 2018) | **29531354** | 10.1038/s41588-018-0058-3 |

> ⚠ **Myocardial infarction (GCST90038610) is on a linear-probability scale**, not a log-odds scale (median instrument se ≈ 0.030× the reference median). Its β/SE must not be interpreted as a log-OR; z and p are unaffected. Screening code: `step9/scripts/s9_power.py` (see `step9/results/scale_check.tsv`).

---

## 4. Reference panel

- **1000 Genomes Phase 3, GRCh38, European ancestry** — used for LD clumping and harmonisation
- Zenodo DOI: **10.5281/zenodo.19068199**
- Local path: `ref/1kg_eur_hg38/` (~1.9 GB, not redistributed)

---

## 5. Analysis software

| Tool | Version | Used for |
|---|---|---|
| PLINK | **1.90b7** | LD clumping (`--clump`, r² < 0.001, 10 Mb) |
| Python | 3.13 | all statistical analysis, figure generation, table generation |
| numpy | 2.5.3 | numerical routines |
| matplotlib | 3.11.2 | figures |
| openpyxl | 3.1.5 | Excel output |

MR estimators (IVW fixed/random, MR-Egger, weighted median, weighted mode, Cochran's Q, MR-PRESSO-equivalent outlier detection, leave-one-out, MVMR, `coloc.abf`) are implemented from scratch in Python (`step3/scripts/mrlib.py`, `step8/scripts/coloc_abf.py`).

---

## 6. Data not redistributed here

| Item | Size | Reason |
|---|---|---|
| Raw GWAS summary statistics | ≈3.6 GB | publicly downloadable from the accessions above |
| 1000 Genomes reference panel | ≈1.9 GB | publicly available (Zenodo DOI above) |
| `step8/results/precheck_raw.json` | 14 MB | per-variant pre-check dump; regenerated by `s8_precheck.py` |
