# Data source

## Dataset used (real data)

- **Title:** OTU Abundance Data: Oral Microbiome in Hypertensive Patients with/without Periodontitis
- **Deposited by:** Ziyin Guo (2025), figshare. DOI: [10.6084/m9.figshare.29897750.v1](https://doi.org/10.6084/m9.figshare.29897750.v1)
- **Direct file URL used by the code:** https://ndownloader.figshare.com/files/57145445 (`data.zip`, 667,945 bytes)
- **License:** CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). The cached zip is redistributed unchanged in `data/raw/` under that license, with this attribution.
- **Accessed:** 2026-09-27
- **Associated paper:** Guo Z. et al. "A comparative analysis of oral microbial communities in hypertensive patients with and without chronic periodontitis." *BMC Oral Health* 26:836 (2026). DOI: 10.1186/s12903-026-08144-6. Raw reads: NCBI SRA BioProject PRJNA1304526.

## What is in the file

| File inside the zip | Content |
|---|---|
| `OTU_Feature_Table.tsv` | Despite the name, a BIOM 1.0 JSON table: 25,540 features ("ASV1"...) x 67 samples, integer read counts, with a 7-rank SILVA lineage per feature |
| `Taxonomic_Annotations.tsv` | JSON list of the taxon names at each rank (not needed; the lineages are in the BIOM file) |
| `Sample_Metadata.csv` | Vendor project report (partly in Chinese) containing a `SamplesID / treat1` table with the group of every sample |

Per the paper: unstimulated saliva, 16S rRNA V3-V4 region (primers 338F
`ACTCCTACGGGAGGCAGCA` / 806R `GGACTACHVGGGTWTCTAAT`), Illumina NovaSeq 6000,
processed with QIIME2. Periodontitis was diagnosed with the 2018 World
Workshop (AAP/EFP) classification.

## Groups

| Code | Meaning | n |
|---|---|---|
| H | Periodontally healthy, normotensive | 16 |
| P | Periodontitis, normotensive | 18 |
| T | Hypertension, periodontally healthy | 17 |
| TP | Hypertension + periodontitis | 16 |

Only group labels are available per sample. Age, sex, smoking, and clinical
measurements (probing depth, attachment loss) are **not** in the public
files, so none of them can be adjusted for.

## Caching and offline use

`oralbiome.data.download_raw()` downloads the zip once to
`data/raw/guo2025_figshare_29897750.zip`; if that file exists it is used
without any network access. Parsed tables are written to `data/processed/`.

## Fallback

If the download fails and no cached zip exists, `load_dataset()` generates a
Dirichlet-multinomial synthetic dataset (`oralbiome.synthetic`), writes a
`data/processed/DATA_IS_SYNTHETIC` marker, prefixes every figure filename and
title with `SYNTHETIC`, and saves the planted ground truth. **The results in
this repository were produced from the real dataset, not the fallback.** The
synthetic generator is still used as a positive control (see `run_all.py` and
`tests/test_synthetic_recovery.py`).
