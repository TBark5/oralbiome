# Data source

## Citations

Source study (the paper):

> Guo Z., Yu X., Liu Y., Hu Q., Zhang Z., Zhang C., Li J. (2026). "A comparative analysis of oral microbial communities in hypertensive patients with and without chronic periodontitis." BMC Oral Health 26(1):836. doi:10.1186/s12903-026-08144-6

- PubMed ID: 41923033
- PMCID: PMC13169644
- DOI: 10.1186/s12903-026-08144-6
- Raw reads: NCBI SRA BioProject PRJNA1304526
- Citation checked against the PubMed record (authors, year, volume, issue, DOI) on 2026-09-27; the article number 836 is the journal's.

Dataset (the files this project downloads), citation as provided by figshare:

> Guo, Ziyin (2025). OTU Abundance Data: Oral Microbiome in Hypertensive Patients with/without Periodontitis. figshare. Dataset. https://doi.org/10.6084/m9.figshare.29897750.v1

The dataset (deposited 2025, one depositor) and the paper (published 2026,
seven authors) are separate works with separate citations. Both are given in
full wherever this repository cites them.

## Dataset record (verified 2026-09-27)

| Field | Value |
|---|---|
| DOI | 10.6084/m9.figshare.29897750.v1 (resolves via doi.org to the figshare record, HTTP 302 -> figshare.com) |
| Version | 1 (the only version listed by the figshare API) |
| Published | 2025-08-13 (record last modified 2025-11-28) |
| License | CC BY 4.0, https://creativecommons.org/licenses/by/4.0/ (read from the figshare API record) |
| File | `data.zip`, 667,945 bytes, MD5 `6a49330fadc5e773059c52a84f6af67c` |
| Direct URL used by the code | https://ndownloader.figshare.com/files/57145445 |
| Cached copy in this repository | `data/raw/guo2025_figshare_29897750.zip`, MD5 `6a49330fadc5e773059c52a84f6af67c` (identical to the figshare file) |
| First accessed / last verified | 2026-09-27 / 2026-09-27 |

**Redistribution.** CC BY 4.0 allows copying and redistribution, including
of unmodified files, provided that appropriate credit is given, a link to the
license is provided, and any changes are indicated. The zip is redistributed
unchanged; attribution, the license link and a no-changes statement are in
[LICENSE-DATA](LICENSE-DATA). The *article* is published under CC BY-NC-ND
4.0; no part of the article or its supplementary files is redistributed
here.

## What is in the zip

| File inside the zip | Content |
|---|---|
| `OTU_Feature_Table.tsv` | Despite the name, a BIOM 1.0 JSON table (written by BIOM-Format 2.1.6, dated 2024-10-31): 25,540 features x 67 samples, integer read counts, with a 7-rank taxonomy string per feature. Feature IDs are labelled "ASV1", "ASV2", ...; the paper describes them as OTUs clustered at 97% similarity with QIIME2 |
| `Taxonomic_Annotations.tsv` | JSON list of the taxon names at each rank (not needed; the lineages are in the BIOM file) |
| `Sample_Metadata.csv` | Vendor project report (Biomarker Technologies pipeline "v3.2", partly in Chinese) with primers, region and a `SamplesID / treat1` table giving the group of every sample |

The deposit contains **no representative sequences** (no FASTA or QIIME2
artefacts), so individual features cannot be re-classified or searched
against a reference from these files.

Per the paper: unstimulated saliva, 16S rRNA V3-V4 region (primers 338F
`ACTCCTACGGGAGGCAGCA` / 806R `GGACTACHVGGGTWTCTAAT`), Illumina NovaSeq 6000,
reads filtered with Trimmomatic v0.33, primers trimmed with Cutadapt v1.9.1,
then denoised, merged, chimera-filtered and clustered into 97% OTUs in
QIIME2. Periodontitis was diagnosed with the 2018 World Workshop (AAP/EFP)
classification. Smokers, people with systemic disease and people who took
antibiotics or probiotics in the previous 3 months were excluded.

## Taxonomy reference database

**Not stated in any source available to this project.** Checked on 2026-09-27:

- the figshare deposit: the BIOM rows carry only a `taxonomy` list, with no
  classifier name, database version or confidence value; the vendor
  metadata file records pipeline version, primers and region only;
- the article full text (PMC13169644): "Taxonomic annotation of OTUs against
  a reference library", with no database named;
- the article's supplementary files (supplementary methods, consent form,
  clinical spreadsheet): no database named.

What the taxonomy strings themselves show: the naming conventions are
characteristic of the SILVA 138 family of releases (for example
`Prevotella_7`, `Rikenellaceae_RC9_gut_group`, `[Eubacterium]_nodatum_group`,
`*_UCG_*` genera, and phyla such as `Actinobacteriota`, `Desulfobacterota`,
`Patescibacteria`), and they are not HOMD/eHOMD names. One phylum is spelled
`Campylobacterota`, whereas SILVA 138/138.1 use `Campilobacterota`, so the
exact release (or any vendor renaming) cannot be confirmed. The
classification confidence threshold is not recorded anywhere. This project
therefore treats the database as "unknown, most likely SILVA 138-era, not
HOMD".

## Groups

| Code | Meaning | n |
|---|---|---|
| H | Periodontally healthy, normotensive | 16 |
| P | Periodontitis, normotensive | 18 |
| T | Hypertension, periodontally healthy | 17 |
| TP | Hypertension + periodontitis | 16 |

The figshare deposit contains only the group of each sample. Per-sample
age, sex, height, weight, blood pressure, probing depth and attachment loss
exist in the article's supplementary spreadsheet (not in the deposit, and not
used by this project). Smoking cannot confound the comparison because
smokers were excluded by the study protocol.

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
