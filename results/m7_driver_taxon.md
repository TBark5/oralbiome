# Driver taxon of the M7 classifier: Unclassified Bacilli

Selected in 100% of L1 folds (mean standardised coefficient +1.411). Lineage in the deposit: kingdom = Bacteria; phylum = Firmicutes; class = Bacilli; order = Unclassified Bacilli; family = Unclassified Bacilli; genus = Unclassified Bacilli; species = Unclassified Bacilli.

## Can it be identified?

No. The figshare deposit contains a feature table, taxonomy strings and sample groups but no representative sequences, so there is nothing to extract to FASTA and nothing to search against a 16S reference. `m7_driver_taxon.fasta` was therefore not created. Identification would require re-processing the raw reads (SRA BioProject PRJNA1304526) to recover the sequences of these features.

## Prevalence and abundance by group

| Group | n | Detected | Median rel. abundance (%) | Max rel. abundance (%) | Reads |
|---|---|---|---|---|---|
| H | 16 | 5/16 | 0.0000 | 0.0102 | 21 |
| P | 18 | 16/18 | 0.0288 | 0.6522 | 1399 |
| T | 17 | 8/17 | 0.0000 | 0.0879 | 210 |
| TP | 16 | 12/16 | 0.0114 | 0.2678 | 447 |

Presence/absence, Fisher's exact test (descriptive, not multiplicity-corrected): P vs H p = 0.0011; TP vs T p = 0.16.

## Which features make up the label

32 features carry this label. The most abundant, ASV8652, holds 881 of 2077 reads and is detected in H 1/16, P 12/18, T 3/17, TP 11/16 samples. The other 31 features contribute the remaining reads; 31 of all features are seen in a single sample only. Full table: `m7_driver_taxon_features.csv`.

## Overlap with flagged and high-richness samples

- Flagged samples where it is detected: H7, P9, P18, T9, TP15 (of H7, P18, P7, P9, T10, T9, TP15, TP7).
- High-richness healthy samples where it is detected: H11, H14, H6, H7, H8 (of H11, H13, H14, H6, H7, H8).

## Network partners (M6, |rho| >= 0.6, q < 0.05)

| Partner | rho | q |
|---|---|---|
| Acidibacter | -0.623 | 0.00175 |
| Candidatus_Solibacter | -0.606 | 0.00257 |
| Unclassified Isosphaeraceae | -0.604 | 0.00268 |
