# Results and discussion

This document checks each prediction in `HYPOTHESIS.md` against what the
pipeline actually found. All numbers come from files in `results/`. Primary
comparison: healthy (H, n = 16) vs periodontitis (P, n = 18). Replication:
hypertension (T, n = 17) vs hypertension + periodontitis (TP, n = 16).

## Before the predictions: a data-quality finding

The QC step found something that affects how every later result should be
read. Of the 25,540 features in the table, 22,781 occur in only one sample.
Six healthy samples (H6, H7, H8, H11, H13, H14) have 1,700-2,300 observed
features, against a median of about 350, and many of the extra features belong
to soil- or gut-associated lineages (Acidobacteriota, Chloroflexi,
Gemmatimonadota, Muribaculaceae) that are not typical of saliva. Eight
samples (four in the primary comparison: H7, P7, P9, P18) have more than half
of their reads in features found in no other sample.

Contamination during extraction or library preparation, or a processing
batch, are the most likely explanations. The public files contain no
negative controls or batch information, so this cannot be confirmed. The
analysis handles it in three ways: genus-level analysis (which absorbs most
of the one-off features; median 99% of reads are kept after filtering),
flagged-sample sensitivity analyses, and a diagnostic table showing which
findings are carried by the high-richness subset.

## P1 - Alpha diversity differs between groups: **not supported**

| Metric | Median H | Median P | p | q (BH, 6 metrics) | rank-biserial r [95% CI] |
|---|---|---|---|---|---|
| Shannon (primary) | 2.884 | 3.086 | 0.796 | 0.796 | +0.06 [-0.35, +0.45] |
| Gini-Simpson | 0.899 | 0.915 | 0.234 | 0.281 | +0.24 [-0.16, +0.62] |
| Observed genera | 179.5 | 140.0 | 0.133 | 0.200 | -0.31 [-0.67, +0.11] |
| Pielou's evenness | 0.561 | 0.640 | 0.037 | 0.200 | +0.42 [+0.03, +0.78] |
| Chao1 | 182.2 | 141.4 | 0.133 | 0.200 | -0.31 [-0.66, +0.11] |
| ACE | 184.2 | 141.4 | 0.133 | 0.200 | -0.31 [-0.66, +0.11] |

Chao1 and ACE were added in a later pass so the results could be compared
with the source study, which reports those two indices; the BH correction
now runs over six metrics instead of four (DECISIONS.md D45). This changed
q-values but no conclusion.

The pre-declared falsification criterion was "Shannon p ≥ 0.05 and the effect
size CI includes 0". Both are true, so P1 is falsified for the primary metric.
Evenness is nominally higher in periodontitis (p = 0.037), but it does not
survive correction for testing six metrics, and it becomes borderline
(p = 0.051) once the flagged samples are removed. The feature-level
sensitivity analysis agrees (Shannon p = 0.904; Chao1 p = 0.479).

**Interpretation.** The literature disagrees on the direction of alpha
diversity change in periodontitis, and subgingival plaque studies often
report *higher* richness in disease. The richness estimators here point the
same way as the source study (lower in periodontitis) but are not
significant; the section "Concordance and discordance with the source study"
compares the two analyses in detail. With 16 vs 18 samples, a small true
difference cannot be detected: the Shannon CI (-0.35 to +0.45) cannot rule
out a moderate difference in either direction.

## P2 - Groups separate in beta-diversity space: **supported, with a caveat**

- Bray-Curtis PERMANOVA: R² = 0.082, pseudo-F = 2.86, p = 0.015.
- PERMDISP (Bray-Curtis): p = 0.065; mean distance to centroid 0.366 (H) vs 0.447 (P).
- Jaccard (presence/absence): R² = 0.058, p = 0.072; PERMDISP p = 0.896.
- Without the four flagged samples: Bray-Curtis R² = 0.119, p = 0.002;
  Jaccard R² = 0.097, p = 0.007.

The falsification criterion (PERMANOVA p ≥ 0.05) is not met for the primary
distance, so P2 survives. The caveat: PERMANOVA is sensitive to differences in
spread as well as location, and PERMDISP is borderline (p = 0.065), so part of
the signal could be that periodontitis samples are more variable rather than
shifted. Group membership explains only about 8% of the between-sample
variation. The signal gets stronger, not weaker, when the atypical samples are
removed, so it is not produced by them.

## P3 - Red-complex species enriched in periodontitis: **supported**

Pre-declared, species-level test (BH across the three species):

| Species | CLR difference (P - H) [95% CI] | r | p | q | Detected in H / P |
|---|---|---|---|---|---|
| *Porphyromonas gingivalis* | +2.98 [1.15, 4.81] | +0.60 | 0.0032 | 0.0048 | 50% / 78% |
| *Tannerella forsythia* | +1.90 [0.51, 3.29] | +0.51 | 0.0124 | 0.0124 | 62% / 78% |
| *Treponema denticola* | +2.63 [0.87, 4.40] | +0.65 | 0.0013 | 0.0038 | 62% / 78% |

All three are more abundant in periodontitis saliva, so none meets its
falsification criterion. In the untargeted scan of 354 genera, none of the
three genera passes FDR on its own: *Treponema* ranks 11th (q = 0.094),
*Tannerella* 34th (q = 0.122) and *Porphyromonas* 108th (q = 0.323).

Two points are worth explaining:

1. **Why the species is significant but the genus *Porphyromonas* is not.**
   The genus also contains commensal species that are often reported in
   healthy saliva (for example *P. pasteri*). Pooling them with *P. gingivalis* dilutes the
   signal. This is a concrete example of why taxonomic resolution matters.
2. **Why a pre-declared test is not cherry-picking.** The three species were
   named in `HYPOTHESIS.md` before any analysis, and the correction is over
   those three tests. Choosing targets *after* seeing the scan and then
   applying a 3-test correction would be cherry-picking. The untargeted scan
   is reported next to it so readers can judge both.

**How much weight the species-level result can bear.** The species labels
come from a taxonomy whose reference database and confidence threshold are
not stated anywhere in the deposit or the paper; the names follow SILVA
138-era conventions, not the curated oral reference eHOMD (DATA_SOURCE.md).
A ~470 bp V3-V4 amplicon often cannot separate close relatives within
*Porphyromonas*, *Tannerella* or *Treponema*, so each "species" here should
be read as "features the unknown classifier assigned to that species name".
The result is consistent with the genus-level evidence (all three genera lean
towards periodontitis and co-occur strongly), but it should be treated as
supporting evidence that needs re-annotation against eHOMD, not as a
species-level identification.

For each species, three of the four lowest periodontitis values come from
flagged samples (P7, P9, P18), whose reads are dominated by non-oral
features. In this dataset *P. pasteri* is almost as abundant as
*P. gingivalis*, which explains the dilution at genus level.

## P4 - Composition predicts disease above chance: **supported, with a caveat**

| Model | AUC [bootstrap 95% CI] | Null mean AUC | Permutation p |
|---|---|---|---|
| L1 logistic regression | 0.92 [0.81, 1.00] | 0.493 | 0.0099 |
| Random forest | 0.81 [0.65, 0.94] | 0.498 | 0.0099 |

Both CIs exclude 0.5 and both models beat every one of 100 label shuffles
(0.0099 is the smallest p possible with 100 shuffles), so P4 is not
falsified. All filtering and transformation happened inside the CV folds, so
the AUC is not inflated by leakage; the permuted-label null confirms that the
pipeline itself does not produce high AUCs (null mean ≈ 0.5).

The caveat is large. One feature, "Unclassified Bacilli" (reads assigned to
the class Bacilli with no order, family or genus), was selected in 100% of
L1 folds. It is detected in 16 of 18 periodontitis samples, 0 of the 10
typical healthy samples and 5 of the 6 high-richness healthy samples (at low
counts). Removing it drops the AUC to 0.64 (L1) and 0.72 (random forest). The
"prediction" therefore rests mainly on one unidentified taxon.

**Could it be identified? No.** The figshare deposit contains no
representative sequences, so no FASTA could be extracted and no search against
a 16S reference was possible (`results/m7_driver_taxon.md`). What the deposit
does show:

- The label pools 32 features. One of them (ASV8652) holds 881 of its 2,077
  reads and is detected in 12/18 P, 11/16 TP, 3/17 T and 1/16 H samples; each
  of the other 31 features occurs in a single sample.
- At label level it is detected in 16/18 P, 12/16 TP, 8/17 T and 5/16 H
  samples (presence P vs H: Fisher p = 0.0011; TP vs T: p = 0.16; descriptive,
  not corrected).
- All five healthy detections are in the high-richness subset (H6, H7, H8,
  H11, H14); it is absent from the other 10 healthy samples.
- Its only network edges are negative correlations (rho -0.60 to -0.62) with
  soil-associated taxa (*Acidibacter*, *Candidatus Solibacter*, an
  unclassified Isosphaeraceae).

A contaminant introduced at random would not be expected to track
periodontitis status in two separate groups of people, so the dominant
feature looks more like a genuine low-abundance organism associated with
periodontitis. However, its presence only in the atypical healthy samples and
the many one-off features pooled under the same label mean a technical
component cannot be excluded. Identifying it requires the raw reads.

## P5 (exploratory) - Red-complex genera co-occur: **supported**

| Pair | Spearman rho | q |
|---|---|---|
| *Porphyromonas* - *Tannerella* | 0.68 | < 0.001 |
| *Porphyromonas* - *Treponema* | 0.85 | < 0.001 |
| *Tannerella* - *Treponema* | 0.81 | < 0.001 |

All three genera are also among the top 10 network hubs (degree 49-58), along
with *Campylobacter*, *Fusobacterium*, *Prevotella*, *Catonella* and
*Selenomonas*, which are anaerobes known from periodontal pockets. This fits
the idea that these organisms occur together as a consortium. Two caveats:
samples of both groups were pooled, so taxa that are all higher in
periodontitis will correlate partly because of the group difference; and a
correlation in saliva across people says nothing about physical or metabolic
interaction inside a pocket.

## Replication in T vs TP: **did not replicate**

| Measure | Primary (H vs P) | Replication (T vs TP) |
|---|---|---|
| Bray-Curtis PERMANOVA | R² = 0.082, p = 0.015 | R² = 0.031, p = 0.402 |
| L1-LR AUC [95% CI] | 0.92 [0.81, 1.00] | 0.52 [0.32, 0.72] |
| L1-LR permutation p | 0.0099 | 0.495 |
| *P. gingivalis* r (q) | +0.60 (0.005) | +0.28 (0.530) |
| *T. forsythia* r (q) | +0.51 (0.012) | +0.07 (0.759) |
| *T. denticola* r (q) | +0.65 (0.004) | +0.12 (0.759) |

In people with hypertension, periodontitis status shows no detectable
community-level difference, and the classifier performs at chance. The
red-complex effects still point the same way but are much smaller, and their
CIs include zero (`figures/15_replication_forest.png`).

Possible explanations, which this data cannot tell apart:

1. **The primary result is partly a false positive or inflated.** Small
   studies overestimate effect sizes, and the primary classifier depends on
   one feature.
2. **Hypertension changes the picture.** The source paper reports that
   hypertension alone shifts the oral microbiome, and antihypertensive drugs
   (some reduce saliva flow) could mask a periodontitis signal.
3. **Batch or contamination.** If the H-vs-P difference is partly driven by
   the unusual healthy samples, it would not be expected to appear in T vs TP.

The honest summary is that the evidence for a *robust* salivary
periodontitis signature in this dataset is weak. The strongest single result
is the pre-declared red-complex enrichment in the primary comparison, which
leans the same way in the replication.

## Concordance and discordance with the source study

The source study analysed the same samples with a different pipeline. Its
findings (quoted from the article, PMC13169644) and how they compare with
this project:

| Source study finding | This project | Status |
|---|---|---|
| Chao1 and ACE significantly higher in healthy controls than in each disease group (P, T, TP), P < 0.05, Tukey's HSD after one-way ANOVA | H vs P: Chao1 182.2 vs 141.4, ACE 184.2 vs 141.4, observed 179.5 vs 140.0 (r = -0.31, p = 0.133, q = 0.200) | Same direction, **not reproduced** at significance |
| No significant Shannon or Simpson differences across groups (P > 0.05) | Shannon p = 0.796, Gini-Simpson p = 0.234 | **Reproduced** |
| Four-group PERMANOVA significant (R² = 0.126, P < 0.01); weighted UniFrac PCoA separates healthy controls from the overlapping disease groups | H vs P Bray-Curtis PERMANOVA R² = 0.082, p = 0.015 (Jaccard p = 0.072) | Consistent in kind (small but significant composition difference); **not directly comparable** (different contrast and distance) |
| Disease groups more similar to each other than to healthy controls | T vs TP: PERMANOVA p = 0.402 | Consistent with it, but not a test of it |
| *P. gingivalis* enriched in P/TP; *T. denticola* and *T. forsythia* enriched in P, T and TP | All three higher in P than H (q = 0.004-0.012); T vs TP small, non-significant effects in the same direction | **Reproduced** for H vs P |
| *Fretibacterium* and *Filifactor* higher, *Haemophilus* lower in P/TP vs HC | *Fretibacterium* q = 0.089, *Filifactor* q = 0.105 (both higher in P); *Haemophilus* lower in P but p = 0.546 | Same direction, **not significant** after FDR |
| Weighted UniFrac, *Streptococcus* sp. I-G5 in T, 4-group comparisons, periodontitis stage, PICRUSt2 function | - | **Not attempted** (no sequences or tree for UniFrac/PICRUSt2; T-specific and stage questions are outside this project) |

**Why richness and diversity indices disagree, and why that is expected.**
Chao1 and ACE estimate how many taxa are present, including unseen ones, and
are driven almost entirely by rare taxa (singletons and doubletons). Shannon
and Simpson weight taxa by abundance, so rare taxa barely move them. A group
that carries many extra rare taxa will score higher on Chao1/ACE with no
change in Shannon. The source study found exactly this pattern (richness
differs, Shannon/Simpson do not), and this project finds the same pattern in
direction.

**Why the richness difference is significant there but not here.** These
methodological differences could contribute; the data available here cannot
say which matters most:

1. **Contrast.** The paper tested four groups with one-way ANOVA and compared
   each disease group with healthy controls (Tukey's HSD). The disease groups
   were compared separately with the control group, not pooled. The primary
   comparison here is healthy vs periodontitis only (H vs P) with a
   Mann-Whitney test.
2. **Feature level.** The paper computed indices on 97% OTUs. The primary
   analysis here uses genera. At feature (OTU) level the direction is the same
   and still not significant (Chao1 403.0 vs 352.2, p = 0.479).
3. **Filtering and rarefaction.** This project rarefies all samples to 44,786
   reads (the minimum depth) with one seeded draw; the paper's rarefaction
   depth is not reported. Prevalence filtering is not applied before alpha
   diversity here, but genus aggregation merges rare features.
4. **Multiple testing.** Here q-values are corrected across six metrics; the
   paper reports per-index P-values.
5. **High-richness samples.** A post-hoc check (not part of any
   pre-registered test) removed the six healthy samples with > 1,000 observed
   features. Without them the healthy group's richness is no longer higher:
   genus-level Chao1 medians 139.8 (H) vs 141.4 (P), p = 0.649; feature-level
   305.8 vs 352.2, p = 0.119
   (`results/m3_richness_high_richness_check.csv`). The paper's own Venn
   analysis reports that the healthy group contained 15,474 OTUs, 11,634 of
   them unique to that group, and describes the healthy samples as more
   dispersed. This project cannot tell whether the same samples influenced the
   published richness result, only that in this dataset the healthy-group
   richness advantage depends on them.

None of this shows that the published analysis is wrong. Both analyses point
the same way on richness, agree that Shannon and Simpson do not differ, and
agree that the red-complex species are enriched in periodontitis. They differ
in whether the richness difference is statistically significant, which
depends on the contrast, the features, the test and a handful of unusual
samples.

**What this means for the replication result.** The source study reports
that hypertension-only patients (T) already show periodontitis-like
dysbiosis, including enrichment of *T. denticola* and *T. forsythia*. If that
is right, T is not a healthy reference, and the T vs TP comparison here is
expected to show a smaller difference than H vs P. The failure to replicate
is therefore consistent with the source study's findings, and it does not by
itself show that the H vs P signal is false.

## Negative and positive controls

- **Differential abundance null:** with shuffled labels, an average of 0.07
  genera pass FDR (200 shuffles), versus 2 on the real labels; 72 genera have
  raw p < 0.05, against 17.7 expected by chance.
- **Classifier null:** mean null AUC 0.493 (L1) and 0.498 (RF).
- **Positive control:** on synthetic data with 15 planted differences, the
  pipeline recovered 10 in the correct direction (sensitivity 67%) with 1
  false positive (false discovery proportion 9%), and PERMANOVA p = 0.002.
  The pipeline can find real signals when they exist, and its false discovery
  rate is close to the nominal 5-10%.

## What I would do next with more data

1. Re-process the raw reads (SRA PRJNA1304526) with DADA2 myself and
   re-annotate the sequences against eHOMD, the curated oral reference, to
   firm up the species-level calls; add decontamination (e.g. the `decontam`
   approach) if negative controls can be obtained from the authors.
2. Use the recovered sequences to identify the "Unclassified Bacilli" features
   (BLAST against eHOMD and NCBI 16S).
3. Use a compositionally aware differential abundance method (ANCOM-BC or
   ALDEx2) to check the M5 results.
4. Meta-analyse several public periodontitis datasets to test whether the
   classifier generalises across studies, which is the real test of a
   signature.
5. Use subgingival plaque instead of saliva and shotgun metagenomics to get
   species/strain resolution and functional genes.
6. Use a longitudinal design (sampling before and after periodontal
   treatment) to get closer to temporal order, although still not causality.
