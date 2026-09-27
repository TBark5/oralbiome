# Hypothesis

Written in Phase 0, before any analysis was run. At that point I had seen only
the dataset's structure (sample groups, read depths, the fact that the three
named pathogens appear in the taxonomy) and the abstract of the source paper
(Guo et al., BMC Oral Health 2026), which reports lower alpha diversity and
enrichment of *Treponema denticola* in disease. That abstract is prior
knowledge, so these predictions are not fully blind.

## Biological question

How does the composition of the oral (salivary) bacterial community differ
between periodontally healthy people and people with periodontitis, and can
community composition alone predict which group a sample came from?

Periodontitis is a chronic inflammatory disease of the tissues supporting the
teeth. The "polymicrobial synergy and dysbiosis" model describes it as a shift
in the whole community rather than a single infection, with the "red complex"
species *Porphyromonas gingivalis*, *Tannerella forsythia* and *Treponema
denticola* (Socransky et al., 1998) repeatedly associated with deep pockets.
Most of that evidence comes from subgingival plaque; this dataset is saliva,
which pools bacteria from all oral surfaces, so any disease signal should be
weaker and more diluted.

## Primary comparison

Healthy controls (group `H`, n = 16) versus periodontitis without hypertension
(group `P`, n = 18). Hypertension is a separate factor in the source study, so
the primary comparison uses only the two normotensive groups to avoid mixing in
a second condition.

## Predictions and falsification criteria

All tests are two-sided unless stated. alpha = 0.05; "FDR" means
Benjamini-Hochberg q-value.

| # | Prediction | Test | Falsified if |
|---|---|---|---|
| P1 | Alpha diversity differs between groups (direction not assumed: subgingival studies often report *higher* richness in disease, the source abstract reports *lower*) | Mann-Whitney U on Shannon (primary), Simpson, observed richness, Pielou's evenness; rank-biserial effect size with bootstrap 95% CI | Shannon p >= 0.05 **and** the effect-size CI includes 0 |
| P2 | Periodontitis samples occupy a different region of beta-diversity space | PERMANOVA on Bray-Curtis (primary) and Jaccard, 999 permutations; PERMDISP to check dispersion | PERMANOVA p >= 0.05. If PERMANOVA p < 0.05 but PERMDISP p < 0.05 too, the prediction is only partly supported, because the difference could be spread rather than location |
| P3 | The red-complex taxa are enriched in periodontitis | (a) Targeted: Mann-Whitney U on CLR abundance of *P. gingivalis*, *T. forsythia*, *T. denticola*, BH-corrected across the 3 tests. (b) Untargeted: genus-wide scan, BH-corrected, checking where *Porphyromonas*, *Tannerella*, *Treponema* rank | For each species separately: CLR median not higher in periodontitis, or q >= 0.05 |
| P4 | Community composition predicts disease state better than chance | L1 logistic regression and random forest, repeated stratified 5-fold CV with all preprocessing inside the folds; permuted-label null | Permutation p >= 0.05, or the 95% CI of AUC includes 0.5 |
| P5 (exploratory) | The red-complex genera co-occur (positively correlated across samples) | Spearman on CLR genus abundances, BH-corrected | No positive, FDR-significant edge between any pair of the three genera |

## Secondary (replication) comparison

The dataset also contains hypertensive patients without (`T`, n = 17) and with
(`TP`, n = 16) periodontitis. Re-running P1-P4 on T vs TP asks whether the
signal found in normotensive people shows up again in an independent set of
participants. This is a partial replication only: same study, same lab, same
sequencing run, different people.

## What this project cannot answer

The design is cross-sectional. Even a perfect association cannot tell whether
the community shift causes periodontitis, results from it (e.g. deeper pockets
creating anaerobic niches), or is driven by a shared third factor (smoking,
oral hygiene, age). Nothing here has clinical or diagnostic value.
