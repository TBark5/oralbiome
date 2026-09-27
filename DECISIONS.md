# Decisions log

Every judgement call made while building the project, with the reason. Numbers
refer to constants in `src/oralbiome/config.py`.

## Phase 0 - Setup

- **D1. Python version.** The machine's existing `.venv` uses Python 3.14, so all
  packages are pinned to versions that install on it (`requirements.txt`). The
  code only uses standard library features available in 3.11+.
- **D2. No scikit-bio.** Diversity indices, Bray-Curtis/Jaccard, PCoA, PERMANOVA
  and PERMDISP are implemented directly with NumPy/SciPy (about 150 lines). This
  avoids a heavy dependency that often lacks wheels for new Python versions, and
  it means every formula in METHODS.md maps to code I can point at.
- **D3. Removed the PyCharm template `main.py`.** It was the default "print_hi"
  sample file and unrelated to the project.
- **D4. Colour palettes.** Groups use Okabe-Ito blue (#0072B2) and vermillion
  (#D55E00); a colour-vision-deficiency check gives a worst-case (protan) colour
  difference of 21.9 (target >= 8). Stacked bars show the top 8 taxa in a
  validated 8-hue categorical order plus a grey "Other" band.
- **D5. Package layout.** Code lives in `src/oralbiome`, installed in editable mode
  (`pip install -e .`) so tests, `run_all.py` and `app.py` import the same code.

## Phase 1 - Data

- **D6. Dataset choice.** Searched for public healthy-vs-periodontitis 16S tables.
  Most studies deposit only raw reads (SRA), which would need a full DADA2 run.
  The Guo et al. figshare deposit has a ready-made ASV table, per-sample group
  labels and a CC BY 4.0 license, and it contains the three red-complex species,
  so it was chosen. Saliva (not subgingival plaque) is a limitation: it dilutes
  the pocket community.
- **D7. Primary comparison H vs P only.** T and TP carry a second condition
  (hypertension). Pooling them would mix two factors, so they are used only as a
  replication set (T vs TP).
- **D8. Atypical-sample flag.** The QC report showed a quality problem: 22,781 of
  25,540 features occur in only one sample, and several samples have 1,700-2,300
  features against a median of about 350, many from soil- or gut-associated
  lineages (Acidobacteriota, Chloroflexi, Gemmatimonadota, Muribaculaceae).
  A sample is flagged when more than 50% of its reads come from features seen in
  no other sample. Flagged: H7, P7, P9, P18 (primary) and T9, T10, TP7, TP15
  (replication). The 50% cut-off is a round number chosen before looking at
  results, not tuned. Flagged samples are **kept** in the primary analysis
  (removing data after looking at it invites bias) and **removed** in a
  sensitivity analysis.
- **D9. Taxonomy cleanup.** Names like "uncultured", "unclassified_X",
  "Unassigned" are replaced with "Unclassified <nearest named higher rank>", so
  unrelated unknown features do not merge into one fake genus.

## Phase 2 - M1 preprocessing and M2 composition

- **D10. Genus is the main unit of analysis.** At feature (ASV) level only a
  median 79% of reads survive a prevalence filter because features are very
  sample-specific. At genus level the median is 99%. Genus-level data are also
  more robust to sequencing error and easier to compare with the literature.
  Feature-level alpha diversity is reported as a sensitivity check.
- **D11. Prevalence filter: >= 10% of samples and mean relative abundance
  >= 0.01%.** With 34 samples, 10% means a genus must appear in at least 4
  people. That removes one-off genera that cannot be tested statistically and
  would only add noise and multiple-testing burden, while keeping low-abundance
  taxa that are consistently present (the red-complex genera are
  low-abundance). Result: 1,648 -> 354 genera, median 99% of reads kept.
- **D12. Filtering is not applied before alpha diversity.** Richness counts rare
  taxa, so removing them would bias it. Alpha diversity uses the unfiltered,
  rarefied genus table.
- **D13. Rarefaction depth = minimum sample depth (44,786 reads).** No sample is
  lost. Rarefaction curves are nearly flat there (median slope 0.1 new features
  per 1,000 reads). One random subsample with seed 42 is used.
- **D14. Zeros: pseudocount 0.5 before the CLR transform.** This is the simplest
  standard choice. It is documented as a limitation because the result depends
  slightly on the pseudocount value.
- **D15. Core microbiome = genus detected (count > 0) in >= 90% of a group's
  samples.** A common convention; detection is on raw (unrarefied) counts.
- **D16. Stacked bars show the top 8 taxa by overall mean plus "Other".** More
  than 8 colours cannot be told apart reliably, especially by colour-blind
  readers.

## Phase 3 - M3 alpha and M4 beta diversity

- **D17. Alpha metrics.** Shannon uses the natural log. "Simpson" is reported as
  Gini-Simpson (1 - sum p^2), so higher means more diverse, like the others.
  Shannon is the pre-declared primary metric (HYPOTHESIS.md). BH correction is
  applied across the four metrics inside each analysis.
- **D18. Effect size for Mann-Whitney = rank-biserial r** (equal to Cliff's
  delta), with a 2,000-resample percentile bootstrap CI resampled within each
  group. Positive r means higher in periodontitis.
- **D19. Alpha sensitivity analyses:** (a) feature (ASV) level instead of genus,
  (b) genus level without the flagged samples. The primary result is not
  replaced by either one.
- **D20. Beta diversity input.** Bray-Curtis uses relative abundances of the 354
  prevalence-filtered genera. Jaccard uses their presence/absence. Bray-Curtis
  is primary (pre-declared).
- **D21. 999 permutations for PERMANOVA and PERMDISP.** The smallest possible
  p-value is 0.001, and both tests together run in about 3 s.
- **D22. PERMDISP p-value** comes from permuting group labels on the
  distance-to-centroid values. This is slightly simpler than vegan's
  residual permutation, but for a two-group, one-factor design the two give
  very similar answers. Negative eigenvalues from Bray-Curtis are handled as in
  Anderson (2006): their axes subtract from squared distances.
- **D23. NMDS** uses scikit-learn non-metric MDS on the precomputed distance
  matrix with 8 random starts (seed 42), and reports Kruskal stress-1.
- **D24. Ellipses** are 95% normal-theory data ellipses (chi-square, 2 df) fitted
  to each group's points. They are for visual guidance only; the tests are
  PERMANOVA and PERMDISP.

## Phase 4 - M5 differential abundance

- **D25. Test = two-sided Mann-Whitney U on CLR abundances, per genus, BH-FDR
  at q < 0.05.** Rank-based, so it makes no normality assumption and needs no
  extra package. Tools such as ANCOM-BC, ALDEx2 or MaAsLin2 are R packages;
  choosing a simple, transparent test was preferred (rule "prefer standard,
  well-documented methods"). The weakness: CLR is still affected by
  compositional effects, which ANCOM-BC corrects more carefully.
- **D26. Effect sizes:** CLR mean difference with a 95% Welch t interval
  (magnitude on a log scale) plus rank-biserial r (probability-based).
- **D27. Null control:** group labels shuffled 200 times and the full
  scan + BH re-run each time. Also reported: the expected number of raw
  p < 0.05 hits by chance (0.05 x number of genera).
- **D28. Red-complex test is pre-declared and separate.** The three species were
  named in HYPOTHESIS.md before analysis. They are tested on a species-level
  CLR table (species prevalence filter, targets always kept), and BH is
  applied across these 3 tests only. The genus-wide scan (354 tests) is
  reported alongside so the reader can see how the genera rank without any
  prior.
- **D29. High-richness check (added after seeing the M3/M4 figures).** Six
  healthy samples have > 1,000 observed features at rarefaction depth. For the
  top 10 genera, `results/m5_high_richness_check.csv` counts detection inside
  vs outside that subset. This is a descriptive diagnostic, not a new
  exclusion rule; it was added because the healthy-enriched hits looked like
  environmental taxa.
- **D30. Volcano labels:** top 8 genera by p plus the three red-complex genera,
  placed by a simple collision-avoiding routine (`style.place_labels`).

## Phase 5 - M6 network, M7 classifier, replication

- **D31. Network input: genera present in >= 30% of samples (217 genera).** With
  34 samples, Spearman correlations between mostly-zero taxa are driven by
  shared zeros, so rarer genera are left out of the network (they stay in M5).
- **D32. Edge rule: |rho| >= 0.6 and BH q < 0.05 over all 23,436 pairs.** Both
  a strength and a significance cut-off, so edges are neither weak nor
  chance-level. Samples of both groups are pooled, so two taxa that are both
  higher in periodontitis will correlate partly *because of* the group
  difference; the network describes co-variation across all people, not
  interactions within a person.
- **D33. Hubs = highest degree** (number of edges), with betweenness
  centrality reported as a tie-breaker. Degree is the simplest, most
  explainable centrality measure.
- **D34. Network figure** shows the largest connected component (167 of 177
  genera, 1,418 of 1,426 edges) with a Kamada-Kawai layout; hubs are numbered
  with a key because names do not fit in the dense core.
- **D35. Classifier input is raw genus counts (all 1,648 genera).** Prevalence
  filtering and CLR run inside each training fold (`PrevalenceCLR`
  transformer), and scaling too, so nothing learned from test samples leaks
  into training.
- **D36. Nested CV for L1 logistic regression.** C is tuned over 8 values
  (0.01-31.6) by 3-fold inner CV on each outer training set (about 27
  samples, so 3 inner folds keep about 9 test samples per inner fold). Class
  weights are balanced. Solver: liblinear.
- **D37. Random forest: 300 trees, sqrt features, balanced class weights, no
  tuning.** Random forests work well at default settings; tuning them on 27
  samples would mostly fit noise and cost runtime.
- **D38. Outer CV: 5-fold stratified, repeated 10 times.** Two AUCs are saved:
  (a) the mean of the 10 per-repeat AUCs, compared with the null; (b) the AUC
  of per-sample probabilities averaged over the repeats, with a 2,000-resample
  stratified bootstrap 95% CI and ROC band. They differ slightly (e.g. RF 0.79
  vs 0.81), and the README says which is which.
- **D39. Permutation null: 100 label shuffles, 3 CV repeats each** (to fit the
  2-minute budget). A 3-repeat mean varies more than a 10-repeat mean, which
  widens the null distribution and makes the test slightly *conservative*.
  With 100 shuffles the smallest possible p-value is 1/101 = 0.0099.
- **D40. Ablation (added after seeing the importance plot).** One genus
  ("Unclassified Bacilli") was selected in every L1 fold, so the classifiers
  were re-run without it (5 repeats). This is a post-hoc diagnostic of how
  much the AUC rests on one feature, not a model-selection step.
- **D41. Replication = T vs TP** through the same code (M1 filtering, alpha,
  Bray-Curtis PERMANOVA/PERMDISP, red-complex test, M7 with its own null).
  Effect sizes for both cohorts are compared in a forest plot.
