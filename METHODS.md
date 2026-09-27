# Methods and formulas

This is a study guide for the formulas the code implements. Each section
names the function that computes it. Notation: a sample has S observed taxa
with counts x_1..x_S, total N = sum x_i, and proportions p_i = x_i / N.

## 0. Where the taxonomy comes from

The taxonomy strings are taken as deposited. The reference database and the
classification confidence threshold are **not stated** in the figshare
deposit, the article or its supplementary files. The names follow SILVA
138-era conventions (for example `Prevotella_7`, `*_UCG_*` genera,
`Actinobacteriota`), not HOMD/eHOMD, but one phylum spelling
(`Campylobacterota`) differs from SILVA 138, so the release cannot be
confirmed. Because the deposit contains no sequences, the features cannot be
re-classified here. Species-level labels should therefore be treated as less
certain than genus-level ones. Details: DATA_SOURCE.md.

## 1. Preprocessing (`preprocessing.py`)

**Relative abundance:** p_i = x_i / N.

**Prevalence filter:** keep taxon j if
(number of samples with x_j > 0) / (number of samples) >= 0.10
and mean over samples of p_j >= 0.0001.

**Rarefaction:** draw exactly d reads from each sample *without replacement*
(multivariate hypergeometric draw), where d = the smallest sample depth. This
gives every sample the same sequencing effort, which matters for richness.

*Why rarefaction is used here, and only here.* Rarefying is contested:
McMurdie and Holmes (2014), "Waste Not, Want Not: Why Rarefying Microbiome
Data Is Inadmissible" (PLoS Computational Biology 10(4):e1003531), showed that
throwing away reads to reach a common depth loses statistical power and
adds random noise, especially for differential abundance testing. This
pipeline therefore uses rarefied counts **only for alpha diversity**, where
unequal sequencing depth directly biases richness (a sample sequenced twice
as deeply will show more taxa). Every other analysis (beta diversity,
differential abundance, network, classifier) uses all reads, either as
relative abundances or as CLR-transformed values. Rarefaction is a single
seeded draw (seed 42) to the minimum depth of 44,786 reads, so no sample is
dropped.

**Expected richness at depth n** (rarefaction curves, `qc.py`, Hurlbert 1971):

E[S_n] = sum over taxa i of [ 1 - C(N - x_i, n) / C(N, n) ]

where C(a, b) is "a choose b". Each term is the probability that taxon i
appears at least once in a random subsample of n reads.

**Centred log-ratio (CLR):** with pseudocount c = 0.5,

clr(x)_i = ln(x_i + c) - (1/D) * sum over j of ln(x_j + c)

where D is the number of taxa. The second term is the log of the geometric
mean, so the CLR values of a sample sum to zero. CLR turns ratios into
differences on a log scale, which is why it is used before correlation,
t-intervals and classifiers.

## 2. Alpha diversity (`alpha.py`)

**Shannon index:** H' = - sum over i of p_i ln(p_i). It increases with both
the number of taxa and how evenly reads are spread. Maximum = ln(S), when all
taxa are equally abundant.

**Simpson:** D = sum over i of p_i^2 is the probability that two reads drawn
at random (with replacement) belong to the same taxon. The code reports the
**Gini-Simpson index 1 - D**, the probability that they belong to *different*
taxa, so that higher = more diverse.

**Observed richness:** S = number of taxa with x_i > 0.

**Pielou's evenness:** J' = H' / ln(S). It ranges from 0 (one taxon dominates)
to 1 (perfectly even).

**Chao1 (bias-corrected):** with F1 = number of taxa seen exactly once and F2
= number seen exactly twice,

Chao1 = S_obs + F1 (F1 - 1) / (2 (F2 + 1)).

The idea: if many taxa are seen only once, many more were probably missed.

**ACE (Abundance-based Coverage Estimator, Chao & Lee 1992):** taxa with at
most 10 reads are "rare". With S_rare rare taxa holding N_rare reads, S_abund
taxa with more than 10 reads, and F_i taxa seen exactly i times:

- sample coverage of rare taxa: C = 1 - F1 / N_rare
- gamma^2 = max( (S_rare / C) * sum_{i=1..10} i (i - 1) F_i / (N_rare (N_rare - 1)) - 1, 0 )
- ACE = S_abund + S_rare / C + (F1 / C) * gamma^2

If every rare read is a singleton, C = 0 and ACE is undefined; the code then
returns Chao1. This never happened in the real data (checked for every sample, at
genus and feature level, in both cohorts).

**Why richness estimators and Shannon can disagree.** Chao1, ACE and observed
richness count taxa and are driven by rare taxa (singletons and doubletons).
Shannon and Simpson weight taxa by abundance, so adding many rare taxa raises
Chao1 a lot and Shannon hardly at all. Disagreement between them is expected,
not contradictory. All six metrics are tested and BH-corrected together.

## 3. Beta diversity (`beta.py`)

For two samples A and B with relative abundances a_i and b_i:

**Bray-Curtis dissimilarity:**

BC(A, B) = sum_i |a_i - b_i| / sum_i (a_i + b_i)

Equivalently, 1 - 2 * sum_i min(a_i, b_i) / sum_i (a_i + b_i). It is 0 for
identical profiles and 1 when no taxa are shared. It is weighted by
abundance, so dominant taxa drive it.

**Jaccard distance** (presence/absence): with A and B as the sets of taxa present,

J(A, B) = 1 - |A ∩ B| / |A ∪ B|

It ignores abundance, so rare taxa count as much as dominant ones.

**PCoA (principal coordinates analysis):** from the n x n distance matrix D,

1. A = -1/2 * D^2 (element-wise square),
2. B = C A C, where C = I - (1/n) 1 1^T is the centring matrix,
3. eigen-decompose B; coordinates on axis k = eigenvector_k * sqrt(eigenvalue_k).

"% variance explained" by axis k = eigenvalue_k / sum of positive eigenvalues.
Bray-Curtis is not Euclidean, so some eigenvalues are negative; they are
dropped for plotting.

**NMDS:** finds 2-D positions whose *rank order* of distances matches the rank
order of the original dissimilarities as closely as possible. Quality is
measured by Kruskal stress-1: below 0.1 is good and below 0.2 is usable.

## 4. PERMANOVA (`beta.permanova`, Anderson 2001)

With N samples in a groups, where group g has n_g samples, and d_ij the
distance between samples i and j:

- Total sum of squares: SS_T = (1/N) * sum over all pairs i<j of d_ij^2
- Within-group sum of squares: SS_W = sum over groups g of (1/n_g) * sum over pairs i<j in g of d_ij^2
- Between-group sum of squares: SS_A = SS_T - SS_W

Pseudo-F statistic:

F = [ SS_A / (a - 1) ] / [ SS_W / (N - a) ]

R² = SS_A / SS_T is the fraction of total variation explained by group.

p-value: shuffle the group labels 999 times, recompute F each time, and

p = (1 + number of permuted F >= observed F) / (1 + 999).

The "+1" counts the observed labelling as one of the possible permutations,
so p can never be 0.

**What PERMANOVA tests:** the null hypothesis is that group centroids are the
same in distance space. PERMANOVA is also sensitive to differences in spread
(dispersion), which is why PERMDISP is run as well.

## 5. PERMDISP (`beta.permdisp`, Anderson 2006)

1. Place all samples in PCoA space (all axes).
2. For each sample, compute its distance z_i to its own group's centroid.
   Axes with negative eigenvalues contribute negatively to the squared
   distance.
3. Run a one-way ANOVA F-test on the z_i values; get the p-value by shuffling
   group labels 999 times.

A significant PERMDISP means the groups differ in spread. If PERMANOVA is
significant *and* PERMDISP is not, the difference is more likely a real
location shift.

## 6. Mann-Whitney U and effect sizes (`stats.py`)

**Mann-Whitney U:** for groups X (n_x) and Y (n_y),
U_x = number of pairs (x, y) with x > y, plus 0.5 per tie. Under the null,
X and Y come from the same distribution. It uses ranks only, so it needs no
normality assumption.

**Rank-biserial correlation:**

r = 2 U_x / (n_x n_y) - 1 = P(X > Y) - P(X < Y)

It ranges from -1 to +1, where 0 means no tendency either way. This is the
effect size reported next to every Mann-Whitney p-value. Its 95% CI is a
percentile bootstrap: resample each group with replacement 2,000 times,
recompute r, and take the 2.5th and 97.5th percentiles.

**Welch CI for a mean difference** (CLR difference in M5):

diff = mean(X) - mean(Y), SE = sqrt(s_x^2/n_x + s_y^2/n_y),
CI = diff ± t_{0.975, df} * SE, with Welch-Satterthwaite df.

## 7. Benjamini-Hochberg FDR (`stats.benjamini_hochberg`)

For m p-values sorted ascending, p_(1) <= ... <= p_(m):

1. Compute p_(i) * m / i for each rank i.
2. Going from the largest rank down, replace each value by the minimum of
   itself and all values above it (this keeps q monotone).
3. Cap at 1. These are the q-values.

Calling everything with q < 0.05 "significant" keeps the *expected* fraction
of false discoveries among the calls at or below 5%. This is different from
Bonferroni (p < 0.05/m), which controls the chance of *any* false positive
and is much stricter.

## 8. Spearman correlation network (`network.py`)

Spearman rho is the Pearson correlation of the ranks. Correlations are
computed on CLR values across samples for genera present in >= 30% of
samples; each correlation has a p-value (t approximation) and BH is applied
over all 23,436 pairs. An edge is kept if |rho| >= 0.6 and q < 0.05.
**Degree** is the number of edges of a node; the highest-degree nodes are
called hubs. **Betweenness** is the fraction of shortest paths between other
nodes that pass through a node.

## 9. Classifiers (`classifier.py`)

**L1 (lasso) logistic regression:** models P(periodontitis) = 1 / (1 + e^(-(b0 + b·x))),
fitted by minimising

-log-likelihood + (1/C) * sum_j |b_j|.

The absolute-value (L1) penalty pushes many coefficients to exactly zero,
which acts as feature selection. C is chosen by 3-fold CV inside each
training set (nested CV).

**Random forest:** 300 decision trees, each grown on a bootstrap sample and
considering sqrt(p) random features per split; the prediction is the average
vote. Importance = mean decrease in Gini impurity.

**How much to trust the AUC confidence interval.** The reported 95% CI
bootstraps the out-of-fold predicted probabilities (averaged over the 10 CV
repeats). That captures sampling variability of the 34 samples but not the
variability of the whole modelling procedure (feature filtering, penalty
tuning, feature selection among hundreds of genera), and repeated CV splits
of 34 samples are strongly correlated. With n = 34 and 354+ candidate
features, bootstrap CIs on a cross-validated AUC therefore tend to be
optimistic (too narrow). The permuted-label null (same full pipeline on
shuffled labels) and the ablation (removing the top feature) are the stronger
evidence, and the README says so beneath the classifier table.

**ROC AUC:** the probability that a randomly chosen periodontitis sample gets
a higher predicted probability than a randomly chosen healthy sample. It is
mathematically the same quantity as the Mann-Whitney U statistic divided by
n_x n_y. 0.5 = chance, 1.0 = perfect ranking.

**Repeated stratified k-fold CV:** split samples into 5 folds with the same
group ratio in each, train on 4 and predict the 5th, rotate, and repeat with
10 different shuffles. Every prediction is made on a sample the model did not
see, and all preprocessing (prevalence filter, CLR, scaling) is re-fitted on
the training folds only.

**Permuted-label null:** shuffle y, rerun the entire CV pipeline, and record
the AUC; repeat 100 times. p = (1 + number of null AUCs >= observed) / 101.

## 10. Dirichlet-multinomial simulation (`synthetic.py`)

For group g with mean composition pi_g, each sample draws p ~ Dirichlet(theta * pi_g),
then counts ~ Multinomial(depth, p). theta (300 here) controls how much
samples vary around the group mean: smaller theta gives more overdispersion,
which is typical of microbiome data. Planted differences multiply selected
entries of pi by 4 (or 1/4) in the disease group before renormalising.
