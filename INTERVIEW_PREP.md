# Interview preparation

Twenty-two questions I should be able to answer about this project, with short,
accurate answers. Numbers refer to files in `results/`.

---

**1. What was the question, and what did you find?**
Whether salivary bacterial composition differs between periodontally healthy
people and people with periodontitis, and whether composition can predict the
group. In 16 healthy vs 18 periodontitis samples, overall composition
differed modestly (PERMANOVA R² = 0.082, p = 0.015), the three red-complex
pathogens were enriched in periodontitis (q ≤ 0.012), and a classifier
reached AUC 0.92. But the classifier depended mostly on one unidentified
taxon, and none of the community-level signals replicated in a second group
of hypertensive participants.

**2. Why is microbiome data compositional, and why does that matter?**
The sequencer produces a roughly fixed number of reads per sample, so counts
only tell you *proportions*, not how many bacteria were actually there. If one
taxon doubles in absolute terms, every other taxon's share goes down even if
they did not change. This creates spurious negative correlations and false
"decreases". I handled it with a centred log-ratio (CLR) transform, which
looks at each taxon relative to the sample's geometric mean, and with
rank-based tests. CLR reduces the problem but does not remove it.

**3. What does PERMANOVA actually test?**
Whether the group centroids differ in the multivariate space defined by a
distance matrix. It splits the total sum of squared distances into
between-group and within-group parts, forms a pseudo-F ratio, and gets a
p-value by shuffling group labels. It is also sensitive to differences in
spread, so I ran PERMDISP alongside it. Here R² = 0.082 means group explains
about 8% of the variation between samples.

**4. What is PERMDISP and why did you run it?**
It tests whether groups differ in dispersion (the average distance of samples
to their own group centroid). If one group is simply more variable,
PERMANOVA can come out significant even with no shift in location. Here
PERMDISP gave p = 0.065 for Bray-Curtis: not significant, but close enough
that I describe the beta-diversity result as "supported with a caveat".

**5. What is the difference between alpha and beta diversity?**
Alpha diversity describes *within* one sample: how many taxa are present and
how evenly reads are spread (richness, Shannon, Simpson, evenness). Beta
diversity describes *between* samples: how different two communities are
(Bray-Curtis, Jaccard). Groups can differ in beta diversity without differing
in alpha diversity, which is what happened here: the communities contain
different members at similar overall diversity.

**6. Why correct for multiple comparisons?**
I tested 354 genera. At p < 0.05 you expect about 18 false positives by
chance alone, even if nothing is different. The Benjamini-Hochberg procedure
controls the false discovery rate: among the taxa I call significant, the
expected fraction of false ones is at most 5%. In this data, 72 genera had
raw p < 0.05 but only 2 survived FDR.

**7. Why BH and not Bonferroni?**
Bonferroni controls the probability of making *any* false positive, which
with 354 tests requires p < 0.00014 and is very conservative. For exploratory
screening, controlling the *proportion* of false discoveries (FDR) is a
better balance between false positives and missed signals.

**8. Why is a permuted-label control necessary for the classifier?**
It shows what AUC the whole pipeline produces when there is, by
construction, no real signal. If a mistake leaks information (for example,
filtering features using all samples before cross-validation), the null AUCs
come out well above 0.5 and reveal the leak. Here the null averaged 0.49-0.50
and the real AUC beat all 100 shuffles, so the signal is real for this
dataset. The ablation then showed *where* that signal came from.

**9. What is data leakage, and how did you avoid it?**
Leakage means information from the test samples influences training, which
inflates performance. I passed raw counts to a scikit-learn Pipeline in which
the prevalence filter, CLR transform and scaling are re-fitted inside every
training fold, and the L1 penalty was tuned in an inner CV loop (nested CV).
A unit test checks that the filter learns its columns from training rows only.

**10. What does an AUC of 0.92 mean here, and why not trust it fully?**
If you pick one periodontitis and one healthy sample at random, the model
ranks the periodontitis sample higher 92% of the time. With 34 samples the
bootstrap CI is 0.81-1.00, but I don't lean on that interval: bootstrapping
out-of-fold predictions from 34 samples and hundreds of candidate features
ignores the variability of feature selection and tuning, so it tends to be
too narrow. The stronger evidence is the permuted-label null (the same
pipeline on shuffled labels averages 0.49) and the ablation: removing one
feature ("Unclassified Bacilli") drops the AUC to 0.64. The model also does
not transfer to the hypertensive cohort (AUC 0.52). So it is a
dataset-specific association, not a general signature. The driver taxon
could not be identified, because the public deposit has no sequences.

**11. What are the limits of 16S sequencing compared with shotgun metagenomics?**
16S amplifies one ~470 bp region of one gene. It usually resolves genus but
not reliably species or strain, it says nothing about function (for example,
virulence genes), it has PCR and primer biases, and copy number varies
between species. Species calls also depend on the reference database: here
the database is not stated (the names look like SILVA 138), whereas eHOMD is
the curated reference for oral bacteria, so I treat species labels as less
certain than genera. Shotgun sequencing reads all DNA, which gives species and
strain resolution plus gene content. It costs more, needs more depth, and in
saliva much of the DNA is human.

**12. Why did *P. gingivalis* come out significant but the genus *Porphyromonas* did not?**
The genus includes commensal species such as *P. pasteri*, which in this
dataset is almost as abundant as *P. gingivalis*. Pooling them at genus level
dilutes the disease-associated species. Taxonomic resolution can change the
answer.

**13. Isn't testing only three species after a big scan cherry-picking?**
It would be if I had picked them after seeing the results. They were named in
`HYPOTHESIS.md` before any analysis, the multiple-testing correction covers
exactly those three pre-declared tests, and the untargeted 354-genus scan is
reported alongside.

**14. What was the data-quality problem, and how did you deal with it?**
22,781 of 25,540 features appear in only one sample, several samples have
five times the usual number of features, and many belong to soil or gut
lineages, which points to likely contamination or batch effects. I analysed
at genus level (99% of reads kept after filtering), flagged samples with more
than half their reads in one-off features, ran sensitivity analyses without
them, and showed which findings were carried by the high-richness subset. I
did not silently drop data.

**15. Why did you rarefy for alpha diversity but use CLR elsewhere? Isn't rarefying frowned upon?**
Rarefying is contested: McMurdie and Holmes (2014, PLoS Computational
Biology, "Waste Not, Want Not: Why Rarefying Microbiome Data Is
Inadmissible") showed it throws away data and power, especially for
differential abundance. So I use it only where unequal depth directly biases
the quantity being measured: richness, since a more deeply sequenced sample
simply shows more taxa. Every other analysis keeps all reads and uses
proportions (Bray-Curtis) or CLR log-ratios. Here the cost is small because
the minimum depth (44,786 reads) is well past where the rarefaction curves
flatten.

**16. Why didn't the result replicate in T vs TP?**
I can't tell which of three explanations is right: (1) the primary effect is
inflated or partly false, because small studies overestimate effects;
(2) hypertension or its medication changes the oral microbiome and masks the
periodontitis signal. The source study itself reports that hypertension-only
patients already show periodontitis-like dysbiosis, so T is not a healthy
reference; (3) the primary difference was partly driven by the
unusual healthy samples. The red-complex effects still point the same way,
just more weakly.

**17. Can you say these bacteria cause periodontitis?**
No. The data are cross-sectional (one sample per person), so the shift could
come before the disease, result from it (deeper pockets create anaerobic
niches), or be caused by a third factor such as oral hygiene, diet or age.
(Smokers were excluded by the study, so smoking is not a confounder here.)
All findings are associations.

**18. What does the co-occurrence network show?**
Genera whose CLR abundances rise and fall together across people (|rho| ≥ 0.6,
FDR < 0.05). The three red-complex genera correlate strongly with each other
(rho 0.68-0.85) and are hubs together with *Campylobacter*, *Fusobacterium*
and *Prevotella*, which fits the idea of an anaerobic consortium. Caveats:
correlation is not interaction, and pooling both groups creates some
correlation just from the group difference.

**19. How do you know your pipeline could find a real effect if one existed?**
I generated Dirichlet-multinomial data with 15 planted differences and ran
the same code. It recovered 10 of 15 (sensitivity 67%) in the right direction
with 1 false positive (9%). That is part of `run_all.py` and a unit test.

**20. What would you do next with more data?**
Re-process the raw reads with DADA2 and a decontamination step, and
re-annotate them against eHOMD; identify the "Unclassified Bacilli"
sequences; confirm differential abundance with ANCOM-BC
or ALDEx2; combine several public periodontitis cohorts and test whether a
classifier trained on one study works on another; move to subgingival plaque
and shotgun metagenomics for species and function; and sample people before
and after periodontal treatment to study change over time.

**21. Why do Chao1 and Shannon disagree here?**
They measure different things. Chao1 (and ACE) estimate how many taxa are
present, including unseen ones, and are driven by rare taxa: singletons and
doubletons. Shannon weights taxa by abundance, so a handful of extra rare
taxa barely moves it. In this data the healthy group's median Chao1 is higher
(182.2 vs 141.4 genera, p = 0.133) while Shannon is essentially equal
(p = 0.796), which means the difference lies in rare taxa, not in the
dominant community. The source study saw the same split (Chao1/ACE differ,
Shannon/Simpson do not). A post-hoc check showed the healthy group's extra
richness comes from six unusually rich healthy samples; without them the
medians are 139.8 vs 141.4.

**22. Your result differs from the published paper's. Explain.**
Less than it first appears. The paper reports higher Chao1/ACE in healthy
controls and no Shannon/Simpson difference; I find no Shannon/Simpson
difference either, and Chao1/ACE lower in periodontitis in the same direction,
just not significant (p = 0.133, q = 0.200 after correcting across six
metrics). We also agree that *P. gingivalis*, *T. forsythia* and
*T. denticola* are enriched in periodontitis. The differences are
methodological: they ran a four-group ANOVA with Tukey post-hoc tests on 97%
OTUs; my primary test is healthy vs periodontitis only, with Mann-Whitney U,
at genus level and with FDR correction. They used weighted UniFrac; I used
Bray-Curtis and Jaccard, because the deposit has no sequences to build a
tree. I don't claim their analysis is wrong. The honest statement is that the
richness difference is sensitive to the contrast, the test and six unusual
healthy samples.
