# Resume bullets

## Short (one line)

- Built a reproducible Python pipeline analysing public 16S rRNA saliva data (34 people) to test whether oral bacterial composition is associated with periodontitis, including permutation-based null controls and a replication cohort.

## Medium (three bullets)

- Designed a hypothesis-driven microbiome analysis of a public 16S rRNA dataset (healthy n=16 vs periodontitis n=18), covering diversity, PERMANOVA, FDR-controlled differential abundance, co-occurrence networks and machine-learning classification in Python.
- Confirmed pre-registered enrichment of three red-complex periodontal pathogens (FDR q ≤ 0.012) and a classifier AUC of 0.92 that beat a permuted-label null (p = 0.0099), then showed with ablation and a second cohort (n=33) that the classifier depended on one taxon and did not replicate.
- Delivered a tested codebase (36 pytest tests, synthetic positive control), a Streamlit dashboard and documentation, with every reported number regenerated from raw data by one command.

## Technical (detailed)

- Implemented Bray-Curtis/Jaccard distances, PCoA, PERMANOVA and PERMDISP from scratch in NumPy/SciPy, with unit tests checking them against hand-calculated cases and classical ANOVA.
- Built leakage-free nested cross-validation for L1 logistic regression and random forest (in-fold prevalence filtering and CLR transform via a custom scikit-learn transformer), with bootstrap AUC confidence intervals, 100-permutation null models and feature ablation.
- Found and documented a contamination or batch signal during QC (22,781 of 25,540 features seen in only one sample) and handled it with genus-level aggregation, pre-specified sample flags and sensitivity analyses instead of silent exclusion.
- Validated the pipeline on Dirichlet-multinomial synthetic data with planted effects (sensitivity 67%, false discovery proportion 9%) and made all outputs byte-for-byte reproducible from a fresh environment.
- Tools: Python, pandas, NumPy, SciPy, scikit-learn, NetworkX, matplotlib, Altair, Streamlit, pytest, git.
