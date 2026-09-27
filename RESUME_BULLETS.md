# Resume bullets

Every variant leads with the replication failure, not the AUC. The negative
result is the stronger and more honest signal: it shows the controls did
their job.

## Short (one line)

- Built a reproducible 16S microbiome pipeline testing pre-registered hypotheses on a public periodontitis cohort (n=67); implemented PERMANOVA/PERMDISP from scratch, nested cross-validation with permuted-label controls, and an internal replication cohort, finding that the disease classifier (AUC 0.92) failed to replicate and depended on a single unidentified taxon.

## Medium (three bullets)

- Showed that an apparently strong microbiome classifier of periodontitis (AUC 0.92 on 34 saliva samples) did not replicate in an independent cohort from the same study (AUC 0.52, n=33) and relied on one unidentified taxon (AUC 0.64 without it), using pre-registered hypotheses, permuted-label controls and ablation.
- Built the full analysis in Python on a public 16S rRNA dataset (n=67): diversity (including Chao1/ACE), PERMANOVA/PERMDISP implemented from scratch, FDR-controlled differential abundance, co-occurrence networks and nested cross-validation; confirmed pre-declared enrichment of three red-complex pathogens (FDR q ≤ 0.012).
- Compared results with the published analysis of the same data, documented where they agree and why they differ, and shipped a tested codebase (pytest suite, synthetic positive control, Streamlit dashboard) whose every reported number regenerates from raw data with one command.

## Technical (detailed)

- Tested a salivary microbiome signature of periodontitis on a public 16S dataset (n=67) with a pre-registered hypothesis file, an internal replication cohort (hypertensive participants with and without periodontitis) and negative controls. The L1 classifier's AUC of 0.92 [bootstrap CI 0.81-1.00, likely optimistic at n=34] beat 100 label permutations (p = 0.0099) but fell to 0.52 in the replication cohort and to 0.64 when one unidentified taxon was removed.
- Implemented Bray-Curtis/Jaccard distances, PCoA, PERMANOVA and PERMDISP from scratch in NumPy/SciPy, with unit tests against hand-calculated cases and classical ANOVA; added Chao1 and ACE to reconcile the results with the source publication.
- Built leakage-free nested cross-validation (in-fold prevalence filtering and CLR transform via a custom scikit-learn transformer) for L1 logistic regression and random forest, with permuted-label nulls and feature ablation.
- Found and documented a contamination or batch signal during QC (22,781 of 25,540 features seen in only one sample) and showed that the healthy group's apparent richness advantage depended on six atypical samples; audited provenance (unstated taxonomy database, data license, citations).
- Validated the pipeline on Dirichlet-multinomial data with planted effects (sensitivity 67%, false discovery proportion 9%); outputs reproduce byte for byte from a fresh environment.
- Tools: Python, pandas, NumPy, SciPy, scikit-learn, NetworkX, matplotlib, Altair, Streamlit, pytest, git.
