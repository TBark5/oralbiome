# Progress

## Done
- Phase 0: venv, pinned requirements, git, HYPOTHESIS.md, DECISIONS.md, plotting style (`src/oralbiome/style.py`).
- Phase 1: real dataset (Guo et al., figshare, CC BY 4.0) downloaded and cached in `data/raw/`; parser, synthetic fallback, QC report (`results/data_quality.json`, `figures/00_data_quality.png`); atypical-sample flag.
- Phase 2: M1 preprocessing, M2 composition + core microbiome.
- Phase 3: M3 alpha diversity (+2 sensitivity analyses), M4 beta diversity (Bray-Curtis/Jaccard, PCoA, NMDS, PERMANOVA, PERMDISP, heatmap).
- Phase 4: M5 differential abundance (genus scan + BH, permutation null, pre-declared red-complex species test, high-richness diagnostic).
- Phase 5: M6 co-occurrence network, M7 classifiers (nested CV, permutation null, ablation), replication cohort T vs TP.

## In progress
- Phase 6: visuals pass + figures/CAPTIONS.md.

## Next
- Phase 7 app, Phase 8 tests, Phase 9 docs.

## Known bugs
- None.

## How to run
```
.venv/Scripts/python -m pip install -r requirements.txt && .venv/Scripts/python -m pip install -e .
.venv/Scripts/python run_all.py
```

## Key files
- `run_all.py` - full pipeline
- `src/oralbiome/` - one module per analysis step
- `DECISIONS.md` - every threshold and judgement call
