# Progress

## Done
- Phase 0: venv, pinned requirements, git, HYPOTHESIS.md, DECISIONS.md, plotting style (`src/oralbiome/style.py`).
- Phase 1: real dataset (Guo et al., figshare, CC BY 4.0) downloaded and cached in `data/raw/`; parser, synthetic fallback, QC report (`results/data_quality.json`, `figures/00_data_quality.png`); atypical-sample flag.
- Phase 2: M1 preprocessing, M2 composition + core microbiome.
- Phase 3: M3 alpha diversity (+2 sensitivity analyses), M4 beta diversity (Bray-Curtis/Jaccard, PCoA, NMDS, PERMANOVA, PERMDISP, heatmap).
- Phase 4: M5 differential abundance (genus scan + BH, permutation null, pre-declared red-complex species test, high-richness diagnostic).
- Phase 5: M6 co-occurrence network, M7 classifiers (nested CV, permutation null, ablation), replication cohort T vs TP.
- Phase 6: visuals pass (16 figures at 300 dpi, colorblind-safe palette, n on every group figure), figures/CAPTIONS.md.
- Phase 7: Streamlit dashboard `app.py` (9 tabs, sidebar filters: groups, flagged samples, taxonomic level, top-N); verified with AppTest and a live launch; screenshots in docs/screenshots/
- Phase 8: 36 pytest tests (formulas, statistics, leakage, synthetic recovery, data layer, dashboard); fresh clone + fresh venv reproduced all outputs byte for byte.

## In progress
- Phase 9: documentation.

## Next
- Final checks and MORNING_REPORT.md.

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
