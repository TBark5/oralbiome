# Progress

## Done
- Phase 0: venv, pinned requirements, git, HYPOTHESIS.md, DECISIONS.md, plotting style (`src/oralbiome/style.py`).
- Phase 1: real dataset (Guo et al., figshare, CC BY 4.0) downloaded and cached in `data/raw/`; parser, synthetic fallback, QC report (`results/data_quality.json`, `figures/00_data_quality.png`); atypical-sample flag.
- Phase 2: M1 preprocessing, M2 composition + core microbiome.
- Phase 3: M3 alpha diversity (+2 sensitivity analyses), M4 beta diversity (Bray-Curtis/Jaccard, PCoA, NMDS, PERMANOVA, PERMDISP, heatmap).
- Phase 4: M5 differential abundance (genus scan + BH, permutation null, pre-declared red-complex species test, high-richness diagnostic).
- Phase 5: M6 co-occurrence network, M7 classifiers (nested CV, permutation null, ablation), replication cohort T vs TP.
- Phase 6: visuals pass (16 figures at 300 dpi, colorblind-safe palette, n on every group figure), figures/CAPTIONS.md.
- Phase 7: Streamlit dashboard `app.py` (9 tabs, sidebar filters: groups, flagged samples, taxonomic level, top-N); verified with AppTest and a live launch; screenshots in docs/screenshots/.
- Phase 8: 36 pytest tests (formulas, statistics, leakage, synthetic recovery, data layer, dashboard); a fresh clone + fresh venv reproduced all outputs byte for byte.
- Phase 9: README.md, RESULTS_DISCUSSION.md, METHODS.md, INTERVIEW_PREP.md, RESUME_BULLETS.md, LINKEDIN_POST.md, MORNING_REPORT.md.
- Final checks: no AI attribution in git log or tracked files; no placeholders; README numbers audited against results/.

## In progress
- Nothing.

## Next (optional extensions, not required)
- Re-process raw reads (SRA PRJNA1304526) with DADA2 + decontamination.
- Cross-check M5 with ANCOM-BC/ALDEx2 (R) and BLAST the "Unclassified Bacilli" features.

## Known bugs
- None. Note: permutation p-values have a floor of 1/101 = 0.0099 (100 shuffles, runtime cap).

## BLOCKED
- Nothing blocked.

## How to run
```
.venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python -m pip install -e .
.venv/Scripts/python run_all.py        # ~3 min, writes results/ and figures/
.venv/Scripts/python -m pytest         # 36 tests
.venv/Scripts/python -m streamlit run app.py
```

## Key files
- `run_all.py` - full pipeline (data -> M1..M7 -> replication -> synthetic positive control)
- `src/oralbiome/` - one module per step; thresholds in `config.py`
- `app.py` - dashboard (reads results/, does not recompute statistics)
- `tests/` - pytest suite
- `DECISIONS.md` - every threshold and judgement call (D1-D44)
- `results/summary.json` - all headline numbers
