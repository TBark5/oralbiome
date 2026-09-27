# Progress

## Done
- Phase 0: venv, pinned requirements, git, HYPOTHESIS.md, DECISIONS.md, plotting style (`src/oralbiome/style.py`).
- Phase 1: real dataset (source study's figshare deposit, CC BY 4.0) downloaded and cached in `data/raw/`; parser, synthetic fallback, QC report (`results/data_quality.json`, `figures/00_data_quality.png`); atypical-sample flag.
- Phase 2: M1 preprocessing, M2 composition + core microbiome.
- Phase 3: M3 alpha diversity (+2 sensitivity analyses), M4 beta diversity (Bray-Curtis/Jaccard, PCoA, NMDS, PERMANOVA, PERMDISP, heatmap).
- Phase 4: M5 differential abundance (genus scan + BH, permutation null, pre-declared red-complex species test, high-richness diagnostic).
- Phase 5: M6 co-occurrence network, M7 classifiers (nested CV, permutation null, ablation), replication cohort T vs TP.
- Phase 6: visuals pass (16 figures at 300 dpi, colorblind-safe palette, n on every group figure), figures/CAPTIONS.md.
- Phase 7: Streamlit dashboard `app.py` (9 tabs, sidebar filters); verified with AppTest and a live launch; screenshots in docs/screenshots/.
- Phase 8: pytest suite; a fresh clone + fresh venv reproduced all outputs byte for byte.
- Phase 9: README.md, RESULTS_DISCUSSION.md, METHODS.md, INTERVIEW_PREP.md, RESUME_BULLETS.md, LINKEDIN_POST.md, MORNING_REPORT.md.
- Correctness/publication pass (T1-T8):
  - canonical citations (PMID 41923033) and a test enforcing them;
  - Chao1 + ACE added to M3 (BH over 6 metrics) and a concordance section comparing with the source study;
  - taxonomy database documented as unstated; driver taxon profiled (`driver_taxon.py`), not identifiable without sequences;
  - LICENSE (MIT) + LICENSE-DATA (CC BY 4.0, verified 2026-09-27);
  - methods defences; resume/LinkedIn rewritten;
  - 48 tests pass; 53 document values audited against results/.
  - Decisions D45-D53.

## In progress
- T9: publish to GitHub as a public repository named `oralbiome` (account TBark5, gh authenticated over https).

## Next
- After pushing: verify the remote log (authors, committers, no attribution strings), add topics, and add the repository URL to MORNING_REPORT.md.
- Optional science: re-process SRA PRJNA1304526 with DADA2, re-annotate against eHOMD, identify the "Unclassified Bacilli" features.

## Known bugs
- None. Permutation p-values have a floor of 1/101 = 0.0099 (100 shuffles, runtime cap).

## BLOCKED
- Nothing blocked.

## How to run
```
.venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python -m pip install -e .
.venv/Scripts/python run_all.py        # ~3.5 min, writes results/ and figures/
.venv/Scripts/python -m pytest         # 48 tests
.venv/Scripts/python -m streamlit run app.py
```

## Key files
- `run_all.py` - full pipeline (data -> M1..M7 -> replication -> driver taxon -> synthetic positive control)
- `src/oralbiome/` - one module per step; thresholds in `config.py`
- `app.py` - dashboard (reads results/, does not recompute statistics)
- `tests/` - pytest suite (48 tests)
- `DECISIONS.md` - every threshold and judgement call (D1-D53)
- `DATA_SOURCE.md`, `LICENSE`, `LICENSE-DATA` - provenance and licenses
- `results/summary.json` - all headline numbers
