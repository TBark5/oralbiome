# Morning report

## Status: complete

All phases (0-9) are done and committed. Every module M1-M7 was built; none
was cut. Real data was used throughout.

## What works

| Item | Status | How it was checked |
|---|---|---|
| `python run_all.py` | Works; about 3 min end to end; slowest steps are M7 (85 s) and replication (82 s), each under the 2-minute cap | Ran from scratch after deleting all outputs, and again in a fresh `git clone` with a brand-new virtual environment; all outputs were byte-identical |
| `python -m pytest` | 36 of 36 tests pass (about 10 s) | Also passed in the fresh environment |
| `streamlit run app.py` | Launches cleanly; 9 tabs; filters for groups, flagged samples, taxonomic level and top-N | Streamlit AppTest across all filters with no exceptions, plus a live launch with headless-Chrome screenshots (`docs/screenshots/`) |
| README numbers | Match `/results` | An automated check of 31 key values against `results/*.json` found no mismatches |
| AI attribution | None in git history or tracked files | Searched `git log` and every tracked file for AI-tool names, co-author trailers and tool footers: no hits. All commits use the machine's existing git identity |
| Placeholders / TODOs | None | `git grep` for TODO, FIXME, TBD and placeholder links |

## What does not work, or is weaker than it looks

- **The findings do not replicate.** In the hypertensive cohort (T vs TP), PERMANOVA
  p = 0.402 and the L1 classifier AUC = 0.52 (permutation p = 0.495). This is a
  scientific result, not a bug, but it is the most important caveat.
- **The classifier rests on one taxon.** Dropping "Unclassified Bacilli" lowers
  the AUC from 0.92 to 0.64 (L1) and from 0.81 to 0.72 (random forest).
- **Likely contamination or batch effects in the source data.** 22,781 of 25,540
  features appear in only one sample; six healthy samples have about 5x the
  usual feature count. The public files have no negative controls, so this
  cannot be resolved.
- **Permutation p-values bottom out at 0.0099** because the null uses 100 shuffles
  (runtime cap). This is stated wherever p = 0.0099 appears.

## Data

**Real data**: Guo et al., *BMC Oral Health* 2026, figshare
10.6084/m9.figshare.29897750 (CC BY 4.0). Saliva, 16S V3-V4, 67 samples:
H 16, P 18, T 17, TP 16. It is cached in `data/raw/` and committed, so the project
runs offline. The synthetic fallback was **not** needed; the synthetic generator
is used only as a positive control.

## Final numbers (primary comparison: healthy n=16 vs periodontitis n=18)

| Analysis | Result |
|---|---|
| Alpha diversity (Shannon) | p = 0.796, r = +0.06 [-0.35, +0.45]; no metric passes BH (P1 not supported) |
| Beta diversity (Bray-Curtis) | PERMANOVA R² = 0.082, p = 0.015; PERMDISP p = 0.065 |
| Differential abundance | 2 of 354 genera pass FDR 5%; 72 with raw p < 0.05 vs 17.7 expected; null average 0.07 |
| Red-complex species (pre-declared) | *P. gingivalis* q = 0.005, *T. forsythia* q = 0.012, *T. denticola* q = 0.004, all higher in periodontitis |
| Network | 177 genera, 1,426 edges; red-complex genera co-occur (rho 0.68-0.85) and are hubs |
| Classifier | L1-LR AUC 0.92 [0.81, 1.00]; RF 0.81 [0.65, 0.94]; permutation p = 0.0099 for both |
| Replication (T n=17 vs TP n=16) | PERMANOVA p = 0.402; L1-LR AUC 0.52; red-complex q = 0.53-0.76 |
| Positive control (synthetic) | Sensitivity 67%, false discovery proportion 9% |

## Exact run commands (Windows, from the project folder)

```
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m pip install -e .
.venv\Scripts\python run_all.py
.venv\Scripts\python -m pytest
.venv\Scripts\python -m streamlit run app.py
```

On macOS/Linux, replace `.venv\Scripts\python` with `.venv/bin/python`.

## Modules cut

None.

## The 5 things you must understand before discussing this with anyone

1. **Compositionality.** Sequencing gives proportions, not amounts. That is why
   the pipeline uses CLR transforms and rank tests, and why "taxon X increased"
   always means "increased *relative to the rest*". See INTERVIEW_PREP Q2 and
   METHODS section 1.
2. **What PERMANOVA + PERMDISP say here.** Group explains about 8% of
   between-sample variation (p = 0.015). The spread difference is borderline
   (p = 0.065), so part of the signal could be that periodontitis samples are
   more variable. Be able to write the pseudo-F formula (METHODS section 4).
3. **Why the classifier result is weaker than AUC 0.92 suggests.** Leakage was
   prevented (all preprocessing in-fold, nested CV) and the permuted-label null
   is about 0.5, so the AUC is real *for this dataset*. But ablation shows one
   unidentified taxon carries most of it, and it does not transfer to the second
   cohort. Lead with this caveat; don't wait to be asked.
4. **Pre-declared vs exploratory.** The red-complex test was written into
   HYPOTHESIS.md before analysis, so correcting over 3 tests is legitimate; the
   354-genus scan is the exploratory counterpart (only 2 FDR hits). Know why
   *P. gingivalis* is significant while the genus *Porphyromonas* is not
   (*P. pasteri* dilutes it).
5. **Association, not causation, and not diagnosis.** Cross-sectional design, a
   small sample, saliva rather than plaque, possible contamination, and no
   covariates (smoking, age, hygiene). The honest headline is "some expected
   associations, a fragile classifier, and no replication". Being able to
   explain *why* is the strongest part of this project.

## Where to look

- `README.md`: overview, results tables, limitations
- `RESULTS_DISCUSSION.md`: each prediction evaluated
- `DECISIONS.md`: 44 logged decisions with reasons
- `METHODS.md`: formulas; `INTERVIEW_PREP.md`: 20 Q&As
- `PROGRESS.md`: handoff state
