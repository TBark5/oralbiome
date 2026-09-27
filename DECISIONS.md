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
