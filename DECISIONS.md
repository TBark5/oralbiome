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
