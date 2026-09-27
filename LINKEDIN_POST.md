# LinkedIn post (draft)

I just finished a personal research project: ORALBIOME, an analysis of the
oral microbiome in periodontitis (gum disease) using public 16S rRNA
sequencing data.

The question: does the bacterial community in saliva differ between people
with healthy gums and people with periodontitis, and can the community alone
predict who has the disease?

What I found, in 16 healthy vs 18 periodontitis samples:
- Overall community composition differed modestly between groups (group
  explained about 8% of the variation).
- Three classic periodontal pathogens (P. gingivalis, T. forsythia,
  T. denticola) were more abundant in periodontitis, as I had predicted
  before running the analysis.
- A classifier separated the groups with AUC 0.92, clearly better than the
  same pipeline on shuffled labels.

And what made me more careful:
- Most of the classifier's performance came from a single unidentified taxon.
  Without it, the AUC fell to 0.64.
- When I ran the same pipeline on a second group of participants (people with
  high blood pressure, with and without periodontitis), the community-level
  signal did not replicate.
- Quality control showed probable contamination in some samples, which I
  flagged and tested instead of quietly removing.

My main takeaway: with small microbiome datasets, controls matter as much as
the headline result. Permuted-label nulls, pre-registered predictions,
sensitivity analyses and replication changed how I would describe this result,
from "a microbial signature of periodontitis" to "some associations worth
testing in larger studies". These are associations only; nothing here is a
diagnostic tool.

Built in Python (pandas, SciPy, scikit-learn, NetworkX, Streamlit) with a
test suite and a one-command reproducible pipeline. Data: Guo et al., BMC
Oral Health 2026, shared under CC BY 4.0.

The code, figures and full write-up (including what did not work) are on my GitHub.

#bioinformatics #microbiome #datascience #python #research
