# LinkedIn post (draft)

My latest project ended with a negative result, and I think that is the most
useful thing it produced.

I analysed public 16S rRNA sequencing data from saliva (67 people) to ask
whether the oral bacterial community differs in periodontitis (gum disease),
and whether the community alone can tell who has it. I wrote the predictions
down before running anything.

In the main comparison (16 healthy people vs 18 with periodontitis), a
classifier separated the groups with an AUC of 0.92, well above the same
pipeline run on shuffled labels. On its own, that looks like a finding.

Then I checked it.

- Removing a single unidentified taxon dropped the AUC to 0.64. The model was
  mostly one feature.
- In a second group of people from the same study (people with high blood
  pressure, with and without periodontitis), the classifier performed at
  chance (AUC 0.52), and the overall community difference was no longer
  detectable.
- Quality control showed probable contamination in some samples, and the
  healthy group's apparent extra richness turned out to come from six unusual
  samples.

Some things did hold up. Three well-known gum-disease bacteria (P. gingivalis,
T. forsythia, T. denticola) were more abundant in periodontitis, as predicted
and as the original publication reports. I also compared my results with
that paper point by point, including where we differ and why.

Why post a result that did not replicate? Small microbiome studies produce
impressive numbers easily. The permuted-label control, the ablation and the
replication cohort are what showed how much this one could actually support.
Without them I would have reported a "microbial signature of periodontitis"
that probably isn't one.

Everything is associational, from a cross-sectional study, and nothing here
is a diagnostic tool. The code is tested and regenerates every number from
the raw data with one command.

Built in Python (pandas, SciPy, scikit-learn, NetworkX, Streamlit). Data
shared under CC BY 4.0 by the authors of: Guo Z., Yu X., Liu Y., Hu Q., Zhang Z., Zhang C., Li J. (2026). "A comparative analysis of oral microbial communities in hypertensive patients with and without chronic periodontitis." BMC Oral Health 26(1):836. doi:10.1186/s12903-026-08144-6

The code, figures and full write-up (including what did not work) are on my GitHub.

#bioinformatics #microbiome #reproducibility #python
