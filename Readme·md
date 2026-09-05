# Hippocampal Proteomics: Differential Abundance & Pathway Enrichment in Alzheimer's Disease

## Purpose

Identify differentially abundant proteins between Alzheimer's disease (AD) and
control hippocampal tissue, and the biological pathways they implicate, using
real published quantitative proteomics data.

This is an independent portfolio project, not connected to my published
research. My co-authored Scientific African paper analyzed hippocampal
**transcriptomics** (mRNA) on Alzheimer's disease; this project applies the
same general differential-analysis logic to real **proteomics** (protein-level)
data instead - a genuinely different data modality, on a different public
dataset, with its own distinct preprocessing challenge (detection-limit-driven
missing values) that mRNA-based methods don't have.

## Background

Alzheimer's disease pathology involves widespread changes in protein
expression in affected brain regions, not just transcriptional changes -
proteins are the functional molecules pathology actually acts through, and
mRNA and protein abundance often diverge. Studying the proteome directly, in
the same brain region already characterized transcriptomically, complements
that earlier analysis rather than repeating it.

## Data

Merrihew et al. (2023), *Scientific Data* - "A peptide-centric quantitative
proteomics dataset for the phenotypic assessment of Alzheimer's disease."
ProteomeXchange PXD034525, DOI: 10.6069/wefm-vv52, CC BY 4.0, hosted publicly
on Panorama Public (no login required).

- Region used: Hippocampus
- 44 individual, clinically characterized postmortem brain samples (Sporadic
  AD, Autosomal-Dominant AD, Control-High-Path, Control-Low-Path), collapsed
  here into a 2-group AD vs. Control comparison
- DIA mass spectrometry, protein-group-level abundance (log2 scale, already
  normalized/batch-adjusted by the original authors)
- 8 pooled QC/reference samples in the raw file are explicitly excluded before
  analysis (not real patient samples)

## Implementation

1. **Data acquisition** (`01_download_data.py`) - download the published
   protein abundance table and sample metadata directly from their public
   hosts (Panorama Public, GitHub).
2. **Preprocessing** (`02_preprocess.py`) - exclude QC/reference samples,
   collapse the 4-way clinical grouping into AD vs. Control, treat
   zero-as-missing (detection limit, not a true zero), filter out proteins
   with excessive missingness, and impute remaining missing values with a
   simple left-censored method appropriate for MNAR mass-spec data.
3. **Differential abundance** (`03_differential_expression.py`) - Welch's
   t-test per protein (AD vs. Control), log2 fold-change directly from the
   already-log2 data.
4. **Multiple-testing correction** - Benjamini-Hochberg FDR, applied in the
   same script.
5. **Visualization** (`04_visualize.py`) - volcano plot, heatmap of top
   differentially abundant proteins, and a per-sample missing-value pattern
   plot (a proteomics-specific diagnostic with no transcriptomics equivalent).
6. **Functional enrichment** (`05_enrichment_analysis.py`) - KEGG/GO
   enrichment on the significant protein list via Enrichr (`gseapy`) - the
   same tool already used in my Scientific African paper.
7. **Optional: PPI network & hub proteins** (`06_ppi_network.py`) - builds a
   protein-protein interaction network on the significant hits via the STRING
   API and ranks hub proteins by degree centrality, mirroring the hub-gene
   approach from my Scientific African paper, now applied to real
   protein-level data.

## Technical Stack

Python - pandas, numpy, scipy, statsmodels (BH-FDR), matplotlib, seaborn,
gseapy (Enrichr), networkx + STRING API (PPI network).

## Usage

```bash
pip install -r requirements.txt
python 01_download_data.py
python 02_preprocess.py
python 03_differential_expression.py
python 04_visualize.py
python 05_enrichment_analysis.py
python 06_ppi_network.py   # optional
```

## File Structure

```
proteomics-de-pipeline/
├── 01_download_data.py
├── 02_preprocess.py
├── 03_differential_expression.py
├── 04_visualize.py
├── 05_enrichment_analysis.py
├── 06_ppi_network.py
├── requirements.txt
├── data/
└── results/
```

## Example Output

Volcano plot, heatmap of top differentially abundant proteins, missingness
pattern plot, enrichment results table, and (optional) PPI network diagram -
generated into `results/` when the scripts are run.

## Limitations

- The t-test + per-protein imputation approach here is simpler than
  field-standard proteomics-specific tools (e.g., limma with moderated
  variance shrinkage, MSstats, or DEP's QRILC/MinProb imputation methods) - a
  defensible first-pass method for a portfolio project, not a claim of full
  statistical rigor.
- The 4 clinical subgroups are collapsed into a simple AD vs. Control
  comparison for clarity; the richer subgroup structure (autosomal-dominant
  vs. sporadic AD, histopathologic burden in cognitively normal controls) is
  present in the metadata but not modeled here.
- This is an independent portfolio project on public data, not connected to
  my published research.
