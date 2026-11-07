# SwapSeq Atlas

**Live page: https://bymankirat.github.io/swapseq/atlas/**

How checkable are public longitudinal microbiome datasets for swapped and mislabelled samples?
SwapSeq was run on 30 public tables (16,336 samples, 21 studies). In each table, a few real
samples' labels were swapped on purpose (hidden test errors, 10 repeats), and we counted how many
SwapSeq caught.

- `data/atlas.csv`, `data/atlas.json` – one row per table: size, re-identification, hidden-test
  catch rates, false flags, strong flags on the labels as published (counts only), and the power
  figure shown by SwapSeq 1.0.1 and 1.1.1.
- `results/` – the per-table run summaries behind those rows.
- `code/` – loaders for each public dataset (`load2.mjs`), the Atlas run (`run.mjs`), the summary
  (`aggregate.py`), the power validation for SwapSeq 1.1 (`val11.mjs`), the extra-algorithm study
  (`ens.mjs`), and the page build (`build_site.py`). Paths in the scripts point to the build
  machine; change them to where you clone the source datasets.

Study-level results only. Sample IDs of flagged samples are not published here; they are shared
with each study's authors first. A low checkability grade means SwapSeq cannot vouch for a table
either way, not that its data are poor.

Data sources: the five SwapSeq validation studies, the GPZ longitudinal collection, the MITRE
example data, the parafac4microbiome data, MicrobeDS and the Borenstein lab curated collection.
See the Method section of the page for the full list of studies.
