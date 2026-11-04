# SwapSeq

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23963612.svg)](https://doi.org/10.5281/zenodo.23963612)

**A sample mix-up checker for longitudinal microbiome studies.** SwapSeq finds swapped tubes,
wrong labels and sample-sheet shifts from each person's microbial fingerprint. It runs entirely
in your browser: files you open are never sent anywhere.

**Use it:** open `index.html`, or the hosted page. Load a feature table (QIIME 2 `.qza`, BIOM 1 or 2,
or a text table) and the sample metadata, choose the column that says who each sample came from,
and press *Check samples*. Try the built-in demo first.

## What it reports

- **Strong flags** – samples that very likely belong to someone else, with the likely swap partner.
- **Check list** – weaker look-alikes to compare with lab records.
- **Evidence** – rare ASVs a flagged sample shares with the person it looks like, household
  look-alikes, plate positions of swap pairs, and possible sample-sheet shifts.
- **Power on your data** – how many swaps it would catch in your study, learnt from errors it
  plants in a copy of your table.
- A QC audit table and a ready-made sentence for your Methods section.

## Method (v3)

1. Fingerprints: soft Jaccard on expected presence at a matched sequencing depth, and Aitchison
   distance (CLR).
2. Single-sample score: nearest own-person sample against nearest other-person sample.
3. Swap score: the gain in fit if two samples' labels were exchanged.
4. Self-calibration: errors are planted in a copy of the user's own table; scores become
   probabilities at two levels (flag and check).
5. Sample-sheet shift scan along plate wells or plating order.

Code: `method/core.mjs` (the same code runs in the page).

## Validation

Tested on 10 public tables from 5 studies (4,732 samples), with errors planted at about a 1% rate,
30–60 repeats each. Full numbers are in `validation/` and on the page.

| Data | Swapped samples flagged | Same pipeline with Bray-Curtis | False flags per 1,000 |
| --- | --- | --- | --- |
| 16S, adults (Student Microbiome Project gut, forehead, palm, tongue; Immigration Microbiome Project) | 43–92% | 2–38% | ≤ 0.4 |
| 16S, before/after (2 samples per person) | 40–53% | 0–3% | 0 |
| Shotgun, daily (Minnesota diet study) | 88–94% | 84–92% | ≤ 0.4 |
| Infants and same-donor mice (negative controls) | 0–3%, reported as low power | 0% | ≤ 0.3 |

In the Minnesota diet study (Johnson et al. 2019), SwapSeq found a real tube swap the authors had
corrected (MCT.f.0086 and MCT.f.0417) without being told where to look.

Limits: one confirmed real error so far; adult 16S results come from two studies; plate patterns
were tested on one real plating order. A flag is a reason to check lab records, not proof.

### Studies used

- Flores GE et al. 2014, *Genome Biology* – Student Microbiome Project
- Vangay P et al. 2018, *Cell* – Immigration Microbiome Project
- Bokulich NA et al. 2016, *Science Translational Medicine* – ECAM
- Johnson AJ et al. 2019, *Cell Host & Microbe* – Minnesota diet study
- Gibson TE et al. 2025, *Nature Microbiology* – MDSINE2 mice

## Rerun the benchmark

Node.js 22, no other software. Clone the source repositories listed in `method/load.mjs` into one
folder, point `SWAPSEQ_DATA` at it, then run for example
`SWAPSEQ_DATA=~/swapseq-data SETS=gut,tongue node method/bench.mjs 30 0.01`.

## Privacy and use

SwapSeq has no server, sets no cookies and uses no analytics. Fonts are served from this site,
not from Google. It is a research quality-control tool, not a medical device, and must not be
used for diagnosis or clinical decisions.

## How to cite

Singh M. (2026). *SwapSeq: a sample mix-up checker for longitudinal microbiome studies* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23963612

## Licence

Copyright (C) 2026 Mankirat Singh. SwapSeq is free software under the GNU General Public License,
version 3 or later (`LICENSE`). You may use, study, share and change it; versions you share must
stay under the same licence and keep this copyright notice. Third-party parts keep their own
licences (`THIRD_PARTY_NOTICES.md`). If you use SwapSeq in published work, please cite it
(`CITATION.cff`).

Made by [Mankirat Singh](https://bymankirat.github.io/myportfolio/).
