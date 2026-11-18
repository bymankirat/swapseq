# Strain markers against species tables

Question: do strain-level marker tables catch swaps in infants, where species tables fail?

Answer: only a little. Strain markers are not a fix for infant studies.

| Study | Who | Samples | Swaps flagged, species | Swaps flagged, markers | With check list, species → markers |
|---|---|---|---|---|---|
| HMP 2019 (IBDMDB) | adults and children | 1,183 | 95% | 95% | 98% → 98% |
| Kostic 2015 | infants | 120 | 70% | 78% | 72% → 90% |
| Vatanen 2016 | infants | 682 | 0% | 0% | 5% → 50% |
| Shao 2019 | newborns | 961 | 0% | 18% | 45% → 25% |
| Bäckhed 2015 | newborns to 12 months | 300 | 0% | 0% | 0% → 0% |

Other flags stayed under 2 per 1,000 samples in every run.

## How it was run

- Data: `curatedMetagenomicData` 3.18.0 (Bioconductor 3.22), MetaPhlAn 3 `relative_abundance` (as counts) and `marker_presence`.
- People with at least 3 samples. The IBDMDB table was cut to 1,200 samples, whole people only. Yassour 2016 had too few repeat samples and was skipped.
- Markers present in at least 2 samples and at most half of them; a random 40,000 where more remained.
- `swapseqr` 0.3.2, `check_swaps(species, subject)` and `check_swaps(markers, subject, type = "markers")`. Figures are the package's own planted-swap power estimate.
- Run once, on 10 November 2026, on a laptop (R 4.5.1). `results.csv` holds the numbers.

## Limits

- One run. The marker subset is random, so figures can move by a few points.
- MetaPhlAn marker presence is a coarse strain signal. SNV-level profiles (for example StrainPhlAn) were not tested and may do better.
- No lab-animal data in this collection, so this says nothing about mice.

## Run it again

```r
source("strain_test.R")
```

`curatedMetagenomicData` 3.18.0 needs `rbiom` below version 3. If it fails to install with
"object 'unifrac' is not exported by 'namespace:rbiom'", install an older `rbiom` first:

```r
install.packages("rbiom", repos = "https://packagemanager.posit.co/cran/2025-06-01")
```
