# Robustness tests (second November update)

Study-level numbers only; no sample IDs. Source tables are fetched from their own repositories
(see `atlas/README.md`); set `ATLAS_SRC`, `SWAPSEQ_DATA`, `YE_DATA` and `MLREPO` to the folders that hold them.

- `small_studies.mjs KEY`: random subsamples of a real table (6 to 40 people, 4 to 6 samples each), one hidden swap,
  40 repeats; strong flags before and now (`OLD_CORE` = an earlier `core.mjs`).
- `full_tables.mjs KEY`: the whole table, published labels, one swap per about 300 samples, and 3% of samples swapped.
- `dup16.mjs KEY`: planted technical duplicates and planted mislabels in count tables.
- `results.json`: what these gave on 19, 34 and 28 tables.
- `cli/test_inputs.mjs` (in the repository root): 60 odd and broken input files.

| One hidden swap, 12 checkable tables | 24 samples | 32 | 40 | 75 | 125 | 240 |
|---|---|---|---|---|---|---|
| Strong flags before | 1% | 15% | 30% | 65% | 73% | 83% |
| Strong flags now | 60% | 70% | 70% | 76% | 78% | 85% |
| Other samples flagged per 1,000, now | 0.4 | 0.7 | 0.0 | 0.3 | 0.8 | 1.2 |

On 34 full tables: strong flags on the published labels 75 -> 77, check list 201 -> 201; with 3% of samples
swapped, 35% -> 43% caught (shotgun tables of Ye et al.: 41% -> 60%), other samples flagged 1.0 -> 1.3 per 1,000.
Duplicates in 28 count tables: 309 of 420 planted duplicates found, no unplanted pair, 5 of 420 planted
mislabels read as a duplicate. Weak where people's samples are nearly alike (vaginal swabs, skin) or a
person's own samples are as close as a replicate (same-stool aliquots).
