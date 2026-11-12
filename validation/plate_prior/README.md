# Plate-aware prior

Idea: a swap between two tubes that sit next to each other is more likely than a swap between two tubes anywhere on
the plate. When two samples are each other's best swap partner and sit in directly neighbouring wells of the same
plate (or next to each other in a plating order), their swap probability is recomputed with the expected error rate
multiplied by 5. Nothing is lowered for other pairs.

Test: 24 public tables (17 from the Atlas, the 7 longitudinal shotgun tables of Ye et al.). The samples of each table
were laid out on 96-well plates, in random order or sorted by person. Reciprocal swaps were then hidden between
direct neighbours, or between random samples, among people with 3 or more samples (10 repeats; 1,000 hidden swapped
samples per row). Strong flags: swaps caught / other samples flagged per 1,000.

| Layout | Hidden swaps | Without plate map | With plate map (x5) |
|---|---|---|---|
| Random | between neighbours | 68% / 4.00 | 76% / 4.03 |
| Random | between random samples | 67% / 3.59 | 67% / 3.59 |
| Random | none (published labels) | 5.06 | 5.06 |
| Sorted by person | between neighbours | 60% / 3.91 | 71% / 3.98 |
| Sorted by person | between random samples | 68% / 3.74 | 68% / 3.76 |
| Sorted by person | none (published labels) | 5.06 | 5.06 |

In the 21 tables SwapSeq can check, neighbour swaps went from 85% to 93% (random layout) and 75% to 86% (sorted).
In the three infant tables almost nothing is caught either way (3% to 6%, 0% to 10%).
Multipliers of 3 and 10 were also run (`results.json`): 74-79% and 66-74% for neighbour swaps; 10 began to add
other flags on published labels (5.06 to 5.11-5.20), so 5 was kept.

Limits:

- **That real swaps are mostly between neighbours is an assumption.** We have no data on where real swaps happen.
  The rule only raises neighbour pairs and lowers nothing, so a plate map can make the check stronger but not weaker.
- The plate layouts were made in the computer; the public tables come without real well positions.
- It needs a well (A01-H12) or plating-order column in which no two samples share a position (per plate).
  Otherwise it does nothing.
- The measured power shown by the tool is for a swap between any two samples and does not include this gain.

Files: `plate_test.mjs KEY` (keys as in `validation/robustness`; `ye:NAME` for the shotgun tables), `results.json`.
