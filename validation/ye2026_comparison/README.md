# SwapSeq vs Find-abnormality (Ye et al. 2026)

Both methods run on the 15 species-level metagenomic tables released with Ye et al. (2026)
(github.com/zhouyonog0713/Find-abnormality, MIT licence; stages 1-2 only). Study-level results only;
sample-level candidates are not published here and go to each study's authors first.

- `compare.py`: runs both methods on the published labels, then hides reciprocal swaps (about 1.5% of samples,
  10 repeats) among samples neither method flagged, and counts what each finds.
- `swapseq_run.mjs`: SwapSeq on one table (relative abundance x 1000 as pseudo-counts).
- `results.json`: per-table numbers.

- `fair.mjs`: the fair head-to-head. Swaps are hidden among ALL samples of people with 3+ samples (nothing excluded),
  10 repeats per table; for each method, the share of hidden swapped samples it flags and how many other samples it flags.
  Results in `fair_results.json` and `fair_log.txt`.
- `duptest2.mjs`: the duplicate check, with a sample re-sequenced under another person's name (abundances jittered,
  rare taxa dropped) and with planted mislabels. Log in `duplicates_log.txt` (`pair .7` is the shipped rule).
- `loadye.mjs`: reads the published tables. Set `YE_DATA` to the folder of profiles and metadata, and `YE_TOOL` to
  Find-abnormality's `find_abnormality.py`.

Fair test, 7 longitudinal tables, 340 hidden swapped samples (share flagged / other samples flagged per 1,000):

| Method | Swaps flagged | Other samples per 1,000 |
|---|---|---|
| Find-abnormality | 96% | 52 |
| SwapSeq 1.1.1, strong flags, 1% expected rate | 54% | 8.2 |
| SwapSeq, November update (reciprocal pairs promoted), 1% | 66% | 8.6 |
| SwapSeq, November update, 2% | 84% | 14.7 |
| SwapSeq, November update, 3% | 91% | 17.7 |
| SwapSeq 1.1.1 with the check list, 1% | 89% | 17.2 |

Duplicates (15 tables, 5 repeats of 3 planted duplicates each): 102 of 105 found in cross-sectional tables,
66 of 105 in longitudinal ones, no unplanted pair reported; 7 of 105 planted mislabels also appeared as a
near-identical pair. An earlier run of `compare.py` hid swaps only among samples neither method flagged; that
favoured both methods and is kept only for the record (`results.json`).

Second November update (small-table rule): with the same hidden swaps at the 1% expected rate, strong flags catch
79% instead of 65% in these seven tables (other samples flagged 9.2 instead of 8.4 per 1,000), and 60% instead of 41%
with 3% of samples swapped; see `validation/robustness`.

## Matched settings (current version)

`settings.mjs` repeats the fair test with the current pipeline at several expected error rates (results in
`settings_results.json`). Same seven tables, same 340 hidden swapped samples. Swaps flagged / other samples per 1,000:

| Method and setting | Strong flags | Strong flags + check list |
|---|---|---|
| Find-abnormality (one setting) | 96% / 52 | |
| SwapSeq, 1% (default) | 75% / 9.3 | 89% / 17.2 |
| SwapSeq, 3% | 91% / 17.7 | 97% / 24.0 |
| SwapSeq, 5% (the "known problems" setting) | 94% / 20.7 | 99% / 29.7 |
| SwapSeq, 10% | 98% / 26.9 | 100% / 39.5 |

So with a net as wide as Find-abnormality's or narrower, SwapSeq lists more of the hidden swaps and fewer other samples.

Read this with care:

- The SwapSeq setting was chosen after seeing the results; Find-abnormality has one setting. That favours SwapSeq.
- It is not a win in every table. At 5% with the check list, Find-abnormality still finds more in one table
  (95% vs 90%, 20 hidden samples) and flags fewer other samples in another (37 vs 61 per 1,000).
- At its default 1% setting SwapSeq finds fewer swaps than Find-abnormality (75% vs 96%), with fewer other samples.
- Seven tables, swaps made in the computer, "other samples" include real errors already in the published labels.
- A wider net costs work everywhere: on the published labels of 23 other public tables (13,265 samples) the 5%
  setting lists 31 samples per 1,000 (8 flags, 23 to check) against 7.5 per 1,000 at 1% (2 flags, 5 to check).
- Eight other fingerprints were tried on these tables, among them Find-abnormality's log Bray-Curtis; none was better
  across the seven tables, so the method is unchanged.
