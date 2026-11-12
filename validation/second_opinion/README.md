# Second opinion: a rank rule on the listed samples

Question: does a second, simpler rule add anything to SwapSeq's lists?

The rule is the rank rule of Find-abnormality (Ye et al. 2026, MIT licence), written again on SwapSeq's own
soft-Jaccard distance. A person is in doubt when two of their samples are far from each other in both directions
(each ranks beyond 5% of all samples in the other's list of neighbours). Within such a person, a sample is named
when no close link (rank under 10 either way) joins it to the person's largest group of samples.

Test: 29 public tables (22 from the Atlas, the 7 longitudinal shotgun tables of Ye et al.), reciprocal swaps hidden
among people with 3 or more samples, 10 repeats, 1,560 hidden swapped samples. Swaps named / other samples per 1,000:

| | Hidden swaps | Other samples per 1,000 |
|---|---|---|
| Rank rule alone | 93% | 85 (287 in the three infant tables) |
| SwapSeq strong flags | 70% | 3.3 |
| ... of which the rank rule agrees | 67% | 2.0 |
| SwapSeq check list (not flagged) | 9% | 7.3 |
| ... of which the rank rule agrees | 8% | 4.0 |

What was decided:

- **Shown as a note, not used to flag.** On listed samples the tool says whether the second check agrees. Where it
  agrees, almost all hidden swaps are kept (95% of flagged, 94% of check-list ones) and 41-46% of the other samples
  drop out, so "agrees" marks the better candidates.
- **Not used to turn check-list samples into flags.** That would have raised strong flags from 70% to 78% of hidden
  swaps but more than doubled the other samples flagged (3.3 to 7.3 per 1,000), most of them in infant, skin and
  tongue tables.
- **Not offered as a list of its own.** Alone the rule names 85 samples per 1,000, and more than one in four in
  infant cohorts.

The same test with Find-abnormality's own distance (Bray-Curtis on log abundance) gave the same picture
(`RD=logbc`), so no second distance is computed.

Files: `rank_rule_test.mjs KEY` (keys as in `validation/robustness`; `ye:NAME` for the shotgun tables),
`results.json`. "Other samples" include real errors already in the published labels.
