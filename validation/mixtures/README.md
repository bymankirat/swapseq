# Mixtures of two people: tried and dropped

Question: can SwapSeq tell when one sample holds material from two people (cross-contamination, or two tubes pooled)?

Answer: not reliably from a taxon table. No mixture check was added to the tool.

Test: in 23 public count tables, a sample was replaced by a mixture of itself and a sample from another person
(50%, 30% or 10% of the reads from the other person; counts drawn again). Study-level numbers only.

| Share of reads from the other person | 50% | 30% | 10% |
|---|---|---|---|
| The tool as it is: mixture put on the flag or check list | 37% | 18% | 3% |
| The tool as it is: strong flag | 20% | 11% | 2% |
| New statistic (personal taxa of another person): median table | 26% | 23% | 8% |
| New statistic: tables that caught half or more | 8 of 23 | 4 of 23 | 3 of 23 |

- When the tool as it is lists a 50% or 30% mixture, it names the right other person (99-100%).
- New statistic: for each sample, the share of another person's personal taxa (kept by that person, rare in others,
  never seen in the labelled person's other samples) that turn up in it; a mixture should show one person far above the
  rest. Counted as caught when above the 4th highest unmixed sample of the table. Thirteen other variants (read-weighted,
  prevalence-weighted, gap to the second person) did no better across tables.
- It worked where people differ a lot (immigrant cohort 95%, cafe study 100%) and failed in most other tables: an ordinary
  sample also picks up taxa that some other person carries, so a mixture and a normal change look alike.

Files: `personal_taxa.mjs KEY` (new statistic), `current_tool.mjs KEY` (the tool as it is), `results.json`.
Data and environment variables as in `validation/robustness`.

Limits: mixtures were made in the computer, not in the lab; 8 to 40 planted mixtures per table and level.
