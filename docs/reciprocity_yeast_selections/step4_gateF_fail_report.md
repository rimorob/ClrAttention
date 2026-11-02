# Yeast reciprocity preregistration: gate F FAILED at protocol step 4, so verdict 1

*2 Nov 2026 (job 0986). No truth has been downloaded or opened, so yeast is unconsumed. The preregistration was fixed in 02cbcfd; deviations 1 (743d24e) and 2 (25a5a8d) were both logged before any truth.*

## Result

| Quantity | Value |
|---|---|
| Mutant source | rebuilt with the deviation-2 rules (384 mutants; the file label "rebuilt_376" is the stage name, not the count) |
| Wild-type matrix | 169 cultures (balanced pairs passing the r > 0 pair check) |
| Matched N | 169 |
| Genes | 6,169 |
| **median R′ mutants** | **0.122** (range over 10 draws 0.093–0.161) |
| **median R′ wild types** | **0.296** (range 0.280–0.316) |
| Q, the share of mutant top-1% pairs also top-1% in the WT network | 0.205 (≤ 0.5, so this condition passes) |
| Genes per mutant outside the WT 0.5–99.5% range, median [q10, q90] | 386 [51, 2,053] |

**Gate rule:** pass if median R′_mut > median R′_wt **and** Q ≤ 0.5. The first condition fails in every draw. **By the fixed preregistration the result is verdict 1:** "Not run: no adequately powered held-out test of reciprocity was possible". When F fails, nothing else runs, H1 is not run, and no truth is consumed.

## Other step-4 outputs (label-free, committed)

- 10 random halves of the 652 published-M mutants.
- 286 mutant clusters: 147 singletons, the largest holding 20 mutants (3%), so no fallback was needed. Replicate N ranges from 277 to 389.

## Caveats a reader should know (none changes the verdict)

1. **Gate F ran on rebuilt mutants, not on the published M** (the implementation check did not reproduce the published M exactly, so the preregistered fallback applied). In sub-series `[hs1991]` the fitted orientation rule reproduced the published M poorly on held-out mutants (median r 0.42, against 0.93 in sub-series `none`). Misoriented hybridizations cancel signal when averaged, which lowers R′_mut. The failure may therefore partly reflect the GEO orientation problem rather than the deleteome network itself. That cannot be checked without changing the gate after seeing its result.
2. **The wild-type network is strongly reproducible** (R′ about 0.30 at N = 169, halves of about 85 cultures). This is consistent with strong shared day-to-day or growth programs in the wild-type controls.
3. Only 50.8% of balanced wild-type pairs passed the orientation pair check, so the wild-type orientation is itself uncertain.

## What the PI can do (any change is post hoc, because gate F's result is now known)

- **(a) Accept verdict 1 (recommended).** Report that the preregistered held-out test could not be run: the deleteome's deposited orientation metadata are inconsistent, and the sanity gate failed. Reciprocity stays a synthetic-only, explained candidate. The yeast truths stay unconsumed, but this file's design cannot be rerun on them without disclosing this failure.
- **(b) Continue anyway** (for example with F recomputed on the published M). This would be a deviation made **after** seeing a failed gate. Any result would have to be labelled exploratory, not confirmatory.

## Not done

Splits for t\*, the yeast and K562 selections, and the HepG2/Jurkat download (GSE264667) and human selections were not run. Disk: data/yeast 2.1 GB; 317 GB free.
