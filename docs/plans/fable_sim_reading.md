# Reading of the planted-module simulation (Fable 5.1, read-only; based on `docs/plans/sim_planted_results.md`)

## 1. What survives of the memo's Q3/Q5 account

**Survives (derivable parts):** the t = 1 identity (R6, 4.3e-15); aCLR's 1/n dilution (slopes -0.77 to -0.95 against -1) and CLR's n^-1/2 penalty (-0.30 to -0.46 against -0.5); the leakage profile of in-module mass at n = 5 and 10; the claim that a row-adaptive exponent matters (fixed-exponent inflation reaches 0.29-0.34 against 0.63-0.70); the link between background spread and gain (R10: sigma_bg halves at keff 40 and the gain over aCLR turns negative). tcaCLR beats aCLR at sizes 20 and 40 in 10 of 10 seeds at every N.

**Refuted:** (a) tcaCLR size-invariance (|slope| < 0.3; observed 1.2-2.0 z per nat). With the measured sigma_bg of 1.4-1.8 nats my own formula gives about 0.6 z per nat from the log n term alone, so the "small slope" claim followed from the assumed 3-6 nats, which was wrong. (b) tcaCLR >= CLR for large n: it loses at 70 and 100. (c) tcaCLR close to aCLR at n = 5: it is far below both. (d) aCLR as anti-DPI on direct edges: aCLR is the best of the three on T4. (e) The satellite inversion and unstable lifts: neither appeared (lifts were the best readouts).

**The memo's R1 failure clause is met** for n >= 70. The size-invariance account of Q3 is therefore wrong as stated. What remains is narrower: a row-adaptive log standardisation helps in moderate-size modules (20-40) and hurts elsewhere, for reasons the memo did not identify.

**Not testable here:** non-linear dependence (the generator is linear-Gaussian, hence |r| wins), ribosome-specific structure, experiment-level variance, the K562 selection, real reachability structure.

## 2. Strong-large-module collapse: hypothesis and diagnostic

**Hypothesis (reachability confinement).** In a strong module (lambda 0.85, n >= 70) the 50 candidates per row lie inside the module and retention at t = 9 is above 0.9. The reachable set of a member row is then essentially the module plus a thin leak. The "background" over which the row z-score is taken is dominated by the module itself, so members cannot stand out from it. This is CLR's sparsity violation, reappearing in the log domain through the reachable-only rule. At lambda 0.55 the module leaks, the reachable background is wide, and scores stay at 5-6.

**Diagnostic (minutes, label-free on one N = 400 compendium).** For member rows of modules of size 70 and 100, by lambda, report: (i) the number of reachable off-diagonal entries; (ii) the share of reachable entries that are in-module (module membership is simulation structure, not a method truth); (iii) the row SD of log P^9 and the member mean z. Then recompute the score with row statistics over all entries with a floor of log(1e-300) (the non-robust a3 form).

- If (i) is about n plus a few and (ii) exceeds 0.5 at lambda 0.85 but not at 0.55: hypothesis confirmed. tcaCLR's reachable-pairs rule fails when modules are closed. The all-entries variant should recover scores. In real data check the same two numbers for ribosomal rows (K562 retention 0.88, 0.28% unreachable overall).
- If (ii) is low and the reachable set is wide: reject. Next check the one-sided Z before the Stouffer step, then the clipping fraction. A low one-sided Z with member log P near the row mean would point to the top-50 truncation, which makes within-module log P nearly constant.

## 3. H2 wording

The simulation supports a lower edge of 20 (10 of 10 seeds at 20 and 40, 0 of 10 at 5 and 10). It does not support "positive for all sizes of 20 or more": sizes 70 and 100 give 3-8 of 10 seeds at N = 400 and 800, with only 2-3 modules per seed. It also contradicts "near zero" for 7-19 and "not positive" understated for 3-6: the simulation gives -0.26 to -0.53 at n = 5 and 10, against net4's -0.003 and -0.004. Reword H2 to: "tcaCLR minus aCLR is expected to be positive in complexes of 20 or more genes pooled, and not positive in complexes of fewer than 20 genes, where it may be markedly negative. No monotone dependence on size above 20 and no advantage over CLR in the largest complexes is predicted. Bins are those of net4; the simulation (linear-Gaussian, sizes 20 and 40) gives the same lower edge."

## 4. Replacement wording for memo section 5.2

- Claim 4, replace with: "tcaCLR's row standardisation is invariant to row-multiplicative and row-power distortions of the walk and equals CLR on candidate sets at t = 1. It converts the 1/n dilution into an additive log n offset divided by the row spread of log P^t; in simulation that spread was 1.4-1.8 nats and the within-module score still fell steeply with size."
- Claim 5, replace with: "aCLR's within-module score falls as about 1/n and CLR's as about n^-0.5 (simulation: slopes -0.8 to -0.95 and -0.3 to -0.46). tcaCLR gains over aCLR in modules of 20-40 genes and loses in modules of 5-10; for 70-100 genes the result is seed-dependent and it can collapse for strong modules."
- Drop: "bounded below uniformly in n", "diffusion CLR removes size dependence", any statement that sigma_bg is several nats, and the section 5.3 size-invariance proposition.
- Add: "Column hub effects, reachability confinement, and linear-Gaussian generator limits were not removed or tested."
