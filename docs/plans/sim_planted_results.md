# Planted-module simulation: results

Specification: `docs/plans/fable_theory_memo.md` section 5.4. Implementation choices and the scoring rules for R1 to R12 were fixed and committed before any run (`docs/plans/sim_planted_choices.md`); the predictions below are quoted unchanged from the memo. Synthetic data only. No method parameter was tuned (keff 10, cap 50, alpha 0.5, t = 9; the keff 40 / cap 200 arm and the depth grid are descriptive). Code: `analysis/sim_planted/`; per-compendium results (not committed): `results/sim_planted/units/`.

Compendia analysed: 120 of 120. Values are mean (SD across seeds) [min, max] over the 10 seeds of the cell; the primary arm is Gaussian noise with identity marginals, shown at N = 150, 400, 800. Compendia per cell:

| N | noise | marginal | compendia |
|---|---|---|---|
| 150 | gaussian | copula | 10 |
| 400 | gaussian | copula | 10 |
| 800 | gaussian | copula | 10 |
| 150 | t4 | copula | 10 |
| 400 | t4 | copula | 10 |
| 800 | t4 | copula | 10 |
| 150 | gaussian | identity | 10 |
| 400 | gaussian | identity | 10 |
| 800 | gaussian | identity | 10 |
| 150 | t4 | identity | 10 |
| 400 | t4 | identity | 10 |
| 800 | t4 | identity | 10 |

## Verdicts

A prediction is judged on clause cells (a clause evaluated at one N); held = every cell true, did not hold = none true, partly held = otherwise (rules in the choices file, section 5).

| Prediction | Verdict | Cells_true | Statement |
|---|---|---|---|
| R1 | partly held | 23 of 39 | (T2, n >= 20): tcaCLR > CLR > aCLR, with the aCLR deficit and the tcaCLR-CLR gap both growing with n; slopes ~ -1 (aCLR), ~ -0.5 (CLR), \|slope\| < 0.3 z/nat (tcaCLR). |
| R2 | partly held | 11 of 24 | (T2, n = 5): CLR >= tcaCLR ~ aCLR; all three below their n = 20 level (forced leakage); n = 10 intermediate. |
| R3 | partly held | 7 of 12 | (T3, loose regulons): aCLR ~ tcaCLR > CLR, gap largest for b = 0.35 and m = 40. |
| R4 | partly held | 3 of 9 | (T4, direct edges): CLR >= tcaCLR >= aCLR; \|r\| competitive. |
| R5 | partly held | 7 of 15 | (heavy-tail arm): CLR's size penalty becomes noisier (sigma_i inflated), tcaCLR's advantage on n >= 20 grows; (copula arm): CLR and tcaCLR unchanged in rank, \|r\| loses. |
| R6 | held | 1 of 1 | (identity): tcaCLR at t = 1 equals Stouffer of within-candidate row z of CLR to 1e-10. |
| R7 | partly held | 9 of 15 | (lemma): Sinkhorn-balanced walk, raw readout, ties aCLR on n >= 20 (within seed SD) and reduces the hub share of the top 1%. |
| R8 | partly held | 14 of 33 | (depth): at t = 36 and 72 all methods' T1 falls toward the hub-score level and the hub share of the top 1% rises; aCLR's rises faster than tcaCLR's. |
| R9 | partly held | 24 of 36 | (satellites): aCLR ranks satellite-member pairs above member-member pairs of the same module for n >= 20; CLR and tcaCLR rank member-member higher. |
| R10 | partly held | 5 of 6 | (wide-row arm, keff 40): sigma_bg falls, tcaCLR's T1 gain over aCLR shrinks or vanishes; CLR unaffected. |
| R11 | partly held | 3 of 18 | (lift): pi log-lift unstable (large seed SD, hub-dominated top 1%); column-mean log-lift between aCLR and tcaCLR on T1. |
| R12 | partly held | 8 of 9 | (inflation control): fixed-exponent inflation does not reach tcaCLR on T1 for both r; Sinkhorn(A) -> walk -> tcaCLR is at or above tcaCLR on hub-stratified T1 (speculation; a failure here is informative). |

## H2 bin-edge rule for `docs/plans/tcaclr_confirm_plan.md`

Rule: in the primary arm (Gaussian, identity), the smallest module size at which tcaCLR minus aCLR (T2 per-size-bin pooled AUPR) is positive in at least 9 of the 10 seeds, evaluated at each N. Strictly positive difference per seed.

| N | module size | seeds with tcaCLR > aCLR | mean difference in AUPR | meets 9 of 10 |
|---|---|---|---|---|
| 150 | 5 | 0 of 10 | -0.37 | no |
| 150 | 10 | 0 of 10 | -0.257 | no |
| 150 | 20 | 10 of 10 | 0.146 | yes |
| 150 | 40 | 10 of 10 | 0.366 | yes |
| 150 | 70 | 10 of 10 | 0.277 | yes |
| 150 | 100 | 9 of 10 | 0.175 | yes |
| 400 | 5 | 0 of 10 | -0.531 | no |
| 400 | 10 | 0 of 10 | -0.476 | no |
| 400 | 20 | 10 of 10 | 0.178 | yes |
| 400 | 40 | 10 of 10 | 0.504 | yes |
| 400 | 70 | 8 of 10 | 0.201 | no |
| 400 | 100 | 3 of 10 | 0.0779 | no |
| 800 | 5 | 0 of 10 | -0.443 | no |
| 800 | 10 | 0 of 10 | -0.431 | no |
| 800 | 20 | 10 of 10 | 0.276 | yes |
| 800 | 40 | 10 of 10 | 0.568 | yes |
| 800 | 70 | 7 of 10 | 0.155 | no |
| 800 | 100 | 8 of 10 | 0.223 | no |

Edge per N (smallest qualifying size): N = 150: 20; N = 400: 20; N = 800: 20. Every larger size also qualifies: N = 150: yes; N = 400: no; N = 800: no.

**Overall (the smallest size that qualifies at every N): 20**; the most favourable N gives 20.

## Plain statements on what held and what did not

The verdicts above come from the pre-written rules; this section only states in words what the tables show. All numbers are seed means in the primary arm (Gaussian noise, identity marginals) at N = 150, 400, 800 unless stated.

**What held.**
- R6 (identity at t = 1): the largest difference over 120 compendia is 4.3e-15.
- R1, partly: at n = 40, tcaCLR > CLR > aCLR at every N (T2 AUPR 0.70, 0.82, 0.84 against CLR 0.50, 0.50, 0.42 and aCLR 0.33, 0.32, 0.27). CLR > aCLR holds at every size from 20 to 100. The size slopes of aCLR (-0.77, -0.93, -0.95 on log score per log n) and CLR (-0.30, -0.41, -0.46) match the predicted -1 and -0.5.
- R10: median sigma_bg falls from 1.41, 1.68, 1.80 nats at keff 10 to 0.68, 0.80, 0.90 at keff 40, and tcaCLR's T1 advantage over aCLR (+0.14, +0.05, +0.13) turns negative (-0.05, -0.06, -0.04). The one failed cell (N = 400) is the 9-of-10 sign rule applied to a noisy keff 10 gap (SD 0.18); the means agree with the prediction.
- R12: fixed-exponent inflation reaches only 0.29 to 0.34 on T1 against tcaCLR's 0.63 to 0.70 for both r at every N. Sinkhorn then tcaCLR is at or above tcaCLR on the hub-stratified T1 (true hub loading sum) at N = 400 and 800 but not at 150 (0.324 against 0.345).
- R7, the hub-share part: the Sinkhorn-balanced walk lowers the hub share of the top 1% (0.043, 0.047, 0.042 against aCLR 0.049, 0.058, 0.058) at every N.
- R9, the CLR and tcaCLR halves: both rank member-member above satellite-member pairs.

**What did not hold.**
- R1, tcaCLR slope and the n >= 20 advantage: tcaCLR's within-module score falls by 1.2, 2.0, 1.8 z per unit of log n for n >= 20 (prediction: |slope| < 0.3). In z units that is flatter than CLR's (-2.8, -4.5, -5.2), but in log-score units tcaCLR (-0.34, -0.69, -0.62) is steeper than CLR (-0.30, -0.41, -0.46). tcaCLR beats CLR at n = 20 only at N = 800 and not at n = 70 or 100 at any N; the tcaCLR-CLR gap does not grow with n. The size-invariance account of memo section 5.3 is therefore not supported for n >= 70 in this simulation; the memo's own failure statement for R1 ("tcaCLR's slope is as steep as CLR's or its n >= 20 advantage is absent") is met at n = 70 and 100, not at n = 20 to 40.
- The failure is concentrated in strong large modules: at N = 400, tcaCLR's per-module mean score is 5.2 to 5.9 for lambda = 0.55 at every size, but 1.9 (n = 70) and 1.5 (n = 100) for lambda = 0.85, and module-level AP is 0.62 for lambda = 0.70 and 0.85 at n = 100 and 0.84 and 0.56 at n = 70 (CLR 0.97 to 0.99 and aCLR 0.96 at n = 100 for lambda >= 0.70). The mechanism was not diagnosed here; only the pattern is reported. Seed variation of tcaCLR's T2 AUPR at n = 70 and 100 is large (range 0.003 to 0.80 at N = 400) because each seed contains only three modules of size 70 and two of size 100.
- sigma_bg is 1.4 to 1.8 nats at t = 9 (CV across rows 0.23 to 0.25), not the "several nats" (3 to 6) assumed in memo section 3.2; this is the step the memo marked as plausible rather than derived.
- R2: at n = 5 tcaCLR is far below both others (T2 AUPR 0.13, 0.08, 0.07 against CLR 0.59, 0.78, 0.87 and aCLR 0.50, 0.61, 0.51); CLR >= tcaCLR holds, tcaCLR ~ aCLR does not. CLR and aCLR do not fall from n = 20 to n = 5 in T2 AUPR (CLR's n = 5 level exceeds its n = 20 level), so forced leakage shows in AUPR for tcaCLR only (bins 5 < 10 < 20 at every N). The in-module mass at t = 9 follows the predicted leakage pattern (0.26, 0.33, 0.36 at n = 5; 0.51, 0.66, 0.72 at n = 10; above 0.9 at n >= 70). Raw per-bin AUPR across bins is confounded with the rising base rate.
- R3: aCLR > tcaCLR > CLR on pooled regulon pairs (T3 AUPR 0.034, 0.11, 0.17 for aCLR; 0.029, 0.092, 0.14 for tcaCLR; 0.017, 0.062, 0.11 for CLR). Both beat CLR, but aCLR and tcaCLR are not equivalent. The largest gain over CLR is at m = 40, b = 0.35 only at N = 800; at N = 150 and 400 it is at m = 40, b = 0.50.
- R4: tcaCLR is not at or above aCLR on direct regulator-target edges; aCLR is the best of the three at every N (T4 AUPR 0.006, 0.027, 0.060, about double CLR at N = 400 and 800). |r| is not competitive (0.004, 0.007, 0.008). In this simulation aCLR is not anti-DPI on direct edges.
- R5: the heavy-tail clauses hold at N = 150 and 400 (seed SD of CLR's slope larger; tcaCLR - CLR gap on n >= 20 larger) but the differences are small, and the gap clause fails at N = 800 (0.147 against 0.148). In the copula arm |r| does not lose: its T1 AUPR (0.72, 0.88, 0.92) is above CLR's (0.49, 0.69, 0.73) and tcaCLR's (0.62, 0.77, 0.73) at every N. The sign of tcaCLR - CLR on T1 changes between the identity and copula arms at N = 400.
- R7, the tie clause: the Sinkhorn-balanced raw readout does not consistently tie aCLR (6 of 12 size and N cells within the stated band; at n = 20 it is 0.33 against 0.49 at N = 150).
- R8: T1 AUPR does not fall at t = 36 or 72 for either method (aCLR rises from 0.57 at t = 9 to 0.80 at t = 36 at N = 150; tcaCLR is flat at 0.63 to 0.75) and does not approach the hub-score level (0.02 to 0.025). The hub share of the top 1% does rise (aCLR 0.05 to 0.06 at t = 9 and 0.12 to 0.16 at t = 72; tcaCLR 0.05 at t = 9 and 0.09 to 0.11 at t = 72), faster for aCLR at N = 400 and 800. Within the grid the walk has not reached the regime t >> t_mix.
- R9, aCLR: aCLR ranks member-member above satellite-member pairs (probability that a satellite-member pair outranks a member-member pair 0.16 to 0.33, below 0.5 at every size and N). aCLR's value is higher than CLR's at n >= 40 and tcaCLR's is lowest, so the ordering across methods is in the predicted direction but the INTS8-type inversion is not reproduced.
- R11: neither lift is unstable or hub-dominated in this simulation. The pi log-lift reaches T1 AUPR 0.76, 0.91, 0.92 and the column-mean lift 0.79, 0.92, 0.92, above CLR, aCLR and tcaCLR at every N; hub shares of the top 1% are about 0.05. Only the prediction that the lifts' seed variability exceeds aCLR's holds.

**Unexpected.**
- |r| has the highest T1 AUPR of all methods at every N in the identity arm (0.82, 0.92, 0.93); the generator is linear-Gaussian, so this is the regime where correlation is the natural estimator, and the simulation says nothing about regimes with non-linear dependence.
- tcaCLR's seed variability is much larger than that of the other methods (T1 SD 0.11 to 0.18 against 0.006 to 0.027 for aCLR and CLR at N = 400 and 800), driven by the large strong modules.
- Sinkhorn balancing of the lazy operator did not reach the 1e-9 tolerance within the 2,000-iteration cap in most compendia (final max column-sum error about 3e-4); the Sinkhorn arm used the last iterate.

**H2 bin edge.** The rule gives size 20 at each of N = 150, 400, 800 (10 of 10 seeds at n = 20 and 40 at every N). The rule counts only the single smallest qualifying size: at N = 400 and 800 the sizes 70 and 100 do not qualify (positive in 3 to 8 seeds out of 10, mean difference still positive), so the positive region is not contiguous above 20 at those N. At n = 5 and 10, tcaCLR exceeds aCLR in 0 of 10 seeds at every N.

## R1

Prediction (memo): (T2, n >= 20): tcaCLR > CLR > aCLR, with the aCLR deficit and the tcaCLR-CLR gap both growing with n; slopes ~ -1 (aCLR), ~ -0.5 (CLR), |slope| < 0.3 z/nat (tcaCLR).

Verdict: **partly held** (23 of 39 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T2 AUPR tcaclr n=20 | 0.631 (0.0241) [0.576, 0.656] | 0.72 (0.0385) [0.648, 0.758] | 0.743 (0.0367) [0.691, 0.789] |
| T2 AUPR clr n=20 | 0.582 (0.0399) [0.524, 0.656] | 0.677 (0.0297) [0.639, 0.706] | 0.605 (0.0246) [0.569, 0.63] |
| T2 AUPR aclr n=20 | 0.485 (0.0252) [0.441, 0.514] | 0.542 (0.0183) [0.52, 0.569] | 0.468 (0.0239) [0.436, 0.498] |
| T2 AUPR tcaclr n=40 | 0.695 (0.0948) [0.533, 0.834] | 0.819 (0.0387) [0.766, 0.877] | 0.842 (0.0514) [0.746, 0.919] |
| T2 AUPR clr n=40 | 0.504 (0.0789) [0.393, 0.622] | 0.502 (0.0422) [0.459, 0.565] | 0.421 (0.0168) [0.399, 0.447] |
| T2 AUPR aclr n=40 | 0.329 (0.0258) [0.292, 0.364] | 0.315 (0.00816) [0.299, 0.325] | 0.274 (0.00918) [0.261, 0.289] |
| T2 AUPR tcaclr n=70 | 0.532 (0.212) [0.256, 0.76] | 0.446 (0.227) [0.0988, 0.824] | 0.371 (0.189) [0.152, 0.804] |
| T2 AUPR clr n=70 | 0.419 (0.0268) [0.38, 0.453] | 0.438 (0.00933) [0.421, 0.45] | 0.38 (0.00841) [0.368, 0.39] |
| T2 AUPR aclr n=70 | 0.255 (0.0145) [0.224, 0.277] | 0.245 (0.00418) [0.238, 0.252] | 0.216 (0.00482) [0.208, 0.222] |
| T2 AUPR tcaclr n=100 | 0.381 (0.142) [0.109, 0.59] | 0.28 (0.277) [0.00275, 0.804] | 0.398 (0.248) [0.00275, 0.73] |
| T2 AUPR clr n=100 | 0.327 (0.126) [0.167, 0.522] | 0.38 (0.0468) [0.332, 0.448] | 0.356 (0.018) [0.335, 0.386] |
| T2 AUPR aclr n=100 | 0.206 (0.0236) [0.16, 0.227] | 0.202 (0.00685) [0.193, 0.211] | 0.175 (0.00642) [0.167, 0.186] |
| slope aCLR (log score on log n, n>=20) | -0.772 (0.0628) [-0.881, -0.673] | -0.925 (0.0193) [-0.956, -0.9] | -0.948 (0.0125) [-0.974, -0.932] |
| slope CLR (log score on log n, n>=20) | -0.301 (0.0436) [-0.375, -0.217] | -0.413 (0.0197) [-0.437, -0.385] | -0.456 (0.0111) [-0.476, -0.445] |
| slope tcaCLR (z per nat, n>=20) | -1.21 (0.722) [-2.16, -0.279] | -1.96 (0.887) [-3.4, -0.394] | -1.84 (0.832) [-3.25, -0.56] |
| memo failure test, z units: \|slope tcaCLR\| - \|slope CLR\| (negative = tcaCLR flatter) | -1.58 (0.637) [-2.72, -0.795] | -2.53 (0.867) [-4.5, -1.46] | -3.37 (0.821) [-4.64, -2.16] |
| memo failure test, log-score units: \|slope tcaCLR\| - \|slope CLR\| (negative = tcaCLR flatter) | 0.0416 (0.271) [-0.279, 0.42] | 0.275 (0.371) [-0.377, 0.88] | 0.168 (0.366) [-0.366, 0.814] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T2 n=20: tcaCLR > CLR | false | false | true |
| T2 n=20: CLR > aCLR | true | true | true |
| T2 n=40: tcaCLR > CLR | true | true | true |
| T2 n=40: CLR > aCLR | true | true | true |
| T2 n=70: tcaCLR > CLR | false | false | false |
| T2 n=70: CLR > aCLR | true | true | true |
| T2 n=100: tcaCLR > CLR | false | false | false |
| T2 n=100: CLR > aCLR | true | true | true |
| aCLR deficit (CLR - aCLR) grows with n | false | false | true |
| tcaCLR - CLR gap grows with n | false | false | false |
| slope aCLR in [-1.25, -0.75] | true | true | true |
| slope CLR in [-0.75, -0.25] | true | true | true |
| \|slope\| tcaCLR < 0.3 z per nat | false | false | false |

Failed clauses: T2 n=20: tcaCLR > CLR (N = 150); T2 n=70: tcaCLR > CLR (N = 150); T2 n=100: tcaCLR > CLR (N = 150); aCLR deficit (CLR - aCLR) grows with n (N = 150); tcaCLR - CLR gap grows with n (N = 150); |slope| tcaCLR < 0.3 z per nat (N = 150); T2 n=20: tcaCLR > CLR (N = 400); T2 n=70: tcaCLR > CLR (N = 400); T2 n=100: tcaCLR > CLR (N = 400); aCLR deficit (CLR - aCLR) grows with n (N = 400); tcaCLR - CLR gap grows with n (N = 400); |slope| tcaCLR < 0.3 z per nat (N = 400); T2 n=70: tcaCLR > CLR (N = 800); T2 n=100: tcaCLR > CLR (N = 800); tcaCLR - CLR gap grows with n (N = 800); |slope| tcaCLR < 0.3 z per nat (N = 800).

## R2

Prediction (memo): (T2, n = 5): CLR >= tcaCLR ~ aCLR; all three below their n = 20 level (forced leakage); n = 10 intermediate.

Verdict: **partly held** (11 of 24 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T2 AUPR tcaclr n=5 | 0.126 (0.0259) [0.092, 0.168] | 0.0802 (0.0136) [0.0649, 0.11] | 0.0662 (0.0108) [0.0558, 0.0898] |
| T2 AUPR tcaclr n=10 | 0.339 (0.0583) [0.269, 0.438] | 0.266 (0.0505) [0.189, 0.355] | 0.267 (0.056) [0.207, 0.389] |
| T2 AUPR clr n=5 | 0.586 (0.0152) [0.561, 0.611] | 0.775 (0.0225) [0.732, 0.808] | 0.869 (0.015) [0.838, 0.887] |
| T2 AUPR clr n=10 | 0.609 (0.0251) [0.583, 0.647] | 0.769 (0.0226) [0.73, 0.808] | 0.818 (0.0152) [0.791, 0.838] |
| T2 AUPR aclr n=5 | 0.496 (0.0309) [0.454, 0.553] | 0.612 (0.031) [0.574, 0.648] | 0.509 (0.0344) [0.451, 0.554] |
| T2 AUPR aclr n=10 | 0.597 (0.0118) [0.584, 0.615] | 0.742 (0.0255) [0.693, 0.77] | 0.698 (0.0392) [0.633, 0.765] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| tcaclr: bin 5 below bin 20 | true | true | true |
| tcaclr: bin 10 intermediate | true | true | true |
| clr: bin 5 below bin 20 | false | false | false |
| clr: bin 10 intermediate | false | true | true |
| aclr: bin 5 below bin 20 | false | false | false |
| aclr: bin 10 intermediate | false | false | false |
| n=5: CLR >= tcaCLR | true | true | true |
| n=5: tcaCLR ~ aCLR | false | false | false |

Failed clauses: clr: bin 5 below bin 20 (N = 150); clr: bin 10 intermediate (N = 150); aclr: bin 5 below bin 20 (N = 150); aclr: bin 10 intermediate (N = 150); n=5: tcaCLR ~ aCLR (N = 150); clr: bin 5 below bin 20 (N = 400); aclr: bin 5 below bin 20 (N = 400); aclr: bin 10 intermediate (N = 400); n=5: tcaCLR ~ aCLR (N = 400); clr: bin 5 below bin 20 (N = 800); aclr: bin 5 below bin 20 (N = 800); aclr: bin 10 intermediate (N = 800); n=5: tcaCLR ~ aCLR (N = 800).

## R3

Prediction (memo): (T3, loose regulons): aCLR ~ tcaCLR > CLR, gap largest for b = 0.35 and m = 40.

Verdict: **partly held** (7 of 12 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T3 AUPR pooled aclr | 0.0335 (0.00261) [0.0289, 0.0379] | 0.112 (0.007) [0.104, 0.126] | 0.17 (0.00474) [0.164, 0.177] |
| T3 AUPR pooled tcaclr | 0.0294 (0.00392) [0.0237, 0.0351] | 0.0919 (0.0167) [0.0658, 0.12] | 0.136 (0.0232) [0.107, 0.184] |
| T3 AUPR pooled clr | 0.0168 (0.000992) [0.0152, 0.0182] | 0.0618 (0.00419) [0.0552, 0.0694] | 0.114 (0.00367) [0.109, 0.119] |
| stratum m5_b0.35: mean(aCLR,tcaCLR) - CLR (seed-mean AUPR gap) | -2.87e-05 (1.4e-05) [-5.91e-05, -1.14e-05] | 0.000178 (0.000151) [1.76e-05, 0.000411] | 0.00074 (0.000197) [0.000416, 0.00096] |
| stratum m5_b0.5: mean(aCLR,tcaCLR) - CLR (seed-mean AUPR gap) | 0.000194 (0.000302) [-0.000123, 0.000937] | 0.00112 (0.000543) [0.0003, 0.00187] | -0.0032 (0.00177) [-0.00586, -0.000592] |
| stratum m15_b0.35: mean(aCLR,tcaCLR) - CLR (seed-mean AUPR gap) | -0.000133 (4.41e-05) [-0.000203, -6.33e-05] | 0.000905 (0.000427) [0.000309, 0.00146] | 0.00422 (0.00189) [0.00236, 0.00836] |
| stratum m15_b0.5: mean(aCLR,tcaCLR) - CLR (seed-mean AUPR gap) | 0.00319 (0.00106) [0.00189, 0.00446] | 0.0101 (0.00275) [0.00573, 0.0145] | -0.00171 (0.00299) [-0.00551, 0.00313] |
| stratum m40_b0.35: mean(aCLR,tcaCLR) - CLR (seed-mean AUPR gap) | -0.000311 (0.000187) [-0.000574, -1.11e-05] | 0.0067 (0.00194) [0.00434, 0.0103] | 0.0272 (0.00292) [0.0235, 0.033] |
| stratum m40_b0.5: mean(aCLR,tcaCLR) - CLR (seed-mean AUPR gap) | 0.0257 (0.00499) [0.0202, 0.0353] | 0.0411 (0.00971) [0.0301, 0.0574] | 0.00157 (0.0104) [-0.00991, 0.0243] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T3 pooled: aCLR ~ tcaCLR | false | false | false |
| T3 pooled: tcaCLR > CLR | true | true | true |
| T3 pooled: aCLR > CLR | true | true | true |
| largest gap (mean(aCLR, tcaCLR) - CLR) at m=40, b=0.35 | false | false | true |

Failed clauses: T3 pooled: aCLR ~ tcaCLR (N = 150); largest gap (mean(aCLR, tcaCLR) - CLR) at m=40, b=0.35 (N = 150); T3 pooled: aCLR ~ tcaCLR (N = 400); largest gap (mean(aCLR, tcaCLR) - CLR) at m=40, b=0.35 (N = 400); T3 pooled: aCLR ~ tcaCLR (N = 800).

## R4

Prediction (memo): (T4, direct edges): CLR >= tcaCLR >= aCLR; |r| competitive.

Verdict: **partly held** (3 of 9 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T4 AUPR clr | 0.00322 (0.000316) [0.00276, 0.00378] | 0.0128 (0.00135) [0.0102, 0.0144] | 0.0256 (0.00113) [0.024, 0.0271] |
| T4 AUPR tcaclr | 0.00361 (0.000569) [0.00272, 0.00457] | 0.0105 (0.0022) [0.00672, 0.0138] | 0.0128 (0.00262) [0.00965, 0.0179] |
| T4 AUPR aclr | 0.00614 (0.00131) [0.00388, 0.00826] | 0.0265 (0.00526) [0.0194, 0.0397] | 0.0595 (0.0121) [0.0484, 0.0902] |
| T4 AUPR absr | 0.00419 (0.000334) [0.0038, 0.00491] | 0.0071 (0.000501) [0.00609, 0.00777] | 0.00833 (0.000248) [0.00798, 0.00861] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T4: CLR >= tcaCLR | true | true | true |
| T4: tcaCLR >= aCLR | false | false | false |
| T4: \|r\| >= 0.8 x best of CLR, aCLR, tcaCLR | false | false | false |

Failed clauses: T4: tcaCLR >= aCLR (N = 150); T4: |r| >= 0.8 x best of CLR, aCLR, tcaCLR (N = 150); T4: tcaCLR >= aCLR (N = 400); T4: |r| >= 0.8 x best of CLR, aCLR, tcaCLR (N = 400); T4: tcaCLR >= aCLR (N = 800); T4: |r| >= 0.8 x best of CLR, aCLR, tcaCLR (N = 800).

## R5

Prediction (memo): (heavy-tail arm): CLR's size penalty becomes noisier (sigma_i inflated), tcaCLR's advantage on n >= 20 grows; (copula arm): CLR and tcaCLR unchanged in rank, |r| loses.

Verdict: **partly held** (7 of 15 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| CLR slope, Gaussian noise | -0.301 (0.0436) [-0.375, -0.217] | -0.413 (0.0197) [-0.437, -0.385] | -0.456 (0.0111) [-0.476, -0.445] |
| CLR slope, Student-t noise | -0.289 (0.0512) [-0.363, -0.191] | -0.375 (0.0216) [-0.399, -0.337] | -0.435 (0.016) [-0.465, -0.418] |
| tcaCLR - CLR T2 gap n>=20 (per seed, mean over bins), Gaussian | 0.102 (0.0821) [-0.0189, 0.179] | 0.067 (0.117) [-0.108, 0.274] | 0.148 (0.096) [-0.00278, 0.313] |
| tcaCLR - CLR T2 gap n>=20 (per seed, mean over bins), Student-t | 0.146 (0.0292) [0.0954, 0.192] | 0.18 (0.0677) [0.0313, 0.246] | 0.147 (0.0476) [0.0865, 0.246] |
| T1 AUPR absr, identity | 0.822 (0.0625) [0.713, 0.909] | 0.915 (0.0231) [0.888, 0.949] | 0.928 (0.015) [0.907, 0.951] |
| T1 AUPR absr, copula | 0.72 (0.0859) [0.586, 0.837] | 0.878 (0.0353) [0.841, 0.93] | 0.915 (0.0191) [0.892, 0.944] |
| T1 AUPR clr, identity | 0.637 (0.0736) [0.523, 0.742] | 0.749 (0.0266) [0.724, 0.788] | 0.727 (0.0109) [0.712, 0.744] |
| T1 AUPR clr, copula | 0.492 (0.072) [0.392, 0.588] | 0.689 (0.0394) [0.65, 0.75] | 0.732 (0.0182) [0.707, 0.763] |
| T1 AUPR aclr, identity | 0.567 (0.0209) [0.523, 0.591] | 0.583 (0.00704) [0.571, 0.594] | 0.544 (0.00646) [0.537, 0.559] |
| T1 AUPR aclr, copula | 0.385 (0.0334) [0.328, 0.434] | 0.554 (0.00904) [0.539, 0.567] | 0.547 (0.00594) [0.539, 0.561] |
| T1 AUPR tcaclr, identity | 0.703 (0.108) [0.521, 0.828] | 0.63 (0.178) [0.349, 0.939] | 0.671 (0.166) [0.369, 0.897] |
| T1 AUPR tcaclr, copula | 0.621 (0.0867) [0.489, 0.751] | 0.772 (0.112) [0.58, 0.93] | 0.729 (0.0806) [0.603, 0.82] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| heavy tail: seed SD of CLR slope larger than Gaussian | true | true | true |
| heavy tail: tcaCLR - CLR gap (T2, n>=20, mean over bins) larger than Gaussian | true | true | false |
| copula: sign of tcaCLR - CLR on T1 equals identity arm | true | false | true |
| copula: \|r\| T1 below CLR | false | false | false |
| copula: \|r\| T1 below tcaCLR | false | false | false |

Failed clauses: copula: |r| T1 below CLR (N = 150); copula: |r| T1 below tcaCLR (N = 150); copula: sign of tcaCLR - CLR on T1 equals identity arm (N = 400); copula: |r| T1 below CLR (N = 400); copula: |r| T1 below tcaCLR (N = 400); heavy tail: tcaCLR - CLR gap (T2, n>=20, mean over bins) larger than Gaussian (N = 800); copula: |r| T1 below CLR (N = 800); copula: |r| T1 below tcaCLR (N = 800).

## R6

Prediction (memo): (identity): tcaCLR at t = 1 equals Stouffer of within-candidate row z of CLR to 1e-10.

Verdict: **held** (1 of 1 clause cells true).

| clause | N = 0 |
|---|---|
| max \|difference\| <= 1e-10 in every compendium | true |

Maximum absolute difference over 120 compendia: 4.33e-15; rows with at most 10 candidates (excluded): 0 to 0 per compendium (range).

## R7

Prediction (memo): (lemma): Sinkhorn-balanced walk, raw readout, ties aCLR on n >= 20 (within seed SD) and reduces the hub share of the top 1%.

Verdict: **partly held** (9 of 15 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T2 AUPR sinkhorn raw n=20 | 0.327 (0.0291) [0.291, 0.38] | 0.476 (0.0272) [0.441, 0.518] | 0.473 (0.0265) [0.424, 0.505] |
| T2 AUPR aCLR n=20 | 0.485 (0.0252) [0.441, 0.514] | 0.542 (0.0183) [0.52, 0.569] | 0.468 (0.0239) [0.436, 0.498] |
| T2 AUPR sinkhorn raw n=40 | 0.276 (0.0395) [0.223, 0.332] | 0.331 (0.0183) [0.305, 0.359] | 0.293 (0.0115) [0.273, 0.307] |
| T2 AUPR aCLR n=40 | 0.329 (0.0258) [0.292, 0.364] | 0.315 (0.00816) [0.299, 0.325] | 0.274 (0.00918) [0.261, 0.289] |
| T2 AUPR sinkhorn raw n=70 | 0.243 (0.0171) [0.209, 0.265] | 0.272 (0.00465) [0.265, 0.279] | 0.246 (0.00412) [0.238, 0.252] |
| T2 AUPR aCLR n=70 | 0.255 (0.0145) [0.224, 0.277] | 0.245 (0.00418) [0.238, 0.252] | 0.216 (0.00482) [0.208, 0.222] |
| T2 AUPR sinkhorn raw n=100 | 0.205 (0.0436) [0.141, 0.265] | 0.232 (0.00793) [0.222, 0.244] | 0.215 (0.0058) [0.205, 0.224] |
| T2 AUPR aCLR n=100 | 0.206 (0.0236) [0.16, 0.227] | 0.202 (0.00685) [0.193, 0.211] | 0.175 (0.00642) [0.167, 0.186] |
| hub share top 1% (high-pi set), aCLR | 0.0493 (0.00679) [0.0371, 0.057] | 0.0576 (0.00536) [0.0477, 0.0657] | 0.0582 (0.00381) [0.0516, 0.0634] |
| hub share top 1% (high-pi set), sinkhorn raw | 0.0425 (0.00549) [0.0334, 0.0496] | 0.0467 (0.00303) [0.0408, 0.0522] | 0.0422 (0.00444) [0.0335, 0.0489] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T2 n=20: sinkhorn raw ~ aCLR | false | false | true |
| T2 n=40: sinkhorn raw ~ aCLR | false | true | true |
| T2 n=70: sinkhorn raw ~ aCLR | true | true | false |
| T2 n=100: sinkhorn raw ~ aCLR | true | false | false |
| hub share of top 1%: sinkhorn raw < aCLR | true | true | true |

Failed clauses: T2 n=20: sinkhorn raw ~ aCLR (N = 150); T2 n=40: sinkhorn raw ~ aCLR (N = 150); T2 n=20: sinkhorn raw ~ aCLR (N = 400); T2 n=100: sinkhorn raw ~ aCLR (N = 400); T2 n=70: sinkhorn raw ~ aCLR (N = 800); T2 n=100: sinkhorn raw ~ aCLR (N = 800).

## R8

Prediction (memo): (depth): at t = 36 and 72 all methods' T1 falls toward the hub-score level and the hub share of the top 1% rises; aCLR's rises faster than tcaCLR's.

Verdict: **partly held** (14 of 33 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T1 AUPR of the hub-score reference (pi_i + pi_j) | 0.0228 (0.00346) [0.0182, 0.0276] | 0.0249 (0.00288) [0.0208, 0.0282] | 0.0207 (0.00232) [0.0163, 0.0233] |
| aclr T1 AUPR t=9 | 0.567 (0.0209) [0.523, 0.591] | 0.583 (0.00704) [0.571, 0.594] | 0.544 (0.00646) [0.537, 0.559] |
| aclr hub share top 1% t=9 | 0.0493 (0.00679) [0.0371, 0.057] | 0.0576 (0.00536) [0.0477, 0.0657] | 0.0582 (0.00381) [0.0516, 0.0634] |
| aclr T1 AUPR t=36 | 0.798 (0.0377) [0.719, 0.834] | 0.781 (0.0153) [0.752, 0.809] | 0.731 (0.018) [0.708, 0.767] |
| aclr hub share top 1% t=36 | 0.0773 (0.0128) [0.0569, 0.0903] | 0.0808 (0.0194) [0.06, 0.122] | 0.101 (0.0324) [0.0554, 0.146] |
| aclr T1 AUPR t=72 | 0.722 (0.0706) [0.576, 0.797] | 0.812 (0.0169) [0.779, 0.839] | 0.759 (0.0199) [0.737, 0.802] |
| aclr hub share top 1% t=72 | 0.155 (0.0483) [0.0783, 0.247] | 0.121 (0.036) [0.077, 0.169] | 0.151 (0.0701) [0.0576, 0.266] |
| tcaclr T1 AUPR t=9 | 0.703 (0.108) [0.521, 0.828] | 0.63 (0.178) [0.349, 0.939] | 0.671 (0.166) [0.369, 0.897] |
| tcaclr hub share top 1% t=9 | 0.0467 (0.00593) [0.0379, 0.0555] | 0.0499 (0.0108) [0.0329, 0.065] | 0.053 (0.00938) [0.0299, 0.0624] |
| tcaclr T1 AUPR t=36 | 0.748 (0.118) [0.581, 0.881] | 0.642 (0.181) [0.346, 0.939] | 0.681 (0.172) [0.362, 0.912] |
| tcaclr hub share top 1% t=36 | 0.074 (0.00898) [0.0588, 0.0853] | 0.0691 (0.0128) [0.0478, 0.09] | 0.0706 (0.0164) [0.0344, 0.0967] |
| tcaclr T1 AUPR t=72 | 0.675 (0.121) [ 0.5, 0.851] | 0.635 (0.181) [0.334, 0.924] | 0.675 (0.175) [0.347, 0.908] |
| tcaclr hub share top 1% t=72 | 0.11 (0.0221) [0.0812, 0.153] | 0.0881 (0.0204) [0.0667, 0.122] | 0.0855 (0.0306) [0.0352, 0.151] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| aclr: T1 at t=36 below t=9 | false | false | false |
| aclr: T1 at t=72 below t=9 | false | false | false |
| aclr: T1(72) closer to hub-score level than T1(9) | false | false | false |
| aclr: hub share at t=36 above t=9 | true | true | true |
| aclr: hub share at t=72 above t=9 | true | true | true |
| tcaclr: T1 at t=36 below t=9 | false | false | false |
| tcaclr: T1 at t=72 below t=9 | false | false | false |
| tcaclr: T1(72) closer to hub-score level than T1(9) | false | false | false |
| tcaclr: hub share at t=36 above t=9 | true | true | true |
| tcaclr: hub share at t=72 above t=9 | true | true | true |
| hub-share rise 9 -> 72: aCLR larger than tcaCLR | false | true | true |

Failed clauses: aclr: T1 at t=36 below t=9 (N = 150); aclr: T1 at t=72 below t=9 (N = 150); aclr: T1(72) closer to hub-score level than T1(9) (N = 150); tcaclr: T1 at t=36 below t=9 (N = 150); tcaclr: T1 at t=72 below t=9 (N = 150); tcaclr: T1(72) closer to hub-score level than T1(9) (N = 150); hub-share rise 9 -> 72: aCLR larger than tcaCLR (N = 150); aclr: T1 at t=36 below t=9 (N = 400); aclr: T1 at t=72 below t=9 (N = 400); aclr: T1(72) closer to hub-score level than T1(9) (N = 400); tcaclr: T1 at t=36 below t=9 (N = 400); tcaclr: T1 at t=72 below t=9 (N = 400); tcaclr: T1(72) closer to hub-score level than T1(9) (N = 400); aclr: T1 at t=36 below t=9 (N = 800); aclr: T1 at t=72 below t=9 (N = 800); aclr: T1(72) closer to hub-score level than T1(9) (N = 800); tcaclr: T1 at t=36 below t=9 (N = 800); tcaclr: T1 at t=72 below t=9 (N = 800); tcaclr: T1(72) closer to hub-score level than T1(9) (N = 800).

## R9

Prediction (memo): (satellites): aCLR ranks satellite-member pairs above member-member pairs of the same module for n >= 20; CLR and tcaCLR rank member-member higher.

Verdict: **partly held** (24 of 36 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| P(satellite-member > member-member), n=20, aCLR | 0.207 (0.0194) [0.175, 0.235] | 0.161 (0.0135) [0.142, 0.181] | 0.16 (0.0172) [0.136, 0.187] |
| P(...), n=20, CLR | 0.228 (0.0139) [0.208, 0.254] | 0.191 (0.0176) [0.164, 0.219] | 0.19 (0.0193) [0.153, 0.214] |
| P(...), n=20, tcaCLR | 0.139 (0.0235) [0.11, 0.175] | 0.0778 (0.0209) [0.0447, 0.11] | 0.0825 (0.0267) [0.0414, 0.137] |
| P(satellite-member > member-member), n=40, aCLR | 0.258 (0.0263) [0.213, 0.309] | 0.26 (0.0171) [0.239, 0.285] | 0.276 (0.016) [0.256, 0.305] |
| P(...), n=40, CLR | 0.23 (0.0202) [0.211, 0.269] | 0.21 (0.0167) [0.174, 0.236] | 0.204 (0.02) [0.168, 0.233] |
| P(...), n=40, tcaCLR | 0.162 (0.0314) [0.12, 0.217] | 0.127 (0.0292) [0.0816, 0.188] | 0.155 (0.0281) [0.108, 0.206] |
| P(satellite-member > member-member), n=70, aCLR | 0.29 (0.0428) [0.222, 0.346] | 0.308 (0.025) [0.283, 0.365] | 0.316 (0.0227) [0.293, 0.351] |
| P(...), n=70, CLR | 0.244 (0.0394) [0.171, 0.294] | 0.219 (0.0288) [0.185, 0.279] | 0.208 (0.0186) [0.162, 0.233] |
| P(...), n=70, tcaCLR | 0.203 (0.0604) [0.122, 0.293] | 0.244 (0.0529) [0.15, 0.334] | 0.262 (0.0496) [0.202, 0.355] |
| P(satellite-member > member-member), n=100, aCLR | 0.306 (0.0316) [0.268, 0.35] | 0.319 (0.0205) [0.275, 0.343] | 0.327 (0.0146) [0.311, 0.351] |
| P(...), n=100, CLR | 0.258 (0.066) [0.176, 0.363] | 0.224 (0.0479) [0.172, 0.289] | 0.222 (0.0306) [0.167, 0.262] |
| P(...), n=100, tcaCLR | 0.221 (0.0552) [0.147, 0.294] | 0.271 (0.0653) [0.169, 0.373] | 0.262 (0.0747) [0.143, 0.365] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| n=20: aCLR ranks satellite-member above member-member (P > 0.5) | false | false | false |
| n=20: CLR ranks member-member higher (P < 0.5) | true | true | true |
| n=20: tcaCLR ranks member-member higher (P < 0.5) | true | true | true |
| n=40: aCLR ranks satellite-member above member-member (P > 0.5) | false | false | false |
| n=40: CLR ranks member-member higher (P < 0.5) | true | true | true |
| n=40: tcaCLR ranks member-member higher (P < 0.5) | true | true | true |
| n=70: aCLR ranks satellite-member above member-member (P > 0.5) | false | false | false |
| n=70: CLR ranks member-member higher (P < 0.5) | true | true | true |
| n=70: tcaCLR ranks member-member higher (P < 0.5) | true | true | true |
| n=100: aCLR ranks satellite-member above member-member (P > 0.5) | false | false | false |
| n=100: CLR ranks member-member higher (P < 0.5) | true | true | true |
| n=100: tcaCLR ranks member-member higher (P < 0.5) | true | true | true |

Failed clauses: n=20: aCLR ranks satellite-member above member-member (P > 0.5) (N = 150); n=40: aCLR ranks satellite-member above member-member (P > 0.5) (N = 150); n=70: aCLR ranks satellite-member above member-member (P > 0.5) (N = 150); n=100: aCLR ranks satellite-member above member-member (P > 0.5) (N = 150); n=20: aCLR ranks satellite-member above member-member (P > 0.5) (N = 400); n=40: aCLR ranks satellite-member above member-member (P > 0.5) (N = 400); n=70: aCLR ranks satellite-member above member-member (P > 0.5) (N = 400); n=100: aCLR ranks satellite-member above member-member (P > 0.5) (N = 400); n=20: aCLR ranks satellite-member above member-member (P > 0.5) (N = 800); n=40: aCLR ranks satellite-member above member-member (P > 0.5) (N = 800); n=70: aCLR ranks satellite-member above member-member (P > 0.5) (N = 800); n=100: aCLR ranks satellite-member above member-member (P > 0.5) (N = 800).

## R10

Prediction (memo): (wide-row arm, keff 40): sigma_bg falls, tcaCLR's T1 gain over aCLR shrinks or vanishes; CLR unaffected.

Verdict: **partly held** (5 of 6 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| median sigma_bg (nats), keff 10 | 1.41 (0.138) [1.21, 1.57] | 1.68 (0.116) [1.51, 1.81] |  1.8 (0.0809) [1.66, 1.92] |
| median sigma_bg (nats), keff 40 | 0.683 (0.0496) [0.612, 0.748] | 0.803 (0.0436) [0.73, 0.846] |  0.9 (0.0297) [0.84, 0.952] |
| T1 gap tcaCLR - aCLR, keff 10 | 0.137 (0.11) [-0.0422, 0.259] | 0.0473 (0.177) [-0.233, 0.363] | 0.126 (0.163) [-0.172, 0.35] |
| T1 gap tcaCLR - aCLR, keff 40 | -0.0529 (0.0231) [-0.0822, -0.0182] | -0.0596 (0.0208) [-0.0974, -0.0306] | -0.0445 (0.0134) [-0.0659, -0.0231] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| median sigma_bg at keff 40 below keff 10 | true | true | true |
| tcaCLR - aCLR T1 gap smaller at keff 40 than at keff 10 | true | false | true |

Failed clauses: tcaCLR - aCLR T1 gap smaller at keff 40 than at keff 10 (N = 400).

## R11

Prediction (memo): (lift): pi log-lift unstable (large seed SD, hub-dominated top 1%); column-mean log-lift between aCLR and tcaCLR on T1.

Verdict: **partly held** (3 of 18 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T1 AUPR pi_lift | 0.757 (0.0822) [0.631, 0.878] | 0.91 (0.0162) [0.889, 0.938] | 0.915 (0.011) [0.895, 0.929] |
| relative SD of T1 AUPR pi_lift (single value across seeds) | 0.109 (  NA) [0.109, 0.109] | 0.0178 (  NA) [0.0178, 0.0178] | 0.012 (  NA) [0.012, 0.012] |
| T1 AUPR colmean_lift | 0.792 (0.0721) [0.668, 0.894] | 0.924 (0.00972) [0.909, 0.945] | 0.919 (0.00789) [0.906, 0.929] |
| relative SD of T1 AUPR colmean_lift (single value across seeds) | 0.091 (  NA) [0.091, 0.091] | 0.0105 (  NA) [0.0105, 0.0105] | 0.00858 (  NA) [0.00858, 0.00858] |
| T1 AUPR aclr | 0.567 (0.0209) [0.523, 0.591] | 0.583 (0.00704) [0.571, 0.594] | 0.544 (0.00646) [0.537, 0.559] |
| relative SD of T1 AUPR aclr (single value across seeds) | 0.0369 (  NA) [0.0369, 0.0369] | 0.0121 (  NA) [0.0121, 0.0121] | 0.0119 (  NA) [0.0119, 0.0119] |
| T1 AUPR tcaclr | 0.703 (0.108) [0.521, 0.828] | 0.63 (0.178) [0.349, 0.939] | 0.671 (0.166) [0.369, 0.897] |
| relative SD of T1 AUPR tcaclr (single value across seeds) | 0.154 (  NA) [0.154, 0.154] | 0.283 (  NA) [0.283, 0.283] | 0.248 (  NA) [0.248, 0.248] |
| T1 AUPR clr | 0.637 (0.0736) [0.523, 0.742] | 0.749 (0.0266) [0.724, 0.788] | 0.727 (0.0109) [0.712, 0.744] |
| relative SD of T1 AUPR clr (single value across seeds) | 0.116 (  NA) [0.116, 0.116] | 0.0355 (  NA) [0.0355, 0.0355] | 0.015 (  NA) [0.015, 0.015] |
| hub share top 1%, pi log-lift | 0.0453 (0.00688) [0.0335, 0.0524] | 0.0513 (0.00466) [0.0425, 0.0602] | 0.0476 (0.00677) [0.0347, 0.0582] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| relative SD of T1 AUPR: pi log-lift larger than clr | false | false | false |
| relative SD of T1 AUPR: pi log-lift larger than aclr | true | true | true |
| relative SD of T1 AUPR: pi log-lift larger than tcaclr | false | false | false |
| hub share of top 1%: pi log-lift above aclr | false | false | false |
| hub share of top 1%: pi log-lift above tcaclr | false | false | false |
| column-mean log-lift T1 strictly between aCLR and tcaCLR | false | false | false |

Failed clauses: relative SD of T1 AUPR: pi log-lift larger than clr (N = 150); relative SD of T1 AUPR: pi log-lift larger than tcaclr (N = 150); hub share of top 1%: pi log-lift above aclr (N = 150); hub share of top 1%: pi log-lift above tcaclr (N = 150); column-mean log-lift T1 strictly between aCLR and tcaCLR (N = 150); relative SD of T1 AUPR: pi log-lift larger than clr (N = 400); relative SD of T1 AUPR: pi log-lift larger than tcaclr (N = 400); hub share of top 1%: pi log-lift above aclr (N = 400); hub share of top 1%: pi log-lift above tcaclr (N = 400); column-mean log-lift T1 strictly between aCLR and tcaCLR (N = 400); relative SD of T1 AUPR: pi log-lift larger than clr (N = 800); relative SD of T1 AUPR: pi log-lift larger than tcaclr (N = 800); hub share of top 1%: pi log-lift above aclr (N = 800); hub share of top 1%: pi log-lift above tcaclr (N = 800); column-mean log-lift T1 strictly between aCLR and tcaCLR (N = 800).

## R12

Prediction (memo): (inflation control): fixed-exponent inflation does not reach tcaCLR on T1 for both r; Sinkhorn(A) -> walk -> tcaCLR is at or above tcaCLR on hub-stratified T1 (speculation; a failure here is informative).

Verdict: **partly held** (8 of 9 clause cells true).

| quantity | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T1 AUPR infl_1.5 | 0.33 (0.0157) [0.309, 0.365] | 0.341 (0.0118) [0.327, 0.355] | 0.316 (0.0133) [0.301, 0.332] |
| T1 AUPR infl_2 | 0.316 (0.0173) [0.293, 0.354] | 0.321 (0.0128) [0.305, 0.338] | 0.294 (0.0139) [0.277, 0.312] |
| hub-stratified T1 (true hub loading sum) sink_tca | 0.324 (0.0432) [0.26, 0.385] | 0.347 (0.0691) [0.216, 0.428] | 0.364 (0.0745) [0.209, 0.436] |
| hub-stratified T1 (MI row mean) sink_tca [reported] | 0.41 (0.0254) [0.381, 0.461] | 0.458 (0.0469) [0.36, 0.514] | 0.458 (0.0448) [0.369, 0.512] |
| hub-stratified T1 (true hub loading sum) tcaclr | 0.345 (0.0472) [0.269, 0.406] | 0.334 (0.0713) [0.205, 0.417] | 0.349 (0.0723) [0.204, 0.43] |
| hub-stratified T1 (MI row mean) tcaclr [reported] | 0.484 (0.0347) [0.433, 0.537] | 0.466 (0.0465) [0.371, 0.532] | 0.453 (0.0472) [0.37, 0.515] |
| T1 AUPR tcaCLR | 0.703 (0.108) [0.521, 0.828] | 0.63 (0.178) [0.349, 0.939] | 0.671 (0.166) [0.369, 0.897] |

| clause | N = 150 | N = 400 | N = 800 |
|---|---|---|---|
| T1: tcaCLR > infl_1.5 (inflation does not reach tcaCLR) | true | true | true |
| T1: tcaCLR > infl_2 (inflation does not reach tcaCLR) | true | true | true |
| hub-stratified T1 (true hub loading sum): sinkhorn -> tcaCLR >= tcaCLR | false | true | true |

Failed clauses: hub-stratified T1 (true hub loading sum): sinkhorn -> tcaCLR >= tcaCLR (N = 150).

## Mechanism curves (primary arm)

Mean within-module pair score by module size (seed means; aCLR and CLR on the raw scale, tcaCLR in z units), with fitted slopes (module-level regression on log n with a lambda fixed effect; n >= 20 primary, all sizes secondary).

N = 150

| method | n=5 | n=10 | n=20 | n=40 | n=70 | n=100 | slope n>=20 (log) | slope n>=20 (raw) | slope all (log) |
|---|---|---|---|---|---|---|---|---|---|
| aclr | 0.0493 | 0.0499 | 0.0353 | 0.0207 | 0.0125 | 0.00889 | -0.772 (0.0628) [-0.881, -0.673] | -0.0175 (0.00171) [-0.0205, -0.0154] | -0.31 (0.0657) [-0.389, -0.182] |
| clr | 14.9 | 12.8 | 10.3 | 8.27 | 6.76 | 5.77 | -0.301 (0.0436) [-0.375, -0.217] | -2.79 (0.398) [-3.52, -2.07] | -0.226 (0.0253) [-0.263, -0.19] |
| tcaclr | 5.59 | 6.12 | 6.69 | 6.56 | 5.21 | 4.29 | -0.342 (0.29) [-0.721, -0.0495] | -1.21 (0.722) [-2.16, -0.279] | 0.00475 (0.0873) [-0.123, 0.108] |
| absr | 0.501 |  0.5 | 0.498 | 0.496 | 0.498 | 0.494 | -0.00466 (0.0345) [-0.0708, 0.0356] | 0.000803 (0.0151) [-0.0279, 0.0206] | -0.00236 (0.00657) [-0.0138, 0.00493] |

N = 400

| method | n=5 | n=10 | n=20 | n=40 | n=70 | n=100 | slope n>=20 (log) | slope n>=20 (raw) | slope all (log) |
|---|---|---|---|---|---|---|---|---|---|
| aclr | 0.066 | 0.0655 | 0.044 | 0.0236 | 0.0138 | 0.0097 | -0.925 (0.0193) [-0.956, -0.9] | -0.0231 (0.000834) [-0.0244, -0.0221] | -0.505 (0.0202) [-0.538, -0.462] |
| clr | 24.1 | 18.6 |   14 | 10.5 | 8.25 | 6.93 | -0.413 (0.0197) [-0.437, -0.385] | -4.49 (0.314) [-4.89, -4.08] | -0.364 (0.0245) [-0.403, -0.324] |
| tcaclr | 5.65 | 6.11 | 6.86 | 6.83 | 4.32 |  3.3 | -0.688 (0.375) [-1.31, -0.0588] | -1.96 (0.887) [-3.4, -0.394] | -0.114 (0.108) [-0.296, 0.0812] |
| absr | 0.502 | 0.501 | 0.504 | 0.499 | 0.506 | 0.497 | -0.00244 (0.012) [-0.0196, 0.0179] | -0.00015 (0.00534) [-0.00778, 0.0107] | 0.00372 (0.00553) [-0.00499, 0.00957] |

N = 800

| method | n=5 | n=10 | n=20 | n=40 | n=70 | n=100 | slope n>=20 (log) | slope n>=20 (raw) | slope all (log) |
|---|---|---|---|---|---|---|---|---|---|
| aclr | 0.0705 | 0.0712 | 0.0454 | 0.0239 | 0.0138 | 0.00974 | -0.948 (0.0125) [-0.974, -0.932] | -0.0241 (0.000546) [-0.0251, -0.0234] | -0.534 (0.00895) [-0.548, -0.52] |
| clr | 29.2 | 21.5 | 15.3 | 11.3 | 8.65 | 7.24 | -0.456 (0.0111) [-0.476, -0.445] | -5.22 (0.195) [-5.55, -5.02] | -0.442 (0.0151) [-0.461, -0.418] |
| tcaclr | 5.56 | 6.09 | 6.87 | 6.87 | 3.87 | 4.09 | -0.624 (0.367) [-1.28, -0.0893] | -1.84 (0.832) [-3.25, -0.56] | -0.0937 (0.11) [-0.272, 0.0743] |
| absr | 0.503 | 0.503 | 0.502 | 0.505 | 0.506 | 0.494 | 0.00223 (0.0093) [-0.0145, 0.0171] | 0.000718 (0.00418) [-0.00488, 0.00625] | 0.00329 (0.00493) [-0.00609, 0.0102] |

## Per-row sigma_bg, retention m(n), hub share versus t (primary arm)

Median per-row sigma_bg (SD of log P^t over reachable off-diagonal entries, nats), coefficient of variation across rows, and unreachable share, by depth t (tcaCLR operator, keff 10):

| N | t | median sigma_bg | CV across rows | unreachable share |
|---|---|---|---|---|
| 150 | 3 | 2.42 (0.118) [2.25, 2.53] | 0.294 (0.0577) [0.216, 0.413] | 0.0796 (0.0205) [0.0478, 0.111] |
| 150 | 6 | 1.66 (0.142) [1.45, 1.82] | 0.256 (0.0715) [0.159, 0.383] | 0.0228 (0.0262) [   0, 0.0569] |
| 150 | 9 | 1.41 (0.138) [1.21, 1.57] | 0.25 (0.0624) [0.159, 0.357] | 0.0228 (0.0262) [   0, 0.0569] |
| 150 | 18 | 1.15 (0.126) [0.959,  1.3] | 0.251 (0.0622) [0.156, 0.357] | 0.0228 (0.0262) [   0, 0.0569] |
| 150 | 36 | 1.06 (0.0953) [0.899, 1.18] | 0.225 (0.0603) [0.141, 0.341] | 0.0228 (0.0262) [   0, 0.0569] |
| 150 | 72 | 1.06 (0.0614) [0.946, 1.14] | 0.188 (0.0614) [0.121, 0.315] | 0.0228 (0.0262) [   0, 0.0569] |
| 400 | 3 |  2.7 (0.102) [2.56, 2.81] | 0.296 (0.0501) [0.229, 0.416] | 0.124 (0.011) [0.112, 0.142] |
| 400 | 6 | 1.93 (0.121) [1.76, 2.08] | 0.247 (0.0653) [0.176, 0.386] | 0.0539 (0.0342) [0.00265, 0.114] |
| 400 | 9 | 1.68 (0.116) [1.51, 1.81] | 0.233 (0.0568) [0.173, 0.344] | 0.0536 (0.0346) [   0, 0.114] |
| 400 | 18 |  1.4 (0.103) [1.25, 1.51] | 0.22 (0.0563) [0.162, 0.328] | 0.0536 (0.0346) [   0, 0.114] |
| 400 | 36 | 1.29 (0.079) [1.17, 1.37] | 0.188 (0.0557) [0.138, 0.305] | 0.0536 (0.0346) [   0, 0.114] |
| 400 | 72 | 1.29 (0.0561) [ 1.2, 1.35] | 0.15 (0.0548) [0.0969, 0.275] | 0.0536 (0.0346) [   0, 0.114] |
| 800 | 3 | 2.87 (0.0795) [2.75,    3] | 0.351 (0.0242) [0.303, 0.38] | 0.148 (0.0048) [0.144, 0.159] |
| 800 | 6 | 2.06 (0.0866) [1.92, 2.19] | 0.276 (0.0497) [0.196, 0.349] | 0.0488 (0.0351) [0.000812, 0.114] |
| 800 | 9 |  1.8 (0.0809) [1.66, 1.92] | 0.252 (0.0382) [0.191, 0.313] | 0.0484 (0.0354) [   0, 0.114] |
| 800 | 18 | 1.51 (0.0668) [1.39, 1.61] | 0.232 (0.0355) [0.176, 0.294] | 0.0484 (0.0354) [   0, 0.114] |
| 800 | 36 | 1.39 (0.0452) [ 1.3, 1.46] | 0.195 (0.0334) [0.148, 0.256] | 0.0484 (0.0354) [   0, 0.114] |
| 800 | 72 |  1.4 (0.03) [1.36, 1.44] | 0.151 (0.031) [0.113, 0.21] | 0.0484 (0.0354) [   0, 0.114] |

sigma_bg at t = 1 and t = 9 for the keff 10, keff 40 and Sinkhorn operators (median across rows; CV across rows):

| N | operator | median sigma_bg | CV | unreachable share |
|---|---|---|---|---|
| 150 | keff 10, t = 1 | 1.39 (0.0125) [1.37, 1.41] | 0.46 (0.0428) [0.409, 0.556] | 0.983 (   0) [0.983, 0.983] |
| 150 | keff 10, t = 9 | 1.41 (0.138) [1.21, 1.57] | 0.25 (0.0624) [0.159, 0.357] | 0.0228 (0.0262) [   0, 0.0569] |
| 150 | keff 40, t = 9 | 0.683 (0.0496) [0.612, 0.748] | 0.273 (0.036) [0.213, 0.318] |    0 (   0) [   0,    0] |
| 150 | Sinkhorn, t = 9 | 1.86 (0.216) [1.55, 2.07] | 0.18 (0.0178) [0.143, 0.198] | 0.0228 (0.0262) [   0, 0.0569] |
| 400 | keff 10, t = 1 | 1.43 (0.0114) [1.41, 1.44] | 0.63 (0.0508) [0.552, 0.709] | 0.983 (   0) [0.983, 0.983] |
| 400 | keff 10, t = 9 | 1.68 (0.116) [1.51, 1.81] | 0.233 (0.0568) [0.173, 0.344] | 0.0536 (0.0346) [   0, 0.114] |
| 400 | keff 40, t = 9 | 0.803 (0.0436) [0.73, 0.846] | 0.273 (0.0188) [0.245,  0.3] |    0 (   0) [   0,    0] |
| 400 | Sinkhorn, t = 9 | 2.21 (0.139) [2.05, 2.39] | 0.214 (0.0254) [0.185, 0.265] | 0.0536 (0.0346) [   0, 0.114] |
| 800 | keff 10, t = 1 | 1.45 (0.0121) [1.43, 1.47] | 0.73 (0.0665) [0.631, 0.87] | 0.983 (   0) [0.983, 0.983] |
| 800 | keff 10, t = 9 |  1.8 (0.0809) [1.66, 1.92] | 0.252 (0.0382) [0.191, 0.313] | 0.0484 (0.0354) [   0, 0.114] |
| 800 | keff 40, t = 9 |  0.9 (0.0297) [0.84, 0.952] | 0.262 (0.0103) [0.249, 0.278] |    0 (   0) [   0,    0] |
| 800 | Sinkhorn, t = 9 | 2.33 (0.0852) [2.23, 2.51] | 0.222 (0.0248) [0.193, 0.274] | 0.0484 (0.0354) [   0, 0.114] |

In-module retention by module size: the share of attention A on the other members (t = 1, A only) and the mass sum_{j in module} P^t_ij (self included) at t = 1 and t = 9 (mean over members and modules of that size, then seeds):

| N | n | A in-module share | m_1(n) | m_9(n) | m_9(n), keff 40 |
|---|---|---|---|---|---|
| 150 | 5 | 0.622 (0.0253) [0.582, 0.661] | 0.811 (0.0126) [0.791, 0.83] | 0.256 (0.00983) [0.243, 0.276] | 0.0783 (0.00268) [0.0751, 0.0837] |
| 150 | 10 | 0.776 (0.0219) [0.746, 0.801] | 0.888 (0.011) [0.873, 0.901] | 0.512 (0.0118) [0.493, 0.532] | 0.146 (0.00359) [0.139, 0.152] |
| 150 | 20 | 0.864 (0.0233) [0.822, 0.892] | 0.932 (0.0116) [0.911, 0.946] | 0.726 (0.0312) [0.675, 0.76] | 0.282 (0.0133) [0.263,  0.3] |
| 150 | 40 | 0.93 (0.0258) [0.882, 0.963] | 0.965 (0.0129) [0.941, 0.982] | 0.851 (0.0527) [0.77, 0.926] | 0.514 (0.0564) [0.43,  0.6] |
| 150 | 70 | 0.952 (0.0123) [0.928, 0.966] | 0.976 (0.00617) [0.964, 0.983] | 0.899 (0.0282) [0.851, 0.934] | 0.687 (0.0297) [0.642, 0.736] |
| 150 | 100 | 0.957 (0.0367) [0.879, 0.998] | 0.979 (0.0183) [0.939, 0.999] | 0.914 (0.0741) [0.761, 0.997] | 0.751 (0.145) [0.535, 0.968] |
| 400 | 5 | 0.759 (0.00645) [0.746, 0.768] | 0.879 (0.00322) [0.873, 0.884] | 0.334 (0.00846) [0.315, 0.345] |  0.1 (0.00219) [0.0954, 0.103] |
| 400 | 10 | 0.902 (0.00606) [0.887, 0.907] | 0.951 (0.00303) [0.943, 0.954] | 0.663 (0.0145) [0.628, 0.679] | 0.184 (0.00358) [0.175, 0.188] |
| 400 | 20 | 0.953 (0.00534) [0.943, 0.96] | 0.976 (0.00267) [0.971, 0.98] | 0.896 (0.0147) [0.875, 0.917] | 0.367 (0.0129) [0.349, 0.391] |
| 400 | 40 | 0.983 (0.00344) [0.978, 0.99] | 0.992 (0.00172) [0.989, 0.995] | 0.967 (0.00718) [0.955, 0.977] | 0.655 (0.0361) [0.61, 0.702] |
| 400 | 70 | 0.992 (0.00265) [0.987, 0.997] | 0.996 (0.00133) [0.994, 0.999] | 0.987 (0.00444) [0.978, 0.993] | 0.853 (0.0129) [0.827, 0.878] |
| 400 | 100 | 0.996 (0.00293) [0.988, 0.998] | 0.998 (0.00146) [0.994, 0.999] | 0.994 (0.00285) [0.99, 0.998] | 0.911 (0.0451) [0.87, 0.976] |
| 800 | 5 | 0.778 (0.00266) [0.774, 0.783] | 0.889 (0.00133) [0.887, 0.892] | 0.355 (0.00446) [0.349, 0.364] | 0.106 (0.00111) [0.104, 0.108] |
| 800 | 10 | 0.925 (0.00285) [0.92, 0.928] | 0.962 (0.00142) [0.96, 0.964] | 0.717 (0.00982) [0.704, 0.73] | 0.197 (0.0025) [0.193,  0.2] |
| 800 | 20 | 0.96 (0.00564) [0.953, 0.97] | 0.98 (0.00282) [0.977, 0.985] | 0.924 (0.00983) [0.91, 0.943] | 0.392 (0.00837) [0.38, 0.404] |
| 800 | 40 | 0.987 (0.00303) [0.981, 0.992] | 0.993 (0.00151) [0.99, 0.996] | 0.979 (0.00563) [0.968, 0.988] | 0.711 (0.0295) [0.674, 0.764] |
| 800 | 70 | 0.995 (0.00216) [0.992, 0.997] | 0.997 (0.00108) [0.996, 0.999] | 0.991 (0.00351) [0.985, 0.995] | 0.888 (0.016) [0.863, 0.921] |
| 800 | 100 | 0.997 (0.00158) [0.994, 0.999] | 0.999 (0.000788) [0.997, 0.999] | 0.996 (0.0017) [0.992, 0.998] | 0.941 (0.0206) [0.914, 0.97] |

Hub share of the top 1% of pairs (pairs containing at least one of the 30 highest-pi genes) and T1 AUPR versus depth; chance level 2.0%, hub-score reference in R8:

| N | method | t | T1 AUPR | hub share top 1% | planted-hub share |
|---|---|---|---|---|---|
| 150 | aclr | 3 | 0.308 (0.0151) [0.277, 0.326] | 0.034 (0.00329) [0.0273, 0.0374] | 0.0289 (0.0013) [0.0269, 0.0306] |
| 150 | aclr | 6 | 0.441 (0.0185) [0.402, 0.461] | 0.0446 (0.00568) [0.0337, 0.051] | 0.0319 (0.00173) [0.0295, 0.0349] |
| 150 | aclr | 9 | 0.567 (0.0209) [0.523, 0.591] | 0.0493 (0.00679) [0.0371, 0.057] | 0.034 (0.00221) [0.0308, 0.0378] |
| 150 | aclr | 18 | 0.762 (0.0245) [0.716, 0.79] | 0.0573 (0.00773) [0.0442, 0.066] | 0.0394 (0.00347) [0.0345, 0.0461] |
| 150 | aclr | 36 | 0.798 (0.0377) [0.719, 0.834] | 0.0773 (0.0128) [0.0569, 0.0903] | 0.0389 (0.0102) [0.0214, 0.0598] |
| 150 | aclr | 72 | 0.722 (0.0706) [0.576, 0.797] | 0.155 (0.0483) [0.0783, 0.247] | 0.0265 (0.0206) [0.0109, 0.082] |
| 150 | tcaclr | 3 | 0.38 (0.03) [0.325, 0.43] | 0.0299 (0.00702) [0.0195, 0.0403] | 0.0368 (0.00286) [0.032, 0.0402] |
| 150 | tcaclr | 6 | 0.666 (0.105) [0.482, 0.79] | 0.0453 (0.00633) [0.0356, 0.0544] | 0.0339 (0.00316) [0.0275, 0.0382] |
| 150 | tcaclr | 9 | 0.703 (0.108) [0.521, 0.828] | 0.0467 (0.00593) [0.0379, 0.0555] | 0.0356 (0.00337) [0.0303, 0.0411] |
| 150 | tcaclr | 18 | 0.753 (0.116) [0.577, 0.879] | 0.0552 (0.00569) [0.047, 0.0624] | 0.0435 (0.00624) [0.0361, 0.0543] |
| 150 | tcaclr | 36 | 0.748 (0.118) [0.581, 0.881] | 0.074 (0.00898) [0.0588, 0.0853] | 0.0623 (0.0143) [0.0376, 0.0816] |
| 150 | tcaclr | 72 | 0.675 (0.121) [ 0.5, 0.851] | 0.11 (0.0221) [0.0812, 0.153] | 0.0712 (0.0228) [0.0307, 0.101] |
| 400 | aclr | 3 | 0.347 (0.00581) [0.336, 0.354] | 0.0378 (0.0022) [0.0339, 0.0406] | 0.0268 (0.000761) [0.0258, 0.0281] |
| 400 | aclr | 6 | 0.477 (0.00665) [0.466, 0.486] | 0.0508 (0.00393) [0.0429, 0.0558] | 0.029 (0.000792) [0.0279, 0.0304] |
| 400 | aclr | 9 | 0.583 (0.00704) [0.571, 0.594] | 0.0576 (0.00536) [0.0477, 0.0657] | 0.0291 (0.000832) [0.028, 0.0305] |
| 400 | aclr | 18 | 0.722 (0.0108) [0.701, 0.743] | 0.0646 (0.00726) [0.054, 0.0774] | 0.0312 (0.00147) [0.0292, 0.0334] |
| 400 | aclr | 36 | 0.781 (0.0153) [0.752, 0.809] | 0.0808 (0.0194) [0.06, 0.122] | 0.0338 (0.00197) [0.0303, 0.0375] |
| 400 | aclr | 72 | 0.812 (0.0169) [0.779, 0.839] | 0.121 (0.036) [0.077, 0.169] | 0.034 (0.0035) [0.0292, 0.0383] |
| 400 | tcaclr | 3 | 0.373 (0.0453) [0.317, 0.452] | 0.0361 (0.00715) [0.0257, 0.0459] | 0.0423 (0.00223) [0.0385, 0.046] |
| 400 | tcaclr | 6 | 0.62 (0.175) [0.346, 0.928] | 0.048 (0.0107) [0.0303, 0.0637] | 0.0353 (0.00501) [0.0299, 0.0437] |
| 400 | tcaclr | 9 | 0.63 (0.178) [0.349, 0.939] | 0.0499 (0.0108) [0.0329, 0.065] | 0.036 (0.00535) [0.0299, 0.0447] |
| 400 | tcaclr | 18 | 0.639 (0.18) [0.35, 0.941] | 0.0575 (0.0109) [0.0398, 0.0742] | 0.0434 (0.00939) [0.0318, 0.0604] |
| 400 | tcaclr | 36 | 0.642 (0.181) [0.346, 0.939] | 0.0691 (0.0128) [0.0478, 0.09] | 0.0519 (0.014) [0.035, 0.0784] |
| 400 | tcaclr | 72 | 0.635 (0.181) [0.334, 0.924] | 0.0881 (0.0204) [0.0667, 0.122] | 0.0619 (0.0189) [0.0399, 0.102] |
| 800 | aclr | 3 | 0.338 (0.0055) [0.328, 0.346] | 0.038 (0.00195) [0.0346, 0.0402] | 0.0264 (0.000332) [0.0258, 0.0269] |
| 800 | aclr | 6 | 0.454 (0.00593) [0.445, 0.464] | 0.0498 (0.00265) [0.0461, 0.0539] | 0.0277 (0.000218) [0.0274, 0.028] |
| 800 | aclr | 9 | 0.544 (0.00646) [0.537, 0.559] | 0.0582 (0.00381) [0.0516, 0.0634] | 0.0277 (0.000242) [0.0274, 0.0281] |
| 800 | aclr | 18 | 0.67 (0.0128) [0.652, 0.696] | 0.072 (0.0121) [0.0536, 0.0926] | 0.0282 (0.000302) [0.0276, 0.0286] |
| 800 | aclr | 36 | 0.731 (0.018) [0.708, 0.767] | 0.101 (0.0324) [0.0554, 0.146] | 0.0299 (0.00143) [0.0282, 0.0325] |
| 800 | aclr | 72 | 0.759 (0.0199) [0.737, 0.802] | 0.151 (0.0701) [0.0576, 0.266] | 0.0309 (0.00226) [0.0283, 0.035] |
| 800 | tcaclr | 3 | 0.268 (0.012) [0.252, 0.291] | 0.0254 (0.00512) [0.0173, 0.0315] | 0.0429 (0.00233) [0.0402, 0.0462] |
| 800 | tcaclr | 6 | 0.66 (0.164) [0.366, 0.882] | 0.0501 (0.00961) [0.0283, 0.0616] | 0.0323 (0.00256) [0.0289, 0.0378] |
| 800 | tcaclr | 9 | 0.671 (0.166) [0.369, 0.897] | 0.053 (0.00938) [0.0299, 0.0624] | 0.0323 (0.00266) [0.0285, 0.0382] |
| 800 | tcaclr | 18 | 0.678 (0.168) [0.368, 0.906] | 0.0603 (0.0107) [0.0323, 0.0702] | 0.034 (0.00387) [0.0288, 0.0429] |
| 800 | tcaclr | 36 | 0.681 (0.172) [0.362, 0.912] | 0.0706 (0.0164) [0.0344, 0.0967] | 0.0373 (0.00593) [0.0289, 0.0476] |
| 800 | tcaclr | 72 | 0.675 (0.175) [0.347, 0.908] | 0.0855 (0.0306) [0.0352, 0.151] | 0.0406 (0.00783) [0.0295, 0.0531] |

## All methods, primary arm (pooled AUPR; base rates T1 0.00539, T3 0.00212, T4 0.000133)

N = 150

| method | T1 | T3 | T4 | T1 hub-strat (true) | T1 hub-strat (MI) | satellite share top 1% |
|---|---|---|---|---|---|---|
| absr | 0.822 (0.0625) [0.713, 0.909] | 0.0263 (0.00136) [0.0245, 0.0285] | 0.00419 (0.000334) [0.0038, 0.00491] | 0.373 (0.0246) [0.33, 0.402] | 0.509 (0.0266) [0.467, 0.545] | 0.0413 (0.00349) [0.037, 0.0469] |
| clr | 0.637 (0.0736) [0.523, 0.742] | 0.0168 (0.000992) [0.0152, 0.0182] | 0.00322 (0.000316) [0.00276, 0.00378] | 0.296 (0.0223) [0.257, 0.327] | 0.384 (0.0217) [0.343, 0.422] | 0.0315 (0.00375) [0.0253, 0.0368] |
| aclr | 0.567 (0.0209) [0.523, 0.591] | 0.0335 (0.00261) [0.0289, 0.0379] | 0.00614 (0.00131) [0.00388, 0.00826] | 0.317 (0.0151) [0.296, 0.336] | 0.424 (0.0206) [0.392, 0.46] | 0.0341 (0.00229) [0.0313, 0.0372] |
| tcaclr | 0.703 (0.108) [0.521, 0.828] | 0.0294 (0.00392) [0.0237, 0.0351] | 0.00361 (0.000569) [0.00272, 0.00457] | 0.345 (0.0472) [0.269, 0.406] | 0.484 (0.0347) [0.433, 0.537] | 0.0404 (0.00443) [0.0319, 0.0463] |
| tcaclr_t1 | 0.144 (0.00976) [0.126, 0.157] | 0.0143 (0.000985) [0.0124, 0.0154] | 0.00589 (0.00107) [0.00431, 0.00711] | 0.12 (0.00733) [0.105, 0.13] | 0.142 (0.0113) [0.12, 0.156] | 0.0192 (0.00167) [0.0167, 0.0215] |
| sink_raw | 0.545 (0.0368) [0.481, 0.597] | 0.0237 (0.00183) [0.0201, 0.0267] | 0.00286 (0.000206) [0.00254, 0.00306] | 0.336 (0.0157) [0.313, 0.357] | 0.399 (0.0191) [0.366, 0.425] | 0.0356 (0.0026) [0.0328, 0.0395] |
| sink_tca | 0.637 (0.105) [0.433, 0.754] | 0.0216 (0.00259) [0.0172, 0.0262] | 0.00242 (0.000304) [0.00196, 0.00293] | 0.324 (0.0432) [0.26, 0.385] | 0.41 (0.0254) [0.381, 0.461] | 0.0379 (0.00427) [0.0289, 0.0433] |
| pi_lift | 0.757 (0.0822) [0.631, 0.878] | 0.0235 (0.00179) [0.0202, 0.0254] | 0.00255 (0.000229) [0.00224, 0.0029] | 0.355 (0.0181) [0.321, 0.376] | 0.429 (0.0246) [0.397, 0.469] | 0.0416 (0.00435) [0.0374, 0.0491] |
| colmean_lift | 0.792 (0.0721) [0.668, 0.894] | 0.0232 (0.00177) [0.0197, 0.0257] | 0.00248 (0.0002) [0.00215, 0.00271] | 0.371 (0.0219) [0.336, 0.395] | 0.462 (0.0296) [0.416, 0.515] | 0.0433 (0.00402) [0.0385, 0.0507] |
| infl_1.5 | 0.33 (0.0157) [0.309, 0.365] | 0.0315 (0.00219) [0.0278, 0.0356] | 0.0096 (0.00218) [0.00616, 0.0132] | 0.248 (0.0125) [0.23, 0.269] | 0.294 (0.0177) [0.265, 0.327] | 0.0264 (0.00169) [0.0243, 0.0286] |
| infl_2 | 0.316 (0.0173) [0.293, 0.354] | 0.0282 (0.00213) [0.0246, 0.0322] | 0.00895 (0.0021) [0.00566, 0.0124] | 0.237 (0.0136) [0.218, 0.262] | 0.28 (0.0172) [0.252, 0.308] | 0.0268 (0.0017) [0.0247, 0.0288] |
| aclr_k40 | 0.821 (0.0433) [0.733, 0.877] | 0.0269 (0.00211) [0.0234, 0.0307] | 0.00383 (0.000561) [0.00285, 0.00476] | 0.389 (0.0221) [0.358, 0.416] | 0.515 (0.0309) [0.462, 0.557] | 0.0448 (0.00367) [0.0399, 0.0512] |
| tcaclr_k40 | 0.768 (0.0572) [0.657, 0.849] | 0.0296 (0.00249) [0.0268, 0.0351] | 0.00526 (0.00133) [0.00381, 0.0079] | 0.377 (0.0198) [0.348, 0.405] | 0.497 (0.0287) [0.453, 0.551] | 0.0449 (0.00375) [0.0399, 0.0518] |

N = 400

| method | T1 | T3 | T4 | T1 hub-strat (true) | T1 hub-strat (MI) | satellite share top 1% |
|---|---|---|---|---|---|---|
| absr | 0.915 (0.0231) [0.888, 0.949] | 0.0518 (0.0033) [0.0456, 0.0579] | 0.0071 (0.000501) [0.00609, 0.00777] | 0.443 (0.0238) [0.406, 0.473] | 0.517 (0.0308) [0.47, 0.545] | 0.0493 (0.00184) [0.0469, 0.0521] |
| clr | 0.749 (0.0266) [0.724, 0.788] | 0.0618 (0.00419) [0.0552, 0.0694] | 0.0128 (0.00135) [0.0102, 0.0144] | 0.374 (0.0137) [0.354, 0.395] | 0.462 (0.0167) [0.436, 0.482] | 0.0438 (0.00212) [0.0411, 0.0479] |
| aclr | 0.583 (0.00704) [0.571, 0.594] | 0.112 (0.007) [0.104, 0.126] | 0.0265 (0.00526) [0.0194, 0.0397] | 0.322 (0.0132) [0.299, 0.339] | 0.429 (0.0177) [0.402, 0.455] | 0.0402 (0.00104) [0.0385, 0.0423] |
| tcaclr | 0.63 (0.178) [0.349, 0.939] | 0.0919 (0.0167) [0.0658, 0.12] | 0.0105 (0.0022) [0.00672, 0.0138] | 0.334 (0.0713) [0.205, 0.417] | 0.466 (0.0465) [0.371, 0.532] | 0.0433 (0.00689) [0.0314, 0.054] |
| tcaclr_t1 | 0.184 (0.00283) [0.181, 0.189] | 0.0472 (0.00227) [0.0446, 0.052] | 0.0264 (0.0038) [0.0195, 0.0322] | 0.152 (0.00214) [0.148, 0.155] | 0.184 (0.00704) [0.174, 0.193] | 0.026 (0.000872) [0.0247, 0.0271] |
| sink_raw | 0.62 (0.00788) [0.611, 0.634] | 0.0822 (0.00605) [0.0745, 0.0928] | 0.00858 (0.000591) [0.00765, 0.00964] | 0.369 (0.0162) [0.342, 0.388] | 0.455 (0.0222) [0.421, 0.482] | 0.0439 (0.00188) [0.0412, 0.0467] |
| sink_tca | 0.629 (0.179) [0.345, 0.934] | 0.0726 (0.0112) [0.0555, 0.0876] | 0.00675 (0.00125) [0.00488, 0.00854] | 0.347 (0.0691) [0.216, 0.428] | 0.458 (0.0469) [0.36, 0.514] | 0.0435 (0.00658) [0.0325, 0.054] |
| pi_lift | 0.91 (0.0162) [0.889, 0.938] | 0.0614 (0.00418) [0.0561, 0.0683] | 0.00576 (0.000345) [0.00517, 0.00646] | 0.445 (0.0263) [0.406, 0.473] | 0.516 (0.0399) [0.451, 0.555] | 0.0529 (0.0016) [0.0497, 0.054] |
| colmean_lift | 0.924 (0.00972) [0.909, 0.945] | 0.0637 (0.00378) [0.0597, 0.0709] | 0.00576 (0.000296) [0.00533, 0.00644] | 0.451 (0.0278) [0.409, 0.477] | 0.529 (0.0439) [0.462, 0.577] | 0.0534 (0.00109) [0.0514, 0.054] |
| infl_1.5 | 0.341 (0.0118) [0.327, 0.355] | 0.0966 (0.00426) [0.0924, 0.104] | 0.046 (0.00994) [0.0362, 0.0688] | 0.255 (0.0101) [0.237, 0.268] | 0.313 (0.013) [0.297, 0.335] | 0.0303 (0.00111) [0.0285, 0.0323] |
| infl_2 | 0.321 (0.0128) [0.305, 0.338] | 0.0898 (0.00405) [0.0857, 0.0956] | 0.0448 (0.0103) [0.0349, 0.0694] | 0.241 (0.0113) [0.222, 0.256] | 0.294 (0.0117) [0.279, 0.316] | 0.0304 (0.00116) [0.0282, 0.0322] |
| aclr_k40 | 0.857 (0.012) [0.844, 0.874] | 0.079 (0.00584) [0.0715, 0.0906] | 0.0149 (0.00317) [0.0113, 0.0226] | 0.414 (0.0213) [0.381, 0.443] | 0.507 (0.0275) [0.466, 0.536] | 0.0538 (0.000168) [0.0535, 0.054] |
| tcaclr_k40 | 0.798 (0.027) [0.758, 0.842] | 0.088 (0.0107) [0.0744, 0.104] | 0.0216 (0.00388) [0.0156, 0.0268] | 0.392 (0.0168) [0.365, 0.415] | 0.479 (0.0268) [0.436, 0.51] | 0.054 (2.35e-05) [0.054, 0.054] |

N = 800

| method | T1 | T3 | T4 | T1 hub-strat (true) | T1 hub-strat (MI) | satellite share top 1% |
|---|---|---|---|---|---|---|
| absr | 0.928 (0.015) [0.907, 0.951] | 0.0664 (0.00205) [0.064, 0.0694] | 0.00833 (0.000248) [0.00798, 0.00861] | 0.452 (0.0283) [0.403, 0.48] | 0.505 (0.0185) [0.478, 0.534] | 0.0524 (0.00043) [0.0519, 0.0531] |
| clr | 0.727 (0.0109) [0.712, 0.744] | 0.114 (0.00367) [0.109, 0.119] | 0.0256 (0.00113) [0.024, 0.0271] | 0.364 (0.017) [0.336, 0.383] | 0.46 (0.0123) [0.445, 0.483] | 0.0486 (0.00122) [0.0472, 0.0511] |
| aclr | 0.544 (0.00646) [0.537, 0.559] | 0.17 (0.00474) [0.164, 0.177] | 0.0595 (0.0121) [0.0484, 0.0902] | 0.299 (0.0111) [0.28, 0.312] | 0.403 (0.0125) [0.381, 0.418] | 0.0387 (0.000759) [0.0369, 0.0395] |
| tcaclr | 0.671 (0.166) [0.369, 0.897] | 0.136 (0.0232) [0.107, 0.184] | 0.0128 (0.00262) [0.00965, 0.0179] | 0.349 (0.0723) [0.204, 0.43] | 0.453 (0.0472) [0.37, 0.515] | 0.0444 (0.00706) [0.0313, 0.054] |
| tcaclr_t1 | 0.191 (0.0019) [0.189, 0.195] | 0.0769 (0.00208) [0.0737, 0.08] | 0.0444 (0.00287) [0.0393, 0.0487] | 0.155 (0.00101) [0.154, 0.156] | 0.19 (0.00713) [0.182, 0.201] | 0.0282 (0.000708) [0.0269, 0.029] |
| sink_raw | 0.598 (0.00556) [0.59, 0.606] | 0.134 (0.00543) [0.124, 0.141] | 0.0124 (0.000468) [0.0119, 0.0134] | 0.355 (0.0155) [0.329, 0.372] | 0.45 (0.0116) [0.435, 0.469] | 0.0442 (0.000938) [0.043, 0.0456] |
| sink_tca | 0.686 (0.171) [0.361, 0.905] | 0.108 (0.0172) [0.0878, 0.146] | 0.00853 (0.00153) [0.00696, 0.0119] | 0.364 (0.0745) [0.209, 0.436] | 0.458 (0.0448) [0.369, 0.512] | 0.0447 (0.00676) [0.0323, 0.054] |
| pi_lift | 0.915 (0.011) [0.895, 0.929] | 0.0911 (0.00271) [0.0858, 0.0942] | 0.00722 (0.000147) [0.007, 0.00744] | 0.445 (0.029) [0.395, 0.47] | 0.507 (0.02) [0.472, 0.535] | 0.0534 (0.00114) [0.0509, 0.054] |
| colmean_lift | 0.919 (0.00789) [0.906, 0.929] | 0.096 (0.00273) [0.0907, 0.0995] | 0.00739 (0.00014) [0.00718, 0.00756] | 0.447 (0.0302) [0.394, 0.474] | 0.505 (0.0218) [0.467, 0.537] | 0.0535 (0.00102) [0.0514, 0.054] |
| infl_1.5 | 0.316 (0.0133) [0.301, 0.332] | 0.136 (0.00445) [0.129, 0.142] | 0.0968 (0.00844) [0.0826, 0.113] | 0.237 (0.01) [0.222, 0.249] | 0.294 (0.012) [0.275, 0.311] | 0.029 (0.00104) [0.0271, 0.0305] |
| infl_2 | 0.294 (0.0139) [0.277, 0.312] | 0.126 (0.00464) [0.118, 0.132] | 0.0968 (0.00831) [0.0828, 0.113] | 0.219 (0.0111) [0.204, 0.234] | 0.271 (0.00949) [0.255, 0.285] | 0.0291 (0.00118) [0.0269, 0.0307] |
| aclr_k40 | 0.821 (0.0125) [0.801, 0.833] | 0.122 (0.00443) [0.116, 0.129] | 0.0414 (0.0127) [0.0299, 0.0745] | 0.393 (0.019) [0.362, 0.42] | 0.481 (0.0137) [0.464, 0.501] | 0.0538 (0.00017) [0.0535, 0.054] |
| tcaclr_k40 | 0.776 (0.0188) [0.754, 0.81] | 0.131 (0.00908) [0.119, 0.142] | 0.0432 (0.0134) [0.0279, 0.0668] | 0.379 (0.017) [0.349, 0.397] | 0.452 (0.0158) [0.426, 0.481] | 0.054 (   0) [0.054, 0.054] |

T2 per-bin AUPR (primary arm, seed means; per-seed values in the unit files):

| N | method | n=5 | n=10 | n=20 | n=40 | n=70 | n=100 |
|---|---|---|---|---|---|---|---|
| 150 | absr | 0.357 | 0.458 | 0.55 | 0.641 | 0.718 | 0.707 |
| 150 | clr | 0.586 | 0.609 | 0.582 | 0.504 | 0.419 | 0.327 |
| 150 | aclr | 0.496 | 0.597 | 0.485 | 0.329 | 0.255 | 0.206 |
| 150 | tcaclr | 0.126 | 0.339 | 0.631 | 0.695 | 0.532 | 0.381 |
| 150 | sink_raw | 0.299 | 0.411 | 0.327 | 0.276 | 0.243 | 0.205 |
| 150 | aclr_k40 | 0.462 | 0.523 | 0.614 | 0.676 | 0.643 | 0.564 |
| 150 | tcaclr_k40 | 0.487 | 0.54 | 0.58 | 0.552 | 0.504 | 0.442 |
| 400 | absr | 0.402 | 0.497 | 0.612 | 0.713 | 0.81 | 0.811 |
| 400 | clr | 0.775 | 0.769 | 0.677 | 0.502 | 0.438 | 0.38 |
| 400 | aclr | 0.612 | 0.742 | 0.542 | 0.315 | 0.245 | 0.202 |
| 400 | tcaclr | 0.0802 | 0.266 | 0.72 | 0.819 | 0.446 | 0.28 |
| 400 | sink_raw | 0.634 | 0.663 | 0.476 | 0.331 | 0.272 | 0.232 |
| 400 | aclr_k40 | 0.532 | 0.586 | 0.726 | 0.752 | 0.647 | 0.534 |
| 400 | tcaclr_k40 | 0.558 | 0.574 | 0.595 | 0.583 | 0.522 | 0.435 |
| 800 | absr | 0.411 | 0.509 | 0.616 | 0.735 | 0.821 | 0.831 |
| 800 | clr | 0.869 | 0.818 | 0.605 | 0.421 | 0.38 | 0.356 |
| 800 | aclr | 0.509 | 0.698 | 0.468 | 0.274 | 0.216 | 0.175 |
| 800 | tcaclr | 0.0662 | 0.267 | 0.743 | 0.842 | 0.371 | 0.398 |
| 800 | sink_raw | 0.764 | 0.76 | 0.473 | 0.293 | 0.246 | 0.215 |
| 800 | aclr_k40 | 0.416 | 0.51 | 0.672 | 0.723 | 0.564 | 0.457 |
| 800 | tcaclr_k40 | 0.46 | 0.477 | 0.489 | 0.562 |  0.5 | 0.412 |

T2 per-bin AUPR in the other arms (tcaCLR minus aCLR and CLR minus aCLR, seed mean):

| N | arm | contrast | n=5 | n=10 | n=20 | n=40 | n=70 | n=100 |
|---|---|---|---|---|---|---|---|---|
| 150 | gaussian identity | tcaclr - aclr | -0.37 | -0.257 | 0.146 | 0.366 | 0.277 | 0.175 |
| 150 | gaussian identity | clr - aclr | 0.0904 | 0.0119 | 0.0972 | 0.175 | 0.163 | 0.121 |
| 150 | t4 identity | tcaclr - aclr | -0.237 | -0.178 |  0.2 | 0.363 | 0.42 | 0.334 |
| 150 | t4 identity | clr - aclr | 0.214 | 0.121 | 0.213 | 0.236 | 0.181 | 0.102 |
| 150 | gaussian copula | tcaclr - aclr | -0.198 | -0.0984 | 0.202 | 0.324 | 0.359 | 0.217 |
| 150 | gaussian copula | clr - aclr | 0.185 | 0.123 | 0.194 | 0.206 | 0.145 | 0.0777 |
| 150 | t4 copula | tcaclr - aclr | -0.00801 | -0.00437 | 0.147 | 0.159 | 0.138 | 0.0805 |
| 150 | t4 copula | clr - aclr | 0.362 | 0.33 | 0.301 | 0.168 | 0.0737 | 0.0371 |
| 400 | gaussian identity | tcaclr - aclr | -0.531 | -0.476 | 0.178 | 0.504 | 0.201 | 0.0779 |
| 400 | gaussian identity | clr - aclr | 0.163 | 0.027 | 0.135 | 0.187 | 0.193 | 0.178 |
| 400 | t4 identity | tcaclr - aclr | -0.39 | -0.34 | 0.197 | 0.453 | 0.398 | 0.377 |
| 400 | t4 identity | clr - aclr | 0.234 | 0.0852 | 0.175 | 0.216 | 0.179 | 0.137 |
| 400 | gaussian copula | tcaclr - aclr | -0.417 | -0.357 | 0.156 | 0.456 | 0.416 | 0.261 |
| 400 | gaussian copula | clr - aclr | 0.169 | 0.0575 | 0.138 | 0.193 | 0.173 | 0.131 |
| 400 | t4 copula | tcaclr - aclr | -0.0221 | 0.054 | 0.29 | 0.317 | 0.239 | 0.114 |
| 400 | t4 copula | clr - aclr | 0.43 | 0.327 | 0.257 | 0.137 | 0.0651 | 0.0504 |
| 800 | gaussian identity | tcaclr - aclr | -0.443 | -0.431 | 0.276 | 0.568 | 0.155 | 0.223 |
| 800 | gaussian identity | clr - aclr | 0.36 | 0.12 | 0.137 | 0.147 | 0.164 | 0.181 |
| 800 | t4 identity | tcaclr - aclr | -0.365 | -0.329 | 0.296 | 0.511 | 0.31 | 0.206 |
| 800 | t4 identity | clr - aclr | 0.346 | 0.143 | 0.178 | 0.184 | 0.189 | 0.183 |
| 800 | gaussian copula | tcaclr - aclr | -0.404 | -0.413 | 0.197 | 0.526 | 0.232 | 0.319 |
| 800 | gaussian copula | clr - aclr | 0.307 | 0.121 | 0.159 | 0.187 | 0.186 | 0.174 |
| 800 | t4 copula | tcaclr - aclr | -0.0437 | 0.0381 | 0.269 | 0.318 |  0.2 | 0.0957 |
| 800 | t4 copula | clr - aclr | 0.389 | 0.224 | 0.156 | 0.0788 | 0.0573 | 0.0573 |

## Run record

Compendia 120; wall per compendium  208 (9.04) [ 183,  235] s (R process peak  2.5 (0.163) [2.03, 2.81] GB).
Sinkhorn balancing of the lazy operator hit the 2,000-iteration cap in 115 of 120 compendia; final max |column sum - 1| 0.000298 (0.000202) [9.91e-10, 0.000563] (tolerance 1e-9 not reached where the cap was hit; the walk used the last iterate).
