# Preregistration: attention width (keff) variants

*Fixed 2026-10-16, before any of these variants were scored against a truth. Designed by the statistics reviewer (Opus, high effort) at the PI's request to test all three "smarter keff" options empirically. Any change after this commit is logged below with its date and reason.*

## Shared settings

Top-50 positive Stouffer CLR candidates, lazy walk α = 0.5, the same empty-row handling, S_t = (P^t + P^t′)/2, the same seeds and lab splits. For every variant, report the distribution of row perplexities and the fractions of rows below 3 and above 30.

## Variants

| Code | Variant | Form |
|---|---|---|
| K3, K5, K25 | keff sweep | per-row perplexity calibration at keff = 3, 5, 25 |
| K10 | reference | keff = 10 (current method) |
| MH | multi-head | S_t averaged over heads keff = 3, 10, 30 |
| O1 | one global temperature | w_ij ∝ exp(z_ij / τ) over the cap-50 candidates; primary: width-matched τ = median over rows of the bisected τ_i at keff = 10, per dataset. Secondary O1-tuned: τ ∈ {0.25, 0.35, 0.5, 0.7, 1, 1.4, 2, 2.8, 4} chosen jointly with t on the held-out criterion |
| O2 | per-gene keff from local FDR | primary: permutation null (CLR on permuted data) pooled over pairs; f and f0 conditioned on z > 0 (zeros get lfdr = 1); π0 by Storey on the permutation p-values; keff_i = clamp(Σ over top 50 of (1 − lfdr_ij), 2, 30). Secondary: per-row Efron central matching (central 50%) on the unclipped row-standardized MI z, same clamp |
| O3-A | keff chosen on held-out data | grid {3, 5, 7, 10, 14, 20, 25, 30, 40} jointly with t; existing criterion (AP against the other half's top-1% CLR); width bias reported |
| O3-B | same, hub-robust criterion (primary O3) | mean over genes of the per-row AP of S_t[i, ·] against gene i's top-20 held-out CLR partners |

## Depth

Primary: each variant's own held-out t* (median over lab splits; joint argmax for O1-tuned and O3). Secondary: t = 9 fixed. Diagnostic only: full-data truth curves for t = 1–32, with the truth-chosen best t labelled as an upper bound and never used for claims.

## Truths and endpoints

- M3D and EcoMAC: operons, TF regulons (operon pairs excluded), sigmulons (σ70 and σS excluded).
- K562 essential: CORUM, knockdown-response modules, DoRothEA A–B, CollecTRI.

**Primary endpoint:** difference in pooled co-membership AUPR from K10, each variant at its held-out t*, on four truths: M3D TF regulons, EcoMAC TF regulons, K562 CollecTRI, K562 CORUM. Paired delete-half jackknife over genes; Holm over 8 variants (K3, K5, K25, MH, O1, O2, O3-A, O3-B) × 4 endpoints = 32 tests, familywise α = 0.05. Smallest effect of interest: ±0.005 AUPR.

**Secondary:** regulon-bootstrap coherence (B = 1,000), the remaining truths, the fine-to-coarse granularity pattern.

Truths are touched once, after every held-out selection is frozen and committed.

## Decision rules

The headline (K10) changes only if a variant beats K10 by a Holm-significant margin ≥ 0.005 on at least 3 of the 4 primary endpoints, loses significantly by more than 0.005 on none, and the gains replicate in both M3D and EcoMAC. If O1 or O2 meets this, it becomes the method. If MH gains on CORUM without losing on regulons, MH becomes the main operator.

The headline does not change for: a win on one truth only or on K562 only; wins that need t chosen on the truths; differences within ±0.005; the sweep's expected fine/coarse trade-off (reported as a resolution result in theme 4); O3 selecting another keff without meeting the rule (reported as sensitivity).

## Compute

E. coli: about 20 operators × 2 directions × about 60 s over t ∈ {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32}, 20 lab splits, 4 workers: about 3–4 h per dataset. K562: 10 splits on a coarse t grid for held-out t*, truth scoring only at each variant's t* and t = 9; a separate job. All queued in the main lane after the multiscale pilot, never concurrently. The pilot itself compares label-free spectra for keff 3, 5, 10, 25 and 30.

## Change log

- 2026-10-16: fixed.
- 2026-10-21 (before any full run; reviewer GO-WITH-CHANGES): implementation choices logged as deviations or clarifications of `analysis/keff_test.R`:
  - O1: the width-matched temperature tau is computed on the full-data CLR (median over rows of the bisected tau_i at keff 10) and reused on the training halves.
  - O2: conditional densities on z > 0 are combined with an unconditional pi0 (Storey, lambda = 0.5); on full M3D, P0(z > 0)/P(z > 0) = 0.92 and the keff distribution is unchanged.
  - K562 uses the same t grid as E. coli for held-out selection; K562 truths are scored only at each variant's t* and at t = 9.
  - Decision rule, literal reading: a win needs Holm p < 0.05 against d = 0 and a point estimate >= +0.005 (a loss mirrors it).
  - "Regulons" in the multi-head rule = M3D TF regulons, EcoMAC TF regulons and K562 CollecTRI.
  - Disclosure: smoke tests of the script on 300-gene subsets (all three datasets) scored truths before this entry; O2 on M3D TF regulons showed +0.016/+0.018 there. No design change follows from them.
  - O2/O2e degeneracy, checked on full M3D without truths: on the CLR scale pi0 = 0.94 and the permutation density exceeds the observed density for CLR 2-6, so lfdr is about 1 and 83% of rows sit at the floor of 2; a permutation null on raw MI gives pi0 = 0.016 and 100% of rows at the ceiling of 30; O2e has 94% of rows at the ceiling. The two-group local-FDR approach cannot choose a width in this data on any scale. O2 runs as preregistered; if it meets the decision rule it is reported with its realised keff distribution as a sharp-attention (about keff 2) result, not as an lfdr-calibrated width.
  - Protocol: held-out selections are computed first for all three datasets (`--stop_after_selection 1`), committed to `docs/keff_selections/`, and the truth phase refuses to run unless its selections match the committed copies.
