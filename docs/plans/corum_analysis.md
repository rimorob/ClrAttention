# CORUM loss of CLR attention on K562: diagnosis, cheap fixes, and the principled operator

*22 October 2026. Expert memo (network inference, latent-variable models, diffusion on graphs, attention; Opus, max effort) for the PI. Not committed.*

> **Status: EXPLORATORY, NOT PREREGISTERED.** About 20 operators and readouts were scored on the same K562 truths; treat every effect below as hypothesis-generating. No width (keff) or multi-head variant was scored against a truth: they are preregistered (`docs/preregistration_keff.md`, "truths are touched once"). No E. coli truth was touched. Any fix adopted from this memo must be preregistered before confirmatory use on E. coli (section 6).

## 1. Bottom line

1. **The softmax creates the loss and the walk partly repairs it**, not the other way round. CORUM pooled AUPR: CLR 6.39× random, softmax alone (t = 1) 2.00×, then 2.63, 2.85, 2.96, 3.20 and 3.51× at t = 2, 3, 4, 6 and 9. The walk recovers about a third of what the softmax loses.
2. **Mass dilution is real, but it is about rows (genes), not about pooled module sizes.** Attention turns *strength* into *share*. Each gene has a budget of 1 spread over about 10 effective partners. A ribosomal protein has about 77 interchangeable strong partners, so 76% of the pairs inside large complexes get no weight at t = 1. Meanwhile, peripheral genes such as INTS8, EIF4B, FBH1 and EIF2S3 spend 71-84% of their budget on 60S members. Under the arithmetic symmetrization (P + Pᵀ)/2, these one-sided negatives outrank RP-RP pairs.
3. **After the walk the problem is normalization (calibration), not order or resolution.** At t = 9, attention orders each ribosomal row nearly as well as CLR (row-local AP 0.47 vs 0.58), and 88% of the walk mass stays inside the complex: there is no leak into a "wider program". But the scale differs from row to row, so the per-complex AP collapses (60S: 0.19 vs CLR 0.66) and so does the pooled AP.
4. **Fix (a) works, in the log domain only: CLR → attention → CLR(log Pᵗ).**
   - Each row of log Pᵗ is standardized against its own background (over reachable pairs), then combined with CLR's Stouffer rule. The result (robust version) beats attention on all four jackknifed truths:
     - CORUM 7.95× (vs attention 3.51×, CLR 6.39×, |r| 7.00×);
     - DoRothEA A-B +0.003 [+0.000, +0.006] over attention;
     - CollecTRI +0.005 [+0.002, +0.007] over attention;
     - knockdown response +0.059 [+0.052, +0.067] over attention.
   - Hub-stratified versions of the non-robust forms (a2, a3) keep the CORUM, DoRothEA and CollecTRI gains.
   - The same re-normalization on raw probabilities fails (2.46×).
5. **The PI's latent-variable picture is the right explanation, but explicit coarse-graining fails the guard.** Collapsing cliques, attending among latent factors and expanding back reaches CORUM 7.62×, but costs DoRothEA A-B −0.013 [−0.017, −0.009] and CollecTRI −0.007 against attention. The walk already does the aggregation softly. What is missing is the "expand without dividing by n" step, and the log-domain CLR supplies it.

## 2. Setup and reproduction

- **Input:** `results/baselines_k562ess_truthgenes_genie3/expression.tsv`, 5,671 genes × 707 knockdown profiles (half A). Truth sets are read verbatim from that run's `regulon_members.csv`, and pair universes are rebuilt with `comembership_pairs()`.
- **Operator:** MI with B-splines (bins "hg") → Stouffer CLR → softmax (keff 10, cap 50) → lazy walk α = 0.5 → S_t = (Pᵗ + Pᵗ′)/2.
- **Reproduction:** CLR 0.0847 (6.39×), |r| 0.0928 (7.00×) and attention t = 9 0.0465 (3.51×) on CORUM. On DoRothEA A-B, CollecTRI and knockdown response, attention scores 0.2003, 0.1402 and 0.1870. All of these are identical to `metrics.csv`. The hub-stratified references also match (CLR 0.0451, attention 0.0275).
- **Intervals:**
  - Paired delete-half gene jackknife, 100 replicates. The splits are the same as in `baselines.R` (seed 20261016; the same RNG use is reproduced).
  - Each interval is d ± 1.96 SD, as in `baselines.R`.
  - DoRothEA A-C has point estimates only.
  - Networks are held fixed, so the intervals do not resample experiments.
- **Compute:** about 42 min in total (nice 19, 4 BLAS/OpenMP threads, ≤ 20 GB memory scope).

## 3. Diagnosis

### 3.1 Where the loss is

- **The gap:** the pooled AP difference between CLR and attention is 0.038. Decomposing it exactly by positive pair (fractional weights for pairs in several complexes) gives:
  - **80%** from the 8 RP-dominated complexes (≥ 30% RPL/RPS members): cytoplasmic ribosome 38%, Nop56p pre-rRNA 24%, 60S 10%, the two 40S entries 6%, TRBP-60S 2%;
  - **90%** from the top 10 complexes (adding mitochondrial 55S 5% and PA700-20S-PA28 proteasome 3.5%).
- **Most complexes carry no signal for any method.** The median within-complex CLR is 0.4 in complexes of 5-24 genes, and the mean per-complex AP of CLR is 0.003 for complexes of 5-9 genes.
- **Without the RP-dominated complexes** (1,177 complexes left, base rate 0.0119), attention ties CLR: 2.05× vs 2.18×, difference −0.0016 [−0.0042, +0.0010]. The CORUM loss is essentially one biological module, counted about six times through overlapping CORUM entries.

### 3.2 Softmax vs walk

| t | CORUM × random | KD response | DoRothEA A-B | CollecTRI |
|---|---|---|---|---|
| CLR | 6.39 | 0.191 | 0.194 | 0.132 |
| 1 (softmax, no walk) | 2.00 | 0.099 | 0.185 | 0.129 |
| 2 | 2.63 | 0.144 | 0.194 | 0.132 |
| 4 | 2.96 | 0.166 | 0.201 | 0.137 |
| 9 | 3.51 | 0.187 | 0.200 | 0.140 |

The softmax costs CORUM 0.058 AP, which is 152% of the final gap. The walk then recovers 0.020. On TF regulons the walk is the active ingredient: the softmax alone is at random level. So **the softmax causes the complex loss and the walk produces the regulon gain.** The all-genes granularity curve shows the same shape: CORUM 2.1× at t = 1, with a maximum of 3.9× at t ≈ 17.

### 3.3 The mass-dilution hypothesis, quantified

- **Budget and truncation.** The share of within-complex pairs that receive any softmax weight in either direction is 5% (complexes of ≤ 9 genes), 6% (10-24), 25% (25-49) and 24% (≥ 50). Ribosomal rows have about 77 strong partners against cap 50 and perplexity 10.
  - The per-pair weight falls with complex size (slope of log mean weight on log size −0.50); CLR does not (+0.09 z per e-fold of size).
  - The per-complex loss tracks the members' number of strong partners (CLR ≥ 3: Spearman 0.38, p 7e-19) more than complex size (0.16).
- **Leakage by the walk: not supported.**
  - Ribosomal rows in complexes of ≥ 50 genes keep 93% (t = 1), 91% (t = 2) and 88% (t = 9) of their non-self mass inside their complex.
  - Only 0.6% → 2.1% goes to genes that share no CORUM complex with them (60S rows: RPS genes, EEF1A1, then NACA, RACK1 and EIF3 at ≤ 0.6%).
  - Per-partner mass at t = 9 is about 0.009: the retained ~0.88 divided among 77-95 partners. That is pure dilution, not leakage.
- **One-sided attention (the pooled-ranking mechanism).** Non-members INTS8 (0.84), EIF4B (0.81), FBH1 (0.76), MCM3AP (0.74), EIF2S3 (0.72) and OSBPL8 (0.71) put most of their attention budget on 60S members. Under (P + Pᵀ)/2, these negatives score above RP-RP positives, whose mass is split.
  - Among the top 5,000 pooled pairs: precision is 0.47 for CLR and 0.25 for attention at t = 9; RP-RP pairs make up 34% vs 14%; pairs with a top-1% π gene 32% vs 40%.
  - RP-complex positives (n = 6,002): median pooled percentile 97.3 (CLR), 49 (t = 1, where most are tied at zero) and 95.2 (t = 9). Their share in the top 1% is 44%, 23% and 29%.
  - The other 27,838 positives lose nothing (median 59.6 vs 60.3).
- **Softmax temperature.** It is not the "many near-equal partners" story. Ribosomal rows have *higher* bisected temperatures (median 0.67 vs 0.40) because their CLR scores are large and spread out. The constraint that bites is the fixed perplexity: 10 effective partners out of 77.

### 3.4 Local vs pooled ranking: the decisive test

| Readout | Row-local AP, RP member rows | Per-complex AP, 60S | Per-complex AP, cytoplasmic ribosome |
|---|---|---|---|
| CLR | 0.578 | 0.661 | 0.661 |
| attention t = 1 | 0.300 | 0.146 | 0.186 |
| attention t = 9 | 0.474 | 0.187 | 0.255 |
| CLR on raw S9 | 0.379 | 0.082 | 0.128 |
| CLR on log attention (a2) | 0.668 | 0.837 | 0.827 |

Row-local AP ranks complex partners against genes sharing no complex, within one row, with no comparison across rows.

- **The within-row order is damaged by the softmax and mostly repaired by the walk** (0.30 → 0.47).
- **The comparison across rows is what collapses** (0.47 → 0.19), and CLR does not show that drop (0.58 → 0.66). This is the signature of a calibration problem. Row-wise re-normalization is the matching fix, and it lifts both columns.
- Across all 499 complexes with ≥ 10 pairs, CLR beats attention in 74%. The log-domain readout beats attention in 82% and CLR in 49%.

## 4. Cheap fixes (exploratory; K562 truth genes, t = 9)

The operator is unchanged (CLR → softmax keff 10 cap 50 → lazy walk α 0.5). Only the readout changes, except in (b). Δ is AUPR minus attention's AUPR, with paired jackknife 95% intervals.

| Variant | CORUM × | Δ CORUM | Δ DoRothEA A-B | Δ CollecTRI | Δ KD response | Verdict |
|---|---|---|---|---|---|---|
| attention S9 (reference) | 3.51 | 0 | 0 | 0 | 0 | |
| CLR (reference) | 6.39 | +0.038 [+0.019, +0.057] | −0.006 [−0.011, −0.001] | −0.009 [−0.013, −0.004] | +0.004 [−0.004, +0.012] | |
| **(a1)** CLR on raw S9 | 2.46 | −0.014 [−0.022, −0.006] | −0.008 | −0.007 | −0.044 | fails |
| **(a2)** CLR on log geometric-mean walk, (log Pᵗ + log Pᵗ′)/2 | 8.28 | +0.063 [+0.039, +0.088] | +0.009 [+0.005, +0.012] | +0.0005 [−0.002, +0.003] | +0.007 [−0.000, +0.013] | passes |
| **(a3)** bilateral CLR on log Pᵗ (each row z, then Stouffer) | 9.38 | +0.078 [+0.050, +0.106] | +0.004 [+0.0005, +0.008] | +0.007 [+0.003, +0.010] | +0.078 [+0.068, +0.088] | passes |
| **(a3r)** same, row statistics over reachable pairs only | 7.95 | +0.059 [+0.033, +0.084] | +0.003 [+0.000, +0.006] | +0.005 [+0.002, +0.007] | +0.059 [+0.052, +0.067] | **passes; recommended form** |
| (a2r) a2 over reachable pairs only | 7.44 | +0.052 [+0.029, +0.075] | +0.006 [+0.002, +0.010] | −0.002 [−0.005, +0.001] | −0.005 [−0.012, +0.001] | passes, weaker |
| (d) log geometric mean, no re-normalization | 4.39 | +0.012 [+0.005, +0.019] | +0.009 [+0.005, +0.012] | +0.003 [+0.000, +0.005] | +0.007 [+0.000, +0.014] | partial |
| (d) diffusion density log(Pᵗ_ij/π_j), symmetrized | 1.35 | −0.029 | +0.001 | −0.013 | −0.102 | fails (π max/median 18,457) |
| (d) CLR on log lift | 4.26 | +0.010 | +0.004 | −0.009 | −0.062 | fails |
| (d) attention on residuals (unit factor removed) | 2.25 | −0.017 | −0.000 | +0.000 | −0.055 | fails |
| (d) same + within-unit pairs first | 6.30 | +0.037 | −0.000 | +0.000 | −0.044 | neutral |
| **(b)** coarse-grain r ≥ 0.6 (138 units, 611 genes) | 7.62 | +0.055 [+0.027, +0.082] | −0.013 [−0.017, −0.009] | −0.007 [−0.011, −0.003] | +0.013 [−0.000, +0.026] | **fails guard** |
| (b) r ≥ 0.4 (426 units, 2,381 genes) | 6.62 | +0.041 | −0.018 | −0.009 | −0.039 | fails |
| (b) r ≥ 0.8 (24 units, 83 genes) | 6.21 | +0.036 | −0.004 [−0.005, −0.002] | −0.002 [−0.004, −0.001] | +0.030 [+0.022, +0.037] | fails guard (small) |
| (b0) units' pairs first, then plain S9 (control) | 6.76 | +0.043 [+0.021, +0.065] | −0.000 | −0.000 | +0.005 [+0.003, +0.007] | ties TF; ad hoc |

Notes on the table:

- **Against CLR** (same intervals):
  - a2: CORUM +0.025 [+0.017, +0.033], DoRothEA A-B +0.015, CollecTRI +0.009, KD +0.003 (ns);
  - a3r: CORUM +0.021 [+0.011, +0.030], DoRothEA A-B +0.010, CollecTRI +0.014, KD +0.055.
- **Coarse-graining (b).** Units come from average linkage on 1 − r; the unit factor is the mean of standardized members. MI, CLR, softmax and the walk run among units plus singletons. When expanding back, within-unit pairs are ranked first (by CLR) and other pairs inherit their units' score. Within-unit pairs at r ≥ 0.6 have CORUM precision 0.74.
- **Hub control** (within-stratum percentile in 10 strata of π_i + π_j, compared under the same adjustment):
  - a2 vs attention: CORUM +0.037 [+0.016, +0.057], DoRothEA A-B +0.007 [+0.004, +0.010], CollecTRI +0.007 [+0.005, +0.009], KD +0.034. a2 vs CLR on CORUM: +0.019 [+0.008, +0.030].
  - a3 vs attention: CORUM +0.044 [+0.020, +0.068], DoRothEA A-B +0.005 [+0.001, +0.008], CollecTRI +0.006 [+0.003, +0.009], KD +0.046 [+0.039, +0.054].
  - The CORUM and TF gains are therefore not hub effects. The KD gain shrinks but stays.
- **CORUM without the RP-dominated complexes:** a2 2.73× and a3r 2.53×, against attention 2.05× and CLR 2.18× (a2 vs attention +0.008 [+0.004, +0.012]). The fix also helps outside the ribosome.
- **Depth and reachability caveat.**
  - The log readout needs a walk that has reached most genes. The share of pairs below the numerical floor (Pᵗ = 0, unreachable) is 84% at t = 2, 2.7% at t = 4 and 0.28% at t = 9; at t = 9 these are mostly columns of genes that nobody attends to.
  - With a fixed 1e-300 floor, a2 is poor at t = 2 (2.28×) and t = 4 (3.92×, with TF losses). The robust form (a3r) excludes unreachable pairs from the row mean and SD, so it does not depend on the floor, and it keeps most of the gain.
  - Depth must be re-selected label-free for the new readout before any confirmatory use.
- **KD response.** The large gain of a3/a3r (+0.06 to +0.08) needs care. This truth is built from held-out knockdowns of the same screen, and 62% of its modules are near-duplicates. Treat it as secondary until the planned near-duplicate cluster bootstrap is run.

Reference code for the recommended readout (a3r; P is the lazy walk operator):
```r
Pt <- P; for (s in 2:9) Pt <- Pt %*% P
L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA               # log attention after t steps; unreachable = NA
mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
Z <- (L - mu) / sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1))   # each gene against its own background
Z[is.na(Z) | Z < 0] <- 0
S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0                             # CLR's Stouffer combination
```

## 5. Conceptual synthesis

**Order of operations, resolution or normalization?** Mainly normalization, with a resolution component at the softmax. Order of operations is not the problem.

- **Normalization (the main issue).** CLR's z-score measures surprise against a gene's own background, and it has no budget: a gene can have 77 partners at z = 10. Softmax attention imposes a budget (rows sum to 1) and a fixed perplexity (10). Its weights are therefore *shares*, not strengths. Ranking shares across genes penalizes every gene whose evidence is spread over many interchangeable partners. The mechanism includes the reverse direction: peripheral genes concentrate their budget on the clique. The pooled metric and the per-complex AP both compare across rows, so both collapse even though the within-row order recovers (section 3.4).
- **Resolution (secondary, at t = 1).** Fixed perplexity 10 and cap 50, against 77 equivalent partners, truncate the clique: 76% of within-complex pairs get no weight. The walk is the repair. Inside a well-connected clique it mixes within a few steps, so Pᵗ restricted to the clique converges to the clique's quasi-stationary distribution. That is a soft lumping of the clique, i.e. a coarse-graining in the sense of Markov-chain lumpability and spectral coarse-graining. This is also why the multi-head and keff questions stay orthogonal and preregistered.
- **Order of operations (the PI's intuition).** The intuition explains *why* rows differ: a complex is one latent factor measured n times, and attention splits its budget among the copies. But the walk already performs the aggregation, and 88% of the mass stays in the complex. What the walk gets wrong is the expansion back to genes: each copy receives 1/n of the factor's affinity. The latent-variable model says each copy should *inherit* it.
  - Hard coarse-graining (b) inherits explicitly and rescues CORUM. It also erases each member's individual relations, needs a clique threshold, and costs the TF regulons.
  - Most of its CORUM gain is reproduced by simply putting within-unit pairs first (b0: 6.76× with no TF cost). The coarse-grained walk adds a little on CORUM (7.62×) and costs the regulons.

**The principled operator: normalize, attend, renormalize, with the second CLR acting on attention logits.** S = CLR(log Pᵗ):

- **What log Pᵗ is.** It is the log-probability of a t-step attention walk from i to j: a multi-hop attention logit. At t = 1 it is the softmax logit z_ij/τ_i minus the row's log-partition function.
- **What row standardization removes.** It removes the row's partition function (its budget) and its temperature, which are exactly the quantities through which mass conservation penalizes genes with many equivalent partners.
- **Why the log domain matters.**
  - In probability space the n clique entries dominate the row's SD, so z ∝ √(G/n) and the dilution survives: CLR on raw S9 fails.
  - In log space the 1/n dilution becomes an additive log n offset, measured against a background of thousands of non-partners. So "inherit, don't divide" is approximately restored without choosing cliques.
- **The bilateral Stouffer combination** keeps CLR's rule that each gene judges the pair against its own background.

So the honest one-line framing would be: *attention conserves mass; CLR is the normalization layer that removes the budget again, on both sides of the attention operator* (CLR → softmax → walk → log → CLR).

Related ideas, to cite after checking:
- CLR itself (Faith et al. 2007) as per-gene background normalization;
- NetMF (Qiu et al., WSDM 2018): DeepWalk implicitly factorizes the *log* of a normalized multi-step walk matrix, so the log is the natural readout of walk affinities;
- token merging with proportional attention (Bolya et al., ICLR 2023): the multiplicity s of merged tokens enters attention additively as log s;
- mutual-rank symmetrization in co-expression networks (Obayashi & Kinoshita 2009) for one-sided rankings.

The diffusion-maps and lift readouts (Pᵗ/π) are the classical alternative. They fail here, as in the earlier hub control, because π is extremely skewed: max/median 18,457 on K562.

## 6. Recommendation

1. **Preregister before any confirmatory use.** Write an addendum to the preregistration and commit it *before* E. coli truths are scored with any new readout.
   - Primary readout: a3r (bilateral CLR on log Pᵗ, row statistics over reachable pairs). Secondary: a2r.
   - Operator unchanged.
   - Depth: t = 9 as transfer, plus a held-out t* re-selected label-free for the new readout (the log readout is degenerate at small t).
   - Endpoints mirror the keff preregistration: M3D and EcoMAC TF regulons (SC, all), K562 CollecTRI, DoRothEA and CORUM; KD response secondary. Raw and hub-stratified, paired jackknife, Holm.
   - Proposed decision rule: adopt as the method's readout only if it is non-inferior on E. coli TF regulons (margin −0.005, both compendia) and superior on K562 CORUM. Otherwise report it as the mechanism of the K562 loss.
2. **Do not touch** the keff and multi-head test. The readout is orthogonal to width, and width or multi-head variants must not be scored on K562 truths before phase 2 (not done here).
3. **Do not adopt coarse-graining.** Keep (b) as an explanatory ablation in the supplement. It shows the latent-variable picture, and it shows that collapsing before attending costs the regulons.
4. **Paper text now** (no new claims): replace "the walk trades pair-level resolution for module recovery" with the measured mechanism.
   - The softmax budget turns strength into share. The walk repairs the within-row order but not the scale.
   - The loss is one module (the cytoplasmic ribosome and its overlapping CORUM entries). Without it, attention ties CLR on CORUM.
5. **Before preregistering**, two label-free checks:
   - the held-out depth curve for a3r on M3D lab splits (no truths);
   - the number of unreachable pairs per dataset at the chosen t.
6. **Housekeeping.** The scripts and outputs are in the session scratchpad: `corum_diag.R`, `corum_followup.R`, `corum_extra.R`, `corum_robust.R` and `corum_posthoc.R`, with logs and CSVs in `corum/`, under `/tmp/claude-1000/-home-boris-Projects-ClrAttention/b33914c9-0cfd-485d-8dd4-cf8a2a22b632/scratchpad/`. /tmp is cleaned at boot, so copy them into `analysis/` before the 16:00 shutdown if they should be kept.

## 7. Questions for the PI

1. **Scope:** Must the method handle protein complexes, or is that a stated boundary? If preregistration confirms CLR(log Pᵗ), should it become the method or a second operator?
2. **Fair truth:** 80% of the CORUM gap is the cytoplasmic ribosome counted about six times. Which complex truth do we pre-specify: pooled, overlap-merged, macro, or RP-excluded?
3. **Latent variables:** Did you mean hard collapse (cliques become nodes) or invariance to how many noisy copies of one factor are measured? The data support the second.
4. **Story:** May we frame the operator as "CLR → attention → CLR(log)", although the readout changes after the keff preregistration?
5. **Replication:** Which E. coli complex truth (EcoCyc complexes, ribosomal-protein operons) do you trust for a preregistered replication, and should knockdown-response modules stay secondary?

Context for each: (1) a yes means rerunning the headline, the keff test and the figures. (2) Without the eight RP-dominated entries: attention 2.05×, CLR 2.18× (difference −0.0016 [−0.0042, +0.0010]). (3) The hard version rescues CORUM (7.62×) but costs DoRothEA A-B 0.013 and needs a clique threshold; the walk already aggregates the clique. (4) The alternative is to keep keff-10 attention as the method and present the normalized readout as the explained failure mode plus a preregistered extension. (5) Without an E. coli complex truth, any complex claim rests on K562 alone; the bilateral readout's +0.06 on knockdown response needs the near-duplicate cluster bootstrap.
