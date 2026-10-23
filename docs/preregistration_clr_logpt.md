# Preregistration: CLR on log attention, CLR(log Pᵗ)

*Fixed 2026-10-22, before the readout was computed on any E. coli or RPE1 data and before any E. coli protein-complex truth was built. Drafted by Claude (Opus, high effort) at the PI's request from the exploratory memo `docs/plans/corum_analysis.md` (section 6 and the PI decisions of 22 Oct), then reviewed for multiplicity, researcher degrees of freedom and leakage. Any change after this commit is logged below with its date and reason; a change made after the first truth-phase run is a deviation and is reported as one.*

## Question

Does re-normalizing the attention walk in the log domain, S = CLR(log Pᵗ), rescue protein complexes (where attention loses to CLR) without giving up the regulon gain that is the paper's headline? If yes, CLR → attention → CLR(log Pᵗ) becomes the main method and plain attention becomes its ablation (PI decision 1). If no, the readout is reported as the explained mechanism of the K562 loss.

## Disclosure: what has already been touched

| Data and truth | Status | Role here |
|---|---|---|
| K562 essential (network half A, 5,671 truth genes): CORUM, CollecTRI, DoRothEA A–B and A–C, knockdown response | About 20 operators and readouts, including this readout, were scored on these truths (exploratory, `analysis/corum_exploratory/`). The readout form (a3r), the RP-dominated-complex diagnosis and the overlap-merge idea all came from them. | **In-sample.** Reported, never confirmatory, never able to rescue an E. coli failure. |
| E. coli M3D and EcoMAC: RegulonDB regulons (SC, all), operons, sigmulons | Scored many times for CLR and attention (the headline). The new readout has never been computed on E. coli data, not even label-free (checked: no log-Pᵗ output exists under `results/`). | **Confirmatory** for the new readout. The comparator scores are known. |
| E. coli protein complexes (EcoCyc) | Never downloaded, built or scored. | **Confirmatory.** |
| RPE1 essential Perturb-seq pseudobulk (`results/perturbseq_bulk_rpe1`): CORUM, CollecTRI, DoRothEA | Scored for \|r\|, MI, CLR variants and attention. The new readout has never been computed on RPE1. The author of this file saw the first lines of that `metrics.csv` (\|r\| and MI rows only). | Key secondary: a human replication in another cell line, but with the same CORUM content (the ribosome again). |
| Depth t = 9, keff 10, cap 50, α = 0.5 | Chosen without labels on M3D (paper, Methods); unchanged by the exploration. | Fixed. |

The keff and multi-head test (`docs/preregistration_keff.md`) is separate and unchanged. Its variants are scored only with the plain readout, and no keff or multi-head variant is scored with this readout in any confirmatory analysis. If the readout is adopted, the paper says that it was found after the keff preregistration, on K562, and confirmed here (PI decision 4).

## Operator (fixed)

**Upstream, unchanged from the method.** B-spline MI (`estimate_mi(bins = "hg", transform = "none")`), Stouffer CLR (`calibrate(method = "normal", combine = "stouffer")`), per-row softmax over the top 50 positive CLR scores at perplexity keff = 10, empty rows given a self-loop, lazy walk P = ½(I + A) (α = 0.5), exactly as `build_fit` and `make_A` / `build_operator` in `analysis/lab_splits.R`.

**Pᵗ.** Dense double precision, `Pt <- P; for (s in 2:t) Pt <- Pt %*% P`. No floor, epsilon or smoothing is added at any step.

**Primary readout R (a3r, bilateral CLR on log Pᵗ, row statistics over reachable pairs):**

1. L = log Pᵗ. Entries with Pᵗ_ij = 0 in double precision (unreachable, or underflowed) and the diagonal are set to missing.
2. For each row i, μ_i and σ_i are the mean and SD (n − 1 denominator) of the non-missing entries of row i.
3. Z_ij = (L_ij − μ_i) / σ_i. Missing, non-finite or negative Z (including rows with fewer than 2 reachable entries or σ_i = 0) are set to 0.
4. S = (Z + Zᵀ) / √2, diagonal 0.

This is `rowz_na()` followed by the Stouffer step in `analysis/corum_exploratory/corum_robust.R` (score `a_clr_bilateral_logP_t9_robust`). The new script must reproduce that script's K562 CORUM AUPR at t = 9 (0.1054351) to 6 significant digits before any E. coli truth is read (implementation check; K562 is already in-sample).

**Secondary readout (descriptive only).** a2r: M = (L + Lᵀ)/2 (missing if either entry is missing), then steps 2–4 on M.

**Mechanism control (descriptive only).** a1: steps 2–4 on raw Pᵗ instead of log Pᵗ (expected to fail; it shows that the log is needed).

**Depth.**

- **Primary: t = 9 on every dataset**, for both the new readout and attention. This is the method's existing transferred depth and the configuration that was explored, so the primary contrast isolates the readout (same Pᵗ, different readout).
- **Secondary: held-out t\*** for the new readout, chosen without labels with the existing criterion A (AP of S for the other half's top-1% CLR pairs) on the keff splits (M3D 20 lab splits, EcoMAC 20 author splits, K562 10 random halves of network half A; RPE1 10 random halves of its half A), grid {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32}, aggregation as in `analysis/keff_test.R` (per split, geometric mean of the two directions; median over splits). Attention's own t\* comes from the same run and must equal the committed K10 selections where they exist (M3D 9, EcoMAC 8). If the new readout's t\* is below 4, it is flagged, because the log readout is degenerate when most pairs are unreachable.
- No truth curve over t is computed for E. coli or RPE1. No other depth is scored.

**Label-free checks reported before truth scoring:** per dataset at t = 9 and at t\*, the share of off-diagonal pairs that are unreachable, the number of rows with fewer than 2 reachable entries, and the distribution of σ_i.

## Comparators

- **Attention (primary comparator):** S = (Pᵗ + Pᵗᵀ)/2 at the same t, diagonal 0, the current method.
- **CLR (key secondary comparator):** the Stouffer CLR matrix the operator is built on.
- Descriptive only: \|r\|, a2r, a1, and the hub-stratified versions below.

## Datasets

| Dataset | Construction | Role |
|---|---|---|
| E. coli M3D | 4,297 genes × 466 replicate-averaged experiments (`analysis/baselines.R`, `--dataset m3d`) | Primary |
| E. coli EcoMAC | 2,198 arrays, ComBat by platform (`baselines.R`, `--dataset ecomac`) | Primary (same organism; may share arrays with M3D, so a robustness replicate, not an independent one) |
| K562 essential | Network half A, 5,671 truth genes (`baselines.R --dataset k562ess --gene_subset truths`, seed 20261015); the exact input of the exploration | In-sample, descriptive |
| RPE1 essential | Pseudobulk, significant knockdowns split into halves A/B as in `analysis/perturbseq_bulk.R` (seed 20261015), network on half A, genes as `perturbseq_bulk.R` (8,747) | Key secondary (human replication) |

## Truths

### E. coli protein complexes (primary complex truth)

- **Source.** EcoCyc protein complexes from the newest EcoCyc release available on the download date: every frame of `proteins.dat` that has COMPONENTS is a complex. If EcoCyc data files cannot be obtained under their licence before scoring, the fallback is the EBI Complex Portal E. coli K-12 set (taxon 83333, `complextab/83333.tsv`) of the current release. The source is chosen on availability alone, before any scoring, and logged here with its release and SHA-256. The unused source may be reported only as a descriptive sensitivity analysis; it never replaces the chosen one. If neither can be obtained, P1 and P2 cannot be run and adoption is impossible.
- **Members.** Subcomplexes are flattened recursively to their polypeptide genes; modified forms resolve to the gene of their unmodified protein; non-protein components (RNAs, metabolites, cofactors) are dropped; genes are mapped to b-numbers with `make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")`. Homomers and any complex with fewer than 3 or more than 200 distinct in-universe genes drop out (the CORUM rule of `corum_sets()`).
- **Overlap merge** (the same rule for every complex truth, fixed without looking at the merged result), applied per dataset after the size filter:
  1. remove exact duplicates and every complex whose in-universe gene set is a proper subset of another retained complex (absorption; this changes no positive pair);
  2. cluster the remaining sets with `hclust(as.dist(1 - Jaccard), "complete")` and `cutree(h = 0.5)`, so every two complexes in a merged cluster have Jaccard ≥ 0.5; a merged complex is the union of its cluster (not re-filtered by size).
- **Pairs.** `comembership_pairs()` with the RegulonDB operon map: the universe is all pairs of genes in at least one merged complex; positives share a merged complex; **same-operon pairs are excluded**, as for the regulon truths (operons are their own truth and are trivially co-expressed).
- **Evaluability** (truth structure only, checked before scoring): the endpoint needs ≥ 200 positive pairs and ≥ 10 merged complexes contributing positives in that dataset's gene universe. If it fails, the endpoint is scored with same-operon pairs included instead; this is decided from the counts alone and logged before scoring.
- **Sanity check (PI): the ribosomal-protein module.** All r-protein genes (symbols rplA–rplY, rpmA–rpmJ, rpsA–rpsV in `ecoli_k12_genes.tsv`) as one module; coherence AUROC of within-module cross-operon pairs against pairs of a member and a non-member (`regulon_coherence()`), descriptive.
- **Sensitivities (descriptive):** same-operon pairs included; minimum size 2 (heterodimers in); no merge; ribosome excluded (merged complexes with ≥ 30% r-protein genes removed, the E. coli analogue of the K562 RP-excluded truth).

### E. coli regulons (no-harm truth)

- **Primary:** the headline truth, RegulonDB regulator-agnostic co-membership at strong/confirmed evidence, same-operon pairs excluded (`baselines.R`, `bench$SC`).
- **Secondary (descriptive):** all evidence (`bench$all`); TF-only regulons, operons and sigmulons with the definitions of `analysis/granularity.R` (as in the keff preregistration).

### Human truths (K562 and RPE1)

- **Overlap-merged CORUM** (PI decision 2): `corum_sets()` (3–200 genes), then the overlap merge above. Scored pooled (all pairs) and per merged complex.
- Note: pooled AUPR counts each pair once, so merging cannot reduce the ribosome's weight in the pooled metric (it comes from the number of pairs, not from duplicate entries). The merge matters for the per-complex metric.
- **Sensitivities:** raw CORUM pooled (the exploratory truth; on K562 this is the reproduction check); RP-excluded CORUM (merged complexes with ≥ 30% RP members, RP = `^RP[LS][0-9]+[A-Z]?[0-9]*$|^RPLP[0-2]$|^RPSA$`, removed).
- CollecTRI and DoRothEA A–B (no-harm), DoRothEA A–C, knockdown-response modules from half B (secondary; 62% near-duplicate modules on K562).

## Endpoints

All differences are pooled co-membership AUPR, new readout minus comparator, at t = 9, unless stated.

**Primary (four tests, all confirmatory, all E. coli):**

| Code | Dataset | Truth | Comparator | Hypothesis |
|---|---|---|---|---|
| P1 | M3D | EcoCyc complexes (merged) | attention | superiority: d > 0 |
| P2 | EcoMAC | EcoCyc complexes (merged) | attention | superiority: d > 0 |
| P3 | M3D | RegulonDB regulons, SC | attention | non-inferiority: d > −0.005 |
| P4 | EcoMAC | RegulonDB regulons, SC | attention | non-inferiority: d > −0.005 |

**Key secondary (interpreted only if the primary rule is met; Holm over five):**

| Code | Dataset | Truth | Comparator | Hypothesis |
|---|---|---|---|---|
| S1 | M3D | regulons, SC | CLR | superiority (the headline survives) |
| S2 | EcoMAC | regulons, SC | CLR | superiority |
| S3 | M3D | EcoCyc complexes | CLR | non-inferiority, margin 0.005 (the gap to CLR is closed) |
| S4 | EcoMAC | EcoCyc complexes | CLR | non-inferiority, margin 0.005 |
| S5 | RPE1 | merged CORUM | attention | superiority |

**Descriptive (intervals, no claims):** everything else, including all K562 results; per-complex AP (`analysis/macro_aupr.R` definitions, complexes with ≥ 10 positive pairs, merged complexes); hub-stratified differences (10 strata of π_i + π_j, π the stationary distribution of the walk, both methods under the same adjustment, as `baselines.R`); attention vs CLR on E. coli complexes (is there a loss to rescue?); held-out t\* results; a2r and a1; the truth sensitivities listed above; RPE1 CollecTRI and DoRothEA A–B.

## Statistics

- **AUPR:** average precision of `pr_summary()` / `ap_fast()` in `analysis/regulondb.R` (tie blocks counted at their end; non-finite scores rank last).
- **Primary inference:** paired delete-half gene jackknife, 100 replicates, as `baselines.R`: in each replicate a random half of the truth's genes (`genes_in`) is kept and AUPR is computed over pairs with both genes kept; the same halves serve every method; se = SD of the 100 paired differences. Halves are drawn per truth, in the order of the endpoint tables, after `set.seed(20261022)`; one seed, never rerun with another.
- **p-values (one-sided z):** superiority p = 1 − Φ(d / se); non-inferiority p = 1 − Φ((d + 0.005) / se).
- **Primary rule as an intersection–union test.** Adoption requires all four primary nulls to be rejected, each at one-sided α = 0.025. Because the rule rejects either all four or none, the familywise error of the adoption claim is at most 0.025 without adjustment (Berger 1982). Any claim about a single primary endpoint outside the adoption rule uses Holm-adjusted p over the four primary tests (one-sided familywise α = 0.025).
- **Key secondary:** Holm over S1–S5 at one-sided familywise α = 0.025, tested only if adoption holds (gatekeeping, so no α leaks from the primary family).
- **Experiment-level robustness (M3D):** paired cluster bootstrap as the cluster mode of `analysis/chips_bootstrap.R` (the 466 experiments drawn with replacement, one replicate array of the 907 per occurrence; draws b = 1–200 with seeds 20260924 + b), MI, CLR, attention and the readout refitted in each draw at t = 9; 2.5–97.5% percentile intervals of the paired differences for P1 and P3. EcoMAC has no such bootstrap.
- **Smallest effect of interest:** 0.005 AUPR, as in the keff preregistration.

## Decision rules

**Adopt** CLR(log Pᵗ) as the main method (attention reported as the ablation) only if all of these hold:

1. P1 and P2: one-sided p < 0.025 and point estimate ≥ +0.005;
2. P3 and P4: non-inferiority p < 0.025;
3. M3D cluster bootstrap: P1 interval lower bound > 0 and P3 interval lower bound > −0.005;
4. the implementation check (K562 reproduction) passed before E. coli scoring.

**Anything else is a failure to confirm**, and the main method stays attention:

| Outcome | Conclusion and reporting |
|---|---|
| Complexes improve, regulon NI not shown on M3D or EcoMAC | Not adopted. Reported as a complex-specific readout with an unresolved regulon cost. |
| Regulon loss on M3D or EcoMAC: d ≤ −0.005 and Holm-adjusted one-sided p (for d < 0) < 0.025 | Not adopted. Reported as a trade-off: the readout harms regulons. |
| Complex superiority not shown on M3D or EcoMAC | Not adopted. The K562 rescue is reported as exploratory and not replicated in E. coli. If attention does not lose to CLR on E. coli complexes (descriptive), this is said too (no loss to rescue), but it does not change the verdict. |
| Gene jackknife passes, cluster bootstrap fails (rule 3) | Not adopted. Reported as supported over genes, not robust to resampling experiments. |
| Implementation check fails | No E. coli scoring until it is fixed; the fix and the cause are logged here first. |

After adoption:

- if S1 or S2 fails, the paper does not claim that the new method beats CLR on regulons; the headline over CLR stays attributed to attention;
- if S5 fails, the paper says the human complex result did not replicate in RPE1, and the complex claim rests on E. coli plus in-sample K562;
- a small but significant regulon loss within the margin (NI shown, d < 0 significant) is reported as such;
- if the K562 merged-CORUM difference vs attention is ≤ 0 (descriptive, in-sample), this is reported prominently.

In every case the paper keeps the measured mechanism (softmax budget turns strength into share; the walk repairs the within-row order, not the scale) as an explanation of the K562 loss.

What does **not** count: a win on K562 alone; a win by a2r or a1 when a3r fails; wins at any depth other than t = 9 (including the committed t\*, which is secondary); wins on a sensitivity truth when the primary truth fails; a win on per-complex AP when pooled fails; a win after swapping the complex source, size rule, merge rule or operon rule; results within ±0.005.

## Protocol and gate

1. Commit this file. The script (`analysis/clr_logpt_test.R`, to be written) refuses to read any truth unless `docs/preregistration_clr_logpt.md` is tracked in git and identical to HEAD (`git diff --quiet HEAD -- docs/preregistration_clr_logpt.md`); it writes the hash of the last commit that touched this file into every output.
2. Implementation check on K562 (in-sample): reads only the K562 raw CORUM truth and must reproduce the exploratory AUPR above. Smoke tests use `--quick` gene subsets and never read E. coli or RPE1 truths. The E. coli and RPE1 truth phases refuse to run unless the check's output exists and passed.
3. Label-free phase: reachability checks and held-out t\* for every dataset; selections committed to `docs/clr_logpt_selections/`. The truth phase refuses to run unless its selections are byte-identical to the committed copies (as the keff gate).
4. Build the E. coli complex truth; record source, release, SHA-256 of the raw file and of the built member list, and the evaluability counts (no scores) in the change log, and commit. The truth phase refuses to run if the member list's SHA-256 differs.
5. Code review of the script against this file (a reviewer who has not seen any E. coli result), then the truth phase, once: E. coli first, then RPE1, then K562. The first complete run is the confirmatory one. A rerun is allowed only after a crash, or for a bug that is shown to violate this file (the line and the rule named); both runs, both verdicts and the reason are reported.
6. The verdict is computed by the script from the rules above and written to `verdict.csv` before anyone interprets the numbers.

## Power (rough)

The paired SE of the new readout minus attention on K562 was 0.0013–0.0015 (TF truths), about 60% of the attention–CLR SE. On M3D the attention–CLR SE is 0.003, so the paired SE is likely about 0.002. With a true regulon difference of 0, non-inferiority is shown with probability Φ(0.005/0.002 − 1.96) ≈ 0.71 per compendium; with +0.003 (the K562 value), ≈ 0.98. Both compendia must pass, so at a true difference of 0 the joint probability is between about 0.5 (independent) and 0.71 (fully correlated). The regulon arms are the main risk of a false "not confirmed"; the margin is nevertheless kept at the keff preregistration's smallest effect of interest rather than widened to reach power.

## Compute

E. coli and K562: a dense G² walk to t = 9 takes seconds to minutes; jackknife and held-out selection about 1–2 h per dataset at 4 workers. M3D cluster bootstrap: 200 refits, about 3 h at 4 workers. RPE1 (8,747 genes): a few dense 8,747² matrices, about 5 GB. Every job goes through the scheduler (`jobs_pc/sched`, with memory and CPU headers).

## Change log

- 2026-10-22: fixed.
