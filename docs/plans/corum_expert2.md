# CORUM loss and CLR-NS z-scores: second expert memo

*3 November 2026. Expert memo for the PI (network inference, dependence statistics, diffusion on graphs, attention; Claude Opus). Everything computed here is **EXPLORATORY**. K562 CORUM, DoRothEA and CollecTRI are in-sample: they were already used in `docs/plans/corum_analysis.md`, the keff test, the CLR(log Pᵗ) test and `analysis/rclr_k562/`. BEELINE truths are in-sample for CLR-NS, which was chosen on them. No held-out or sealed truth was opened: no DREAM5 gold standard, no yeast truth, no E. coli complex truth, no RPE1 file, no scFM K562 truth. Labels: **[V]** I re-derived or re-read it; **[R]** recalled from a cited doc and not re-derived; **[S]** speculation.*

> **QC (3 Nov, `docs/plans/corum_expert2_qc.md`): PASS WITH CORRECTIONS.** All numbers reproduce. Read the QC discrepancy table before quoting: (1) the 'real rows are near-symmetric' statement holds for |r| only (B-spline MI row skewness is 1.19 in BEELINE vs 0.70 in the null); (2) 'ranking is not affected' is true only for the ties artefact, row calibration does change the pooled top-1% (Jaccard 0.75 for median/MAD); (3) 'MI is more skewed under the null' holds for K562 only; (4) the false 'standard-normal null' sentence also stood in `methods.md` and `R/calibrate.R` (all three fixed 3 Nov).

## Verdicts

**Task A: no untested reason for the CORUM loss remains that attention can address and that would survive a held-out test.**

- The cap does not bind: cap 100 or 200 at keff 10 changes CORUM AUPR by −0.001 and −0.002.
- Widening past the keff grid lifts raw CORUM above CLR (uniform attention over the top 50 partners: 7.90× vs 6.39×), but the whole gain is the ribosome: merged CORUM without RP complexes changes by −0.0001, and DoRothEA (−0.0035) and small complexes lose.
- CLR's CORUM lead depends on comparing raw scores across rows. Reading CLR by within-row rank cuts its CORUM AUPR from 0.0847 to 0.0445; with rank readouts for both, attention is at or above CLR on all seven K562 truths (point estimates).
- Outside the ribosome a broad but tiny per-complex deficit in 3–9-gene complexes persists under every width and readout. Its mechanism is known (the keff-10 perplexity floor makes small cliques leak); the only candidate fix is reciprocity, which has a held-out design.

**Task B: the concern is right about the paper's wording and wrong about the method.**

- The Methods claim that Stouffer's combination "has a standard-normal null" is false: under a within-gene permutation null, S > 3 occurs in 3.8% of BEELINE pairs instead of 0.13%, because z_ij and z_ji are essentially the same statistic (correlation 0.998), both are clipped, and the null |r| row is half-normal.
- Real rows are not shaped like that null: they are near-symmetric distributions of background dependence (median skewness 0.25) with an SD 3.3× the null's. Real CLR tails are therefore *thinner* than permutation tails in 14 of 14 cases, so a pooled permutation null is very conservative (the known O2 failure, `background_width_exploration.md` §1).
- Ranking is not affected: ties keep the null variance at exactly 1/(n−1), and per-gene tails do not depend on sparsity. None of eight label-free calibrations improves CLR-NS on BEELINE (exploratory; Fisher z is the same ranking, median/MAD +0.021 log2 in 3 of 7 datasets).
- Recommendation: correct the wording, keep the method; no held-out test is warranted.

## 1. What was run

| Job | Script | Output | Compute |
|---|---|---|---|
| 1100 | `analysis/corum_expert2/corum_width_readout.R` (helpers `null_diag.R`) | `results/corum_expert2/k562/` | 34 min, peak 17.7 GB, 8 threads, nice 19 |
| 1101 | `analysis/corum_expert2/clrns_null_beeline.R` | `results/corum_expert2/beeline/` | 9 min, 4 threads, nice 19 |

- K562: half A, 5,671 truth genes × 707 profiles; intervals d ± 1.96 SD over a paired delete-half gene jackknife (100 reps, seed 20261016), networks fixed. BEELINE: 14 cases; stage 1 (label-free) was written before stage 2 read any truth.
- **Disclosure.** Operator and readout lists were fixed before the full runs. A 600-gene smoke test of the K562 script did score truths (in-sample). After it, the degree for the degree-scaled readout was changed from "permutation-BH partners" (empty set) to "1 + top-1% CLR partners". The two BEELINE input variants (detection-residualized, PC1 removed) and the background diagnostics were added after the label-free mDC smoke test showed 99% positive r; no BEELINE truth had been read.
- **Reproduction [V].** K562: CLR 0.08471, K10 0.04654, a3r 0.10544 on CORUM, hub-stratified CLR 0.0451 and attention 0.0275, all as documented. BEELINE: max |Δ| AUPRC ratio vs `results/beeline_baselines/metrics.csv` 4e-14 (CLR-NS), 5e-14 (CLR).

## 2. Task A: candidates

Each candidate has a mechanism, an existing result (with source) and a verdict. "New" marks numbers from job 1100 at t = 9.

| # | Candidate mechanism | Already tested? | New (exploratory) | Verdict | Cheapest decisive test |
|---|---|---|---|---|---|
| A1 | **Cap 50 truncates an ~80-gene clique** | No: cap was 50 in every run | K10 cap 100: −0.0010 [−0.0017, −0.0004]; cap 200: −0.0017 on CORUM raw. TF truths ±0.0004 | **Ruled out** (the extra candidates receive negligible weight at keff 10 [S]) | done |
| A2 | **keff 10 cannot cover the clique (width)** | Yes. keff test on K562: K25 +0.020, O3A (keff 40) +0.039, O1t (near-uniform) +0.047 on CORUM (`results/keff_test/k562ess/metrics.csv`); rule failed (needs 3/4 endpoints and both E. coli; `docs/keff_verdict/verdict.csv`) [V]. Mass-conservation bound m/(n−1) for any width (`asymmetry_expert.md` §3.1) [R] | Beyond the grid: K80c100 7.36×, U50 7.90× (+0.058 [+0.037, +0.080] vs K10; +0.020 [+0.012, +0.028] vs CLR). Merged noRP: −0.001 / −0.000; small 3–6: −0.0012; DoRothEA −0.0039 / −0.0035; per-complex AP for 3–9 genes 0.0015 vs K10 0.0024 | **Explained; not a fix.** Width spreads every row thin. That makes shares comparable across rows, which matters mainly for genes with many strong partners (the ribosome) [S], and costs regulons and small complexes. This is the trade-off the keff rule excludes. | none needed |
| A3 | **Mass leaks into the wider translation program** | Yes: 88% of RP-row mass stays in the complex at t = 9 (`corum_analysis.md` §3.3) [R] | – | **Ruled out**: dilution, not leakage | – |
| A4 | **Depth vs clique diameter** | Yes: diameter 1; granularity maximum 3.9× at t ≈ 17 (`corum_analysis.md` §3.2) [R] | K10 3.51× (t 9) → 3.81× (t 16). Wide operators peak earlier: U100 7.39× at t 4 | **Ruled out** as a lever; it interacts with width | – |
| A5 | **Hubs and one-sided attention in (P + Pᵀ)/2** | Partly: hub-stratified gap persists; INTS8/EIF4B one-sided negatives (`corum_analysis.md` §3.3) [R] | min(P, Pᵀ): +0.0043 [−0.0015, +0.0102] (ns). Degree-scaled (dᵢPᵢⱼ + dⱼPⱼᵢ)/2, with d = 1 + top-1% CLR partners: +0.0094 [+0.0027, +0.0161]. Both stay far below CLR (−0.034, −0.029) | **Minor**: one-sidedness is real but second-order | – |
| A6 | **The pooled metric rewards cross-row calibration that diffusion destroys** | Partly: row-local AP 0.47 vs 0.58; per-complex AP 60S 0.19 vs 0.66 (`corum_analysis.md` §3.4) [R]; coherence CLR 0.574 vs attention 0.538 (`results/baselines_k562ess_truthgenes_genie3/metrics.csv`) [V] | **Row-percentile readout.** CLR 6.39× → 3.36× (−0.040 [−0.058, −0.022]). K10 3.51× → 5.13× (+0.0215 [+0.0115, +0.0315]). K10-rank vs CLR-rank: raw 5.13 vs 3.36; merged noRP 2.26 vs 1.90; small 2.01 vs 1.86; large 2.29 vs 2.18; DoRothEA 1.13 vs 1.08; CollecTRI 1.15 vs 1.04 (point estimates). Per-complex, CLR beats K10 in 74% of complexes (72% of 3–9-gene ones), and CLR-rank still leads K10-rank there (0.0037 vs 0.0026) | **Partly confirmed.** CLR's pooled CORUM lead is a cross-row scale effect, and under a scale-free readout the ordering reverses. The rank readout is a member of the a3r family (row renormalization after attention), which failed confirmation, so it is not a new mechanism. | see §4, item 3 |
| A7 | **Evaluation artefacts (redundancy, overlap, ribosome dominance)** | Yes. Merged CORUM: attention 0.0464 vs CLR 0.0843. Merged noRP: −0.0015 (SE 0.0015) (`results/clr_logpt_test/truth/k562ess/metrics.csv`) [V] | Same values recomputed. Small 3–6: −0.0005 [−0.0038, +0.0028]. Large ≥ 20: −0.046 [−0.066, −0.026] | **Settled**: the pooled loss is the ribosome. Merging does not remove it, because merged RP complexes are still one large clique. | – |
| A8 | **CLR diagonal handling and clipping** | Diagonal 0 is kept in μ and σ (1 of G entries; `src/core/clr_packed.cpp`) [V]. Self-loops have no effect (`self_loop_trace.md`) [R] | 59% of MI row entries are clipped (median), so a typical row keeps ~2,300 positive entries, far above the top-50 candidate set [V, derived] | **Ruled out by mechanism** | – |
| A9 | **Small-complex leakage (perplexity floor) and partner-less junk rows** | Synthetic only: reciprocity raises small complexes ×1.86–6.24 (`background_width_exploration.md` §8) [R]; held-out design exists (`reciprocity_heldout_design.md`) [R] | Per-complex AP, 3–9 genes: CLR 0.0033, K10 0.0024, best attention readout (min) 0.0027; wider is worse (U50 0.0015) | **Open.** It is the only remaining lever, and only for small complexes | reciprocity on K562 CORUM 3–6 (in-sample, ~10 min), then held-out (§4) |

**Interpretation.**
- Widening equalizes shares only for genes whose strong partners fill the candidate set. Row renormalization (a3r, rank) restores comparability across rows after the walk. Neither touches what is left outside the ribosome: the floor that forces small cliques to leak.
- "Spreading or truncation?" Neither. Truncation does not bind (A1) and spreading into the translation program is small (A3). The loss is mass dilution inside the clique combined with comparing shares across rows.

## 3. Task B: CLR-NS and Gaussian z-scores

### 3.1 What the row z-score is used for

| Use | Needs a Gaussian row? | Assessment |
|---|---|---|
| Within-row ranking | No | Any increasing per-row map preserves it [V] |
| Cross-row comparison (pooled ranking) | No, but needs comparable location and scale | Real rows share no null shape. Location and scale are set by background dependence (§3.3) [V] |
| Softmax logits | No | The per-gene temperature absorbs the row's own location and scale. Only the partner term z_ji and the clipping carry the partner's background [V, from `.softmax_rows`] |
| "Stouffer … has a standard-normal null" (paper.md Methods) | Yes, plus independence of z_ij and z_ji | **False.** Under the permutation null, corr(z_ij, z_ji) is 0.998 (|r|) and 0.971 (MI), and both terms are clipped. Null S has mean 0.57 and SD 0.98. P(S > 3) is 0.0377 (BEELINE median, range 0.0371–0.0379) instead of 0.00135 [V] |
| Permutation null for FDR selection | No | Valid as a test of complete independence, but very conservative for real CLR scores (§3.3). K562: with one permutation replicate, 101 of 16.1 million pairs pass BH q = 0.05 (resolution-limited: the smallest attainable p is 1/(1+M)) [V] |

### 3.2 Null distributions

- **|r| of Blom scores.** The permutation variance of r is exactly 1/(n−1) for any fixed vectors, tied zeros included (measured (n−1)·mean(r²) 1.000, range 0.996–1.003) [V]. The null row is half-normal: theory gives skewness 0.995, excess kurtosis 0.869, clipped fraction 0.575 and P(z > 3) = 0.0092, 6.8× the Gaussian value; measured 0.98, 0.77 and 0.0090 [V]. Tied zeros do not make per-gene tails heterogeneous: Spearman(null exceedance, fraction of zeros) lies in −0.08 to +0.06, and the q10–q90 spread (0.0058–0.0120) is binomial noise at ~1,000 partners [V].
- **B-spline MI** is more skewed under the null (skewness 0.71 BEELINE, 1.25 K562; P(z > 3) 0.0074 and 0.0137) [V].
- **The premise holds for the null, not for the data.** Real |r| rows have median skewness 0.25 (−0.16 to 0.80) and excess kurtosis −0.34 (−0.97 to 0.58) [V]: a row is the distribution of a gene's dependence on ~1,000 partners, mostly background structure. In mDC, 99.4% of signed r are positive (median 0.32), and gene scores correlate with log detected-gene count (median |ρ| 0.49) [V].

### 3.3 Contaminated null scale

- Real row mean / null mean: 5.4× (2.7–8.3) for |r|; real row SD / null SD: 3.3× (2.3–5.6) for |r| and 7.6× for MI (K562: 3.6× and 4.3×) [V].
- So observed CLR scores are compressed relative to permuted ones: real P(S > 3) < null in 14/14 cases for both measures (BEELINE |r| 0.0092 vs 0.0377; K562 MI 0.0259 vs 0.0353) [V]. The compression is CLR's intended hub penalty, but it makes the permutation null the wrong reference for calibration, as diagnosed for O2 (`background_width_exploration.md` §1) [R]. An Efron empirical null is an open TODO in `select_threshold` (`R/clr_attention.R`) [V].

### 3.4 The original 'rayleigh' method (`data/CLR versions/CLR_1.1.zip`, `clr.m`) [V]

Per row a Rayleigh MLE (`raylfit`) and its CDF; directions combined as A + Aᵀ − A∘Aᵀ. The Rayleigh is χ₂, while null |r| is χ₁ (half-normal) and MI is χ²-like; it fits scale but no location, so the background shift stays; and the CDF saturates (top 1% of Rayleigh-MI pairs: median 1.6%, up to 5.8%, tied at the maximum). Exploratory: Rayleigh-MI −0.057 vs normal CLR-MI −0.062 log2 relative to CLR-NS; Rayleigh-|r| −0.035.

### 3.5 Label-free alternatives (BEELINE; exploratory and in-sample)

Mean over 7 datasets of the per-dataset mean log2(AUPRC ratio / CLR-NS) on the primary truths (44 cells). Ranking similarity is the median over 14 cases.

| Alternative | Spearman vs CLR-NS | Top-1% Jaccard | log2 vs CLR-NS | Datasets better | Exact Wilcoxon p |
|---|---|---|---|---|---|
| Fisher z, then \|·\|, then CLR | 0.999 | 0.93 | +0.0004 | 4/7 | 0.81 |
| Row median/MAD | 0.992 | 0.75 | +0.021 | 3/7 | 1.00 |
| Half-normal row fit (χ₁ analogue of Rayleigh) | 0.974 | 0.67 | −0.004 | 1/7 | 0.30 |
| Rayleigh on \|r\| (clr.m) | 0.964 | 0.61 | −0.035 | 1/7 | 0.22 |
| Per-row rank to normal | 0.991 | 0.70 | −0.129 | 1/7 | 0.078 |
| Efron central matching on signed r | 0.619 | 0.29 | −0.420 | 2/7 | 0.078 |
| Blom residualized on log detection, then CLR-NS | 0.939 | 0.76 | −0.005 | 2/7 | 0.81 |
| First PC removed, then CLR-NS | 0.166 | 0.36 | −0.083 | 1/7 | 0.031 |

Fisher z is the same ranking (monotone, small r); median/MAD is a coin flip; Efron's two-sided null on signed r flags pairs *less* correlated than a shifted background; removing the global factor does not help (CLR's row standardization absorbs most of it; detection residualization helps only mDC +0.039 and mHSC-E +0.051). Nothing justifies a preregistered test.

**Wording fixes (paper.md):**
- Methods, "Context normalization": replace "which has a standard-normal null" with: "a ranking statistic: the two z-scores are functions of the same dependence value and are clipped at zero, so it has no closed-form null; significance, where used, is assessed by permutation."
- Discussion, "each CLR score is a standardized deviation … with a permutation null": add that this null (complete independence) is conservative for CLR scores, because real rows are dominated by background dependence.

**Not tested here:** reciprocity on real K562 (A9); a3r or rank readouts at a label-free t*; more than one permutation replicate for the K562 BH count.

## 4. What to test next, ranked

| Rank | Test | Cost | Confirmatory dataset | Note |
|---|---|---|---|---|
| 1 | Paper wording (§3.5); no method change for Task B | 30 min | none | Not a test; it removes a false claim |
| 2 | Reciprocal attention, small-complex arm (A9): in-sample K562 sanity check on 3–6-gene merged CORUM and per-complex AP, then the held-out design | 10 min in-sample; held-out per `reciprocity_heldout_design.md` | Yeast Complex Portal (taxon 559292) complexes of 3–6 genes, on an untouched yeast compendium (DREAM5 net4 expression, 5,950 × 536 chips, is a candidate [R]; DREAM5's own gold standard is TF–target, not complexes) [S] | The only remaining complex lever. One preregistration must cover all complex readouts, so that the yeast complex truth is touched once. The PI's 31 Oct rule defers it until after the Collins readout |
| 3 | Rank readout as a robustness twin of a3r (A6) | ~1 h label-free t* on M3D lab splits; no truth | Same yeast complex truth, in the same preregistration as item 2. E. coli complexes and RPE1 CORUM are consumed | Serial-testing risk (the third complex readout). Only worth it if the paper wants the "scale-free readout" framing |
| 4 | Width beyond keff 40 / uniform top-k (A2) | none | – | Do not pursue. Its gain is the ribosome, at a regulon cost |
| 5 | Efron empirical null for CLR FDR (Task B) | 1 day | Synthetic plus permutation checks only (calibration is label-free) | Only if the paper keeps an FDR-selection claim |

## 5. Verification table

| Claim | Source file or script | Value | Mark |
|---|---|---|---|
| CLR / K10 / a3r CORUM raw (5,671 genes) | `results/corum_expert2/k562/repro_check.csv` | 0.08471 / 0.04654 / 0.10544 | V |
| K10 cap 100, cap 200 vs K10 (CORUM raw) | `results/corum_expert2/k562/t9_jackknife_EXPLORATORY.csv` | −0.0010 [−0.0017, −0.0004]; −0.0017 [−0.0027, −0.0006] | V |
| U50 CORUM raw enrichment; vs K10; vs CLR | same | 7.90×; +0.0582 [+0.0365, +0.0799]; +0.0200 [+0.0123, +0.0278] | V |
| K80c100 CORUM raw enrichment | same | 7.36× | V |
| U50 / K80c100 vs K10: merged noRP | same | −0.0001 / −0.0010 | V |
| U50 / K80c100 vs K10: DoRothEA A-B | same | −0.0035 / −0.0039 | V |
| U50 / K80c100 vs K10: small 3–6 | same | −0.0012 / −0.0012 | V |
| min readout vs K10 (CORUM raw) | same | +0.0043 [−0.0015, +0.0102] | V |
| degree-scaled readout vs K10 (CORUM raw) | same | +0.0094 [+0.0027, +0.0161] | V |
| K10 row-percentile vs K10 (CORUM raw) | same | +0.0215 [+0.0115, +0.0315] | V |
| CLR row-percentile vs CLR (CORUM raw) | same | 0.0445 vs 0.0847; −0.0402 [−0.0582, −0.0222] | V |
| Rank-readout enrichments, K10 vs CLR (7 truths) | same | 5.13/3.36, 2.26/1.90, 2.01/1.86, 2.29/2.18, 1.13/1.08, 1.15/1.04 | V |
| K10 vs CLR: small 3–6; large ≥ 20 | same | −0.0005 [−0.0038, +0.0028]; −0.0460 [−0.0661, −0.0259] | V |
| U100 at t = 4; K10 at t = 16 | `results/corum_expert2/k562/pooled_all_t.csv` | 7.39×; 3.81× | V |
| Per-complex AP 3–9 genes, CLR-rank / K10-rank | `results/corum_expert2/k562/t9_macro_summary_EXPLORATORY.csv` | 0.0037 / 0.0026 | V |
| Per-complex AP, 3–9 genes: CLR / K10 / min / U50 | `results/corum_expert2/k562/t9_macro_summary_EXPLORATORY.csv` | 0.0033 / 0.0024 / 0.0027 / 0.0015 | V |
| CLR beats K10 per complex: all / 3–9 genes | same (frac_beats_K10, CLR rows) | 0.737 / 0.719 | V |
| Hub-stratified CLR / attention (repro) | `results/corum_expert2/k562/t9_hubstrat_EXPLORATORY.csv` | 0.0451 / 0.0275 | V |
| Coherence CLR / attention, K562 CORUM | `results/baselines_k562ess_truthgenes_genie3/metrics.csv` | 0.574 / 0.538 | V |
| Merged CORUM attention / CLR; noRP difference | `results/clr_logpt_test/truth/k562ess/metrics.csv` | 0.0464 / 0.0843; −0.0015 (SE 0.0015) | V |
| keff test K562 CORUM wins: K25, O3A, O1t | `results/keff_test/k562ess/metrics.csv`, `docs/keff_verdict/primary_holm.csv` | +0.020, +0.039, +0.047 | V |
| 88% of RP mass stays in the complex | `docs/plans/corum_analysis.md` §3.3 | 88% | R |
| Reciprocity small-complex synthetic gain | `docs/plans/background_width_exploration.md` §8 | ×1.86–6.24 | R |
| BEELINE reproduction, CLR-NS / CLR-MI | `results/corum_expert2/beeline/repro_check.csv` | max \|Δ\| 4.3e-14 / 5.0e-14 | V |
| Null corr(z_ij, z_ji): \|r\| / MI (BEELINE median) | `results/corum_expert2/beeline/diagnostics.csv` | 0.998 / 0.971 | V |
| Null P(S > 3), BEELINE \|r\| | same | median 0.0377 (0.0371–0.0379) vs 0.00135 | V |
| Null S mean / SD (\|r\|) | same | 0.57 / 0.977 | V |
| Null P(z > 3): \|r\| / MI; half-normal theory | same; R `2*pnorm(-(sqrt(2/pi)+3*sqrt(1-2/pi)))` | 0.0090 / 0.0074; 0.00915 | V |
| (n−1)·mean null r² | same | 1.000 (0.996–1.003) | V |
| Null exceedance vs sparsity, Spearman range | same | −0.078 to +0.059 | V |
| Real \|r\| row skewness / excess kurtosis (median, range) | same | 0.25 (−0.16, 0.80) / −0.34 (−0.97, 0.58) | V |
| SD inflation, \|r\| / MI (BEELINE median) | same | 3.3× (2.3–5.6) / 7.6× | V |
| Real < null P(S > 3), cases | same (computed in-session) | 14/14 for both measures | V |
| Real vs null P(S > 3), BEELINE \|r\| | same | 0.0092 vs 0.0377 | V |
| mDC: fraction of r > 0; median r; detection \|ρ\| | same (row mDC_500) | 0.994; 0.315; 0.489 | V |
| K562 MI: null P(z > 3); real vs null P(S > 3); SD inflation MI / \|r\| | `results/corum_expert2/k562/taskB_k562_diagnostics.csv` | 0.0137; 0.0259 vs 0.0353; 4.3× / 3.6× | V |
| K562 BH q = 0.05 pairs (one replicate) | `results/corum_expert2/k562/taskB_k562_perm_bh.csv` | 101 | V |
| Rayleigh-MI top-1% ties at maximum | `results/corum_expert2/beeline/rank_change.csv` | median 1.6%, max 5.8% | V |
| Alternatives, log2 vs CLR-NS, datasets better, Wilcoxon p (§3.5) | `results/corum_expert2/beeline/summary_EXPLORATORY.csv` | as tabulated | V |
| Alternatives, Spearman and top-1% Jaccard vs CLR-NS (§3.5) | `results/corum_expert2/beeline/rank_change.csv` (median over 14) | as tabulated | V |
| Per-dataset detection-residualized gains (mDC, mHSC-E) | `results/corum_expert2/beeline/metrics_EXPLORATORY.csv` (per-dataset mean log2) | +0.039, +0.051 | V |
| clr.m 'rayleigh' = raylfit, raylcdf, A + A′ − A∘A′ | `data/CLR versions/CLR_1.1.zip` → `clr.m` | read | V |
| Efron TODO in `select_threshold` | `R/clr_attention.R` | present | V |
| Yeast Complex Portal plus DREAM5 net4 as a suitable held-out complex test | – | – | S |
