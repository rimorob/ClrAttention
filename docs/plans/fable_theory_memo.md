# Theory memo: what CLR, aCLR and tcaCLR are, and what can honestly be claimed about them

*Reviewer model: Claude Fable 5.1 (model id `claude-fable-5-1`). Independent review for the PI, 16 November 2026. Read-only: no compute was run, no truth was opened, no existing file was edited. Sources read: `R/clr_attention.R` (`build_operator`, `.softmax_rows`, `.topk_rows`), `R/calibrate.R` (header), `docs/jim_short_preread_src/jim_short_preread.md` (appendix E pseudocode), `docs/plans/corum_analysis.md` (variants table and a3r reference code), `docs/plans/corum_expert2.md`, `docs/plans/asymmetry_expert.md` (section 3), `docs/plans/explore_stack_plan.md`, `docs/explore_stack_results.md`, `results/explore_stack/{k562,m3d}/diagnostics_layers.csv`, `docs/plans/clr_logpt_interpretation.md`, `docs/STATUS.md`.*

**Tags.** Every claim below carries one of: **[derivable]** (follows from the definitions by a short argument given or sketched here), **[plausible]** (a mechanism consistent with the definitions and the project's measurements, not derived), **[speculation]** (a conjecture), **[measured]** (a number from the project's own files, cited). References carry **[confident]** (I am sure the work exists as cited) or **[verify]** (I believe it exists but the details should be checked before citing).

**Notation.** G genes, N samples. D is the MI matrix. Row statistics μ_i, σ_i are the mean and SD of row i of D (the diagonal 0 included, per `corum_expert2.md` A8). z_ij = (D_ij − μ_i)/σ_i, clipped at 0. S_ij = (z_ij + z_ji)/√2 is the CLR pair score. C_i is gene i's candidate set (up to 50 largest positive S_ij). A is the row-stochastic attention matrix, A_ij = exp(S_ij/τ_i)/Σ_{k∈C_i} exp(S_ik/τ_i), with τ_i set so that exp(H(A_i·)) = 10. P = (I + A)/2. M = P^t, t = 9. aCLR readout: (M + Mᵀ)/2. tcaCLR readout: L = log M (unreachable = NA, diagonal NA), Z = row z-score of L over reachable entries, clipped at 0, S' = (Z + Zᵀ)/√2.

---

## Q1. CLR: what the per-gene z-score estimates, and what guarantees exist

### 1.1 What CLR is, formally

**[derivable]** CLR is a per-gene empirical-null studentisation of a dependence estimate. For gene i, the row of D is treated as a sample from gene i's own distribution of dependence values, and z_ij measures how far D_ij sits above that distribution's location, in units of its spread. This is exactly the "z-value with an estimated empirical null" of Efron (2004) [confident], applied with the gene as the unit of analysis, with the null location and scale estimated by the raw row mean and SD (not by Efron's central matching), and without the later density-estimation step that produces a local false-discovery rate. CLR therefore stops at the z-value stage of Efron's pipeline; it never produces a calibrated probability.

**[derivable]** Under the sparse-signal assumption (row i contains a fraction π_i of true partners with mean excess Δ_i over background, the rest i.i.d. with mean μ_i⁰ and SD σ_i⁰), μ̂_i = μ_i⁰ + π_i Δ_i + O(G^{-1/2}) and σ̂_i² = σ_i⁰² + π_i(1 − π_i)Δ_i² + O(G^{-1/2}). So z_ij is a consistent estimator of the background-standardised score (D_ij − μ_i⁰)/σ_i⁰ as G → ∞ and π_i → 0, with location bias π_i Δ_i/σ_i⁰ and relative scale bias of order π_i Δ_i²/σ_i⁰² otherwise. That is the whole asymptotic content of CLR: consistency of a *gene-specific* background standardisation under sparsity. There is no consistency statement about recovering a regulatory graph, and none should be claimed.

**[derivable]** The guarantee worth stating is an invariance: CLR's row z-scores are invariant to any affine distortion of row i of D, D_ij → a_i + b_i D_ij with b_i > 0. The finite-sample bias of histogram and B-spline MI estimators depends on the marginal of each gene (bin occupancy, entropy; Daub et al. 2004 [confident], Steuer et al. 2002 [confident]) and on N, and to first order such a bias is gene-specific and additive or multiplicative in the row. CLR removes exactly that class of distortion. This is the principled reason the 2007 method standardises per gene rather than pooling, and it is a statement about the estimator, not about biology.

### 1.2 What CLR does not do

**[derivable]** CLR corrects the query's baseline, not the key's. If the dependence matrix carries additive gene effects, D_ij = h_i + h_j + ε_ij (a "promiscuity" or hub effect), then row standardisation removes h_i but leaves z_ij ∝ (h_j − h̄)/σ_i. The bilateral Stouffer score is then S_ij ∝ (h_i − h̄) + (h_j − h̄): the two hub effects are *reintroduced* by the combination, one from each direction. The Euclidean combination does the same. So CLR is not hub-neutral in the pair score, which is why the project's hub-stratified evaluations are necessary and why no "CLR removes hub effects" sentence should appear in the paper. A double-centred score (D_ij − μ_i − μ_j + μ̄, the two-way ANOVA residual) would remove additive hub effects; CLR does not compute it. (Whether that would help is a separate question; the project's finding that CLR's CORUM lead is a cross-row scale effect, `corum_expert2.md` A6, suggests the row scales carry useful information that double centring would keep but a rank readout would destroy.)

**[derivable]** CLR penalises large modules by a known factor. In the noise-free limit with row i holding n − 1 entries at a and G − n at b, z_ij = √((G − n)/(n − 1)) for a member pair, independent of a − b (`asymmetry_expert.md` section 3.2, which matches its simulation to within 2%). So CLR's within-module score falls as roughly √(G/n): a 4-gene module scores 26 and an 80-gene module 4.8 at G = 1,000. This is the sparsity assumption failing gracefully: the module contaminates its own members' null. It is also why |r| beats CLR on complexes (K562 7.00× vs 6.39× [measured]). This penalty is milder than aCLR's 1/n (Q2) and is the baseline tcaCLR must be compared with (Q3).

**[derivable]** CLR is not calibrated. Under a within-gene permutation null the row of MI is skewed (χ²-like), z_ij and z_ji are nearly the same statistic (corr 0.998 for |r|, 0.971 for MI [measured, `corum_expert2.md` section 3.1]), and both are clipped, so S has no standard-normal null. The Stouffer form (Stouffer et al. 1949 [confident]) is the right combination only for independent standard-normal z's; here it is a ranking statistic, the mean of two one-sided standardised scores times √2. The project has already corrected the paper's wording on this; the memo endorses that correction.

**[derivable]** Stouffer versus Euclidean is a choice about reciprocity. For fixed sum s = z_ij + z_ji, the Euclidean combination √(z_ij² + z_ji²) = √(s² − 2 z_ij z_ji) is maximised when one side is 0 and minimised when the sides are equal. So the 2007 Euclidean form ranks one-sided pairs above balanced pairs of the same total evidence, and the Stouffer (L1) form ranks them equally. Given the project's finding that one-sided pairs (satellites such as INTS8 attending the ribosome) are a false-positive mechanism, the Stouffer form is the more defensible of the two; min(z_ij, z_ji) would be the most reciprocal. This is a small but honest thing a theory section can say about why the project's CLR differs from the 2007 release, and it matches the mutual-rank idea (Obayashi and Kinoshita 2009 [confident]).

### 1.3 Where it breaks

- **Heavy tails [derivable].** The row SD is a non-robust scale estimate; a few extreme partners inflate σ_i and compress all other z's in that row. The project's robust variants (median/MAD, median/IQR) had no consistent effect (`docs/robust_z_results.md` [measured]), which says the tails in real rows are not the binding problem. The real rows are near-symmetric with an SD 3–8× the permutation null's (`corum_expert2.md` section 3.3): the "null" CLR estimates is each gene's typical dependence, not independence. That is the correct description of the estimand and should be stated as such.
- **Hubs [derivable]:** section 1.2, the h_i + h_j reintroduction.
- **Modules [derivable]:** the √(G/n) penalty. It is a bias that the pooled AUPR rewards or punishes depending on the size mix of the truth.
- **Dense dependence [plausible].** When most pairs are dependent (single-cell data with a dominant detection factor, mDC 99.4% positive r [measured]), the sparsity assumption fails for every row at once; CLR still ranks within rows but the cross-row scales become set by the global factor. The project's PC-removal and detection-residualisation results show this is second-order in practice.

### 1.4 Relation to ARACNE and the data-processing inequality

**[derivable]** ARACNE (Margolin et al. 2006 [confident]) prunes, in each triangle, the edge with the smallest MI, using the DPI: for a Markov chain X–Y–Z, I(X;Z) ≤ min(I(X;Y), I(Y;Z)). Its guarantee is structural: with exact MI and a tree-structured true network, DPI pruning recovers the tree (the Chow–Liu 1968 [confident] setting). CLR has no structural guarantee; it has the invariance of section 1.1. The two are complementary post-processors of the same MI matrix. The walk in aCLR goes in the opposite direction from DPI: P^t sums over paths, so it *rewards* indirect dependence. aCLR is in that sense an anti-ARACNE, which is consistent with its loss on direct TF–target edges (M3D AUPR ratio 5.17 vs CLR 6.75 [measured]) and its gain on regulon co-membership. A cheap check of this account is given in Q3 (prediction P8).

---

## Q2. aCLR and tcaCLR: clean descriptions and asymptotics

### 2.1 aCLR is a lazy random walk on a perplexity-calibrated attention graph

**[derivable]** A is, up to the choice of kernel, the conditional-probability matrix p_{j|i} of t-SNE (van der Maaten and Hinton 2008 [confident]): a per-row softmax whose per-row bandwidth is set by bisection so that the row's perplexity equals a fixed number (10 here, over a candidate set of 50). t-SNE uses exp(−d²/2σ_i²) on distances; aCLR uses exp(S_ij/τ_i) on CLR scores. UMAP's fuzzy kNN graph uses the analogous log₂k rule (McInnes et al. 2018 [confident]). So the operator is a perplexity-calibrated directed kNN affinity graph, a construction with a 15-year track record in embedding, and the reader will recognise it. Two consequences follow directly:

- The temperature bisection makes A_i· invariant to affine transformations of row i of S. CLR's row-i standardisation therefore reaches the operator only through (i) the clipping at 0, which removes candidates, and (ii) the partner term z_ji inside S_ij, which demotes hubs as candidates in every row (a hub j has large μ_j, σ_j and hence small z_ji with everyone). This was noted in `corum_expert2.md` section 3.1 and is worth stating in the paper, because it says what "CLR normalisation is essential" (the ablation in appendix A.2 of the pre-read) actually means: it is the partner-side normalisation and the clipping that matter, not the query-side scale.
- Softmax attention turns strength into share. Under exchangeability of a module's members (the lemma of `asymmetry_expert.md` section 3.1, which I re-derived and accept), any row-stochastic operator gives each within-module pair exactly m_i/(n − 1) where m_i ≤ 1 is the row's retained in-module mass, and any row, column or Sinkhorn scaling preserves that, because the scaling is constant on the module by the same symmetry. The 1/n dilution is therefore not a defect of the softmax that a different normalisation of the operator could remove; it is mass conservation.

**[derivable]** P = (I + A)/2 = I − L_rw/2, where L_rw = I − A is the random-walk Laplacian of the directed weighted graph; P^t is a polynomial low-pass filter of L_rw, which is the propagation step of SGC and APPNP without features (Wu et al. 2019 [confident]; Klicpera et al. 2019 [confident]). Expanding, P^t = 2^{−t} Σ_k C(t,k) A^k: the lazy walk at depth t is a binomial mixture over 0..t attention hops. At t = 9 the weights on k = 0..9 are (1, 9, 36, 84, 126, 126, 84, 36, 9, 1)/512, mean 4.5 hops, SD 1.5. This places aCLR in the network-propagation family (random walk with restart, Köhler et al. 2008 [confident]; heat kernel, Kondor and Lafferty 2002 [confident]; personalised PageRank; the review of Cowen et al. 2017 [confident]) as the member with a binomial rather than geometric or Poisson depth profile. The project's finding that PPR (+0.003) and the heat kernel (−0.001) are indistinguishable from P^9 on M3D (`docs/explore_idea3_results.md` [measured]) is what this view predicts: any smooth mixture centred near 4–5 hops is the same operator to within noise.

**[derivable]** The readout (M + Mᵀ)/2 is not a diffusion distance. Diffusion maps (Coifman and Lafon 2006 [confident]; Nadler et al. 2006 [confident]) need a reversible chain (symmetric kernel, P = D⁻¹K), for which P^t_ij/π_j is symmetric and the diffusion distance is defined. A is a directed kNN softmax; P is not reversible; π is not proportional to a degree and has no closed form. The symmetrised transition probability is a heuristic readout with no metric interpretation. tcaCLR's readout is also a heuristic, but with the invariance property of Q3.

**[derivable]** Spectral clustering, label propagation and lumpability are the right neighbours, not identities. P^t applied to a one-hot seed is unclamped label propagation (Zhou et al. 2004 [confident]) with a row-stochastic rather than symmetric-normalised operator. The random-walk view of normalised cut (Meilă and Shi 2001 [confident]) says that a good partition is one the walk rarely crosses, which is what the planted-module asymptotics below formalise. Nothing in aCLR extracts eigenvectors, so spectral-clustering consistency results (Rohe, Chatterjee and Yu 2011 [confident]; Lei and Rinaldo 2015 [confident]) do not transfer; they would need a different method.

### 2.2 Asymptotics of P^t on a graph with planted modules

**[derivable]** Two-timescale (nearly completely decomposable) behaviour. Let module C of n genes have per-step escape probability ε_C (the share of a member row's attention outside C) and let the lazy chain restricted to C have second eigenvalue λ₂(C). This is the setting of Simon and Ando (1961) [confident] and Courtois (1977) [confident], and of lumpability in Kemeny and Snell (1976) [confident]. For i, j ∈ C,

  P^t_ij = (1 − ε_C)^t ν_C(j) · (1 + O(λ₂(C)^t)),

where ν_C is the quasi-stationary distribution of the sub-chain, ν_C ≈ 1/(n − 1) off the diagonal for exchangeable members (the lemma). For a clique with near-uniform attention over its n − 1 members, the sub-chain's non-trivial eigenvalue is −1/(n − 1), lazy (1 − 1/(n − 1))/2 ≈ 1/2, so within-module mixing is complete to 0.2% by t = 9. For j ∉ C, P^t_ij is the leaked mass, Σ_{s<t} (1 − ε_C)^s ε_C (exit distribution · P^{t−s−1})_j, which is spread over many genes and is small per gene. As t → ∞ with the whole chain irreducible and aperiodic (laziness gives aperiodicity; irreducibility fails for the 0.28% unreachable pairs [measured]), P^t_ij → π_j for every i (Perron–Frobenius; Levin, Peres and Wilmer 2009 [confident]).

**[derivable]** Three regimes follow and the project has measured all three:
1. **Large modules (n ≥ keff): dilution, not leakage.** ε_C ≈ 0 (the perplexity constraint is satisfied inside C), so P^t_ij ≈ 1/(n − 1) for members: 88% of ribosomal-row mass stays in the complex at t = 9 and per-partner mass is 0.009 [measured]. The readout scale is set by n. Any pooled metric that compares across rows then under-ranks large modules.
2. **Small modules (n − 1 < keff): forced leakage.** A row with fewer than 10 strong candidates must put mass outside C to reach perplexity 10, so ε_C > 0 by construction (in-module mass 0.72 at n = 4, 0.23–0.26 retained at t = 9 [measured, simulation]). This is the "perplexity floor" and no readout can repair it, because the information has left the row.
3. **Deep walks: hub collapse.** Both readouts converge to a hub score: (M + Mᵀ)/2 → (π_i + π_j)/2 and (Q3) tcaCLR → z(log π_i) + z(log π_j). The rate is governed by the global λ₂, and the K562 π is extremely skewed (max/median 18,457 [measured]) because in-weight is concentrated (top 1% of genes hold 13.6% of in-weight at layer 1 [measured, `diagnostics_layers.csv`]) and the walk compounds it. The hub-stratified evaluation is therefore not optional; it is the correction for a known limit.

**[derivable]** Degree-corrected SBM and configuration-model nulls. On an undirected graph with expected weights w_ij ∝ θ_i θ_j (Chung and Lu 2002 [confident]; the degree-correction parameters of Karrer and Newman 2011 [confident]), the annealed transition probability is P_ij ≈ θ_j/Σθ, independent of i, and so is P^t_ij for every t ≥ 1. The multi-hop affinity under a pure hub null is a column effect: a function of the target only. This is the fact that matters for Q3.

### 2.3 Why aCLR gains on regulons and loses on edges [plausible]

Regulons are loose (target–target correlation small, median size 11 on RegulonDB SC [measured]); a 4–5-hop walk pools weak pairwise evidence along paths, which a pairwise statistic cannot do, so co-membership improves; the same pooling blurs direct edges (anti-DPI, section 1.4). This is plausible and matches every measured contrast, but it is not derived: a derivation would require a model of how regulon members' pairwise dependences are distributed, and the Shapley attribution (depth +0.018, CLR normalisation +0.013 [measured]) is the closest the project has to one.

---

## Q3. Why renormalise in tcaCLR, and what the second CLR is

### 3.1 What the row z-score of log P^t removes [derivable]

Write M_ij = P^t_ij. The row z-score of L = log M is invariant to multiplying row i of M by any positive constant c_i (log c_i is an additive row constant and is removed) and to raising row i to any positive power (a multiplicative constant in L, removed by the SD). So tcaCLR is exactly invariant to per-row *multiplicative* and *power* distortions of the walk. The row's budget (its partition function Z_i, its temperature τ_i, the retention (1 − ε_C)^t) are all row-multiplicative, and the softmax logits' scale is a row power. These are the quantities through which mass conservation makes module members' scores depend on n, and tcaCLR removes them. That is the derivable part of the "normalisation layer on both sides" picture.

It does **not** remove column effects. Under the configuration-model null of section 2.2, M_ij ≈ θ_j/Σθ, so L_ij = log θ_j + const is the same vector in every row, Z_ij = z(log θ_j) is the standardised log-hubness of the target in every row, and S'_ij = (z(log θ_i) + z(log θ_j))/√2. **Under a pure hub null, tcaCLR ranks pairs by the sum of standardised log-hubnesses.** It is therefore *not* a z-score against a configuration-model null. What it is, is a z-score against each gene's own background of reach probabilities, exactly as CLR is a z-score against each gene's own background of MI. "Diffusion CLR" is an accurate name; "hub-aware null" is not. The hub awareness tcaCLR does have is partial and comes from two places: (a) the log compresses the π skew, so a hub's z is bounded (a hub pair cannot dominate the ranking the way π_i + π_j lets it); (b) the Stouffer step adds the reverse direction, and a hub j that everyone reaches is itself spread over many targets, so Z_ji is unremarkable; a true module pair is high in both directions. The project's measurement that the hub-stratified gain "shrinks but stays" (`corum_analysis.md` section 4, hub control) is consistent with partial, not full, hub correction. Prediction P4 below tests this directly.

### 3.2 Why the log [derivable, with one plausible step]

Three independent reasons, each derivable:

1. **The 1/n becomes additive.** For i, j ∈ C, L_ij ≈ t log(1 − ε_C) − log(n − 1). Against a row background of reach log-probabilities with mean μ_bg,i and SD σ_bg,i (over the thousands of non-members i reaches), Z_ij ≈ (−log(n − 1) − μ_bg,i)/σ_bg,i. In the linear domain, by contrast, the n member entries dominate the row's SD and the member z-score is √((G − n)/(n − 1)) again (same algebra as CLR's), i.e. the size penalty survives. `asymmetry_expert.md` section 3.2 derived both and its simulation confirmed them (a3r pair score 5.4–8.3 with no trend in n from 4 to 40; the linear-domain CLR on P^9 fails at 2.46× [measured]). The size dependence of tcaCLR is log(n)/σ_bg per nat of size, which is small if σ_bg is several nats. **The one step that is plausible rather than derived is that σ_bg is several nats and roughly homogeneous across rows.** Prediction P2 checks it.
2. **Products become sums.** P^t_ij is a sum over t-step paths of products of transition probabilities. For background pairs the dominant contributions are products of several small numbers, so the row distribution of M_ij spans orders of magnitude and is extremely right-skewed; its log is far closer to symmetric (a log-domain central-limit argument for multiplicative cascades; no specific reference, the argument is standard). Moment estimators (mean, SD) are only usable on the log scale. Prediction P3.
3. **At t = 1 the log readout is CLR itself.** M_ij = A_ij/2 for j ≠ i, and log A_ij = S_ij/τ_i − log Z_i for j ∈ C_i, so the row z-score over reachable j is Z_ij = (S_ij − mean_{C_i} S_i·)/sd_{C_i}(S_i·): tcaCLR at t = 1 is exactly the Stouffer combination of CLR scores re-standardised within each gene's candidate set (rows with ≤ 10 candidates get uniform A and an undefined z; these are rare, median candidate count 50 [measured]). So tcaCLR is a continuous deformation of CLR in t, with the walk inserted between the two normalisations. Prediction P1 is the identity check.

A further derivable identity places tcaCLR next to MCL. Un-logging, exp(Z_ij) = (M_ij / g_i)^{1/σ_bg,i} with g_i the row geometric mean: tcaCLR is a **row-adaptive power transform (inflation) of P^t**, with the exponent set by each row's background spread, without MCL's re-normalisation (van Dongen 2000 [confident]; Enright, van Dongen and Ouzounis 2002 [confident]). The project's soft-MCL variant with a fixed exponent chose r = 1 label-free and all inflated readouts lost (`docs/explore_ideas47_results.md` [measured]); the per-row exponent is what differs. [plausible] That the adaptive exponent is the useful ingredient is a conjecture, testable with a fixed-exponent control inside the simulation of Q5.

### 3.3 Why a second CLR and not another normalisation

- **Sinkhorn / column scaling of P^t [derivable]:** cannot touch dilution, by the lemma (scalings are constant on an exchangeable module). What Sinkhorn balancing of A *would* do is make the stationary distribution uniform and so remove hub drift entirely. So Sinkhorn fixes the hub limit and not the dilution; tcaCLR fixes the dilution (approximately) and not the column hub effect. They are complementary. [speculation] Sinkhorn(A) → lazy walk → tcaCLR might combine both. It is not among the ~20 K562 readouts listed in `corum_analysis.md`; the asymmetry memo simulated a doubly-stochastic walk with the raw readout only. Caveat: Sinkhorn on a sparse non-symmetric A needs total support, which the lazy diagonal provides (Sinkhorn and Knopp 1967 [confident]). If tried, it is in-sample on K562 and must be labelled so.
- **Lift P^t_ij/π_j and log-lift [derivable why it fails]:** this is the diffusion-map / PMI normalisation (NetMF shows DeepWalk factorises log of the depth-averaged walk divided by the target's stationary mass; Qiu et al. 2018 [confident]; Levy and Goldberg 2014 [confident]). Within an exchangeable module the lift is n-invariant (π_j ∝ 1/n cancels the 1/n), so it is the textbook fix. It failed here (1.35× and 4.26× [measured]) for two reasons that follow from section 2.2: π is the t → ∞ limit and at t = 9 the chain is far from mixed (0.28% of pairs unreachable), so dividing by π_j injects the irrelevant global limit into a local quantity; and the chain is non-reversible, so P^t_ij/π_j is not symmetric and has no kernel interpretation. On the small-π side the ratio explodes for genes nobody attends to. [plausible] A t-matched column normaliser, the column mean of P^t itself (mean_k P^t_kj, the average in-flow to j at depth t), would avoid the first problem. Prediction P5. [speculation] On a symmetrised (reversible) operator the log-lift would coincide with PMI and might work; that is a different method.
- **Rank readout [measured, `corum_expert2.md` A6]:** a within-row percentile is the scale-free member of the same family and it also helps aCLR on CORUM (3.51× → 5.13×) but hurts CLR (6.39× → 3.36×). This is the cleanest evidence that the a3r gain is cross-row calibration and that the *scale* information in log P^t (which the rank discards) is worth keeping: a3r 7.95× > rank 5.13×.
- **Degree-corrected readout (d_i P_ij + d_j P_ji)/2 [measured, A5]:** +0.009 against aCLR, far below a3r. Consistent with section 3.1: an additive degree correction does not remove a multiplicative row factor.

### 3.4 Cheap falsifiable predictions for the plausible parts

The saved caches (`results/explore_stack/*/cache/readout__*.rds`, `results/robust_z/k562/cache/V0__*.rds`) hold scored pairs, not walk matrices, so each check needs P^9 recomputed on K562 half A (`results/baselines_k562ess_truthgenes_genie3/expression.tsv`; MI, CLR, softmax, nine sparse-dense products; the stacking plan quotes about a minute for the softmax and seconds for the walk). Every check is label-free except P5 and P6, which are in-sample K562.

| # | Claim tested | Prediction | Refutes the claim if |
|---|---|---|---|
| P1 | tcaCLR(t = 1) = CLR re-standardised within candidate sets (3.2, item 3) [derivable; the check is for the implementation] | max abs difference between a3r at t = 1 and `stouffer_sym(rowz over C_i of S)` ≤ 1e-10 after clipping, excluding rows with ≤ 10 candidates | any difference above tolerance |
| P2 | σ_bg is several nats and homogeneous (3.2, item 1) [plausible] | per-row SD of log P^9 over reachable off-diagonal entries: median 3–6 nats, coefficient of variation across rows < 0.3; ribosomal-row member entries sit ≥ 2 row-SDs above the row mean | median σ_bg < 2 nats, or CV > 0.5, or member entries < 1 SD above the mean |
| P3 | the log makes moments usable (3.2, item 2) [plausible] | median row skewness of log P^9 (reachable entries) with absolute value < 1.5; median row skewness of P^9 itself > 10 | log-domain skewness comparable to linear |
| P4 | tcaCLR does not remove the column hub effect (3.1) [derivable; the check is of the size of the effect] | over pairs sharing no CORUM complex, Spearman(Z_ij, log π_j) > 0.5 within rows and Spearman(S'_ij, z(log π_i) + z(log π_j)) > 0.4 pooled; the same correlations for aCLR's readout against π_i + π_j are higher | correlation < 0.2 (then the readout is more hub-neutral than derived, and section 3.1's account of the Stouffer step is incomplete) |
| P5 | a t-matched column normaliser helps where π does not (3.3) [plausible] | CLR(log(P^9_ij / mean_k P^9_kj)) scores between a3r and the π log-lift on raw K562 CORUM and at or above a3r hub-stratified | it scores at the π log-lift level (4.3×) or below |
| P6 | the Stouffer step adds reciprocity information (3.1 b) [plausible] | row-local AP on RP member rows of the one-directional Z (before Stouffer) equals that of log P^9 (identical ranks within rows, by construction) and is below the 0.668 of the symmetrised a2/a3r readout | one-directional Z already reaches 0.67 (then the gain is entirely cross-row scale and reciprocity contributes nothing) |
| P7 | size dependence of within-module scores (3.2, item 1) [derivable in the exchangeable limit] | regress the mean within-complex pair score on log complex size over merged CORUM complexes with ≥ 10 pairs: slope ≈ −1 for aCLR on log score, ≈ −0.5 for CLR on log score, and |slope| < 0.3 in z-units per nat for tcaCLR | tcaCLR's slope is as steep as CLR's |
| P8 | aCLR is anti-DPI (1.4) [derivable in direction; size unknown] | among MI triangles in the top 5% of CLR, the DPI-weakest edge is in the method's top 1% more often for aCLR than for CLR, and for tcaCLR in between | aCLR keeps weakest edges less often than CLR |

---

## Q4. Why stacking (a3r2, a3r3) did not help

The results are in `docs/explore_stack_results.md` and `diagnostics_layers.csv` [measured]: K562 CORUM-all 7.95× → 8.19× → 8.15× (within noise of 7.95× and of one layer at matched depth 8.00×/7.86×); every non-clique truth falls (CollecTRI −0.006/−0.011, KD −0.044/−0.079, DoRothEA −0.002/−0.008); M3D complexes rise (+0.036/+0.047) while regulons fall; unreachable pairs grow 0.28% → 5.7% → 12.3% → 17.2%; hub drift grows (Spearman of π_k with π_1: 0.74, 0.52, 0.40); the pure repeat c2 is worse, not equal (5.85×, 6.83×). Four accounts, in order of how much they are supported:

1. **Depth is not the limit; re-truncation is [derivable from the data].** The matched-depth control (one layer, t = 18 or 27) loses nothing on KD (+0.024, +0.032) and CollecTRI (+0.004, +0.005), while the stacked readouts at the same total depth lose heavily on both. Stacking = depth + (clip, top-50 re-selection, re-softmax, re-standardise). Depth alone is neutral-to-positive, so the second operator construction is what hurts. This is the cleanest conclusion the experiment supports.

2. **Community closure and over-smoothing [derivable mechanism, plausible as the cause].** After nine steps, the top 50 entries of S₁ in a row are the row's diffusion neighbourhood (its community), not its direct partners. An operator built from those is more modular than A: candidate sets close within communities, cross-community conductance falls, components appear, and the walk on it mixes less across communities. The measured growth of unreachable pairs (0.28% → 17%) and the fall of π_max/median at deeper layers (18k → 24k → 13k → 4.8k → 1.5k, as the graph fragments into components each with its own stationary distribution) are the signatures. Truths that depend on cross-community paths (KD modules, CollecTRI) lose; cliques do not gain because they were already closed at layer 1. This is the over-smoothing of repeated graph propagation (Li, Han and Wu 2018 [confident]; Oono and Suzuki 2020 [confident]) in a sharpened form: each layer both smooths within communities and cuts between them, the iterated-kNN-sharpening behaviour of MCL-like schemes.

3. **Saturation of the clique signal [derivable from the data].** On K562 the ribosomal rows are already ordered nearly perfectly after one layer (row-local AP 0.67, per-complex AP 0.83 [measured]); the residual pooled error sits in small complexes (perplexity floor, unfixable by any readout, Q2 regime 2) and in cross-complex calibration, which one standardisation already does. A second calibration cannot add information, and a second walk at matched depth adds none (control c1). The ceiling was reached at layer 1.

4. **Idempotence [derivable, and it is the wrong explanation].** The row z-score is idempotent on an unclipped matrix (z of z is the identity), but the full readout is not: clipping is nonlinear (z(clip(z(x))) ≠ clip(z(x))), the Stouffer symmetrisation gives rows of a symmetric matrix different means, and the top-50 truncation discards everything else. The pure repeat c2 therefore changes the ranking and, because each row's zero fraction differs, it changes it for the worse (5.85×). So stacking does not fail by doing nothing; it fails by doing something harmful (accounts 1–2) to a signal that was already saturated (account 3).

[plausible] M3D's complex rise with stacking (+0.047 at three layers, hub-stratified about half) is what account 2 predicts for a truth dominated by one large closed module (60% of positives in r-protein complexes [measured]): sharpening the ribosome's community at the cost of everything else. It is not evidence of a lever.

---

## Q5. Verdict, honest theory section, the one result worth proving, and a simulation

### 5.1 Verdict

This is a heuristic with a clean description and a handful of derivable properties, sitting inside recognisable frameworks; it is not near a theorem of optimality or consistency, and the paper should not imply one. The recognisable frameworks are:

- **CLR:** per-gene empirical-null studentisation (Efron 2004), invariant to per-gene affine distortion of the dependence estimate, with a derivable √(G/n) module-size penalty and a reciprocity-favouring (Stouffer) combination.
- **aCLR:** a lazy random walk on a perplexity-calibrated (t-SNE-style) directed kNN attention graph built from CLR scores; a binomial mixture over 0..t attention hops, i.e. a member of the network-propagation family; nearly-completely-decomposable chain behaviour on planted modules, with an exact 1/(n − 1) dilution in large modules (mass conservation), forced leakage below keff, and convergence to a hub score with depth.
- **tcaCLR:** per-row empirical-null studentisation of log t-step reach probabilities; exactly invariant to row-multiplicative and row-power distortions of the walk; turns the 1/n dilution into an additive −log n against a background of several nats; reduces to CLR on candidate sets at t = 1; equivalent to a row-adaptive MCL inflation; does not remove column (hub) effects, so under a hub null it is a sum of standardised log-hubnesses.
- **Stacking:** re-truncation, not depth, is the harm; one calibration suffices; cliques saturate at one layer.

The transformer framing ("attention conserves mass; CLR is the normalisation layer on both sides", pre-norm and post-norm as in Xiong et al. 2020 [confident]) is a fair informal description, not a result.

### 5.2 What a theory section could claim honestly

1. CLR's row z-score is a gene-specific empirical-null z-value (Efron), consistent under sparsity, invariant to per-gene affine distortions of the MI estimator, not calibrated, not hub-neutral in the pair score, and penalising large modules as √(G/n). [all derivable]
2. The attention operator is a perplexity-calibrated kNN affinity (t-SNE's p_{j|i}) on CLR scores; the lazy walk at t = 9 is a binomial average over 0–9 hops centred at 4.5; the readout tends to the hub score π_i + π_j with depth. [derivable]
3. Within an exchangeable module of n genes, any row-stochastic operator (and any row/column/Sinkhorn rescaling of it) gives each pair exactly m/(n − 1); the readout's scale is set by module size. [derivable; the project's lemma]
4. tcaCLR removes exactly the row-multiplicative and row-power factors of the walk, converts 1/n into an additive log n offset, equals CLR on candidate sets at t = 1, and leaves the column hub effect, which the bilateral combination corrects only partly. [derivable, with σ_bg's size as the one measured input]
5. The size-regime prediction: aCLR's within-module score falls as 1/n, CLR's as n^{−1/2}, tcaCLR's by (log n)/σ_bg; below keff all three lose to forced leakage; above the mixing time all collapse to hub scores. [derivable in the exchangeable limit]
6. Stacking cannot add to a calibration fix and harms by re-truncation. [measured, with a derivable mechanism]

Claims to avoid: calibration or p-value interpretations of any score; "CLR removes hub effects"; "diffusion distance"; "configuration-model null"; any optimality of t = 9, keff = 10 or the Stouffer form; any suggestion that the size-regime prediction has been confirmed out of sample (the yeast size stratification matches it, +0.137 for ≥ 20 genes against −0.004 and −0.003 for 3–6 and 7–19 [measured], but that study was exploratory and the prediction was not registered before it; RPE1 merged CORUM, where a3r is below CLR, 6.14 vs 6.77 [measured], is not explained by anything here and must be reported as such).

### 5.3 The one result that would most strengthen the paper

**A size-invariance proposition with its regime boundaries.** Under a planted-module model (modules of size n exchangeable under the operator, in-module retention m(n) per step, background reach log-probabilities with spread σ_bg bounded below), the expected within-module pair score satisfies: aCLR ∝ m(n)^t/(n − 1); CLR ∝ √((G − n)/(n − 1)); tcaCLR = (−log(n − 1) + t log m(n) − μ_bg)/σ_bg, i.e. bounded below uniformly in n for n ≥ keff, with the two failure boundaries n < keff (m(n) < 1 by the perplexity constraint) and t ≫ t_mix (all three → hub scores). The exchangeable case is a short derivation (it is already three-quarters written across `asymmetry_expert.md` sections 3.1–3.2 and section 3.2 of this memo); the general case with noise is the simulation below. Its value is that it converts "tcaCLR won on K562 and in a yeast estimate" into a mechanistic prediction about *where* it wins (modules ≥ 20, growing with n), *where it ties* (7–19), *where it loses* (≤ 6, direct edges, t too deep). That prediction can be registered before any untouched complex truth is opened and tested without a new selection step, which is exactly what the project's a3r evidence currently lacks.

### 5.4 Simulation specification (for a worker; do not run under this memo)

**Purpose.** Demonstrate the mechanism of 5.3 and map the regime where tcaCLR should beat aCLR and CLR, and where it should not. Known truth; no real data; no tuning of any method parameter (keff 10, cap 50, α 0.5, t = 9 fixed; depth grid descriptive only).

**Generator (Gaussian factor model, one draw = one "compendium").**
- G = 3,000 genes, N ∈ {150, 400, 800} samples (bracketing M3D 466 and K562 707).
- Module factors: 40 modules of size 5, 20 of 10, 10 of 20, 5 of 40, 3 of 70, 2 of 100 (1,210 member genes; positive pairs 400, 900, 1,900, 3,900, 7,245, 9,900). Each module has a loading λ ∈ {0.55, 0.70, 0.85} (member–member r = λ² = 0.30, 0.49, 0.72), assigned round-robin so every size has every λ. Members load λ on their module factor only.
- Regulons: 30 regulators, each with m ∈ {5, 15, 40} targets drawn from non-member genes (a gene may be in at most two regulons). Regulator mRNA loads λ_R = 0.6 on its activity factor; targets load b ∈ {0.35, 0.50} (target–target r 0.12–0.25, regulator–target r 0.21–0.30).
- Hubs: 50 non-member genes, each loading 0.25 on a random 30% of all module and regulon factors.
- Satellites: for each module of size ≥ 20, three non-member genes loading 0.8λ on that module's factor (truth negatives, the INTS8 case).
- Global factors: three, every gene loads N(0, 0.15²) on each (dense background dependence, so rows look like real rows: near-symmetric, inflated SD).
- Unique noise scaled so every gene has unit variance. Noise arms: Gaussian (primary); Student-t with 4 df (heavy-tail arm).
- Marginal arm: identity (primary); Gaussian copula to log-normal marginals with gene-specific log-scale σ_g ~ U(0.3, 1.0) (count-like marginals, so the B-spline MI bias varies by gene and CLR's affine invariance is exercised).
- Cells: N (3) × noise (2) × marginal (2) = 12; 10 seeds per cell; 120 compendia.

**Methods (all from the package; nothing tuned).** |r|; CLR (B-spline MI, bins "hg", Stouffer); aCLR (keff 10, cap 50, α 0.5, t = 9); tcaCLR (t = 9). Diagnostic variants, same t: tcaCLR at t = 1 (identity check P1); Sinkhorn-balanced A then walk, raw readout (lemma check); π log-lift and column-mean log-lift (section 3.3); fixed-exponent inflation r ∈ {1.5, 2} of P^9 with row z (the MCL control of 3.2); Sinkhorn(A) then walk then tcaCLR (the speculation of 3.3). Depth grid t ∈ {3, 6, 9, 18, 36, 72} descriptive for CLR-free methods. One extra arm with keff 40, cap 200 (wide rows) to produce the "small σ_bg" failure regime.

**Truths and negatives.** T1 module co-membership, pooled over sizes; T2 module co-membership by size bin (5, 10, 20, 40, 70, 100); T3 regulon co-membership (regulator plus targets); T4 direct regulator–target edges; negatives are all other pairs among evaluated genes, with satellite–member pairs kept as negatives and also reported separately.

**Metrics.** Pooled AUPR over base rate for T1–T4; per-bin AUPR for T2; macro per-module AP; hub-stratified AUPR (10 strata of true hub loading sum, and of MI row mean, as the project does); satellite false-positive share (fraction of the top 1% of pairs that are satellite–member); mechanism curves: mean within-module pair score against log n with a fitted slope (log score for aCLR and CLR; z-units for tcaCLR), per-row σ_bg (SD of log P^9 over reachable entries), in-module retention m(n) at t = 1 and t = 9, hub share of the top 1% as a function of t. Intervals: across the 10 seeds per cell (SD and range); no p-values needed.

**Pre-specified predictions (written before running).**
- R1 (T2, n ≥ 20): tcaCLR > CLR > aCLR, with the aCLR deficit and the tcaCLR–CLR gap both growing with n; slopes ≈ −1 (aCLR), ≈ −0.5 (CLR), |slope| < 0.3 z/nat (tcaCLR).
- R2 (T2, n = 5): CLR ≥ tcaCLR ≈ aCLR; all three below their n = 20 level (forced leakage); n = 10 intermediate.
- R3 (T3, loose regulons): aCLR ≈ tcaCLR > CLR, gap largest for b = 0.35 and m = 40.
- R4 (T4, direct edges): CLR ≥ tcaCLR ≥ aCLR; |r| competitive.
- R5 (heavy-tail arm): CLR's size penalty becomes noisier (σ_i inflated), tcaCLR's advantage on n ≥ 20 grows; (copula arm): CLR and tcaCLR unchanged in rank, |r| loses.
- R6 (identity): tcaCLR at t = 1 equals Stouffer of within-candidate row z of CLR to 1e-10.
- R7 (lemma): Sinkhorn-balanced walk, raw readout, ties aCLR on n ≥ 20 (within seed SD) and reduces the hub share of the top 1%.
- R8 (depth): at t = 36 and 72 all methods' T1 falls toward the hub-score level and the hub share of the top 1% rises; aCLR's rises faster than tcaCLR's.
- R9 (satellites): aCLR ranks satellite–member pairs above member–member pairs of the same module for n ≥ 20; CLR and tcaCLR rank member–member higher.
- R10 (wide-row arm, keff 40): σ_bg falls, tcaCLR's T1 gain over aCLR shrinks or vanishes; CLR unaffected.
- R11 (lift): π log-lift unstable (large seed SD, hub-dominated top 1%); column-mean log-lift between aCLR and tcaCLR on T1.
- R12 (inflation control): fixed-exponent inflation does not reach tcaCLR on T1 for both r; Sinkhorn(A) → walk → tcaCLR is at or above tcaCLR on hub-stratified T1 (speculation; a failure here is informative).

**Compute.** Per compendium: MI on 3,000 × 800 with the package's C++ a few minutes at 8 threads; softmax about 30 s; walks and scoring seconds. 120 compendia plus variants: roughly 6–10 h at 8 threads, well under the project's scheduler limits; declare memory from the first cell's measured peak.

**What a failure of R1 would mean.** If tcaCLR's slope is as steep as CLR's or its n ≥ 20 advantage is absent in the Gaussian arm, then the size-invariance account of Q3 is wrong, the K562 and yeast gains are something else (most likely ribosome-specific), and the paper should present tcaCLR only as the measured mechanism of the K562 loss, which is the PI's current position.

---

## References

Confidence flags: [confident] = I am sure the work exists as cited; [verify] = check the details before citing.

- Bolya D et al. (2023). Token Merging: Your ViT but faster. ICLR. [confident]
- Chow CK, Liu CN (1968). Approximating discrete probability distributions with dependence trees. IEEE Trans Inf Theory 14:462–467. [confident]
- Chung F, Lu L (2002). The average distances in random graphs with given expected degrees. PNAS 99:15879–15882. [confident]
- Coifman RR, Lafon S (2006). Diffusion maps. Appl Comput Harmon Anal 21:5–30. [confident]
- Courtois PJ (1977). Decomposability: Queueing and Computer System Applications. Academic Press. [confident]
- Cowen L, Ideker T, Raphael BJ, Sharan R (2017). Network propagation: a universal amplifier of genetic associations. Nat Rev Genet 18:551–562. [confident]
- Daub CO, Steuer R, Selbig J, Kloska S (2004). Estimating mutual information using B-spline functions. BMC Bioinformatics 5:118. [confident]
- Efron B (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. JASA 99:96–104. [confident]
- Efron B (2010). Large-Scale Inference. Cambridge University Press. [confident]
- Enright AJ, van Dongen S, Ouzounis CA (2002). An efficient algorithm for large-scale detection of protein families. Nucleic Acids Res 30:1575–1584. [confident]
- Faith JJ et al. (2007). Large-scale mapping and validation of E. coli transcriptional regulation from a compendium of expression profiles. PLoS Biol 5:e8. [confident]
- Karrer B, Newman MEJ (2011). Stochastic blockmodels and community structure in networks. Phys Rev E 83:016107. [confident]
- Kemeny JG, Snell JL (1976). Finite Markov Chains. Springer. [confident]
- Klicpera J, Bojchevski A, Günnemann S (2019). Predict then propagate: graph neural networks meet personalized PageRank. ICLR. [confident]
- Köhler S, Bauer S, Horn D, Robinson PN (2008). Walking the interactome for prioritization of candidate disease genes. Am J Hum Genet 82:949–958. [confident]
- Kondor RI, Lafferty J (2002). Diffusion kernels on graphs and other discrete structures. ICML. [confident]
- Lei J, Rinaldo A (2015). Consistency of spectral clustering in stochastic block models. Ann Stat 43:215–237. [confident]
- Levin DA, Peres Y, Wilmer EL (2009). Markov Chains and Mixing Times. AMS. [confident]
- Levy O, Goldberg Y (2014). Neural word embedding as implicit matrix factorization. NeurIPS. [confident]
- Li Q, Han Z, Wu X-M (2018). Deeper insights into graph convolutional networks for semi-supervised learning. AAAI. [confident]
- Margolin AA et al. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics 7(Suppl 1):S7. [confident]
- McInnes L, Healy J, Melville J (2018). UMAP: uniform manifold approximation and projection for dimension reduction. arXiv:1802.03426. [confident]
- Meilă M, Shi J (2001). A random walks view of spectral segmentation. AISTATS. [confident]
- Nadler B, Lafon S, Coifman RR, Kevrekidis IG (2006). Diffusion maps, spectral clustering and reaction coordinates of dynamical systems. Appl Comput Harmon Anal 21:113–127. [confident]
- Obayashi T, Kinoshita K (2009). Rank of correlation coefficient as a comparable measure for biological significance of gene coexpression. DNA Res 16:249–260. [confident]
- Oono K, Suzuki T (2020). Graph neural networks exponentially lose expressive power for node classification. ICLR. [confident]
- Qiu J et al. (2018). Network embedding as matrix factorization: unifying DeepWalk, LINE, PTE, and node2vec. WSDM. [confident]
- Rohe K, Chatterjee S, Yu B (2011). Spectral clustering and the high-dimensional stochastic blockmodel. Ann Stat 39:1878–1915. [confident]
- Simon HA, Ando A (1961). Aggregation of variables in dynamic systems. Econometrica 29:111–138. [confident]
- Sinkhorn R, Knopp P (1967). Concerning nonnegative matrices and doubly stochastic matrices. Pacific J Math 21:343–348. [confident]
- Steuer R, Kurths J, Daub CO, Weise J, Selbig J (2002). The mutual information: detecting and evaluating dependencies between variables. Bioinformatics 18(Suppl 2):S231–S240. [confident]
- Stouffer SA et al. (1949). The American Soldier, vol. 1. Princeton University Press. [confident]
- van der Maaten L, Hinton G (2008). Visualizing data using t-SNE. JMLR 9:2579–2605. [confident]
- van Dongen S (2000). Graph clustering by flow simulation. PhD thesis, University of Utrecht. [confident]
- Wu F et al. (2019). Simplifying graph convolutional networks. ICML. [confident]
- Xiong R et al. (2020). On layer normalization in the transformer architecture. ICML. [confident]
- Zhou D, Bousquet O, Lal TN, Weston J, Schölkopf B (2004). Learning with local and global consistency. NeurIPS. [confident]
- No reference is given for the log-domain central-limit argument of section 3.2 item 2; it is a standard argument and I would not attach a citation to it. The "multiplicative cascade" phrase should not be cited to any specific paper without checking.
