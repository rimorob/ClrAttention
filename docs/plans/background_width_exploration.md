# Background-adaptive attention width: exploratory design memo

*30 October 2026. Exploratory design memo (statistics expert, Claude Opus, high effort) for the PI. **Not a preregistration.** **No truth set was read or scored**: no CORUM, RegulonDB, DoRothEA, CollecTRI, Complex Portal or operon file was opened, and no network was scored against a real truth. What was run:*

- *synthetic simulations: `analysis/synthetic_width/` (G = 2,000 genes, N = 500 samples, 2 scenarios × 4 seeds, plus two follow-ups);*
- *a label-free diagnostic on real expression only: M3D (4,297 × 466) and K562 essential, half A (5,671 × 707), read from `results/baselines_k562ess_truthgenes_genie3/expression.tsv` (the truth files in that folder were not opened).*

*Everything ran at nice 19, one thread, in a 4 GB memory scope, about one hour of single-thread CPU in total. Nothing in `results/`, the paper, `docs/STATUS.md` or the preregistrations was changed. Each claim below is marked **[ran]** (a number from these runs), **[derived]** (algebra, shown) or **[reasoned]** (argument only, not tested).*

**Glossary.** keff is the number of effective partners per gene (the softmax perplexity). K10 is the method as it stands: softmax over the top 50 positive Stouffer CLR candidates at keff 10, lazy walk with α = 0.5, t = 9, readout (Pᵗ + Pᵗ′)/2. "×" is pooled AUPR divided by the base rate. Ratios are geometric means over 4 paired seeds, against K10, with the number of seeds in which the variant is better in brackets. A **partner-less gene** is one that loads on no module.

## 1. Verdict

**The PI's idea contains one genuinely new and correct insight, but it is attached to the wrong lever.**

- **What is new relative to O2.** O2 took its background from a permutation null, pooled over all pairs, on the CLR scale. The PI's version estimates the background from the observed data themselves, does not assume it is Gaussian, and conditions it on both genes of a pair (the X-by-Y quantity). That removes O2's main defect.
  - O2's permutation null is the null of no dependence at all. On the CLR scale it is not pivotal: real rows have their SD inflated by real dependence, so observed z-scores are compressed relative to the permuted ones. This gives lfdr ≈ 1 and keff at the floor of 2.
  - The same failure reproduces in every synthetic run: π0 = 1.000 in 8/8 runs, the median partner-less gene at keff 2, and K10 ratio ×0.58–0.59 on complexes. **[ran]**
- **The weak link is "width = how many partners rise above background".** Three results argue against it.
  1. **Even a perfect width does not fix large complexes. [derived, ran]**
     - By mass conservation, a pair inside a mixed n-gene clique carries m/(n − 1) at any width.
     - Giving every complex member its true width n − 1 (an oracle) barely moves the 80-gene ribosome in the messy world (×1.04) and helps it modestly in the clean one (×1.60). It costs 20-gene complexes ×0.29 / ×0.53, because a softmax can only reach perplexity n − 1 over n − 1 unequal partners by leaking mass to non-members.
     - It helps only tiny complexes (n = 5: ×7.7 / ×2.9).
  2. **The count is not identifiable in realistic data. [ran, reasoned]**
     - Dependence induced by batches, labs or global programs is real dependence: two genes that share a lab profile form a module as far as MI can tell. So no function of the MI distribution can separate "batch partner" from "biological partner" without covariates.
     - The graphon (X-by-Y) background, the literal form of the idea, puts every gene at roughly 150–540 "above-background" partners: synthetic messy 403–451, synthetic clean 153–231, M3D median 538, K562 median 212. That is 100% at any sensible ceiling, which reproduces O2e.
     - A per-row extreme-value tail model finds about 11 partners for partner-less genes and 12 for genes in small complexes. It cannot tell them apart.
  3. **A count of "partners above background" measures statistical power, not biology. [derived]** The smallest detectable correlation is about tanh(3.5/√(N − 3)): 0.31 at N = 125 and 0.13 at N = 707. A width defined by detection therefore grows with sample size.
- **What survives is the same intuition placed on the row, not on the width.** The lever is what genes *without* real partners do with their budget.
  - Under fixed-budget attention, every partner-less gene spends weight 1 on about 10 junk partners. In the pooled ranking these junk pairs (about 1/10 each) outrank ribosome pairs (about 1/79 each).
  - Silencing partner-less genes in both directions (oracle "symmetric sink"), with nothing else changed, gives:

    | Endpoint | messy | clean |
    |---|---|---|
    | all complexes | ×1.36 | ×1.10 |
    | ribosome | ×1.54 | ×1.11 |
    | small complexes (3–6 genes) | ×2.61 | ×1.79 |
    | regulons | ×3.23 | ×2.25 |

    All 4/4 seeds in both scenarios. **[ran]**
  - Spreading partner-less rows thin instead of silencing them (oracle keff 50 for them only) gives ×1.56 on complexes and ×1.34 on regulons (messy).
  - This also explains the K25 result. K25 spreads *every* row thin: ribosome ×3.73, regulons ×1.05, but small complexes ×0.23. That matches the real keff verdict: K25 won on ribosome-dominated K562 CORUM and was non-negative on regulons.

So the useful question is not "how wide should each gene attend?" but **"does this gene have any partner at all, and if not, why is it spending attention?"**.

- A realizable, label-free detector for that question is the open problem.
  - A cross-half reproducibility gate with a threshold fixed by a chance-overlap rule, combined with reciprocity, recovers part of the oracle gain in both scenarios: all complexes ×1.08 / ×1.04, ribosome ×1.31 / ×1.27, small complexes ×1.71 / ×1.19 and regulons ×1.19 / ×1.40 (messy / clean), at a cost on 20-gene complexes (×0.82 / ×0.87). **[ran]**
- Reciprocity (attention only where it is returned) is a cheap partial answer that already works for small complexes: small ×2.39 messy, ×1.50 clean.
- The real-data diagnostic says how much is at stake: about 10% of K562 genes, but only 0.3% of M3D genes, have no reproducible partner across data halves.

## 2. Wild ideas (generation phase)

Probably-wrong ideas are marked ✗. **[ran]**, **[derived]** and **[reasoned]** say how far each one was taken.

1. **Graphon background (PI, literal).** Bin genes by their background level θᵢ (the row median of MI^⅓). The background of pair (i, j) is the empirical MI distribution of all pairs in cell (bin_i, bin_j), and q_ij = 1 − lfdr. Width = Σⱼ q_ij. **[ran]** ✗: everyone gets hundreds of partners (section 1).
2. **Wilson–Hilferty plus Tukey.** Under independence, 2N·MI is roughly chi-square, so MI^⅓ is roughly Gaussian. Fit a symmetric two-way median polish, MI^⅓_ij = μ + aᵢ + aⱼ + r_ij, and standardize by each gene's MAD. This is a robust CLR, "rCLR". **[ran]** As a score it trades small complexes for the ribosome. As attention input it is mixed (section 4.3).
3. **APC from protein coevolution (Dunn et al. 2008).** MI_ij − MIᵢ·MI·ⱼ/MI··. "Integrate out one gene, then the other" is literally a product (rank-one) background. **[reasoned]** It is a bias correction for the score; it gives no width.
4. **Robust PCA on the MI matrix** (low rank + sparse): background = low-rank part. **[reasoned]** ✗: an 80-gene clique *is* a rank-one block and would be absorbed into the background.
5. **Extreme-value row tails.** Fit an exponential (GPD) tail to the 90th–95th percentiles of each row; tail FDR for the top 200; width = Σ(1 − FDR). **[ran]** ✗ as a width: about 11 for partner-less genes and 12 for small complexes (messy). It is stable in N on M3D (median 1.9, 2.8 and 2.6 at N = 116, 233 and 466).
6. **Rényi spacings / plateau detection.** On a sorted CLR row, a module shows as a plateau followed by a large normalized gap k·(x₍ₖ₎ − x₍ₖ₊₁₎). Width = the position of the largest gap. **[reasoned]** ✗: the oracle shows that even a perfect width does not help large complexes.
7. **Random matrix theory.** Take spikes above the Marchenko–Pastur edge (1 + √(G/N))², with width = the inverse participation ratio of the localized spike a gene sits in. **[ran, derived]** ✗.
   - BBP threshold: a clique is visible only if (n − 1)λ² > √(G/N). In K562 (G/N = 8) this needs λ > 0.84 at n = 5, 0.39 at n = 20 and 0.19 at n = 80.
   - Only 2–7 localized spikes per run; complexes ×0.45 (messy).
8. **Persistent homology / HDBSCAN stability.** Width = size of the most persistent cluster containing the gene. **[reasoned]** ✗ for regulons: hard units cost regulons ×0.30 in the asymmetry memo (C4b).
9. **Percolation.** Threshold the CLR graph at decreasing τ; a gene's width = the size of its component just before it joins the giant component. **[reasoned]** ✗, same failure as 8.
10. **Attention sink / "softmax off by one"** (Miller 2023; attention sinks, Xiao et al. 2023; registers, Darcet et al. 2024). Add a background token with logit bᵢ: wᵢⱼ = exp(zᵢⱼ/τ)/(exp(bᵢ/τ) + Σₖ exp(zᵢₖ/τ)). Mass on the sink stays home. **[ran as an oracle]** The one-sided version makes partner-less genes absorbing (small ×0.38). The symmetric version works (section 1). **Survivor S1.**
11. **Unbalanced optimal transport** (Chizat et al. 2018). Doubly-stochastic attention fails by mass conservation (asymmetry memo, C2), but KL-relaxed marginals let rows send less than 1: a principled, soft symmetric sink. **[reasoned]** A smooth version of S1.
12. **Reproducibility as background** (IDR; Li et al. 2011). What is not reproduced across lab halves is background. **[ran]**
    - As a candidate filter it fails: halving N loses the weak small complexes (small ×0.22).
    - As a row gate: see section 4.4.
13. **Reciprocity / mutual kNN** (hubness reduction; Radovanović et al. 2010). Keep i → j only if j → i; unreturned mass stays home. **[ran]** **Survivor S2.**
14. **Profile cosine readout.** cos(rowᵢ(Pᵗ), rowⱼ(Pᵗ)) is bounded and → 1 inside any mixed clique, whatever n (a topological-overlap cousin; Zhang & Horvath 2005). **[ran]** It works only once junk rows are silenced: ×0.60 on K10, but regulons ×4.26 on the oracle symmetric sink (messy).
15. **Width-rescaled readout** Pᵗᵢⱼ·keffᵢ ("undo the budget"). **[reasoned]** Scale-free only if keff is right, so it inherits the identifiability problem of point 2 in section 1. ✗
16. **LLR-calibrated temperature.** If CLR z were a log-likelihood ratio, the Bayes attention would have temperature 1 on the LLR scale, with no free keff. **[reasoned]** ✗: this is O1 in disguise, and O1 lost on all four endpoints.
17. **Ramanujan numerology.** A k-regular expander mixes an n-clique in about log n / log(k/(2√(k − 1))) steps: 8.6 for k = 10, n = 80, about the method's t = 9. It predicts that a wider keff needs a shallower t. **[derived]** ✗: the verdict shows K25 at t* = 8. Numerology.
18. **Effect-size deconvolution (Efron g-modeling).** Estimate each row's distribution g of true dependence strengths; width = g-mass above a *biological* effect size (|ρ| > 0.3). This makes the width N-invariant. **[reasoned]** **Survivor S5 (conceptual).**
19. **Copula-exact MI null per pair** (the B-spline bias depends on both genes' bin counts). **[reasoned]** A score correction; no width.
20. **"Background" as stationary hubness:** divide by πⱼ (the lift). **[reasoned]** ✗: already failed as a readout (K562 memo: 4.26×) and removes regulator signal (asymmetry memo).

## 3. Survivors, ranked

### S1. Symmetric background sink: partner-less genes neither send nor receive attention

- **Algorithm.** Given a row gate ρᵢ ∈ [0, 1] (the probability that gene i has a real partner):
  - B_ij = A_ij · ρᵢ · ρⱼ for i ≠ j, B_ii = 1 − Σⱼ B_ij, where A is the K10 softmax;
  - then the usual lazy walk and readout.
  - The soft version is unbalanced OT (idea 11).
- **Mechanism. [derived]** In the pooled readout, a within-clique pair scores about m/(n − 1). A partner-less row emits its full budget onto about 10 junk partners (about 1/10 each at t = 1). So roughly G_bg × 10 junk pairs compete with every complex pair. The sink removes them.
  - The one-sided sink (send nothing, still receive) fails: partner-less genes become absorbing states, and complex members that leak into them build large Pᵗ toward them (small complexes ×0.38 messy, ×0.53 clean).
  - Receiving must also be blocked: "X-by-Y" in the PI's sense.
- **Predictions (oracle gate, t = 9). [ran]**

  | Endpoint | messy | clean |
  |---|---|---|
  | n = 5 | ×3.41 | ×2.33 |
  | n = 20 | ×1.14 | ×0.96 |
  | n = 80 | ×1.54 | ×1.11 |
  | regulons | ×3.23 | ×2.25 |

  - The gain scales with the fraction of partner-less genes: 51% in this generator, by construction.
  - With the profile-cosine readout, regulons rise to ×4.26 messy, ×2.49 clean.
- **Realizable version. [ran]** A cross-half reproducibility gate (section 4.4), plus reciprocity, gives:

  | Endpoint | messy | clean |
  |---|---|---|
  | all complexes | ×1.08 | ×1.04 |
  | ribosome | ×1.31 | ×1.27 |
  | small complexes | ×1.71 | ×1.19 |
  | regulons | ×1.19 | ×1.40 |
  | n = 20 | ×0.82 | ×0.87 |

  Pooled complexes and regulons are 4/4 in both scenarios.
- **Failure modes.**
  - The detector, which no tested label-free method provides cleanly: the gate leaves 29% of partner-less genes open in messy and closes 22–31% of small-complex genes (section 4.4).
  - Real data with few partner-less genes. M3D has 0.3% of genes with no reproducible partner, so the expected E. coli gain is small. K562 has 10.1%.
  - Silencing a weakly detectable small complex removes its pairs entirely.
- **Cheap falsifier.** Vary the partner-less fraction (5%, 10%, 25%, 50%) with a realizable gate. If the gain at the fraction measured on K562 (about 10%) is under ×1.05 on complexes, drop S1.

### S2. Reciprocity: attention only where it is returned

- **Algorithm.** Keep A_ij only if A_ji > 0 (both genes are in each other's cap-50 softmax support). Removed mass stays on i. Then the same walk and readout. No parameter.
- **Predictions. [ran]**

  | Endpoint | messy | clean |
  |---|---|---|
  | n = 5 | ×3.45 (4/4) | ×2.19 (4/4) |
  | small complexes (3–6) | ×2.39 (4/4) | ×1.50 (4/4) |
  | n = 20 | ×0.77 (0/4) | ×0.86 (0/4) |
  | n = 80 | ×0.82 (0/4) | ×0.90 (0/4) |
  | regulons | ×1.04 (4/4) | ×0.95 (0/4) |
  | all complexes | ×0.93 | ×0.98 |

- **Mechanism.** Small cliques are mutual; junk attention from partner-less rows and satellites is mostly unreturned. It removes a little of what S1 removes, without needing a detector.
- **Failure modes.**
  - Regulators attend to their top targets only, so target → regulator edges outside that set are cut: −5% regulons in the clean world.
  - It does not help the ribosome.
- **Cheap falsifier.** A regulon-heavy generator (regulators with 150–500 targets, loose loadings) with the TF–target edge metric from the asymmetry memo. If edges or regulons drop by more than 5%, S2 can only be a complex-specific readout.

### S3. Scale-free readouts for the large-complex problem (a3r already exists; profile cosine is new)

- **Why.** Mass conservation means no width or budget change can lift the per-pair mass of an 80-gene complex above about 1/79. Only a readout whose scale does not depend on n can.
- **What the runs show. [ran]**
  - a3r on K10: all complexes ×2.03 messy and ×1.86 clean; ribosome ×3.88 and ×6.24; regulons ×1.14 and ×1.12 (all 4/4).
  - The profile cosine helps only after junk rows are silenced (S1).
  - a3r is fragile to its input scale. On rCLR input in the clean world it collapses on the ribosome (×0.02, 4/4).
- **Status.** a3r was preregistered and not adopted (E. coli complexes not confirmed), so it cannot be recycled as a "new" idea. The cosine is new but conditional on S1.
- **Cheap falsifier for the cosine.** Profile cosine on S1 with a realizable gate: if it does not beat a3r on the same operator in both scenarios, drop it.

### S4. Narrow only the rows that are in small modules

- **Oracle. [ran]** keff = min(10, n − 1) for complex members gives small ×3.72 (messy) and ×1.87 (clean), with the ribosome and regulons unchanged (0.88–1.01).
- **Assessment.** This is the only place where width itself pays.
  - A realizable version needs a "small module" detector, and small modules are exactly the weakest signals.
  - Reciprocity (S2) obtains about two thirds to four fifths of the oracle's small-complex gain (messy ×2.39 of ×3.72; clean ×1.50 of ×1.87) without a detector, so S4 adds little beyond S2.

### S5. N-invariant width by effect-size deconvolution (reasoned only)

- **Algorithm.**
  - For each row, deconvolve the observed Fisher-z of normal-score correlations into a prior g over true |ρ|, as in Efron's g-modeling on a log-spline grid, pooled over genes with a per-gene scale.
  - Width_i = (G − 1) · g_i(|ρ| > ρ_bio), with ρ_bio fixed in advance (for example 0.3).
- **Prediction.** It tracks module size, not N. In the synthetic data, the detection-based widths that should track module size grow with N for real modules (EVT width of ≥ 20-gene complexes 15, 24 and 28 at N = 125, 250 and 500) and not for partner-less genes (about 11 at every N). **[ran]**
- **Failure modes.**
  - Batch dependence above ρ_bio still counts.
  - By section 1 point 1, width does not fix large complexes anyway.
- **Use.** A principled descriptive statistic ("how many partners does this gene have?"). It is not an operator change.

## 4. Synthetic experiments that were run

### 4.1 Generator (`analysis/synthetic_width/sim_width.R`)

- **Size.** G = 2,000 genes, N = 500 samples, Gaussian factor model, B-spline MI (the package's estimator, bins "hg"), the package's Stouffer CLR and `.softmax_rows`. A per-gene keff vector uses the same bisection; it was checked identical at keff 10 (max abs difference < 10⁻⁸).
- **Complexes.** A CORUM-like size spectrum with loadings U(0.35, 0.85):
  - sizes 3 (×10), 4 (×8), 5 (×8), 6 (×5), 8 (×5), 10 (×4), 12 (×3), 15 (×2), 20 (×3), 30 (×2), 40 (×1);
  - one 80-gene ribosome (common factor 0.75 plus a 60S/40S subunit factor 0.35), holding about half of all complex pairs, as the RP complexes do in K562 CORUM.
- **Other structure.**
  - 25 satellites at 0.55 on the ribosome factor;
  - 12 regulators with 10–150 overlapping targets (target loadings U(0.2, 0.5));
  - 30 hubs (0.3 on 5 module factors);
  - a weak global factor;
  - about 1,016 partner-less genes.
- **messy scenario.** Adds:
  - 10 labs with gene-specific lab-effect SD U(0, 0.8);
  - 2 global programs with loadings N(0, 0.25²);
  - t₄ noise;
  - zero-inflation in 20% of genes.
- **clean scenario.** Gaussian noise only.
- **Fidelity check.** The generator reproduces the real pattern that motivates the question:
  - attention trails CLR on complexes (K10 104× vs CLR 127× messy; 147× vs 199× clean);
  - attention beats CLR on regulons (4.9× vs 3.0×; 13.3× vs 8.0×);
  - K25 beats K10 on complexes and is non-negative on regulons.
- **Order of work (disclosure).** One messy seed was run first with fewer operators. After inspecting it, the following were added before the 4-seed runs: the oracle decompositions, EVT, reciprocity, the reproducibility variants and the symmetric sink; the profile-cosine readout was also corrected. The rCLR and gate follow-ups were designed after seeing the 4-seed tables. All of this is synthetic, but the variants were not fixed in advance.
- **Seeds.** Simulation 20261030 + s, s = 1–4; permutation (O2) seed + 1000. Lab halves: labs 1–5 vs 6–10. N-subsampling (seed 1): column seed 7, permutation seed 8.

### 4.2 Main results (t = 9, attention readout unless stated; ratio to K10, seeds better)

| Variant | Complexes | n = 5 | n = 20 | n = 80 | small 3–6 | regulons |
|---|---|---|---|---|---|---|
| *messy* | | | | | | |
| CLR (pairwise, for reference) | 1.22 (4/4) | 4.72 | 1.62 | 1.72 | 2.56 | 0.61 (0/4) |
| K25 | 1.76 (4/4) | 0.22 | 1.68 | 3.73 | 0.23 | 1.05 (4/4) |
| O2 (permutation lfdr width) | 0.59 (0/4) | 0.50 | 0.23 | 0.73 | 0.62 | 0.97 (2/4) |
| EB graphon width (PI, literal) | 1.42 (4/4) | 0.00 | 0.00 | 8.26 | 0.00 | 0.48 (0/4) |
| EVT row-tail width | 0.76 (0/4) | 0.38 | 0.25 | 1.30 | 0.47 | 0.96 (0/4) |
| RMT (IPR) width | 0.45 (0/4) | 1.37 | 0.09 | 0.69 | 1.27 | 0.90 (0/4) |
| q-attention (weights = 1 − lfdr) | 0.04 (0/4) | 0.00 | 0.00 | 0.11 | 0.00 | 0.33 (0/4) |
| ORACLE width = n − 1 | 0.73 (0/4) | 7.66 | 0.29 | 1.04 | 3.54 | 0.96 (0/4) |
| ORACLE narrow small only | 0.96 (0/4) | 8.40 | 0.73 | 0.90 | 3.72 | 1.01 (4/4) |
| ORACLE partner-less rows keff 50 | 1.56 (4/4) | 2.15 | 1.98 | 2.23 | 1.80 | 1.34 (4/4) |
| ORACLE one-sided sink | 1.10 (4/4) | 0.37 | 0.71 | 1.44 | 0.38 | 1.84 (4/4) |
| **ORACLE symmetric sink (S1)** | **1.36 (4/4)** | **3.41** | 1.14 | **1.54** | **2.61** | **3.23 (4/4)** |
| **Reciprocity (S2)** | 0.93 (0/4) | **3.45** | 0.77 | 0.82 | **2.39** | 1.04 (4/4) |
| REPRO10 (reproducible candidates only) | 0.72 (0/4) | 0.26 | 0.58 | 0.63 | 0.22 | 0.97 (1/4) |
| K10, a3r readout | 2.03 (4/4) | 2.19 | 3.46 | 3.88 | 1.97 | 1.14 (4/4) |
| K10, profile-cosine readout | 0.60 (0/4) | 0.73 | 0.83 | 0.39 | 0.75 | 0.60 (0/4) |
| S1 oracle, profile-cosine readout | 1.44 (4/4) | 4.21 | 1.76 | 1.92 | 2.72 | 4.26 (4/4) |
| *clean* | | | | | | |
| CLR | 1.35 (4/4) | 2.26 | 1.43 | 2.24 | 1.46 | 0.60 (0/4) |
| K25 | 1.44 (4/4) | 0.20 | 2.16 | 3.52 | 0.28 | 1.64 (4/4) |
| O2 | 0.58 (0/4) | 0.20 | 0.26 | 0.63 | 0.34 | 1.19 (4/4) |
| EB graphon width | 1.17 (4/4) | 0.00 | 0.01 | 8.80 | 0.00 | 2.31 (4/4) |
| ORACLE width = n − 1 | 1.02 (3/4) | 2.86 | 0.53 | 1.60 | 1.84 | 0.93 (0/4) |
| ORACLE partner-less rows keff 50 | 1.15 (4/4) | 1.04 | 1.24 | 1.42 | 1.03 | 1.39 (4/4) |
| ORACLE one-sided sink | 1.00 (2/4) | 0.41 | 0.87 | 1.12 | 0.53 | 1.63 (4/4) |
| **ORACLE symmetric sink (S1)** | **1.10 (4/4)** | **2.33** | 0.96 | 1.11 | **1.79** | **2.25 (4/4)** |
| **Reciprocity (S2)** | 0.98 (0/4) | **2.19** | 0.86 | 0.90 | **1.50** | 0.95 (0/4) |
| K10, a3r readout | 1.86 (4/4) | 2.07 | 2.66 | 6.24 | 1.46 | 1.12 (4/4) |

Absolute values (mean (SD) over seeds) for every variant and readout are in `out/summary_<scenario>.csv`; per-seed values in `out/<scenario>/metrics_by_seed.csv`. The "n = 5/20/80" columns use the within-complex pairs of complexes of exactly that size as positives and all non-complex pairs as negatives: 8, 3 and 1 complexes per seed. So the n = 5 column is noisy (between-seed SD up to 80% of the mean), while the 4/4 counts are paired and more reliable.

**Mechanism check (messy, t = 9; `out/messy/mechanism_by_size.csv`). [ran]**

| Operator | per-pair mass n = 5 | n = 20 | n = 80 | in-complex mass n = 5 | n = 20 | n = 80 |
|---|---|---|---|---|---|---|
| K10 | 0.030 | 0.029 | 0.0114 | 0.13 | 0.57 | 0.94 |
| K25 | 0.008 | 0.014 | 0.0111 | 0.03 | 0.27 | 0.90 |
| ORACLE width = n − 1 | 0.071 | 0.018 | 0.0083 | 0.33 | 0.36 | 0.66 |
| ORACLE symmetric sink | 0.072 | 0.036 | 0.0116 | 0.34 | 0.74 | 0.96 |

- For n = 80, every operator gives m/(n − 1) to two digits (0.94/79 = 0.0119).
- The ribosome's gain under K25 or the sink comes entirely from the *competing* pairs, not from its own mass. This is the derivation in S1, measured.
- **Entropy bound. [derived]** With perplexity k over a cap of K = 50 candidates, a row whose module has n − 1 < k partners can keep at most m*(n, k) of its mass in the module, where h(m) + m·log(n − 1) + (1 − m)·log(K − n + 1) = log k:

  | Module size n | m* at k = 10 | m* at k = 25 |
  |---|---|---|
  | 3 | 0.69 | 0.42 |
  | 5 | 0.82 | 0.53 |
  | 8 | 0.94 | 0.65 |
  | 10 | 0.98 | 0.72 |

  This is the forced leak that narrow-width oracles remove, and why K25 loses small complexes.

**Widths the background models actually produce (messy, median per role; `out/messy/widths_by_role.csv`). [ran]**

| Role | O2 keff | graphon Σq | additive Σq | EVT Σq | reproducible top-50 candidates | oracle |
|---|---|---|---|---|---|---|
| partner-less | 2.0 | 421 | 393 | 10.7 | 2 (13% have 0) | – |
| complex 3–6 | 2.0 | 405 | 393 | 11.6 | 4 | 2–5 |
| complex 20–40 | 9.1 | 409 | 398 | 26.3 | 17 | 19–39 |
| ribosome | 10.8 | 441 | 429 | 41.2 | 36 | 79 |
| regulon target | 2.0 | 427 | 409 | 11.7 | 5 | – |

- Spearman correlation of width with a partner-less gene's lab-effect SD: graphon 0.11, additive two-way 0.29, EVT −0.01, O2 0.26.
- The additive background partly converts batch exposure into "partners".
- **Width against sample size (messy, seed 1; `out/messy/width_vs_N.csv`). [ran]** At N = 125, 250 and 500, EVT Σq is about 11 for partner-less genes at all N, but grows for ≥ 20-gene complexes (15, 24, 28). O2 keff stays at the floor of 2 except for the large complexes at N = 500 (7.4).

### 4.3 Follow-up: the background inside CLR itself (rCLR; `sim_rclr.R`, same seeds). [ran]

Building K10 on the robust two-way CLR (idea 2) instead of CLR:

| Endpoint | messy | clean |
|---|---|---|
| all complexes | ×1.02 (3/4) | ×0.88 (0/4) |
| n = 20 | ×1.39 | ×0.86 |
| regulons | ×1.08 (4/4) | ×0.84 (0/4) |

- Adding reciprocity in messy gives small ×1.96 and regulons ×1.09, but loses in clean.
- As a pairwise score, rCLR vs CLR:

  | Endpoint | messy | clean |
  |---|---|---|
  | all complexes | 229× vs 127× | 277× vs 199× |
  | ribosome | ×4.6 | ×3.9 |
  | small complexes | ×0.57 | ×0.55 |

  The mechanism: CLR's row SD is inflated by a gene's own module, which compresses large-module z-scores; a bulk-only (MAD) scale removes that. But this is a change to CLR, the baseline itself, and it is outside this question.
- **Not a survivor for attention** (inconsistent sign between scenarios).
- a3r on rCLR collapses on the ribosome in clean (×0.02, 4/4), a fragility worth knowing.

### 4.4 Follow-up: a realizable gate for S1 (`sim_gate.R`, same seeds). [ran]

**Gate.** ρᵢ = 1 if repro_nᵢ ≥ k₀, else 0, where:

- repro_nᵢ = the number of gene i's top-50 positive CLR candidates shared by the two lab halves;
- k₀ = the smallest k with P(Poisson(λ₀) ≥ k) < 0.05, λ₀ = 50²/(G − 1) the chance overlap; k₀ = 4 here, fixed by rule and not tuned.

Unlike REPRO10, an open row attends exactly as K10.

**How well the gate detects partner-less genes**

| Share of genes with the gate open | messy | clean |
|---|---|---|
| partner-less | 29% | 5% |
| complexes < 20 genes | 69% | 78% |
| complexes ≥ 20 | 91% | 96% |
| regulon targets | 63% | 68% |
| regulators | 94% | 96% |

So the gate is imperfect in both directions, and lab-driven spurious reproducibility keeps 29% of partner-less genes open in messy.

**Ratio to K10 (t = 9, attention readout)**

| Variant | Complexes | n = 5 | n = 20 | n = 80 | small 3–6 | regulons |
|---|---|---|---|---|---|---|
| *messy* | | | | | | |
| GATE symmetric sink | 1.08 (4/4) | 1.60 | 0.85 (1/4) | 1.31 (4/4) | 1.17 (3/4) | 1.10 (3/4) |
| GATE symmetric sink + reciprocity | 1.08 (4/4) | 2.55 | 0.82 (1/4) | 1.31 (4/4) | 1.71 (4/4) | 1.19 (4/4) |
| ORACLE symmetric sink | 1.36 (4/4) | 3.41 | 1.14 | 1.54 | 2.61 | 3.23 (4/4) |
| *clean* | | | | | | |
| GATE symmetric sink | 1.05 (4/4) | 1.80 | 0.92 (1/4) | 1.28 (4/4) | 1.17 (4/4) | 1.37 (4/4) |
| GATE symmetric sink + reciprocity | 1.04 (4/4) | 1.88 | 0.87 (0/4) | 1.27 (4/4) | 1.19 (4/4) | 1.40 (4/4) |
| ORACLE symmetric sink | 1.10 (4/4) | 2.33 | 0.96 | 1.11 | 1.79 | 2.25 (4/4) |

**Reading.**

- The realizable gate recovers about a quarter (messy) to a half (clean) of the oracle's pooled-complex gain, and a small part of its regulon gain.
- Combined with reciprocity it is the only label-free variant here that improves all complexes, the ribosome, small complexes and regulons together, in both scenarios, with 4/4 seeds on pooled complexes and regulons.
- Its consistent cost is 20-gene complexes (×0.82–0.92). The likely cause is that some 20-gene members whose partners are not reproducible across halves get silenced.
- On top of the a3r readout the gate adds nothing (messy all complexes ×1.99 vs a3r alone ×2.03). With reciprocity plus a3r it collapses small complexes in clean (×0.24, 3/4 worse). So the readout and the budget fixes are partly redundant, and a3r is fragile.
- **Caveat.** The generator has 51% partner-less genes by construction; K562 has about 10% by the label-free count, so real gains should be smaller.

### 4.5 Label-free diagnostic on real expression (`m3d_label_free.R`; no truths). [ran]

| Quantity per gene | M3D (G 4,297, N 466, 39 labs) | K562 essential half A (G 5,671, N 707) |
|---|---|---|
| graphon Σq, median [q10, q90] | 538 [379, 749]; 100% > 30 | 212 [83, 589]; 99.7% > 30 |
| EVT Σq, median [q10, q90] | 2.6 [0, 25]; 38% < 1 | 0 [0, 11]; 72% < 1 |
| EVT Σq median at 25% / 50% / 100% of N | 1.9 / 2.8 / 2.6 | 0 / 0 / 0 (q90: 6.8 / 9.1 / 11.1) |
| reproducible top-50 candidates across halves, median [q10, q90] | 14 [5, 29] | 6 [0, 23] |
| genes with 0 reproducible candidates | 0.3% | 10.1% |
| genes with ≤ 2 (chance overlap: 0.58 M3D, 0.44 K562) | 2.9% | 31.2% |

- M3D halves are one balanced lab split (seed 20261030). K562 has no labs, so its halves are a random split of the knockdown profiles (same seed); that is a weaker test of reproducibility.
- **What this says.**
  1. On both real datasets the PI's X-by-Y background declares hundreds of partners per gene, so no width can be read from it.
  2. Partner-less genes, the target of S1, are rare in M3D but common in K562 (about 10–30% depending on the threshold). So S1 has something to act on exactly where the complex loss was seen, and little in E. coli.
- This is a statement about co-expression structure only. Nothing was scored.

## 5. Recommended path

1. **Do not pursue "width from background counts" further.**
   - It is dominated, under mass conservation, by the budget lever.
   - The background counts are not identifiable under batch structure (synthetic and both real datasets).
   - O2, O2e and the graphon version all fail in the same way.
   - Record this in the supplement as an explained negative, together with O2's verdict.
2. **Next synthetic round (cheap, about 1 h single-thread).**
   - **S1 robustness:** partner-less fraction 5–50%; gate thresholds fixed by the Poisson chance-overlap rule, with no tuning; lab effects ×0, ×1, ×2.
   - **S2 on a regulon-heavy generator** with the TF–target edge metric.
   - **Profile cosine vs a3r** on S1.
   - **Decision rule written beforehand:** carry forward only variants that beat K10 on complexes *and* are at least ×0.95 on regulons in all three generators at the realistic partner-less fraction (about 10%).
3. **Label-free selection on real data (no truths).**
   - For whatever survives step 2, choose between {K10, S1-gate, S2, S1 + S2} by the existing held-out lab-split criteria (O3A/O3B) on M3D, EcoMAC and K562.
   - Commit the selections before any truth is touched, as was done for keff.
   - Expect the criteria to favour CLR-like operators; report which one is chosen.
4. **Confirmatory test: preregister, then touch a truth once.**
   - **Primary.** Pooled AUPR on complexes and on regulons (non-inferiority, margin 0.005) against K10 at t = 9.
     - Complex endpoint: a ribosome-capped truth (the largest complex at most 20% of positive pairs), with per-complex macro AP as key secondary.
     - Reason: K562 CORUM and E. coli complexes are both ribosome-dominated, and S1/S2 predict gains for small complexes, not the ribosome.
   - **Power.** Project the complex SE from truth structure alone, as the CLR(log Pᵗ) interpretation memo requires. The a3r test showed that about 54 E. coli complexes give SE ≈ 0.018, too few.
   - **Serial-testing disclosure.** This would be the third complex-oriented operator test after keff and CLR(log Pᵗ). State the order and treat it as the last.
5. **Held-out data that could make it honest.** There is no untouched human dataset, so in order of preference:
   - **(a)** A yeast compendium (M3D also distributes *S. cerevisiae*) with CYC2008 or Complex Portal yeast complexes and YEASTRACT regulons. New organism, many small complexes, never scored here. If the a3r follow-up also wants yeast, the two must share one preregistration so the truth is touched once.
   - **(b)** Other essential-gene Perturb-seq screens in new cell lines (HepG2, Jurkat; Nadig et al. 2025; check availability and licence) with CORUM. The expression is new, but CORUM's structure has been studied here, so use the ribosome-capped endpoint.
   - **(c)** RPE1 is not suitable as confirmatory: it has been scored descriptively for complexes (S5 of the CLR(log Pᵗ) test) and shares the K562 screen.
   - **(d)** K562 genome-wide (gwps) is not suitable: same cells and same truths.
   - E. coli (M3D, EcoMAC) can serve only as a no-harm regulon check: its truths have been scored repeatedly, and the label-free diagnostic predicts little for S1 there (0.3% partner-less genes).
6. **Things to specify for later, not to run now.** One descriptive real-data check would test the mechanism directly, but it touches CORUM, so it must wait for a preregistration or be declared descriptive: K25 vs K10 on *ribosome-excluded* K562 CORUM. The prediction from the synthetic data is that K25 loses there (small complexes ×0.23 here).

## 6. Open questions for the PI

1. **Which target matters for the paper?** Is it the ribosome (large complexes, where only a readout can help, and a3r was not confirmed) or small complexes and regulons (where S1 and S2 act)? The synthetic data say no operator in this family fixes both without a readout change.
2. **Is a change to the operator's budget acceptable?** S1 and S2 let genes keep their own mass, so the attention matrix is no longer a full softmax over partners. This weakens the "it is a transformer" framing slightly; the attention-sink literature makes it defensible.
3. **K562 split variable.** What is the closest analogue of a "lab" in the Perturb-seq pseudobulk: gem group, sequencing batch, knockdown programme? A random split of knockdowns (used here) is weaker for detecting reproducibility.
4. **How many more confirmatory shots?** After keff (no change) and CLR(log Pᵗ) (not adopted), a third complex-oriented test needs a stopping rule written in advance.
5. **Should rCLR be followed up?** It is a background-robust CLR (bulk-only scale, two-way) that lifts the ribosome ×3.9–4.6 as a pairwise score but costs small complexes. It is a change to the baseline, not to attention, and a separate question.

## 7. Files

- `analysis/synthetic_width/sim_width.R`: generator, background models (graphon, additive two-way, EVT, O2, RMT), operators, metrics. `Rscript sim_width.R 1:4 messy 1` (the last argument adds the N-scaling check).
- `analysis/synthetic_width/sim_rclr.R`: K10 built on rCLR.
- `analysis/synthetic_width/sim_gate.R`: realizable reproducibility gate for S1.
- `analysis/synthetic_width/summarize.R`: summaries and paired ratios.
- `analysis/synthetic_width/m3d_label_free.R [m3d|k562]`: the label-free real-data diagnostic.
- `analysis/synthetic_width/out/`: `summary_*.csv`, `paired_vs_K10_*.csv`, `<scenario>/{metrics_by_seed, widths_by_role, mechanism_by_size, width_vs_N, rclr_vs_K10, gate_vs_K10, gate_open_by_role}.csv`, `{m3d,k562}_label_free_*.csv` and run logs.
- Run with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19`, under `systemd-run --user --scope -p MemoryMax=4G`. Each sim_width seed takes about 170 s.
