# Asymmetric attention (columns as regulators or latents, rows as observed genes): expert assessment of the PI's proposal

*26 October 2026. Expert memo (network inference, latent-variable models, diffusion on graphs, attention; Claude Opus) for the PI. **Synthetic data only.** No real expression data, truth set or result file was read or scored (no K562, RPE1, E. coli, CORUM, EcoCyc or Complex Portal; nothing under `results/`). Scripts and small outputs: `analysis/synthetic_asymmetry/`. Compute: 12 simulated datasets, about 3.5 min in total, nice 19, one thread, in a 6 GB memory scope.*

## 1. Verdict

**Partial pass.** The proposal as stated fails: asymmetry does not reduce the tension that makes tight cliques fail. One component passes for a different target: aggregating columns into latents gives a principled operator that helps direct regulator–target edges.

1. **Fails as a fix for complexes.** Three findings, from theory and simulation:
   - The clique loss comes from mass conservation. Any stochastic operator, whether it normalizes rows, columns or both (Sinkhorn), spreads a clique's mass over its n exchangeable members, about m/(n−1) per pair. Column normalization leaves this unchanged (section 3.1). In the simulation, the per-pair mass of the doubly-stochastic walk at t = 9 is 0.0115 in the 80-gene clique, against 0.0122 for the plain walk.
   - Column-side normalization in the CLR (Faith 2007) form (C1) ties a3r on complexes (paired ratio 1.01 in the base scenario, 0.97 with global programs). It costs regulons: 0.94 and 0.86 of a3r, worse in 8/8 paired seeds. The reason is that the column z-score is a standardized log lift, log(Pᵗ/π), and it removes the hub signal that regulators carry.
   - Doubly-stochastic attention (C2) costs regulons more: 0.74 to 0.85 of a3r.
2. **Column normalization does not separate hubs from cliques.** The hub AUROC is the probability that a within-module pair outranks a pair made of a hub and a module member. It is 0.781 for a3r, 0.770 for C1 and 0.740 for C2 (base). With global programs it is 0.749, 0.714 and 0.695.
3. **Hard order-of-operations variants fail for a mechanistic reason.** Residual attention after removing an estimated clique factor (C4a) loses most of the signal, and ranking within-unit pairs first (C4b) keeps 0.89 of a3r on complexes but 0.30 on regulons. An estimated clique factor is not orthogonal to the programs the clique sits in: units at r ≥ 0.6 also capture targets that share two regulators. Regressing those factors out deletes the regulon programs (within-regulon |r| 0.20 → 0.07).
4. **Direction is not identifiable from co-expression.** The asymmetric scores orient an edge toward the gene that is the best single proxy of the shared latent (centrality or purity), not toward the regulator. When the regulator's mRNA loading is drawn like its targets' (control scenario), CLR, attention and a3r point the wrong way (target→TF read 0.14 to 0.23), while C3 still reads 0.70 to 0.79 because, in this generator, regulators are single-program genes. This confirms the first reaction that imposing direction on undirected latent structure is not meaningful.
5. **Passes, for another endpoint: latent-column CLR (C3).** Each column j becomes a latent: the attention-weighted mean of the genes that attend to j, leaving out the gene being scored. Gene i is then scored against that latent with CLR's row and column z-scores (section 2).
   - Each member inherits the latent's dependence instead of dividing it: the gene–latent correlation is flat in clique size, at about λ rather than λ².
   - With aggregation through P⁴, C3 beats a3r on direct regulator–target edges in 12/12 paired seeds, by ×1.15, ×1.12 and ×1.25 in the base, global and matched-loading scenarios. It beats CLR by ×1.06, ×1.12 and ×1.32. The gain is largest when the regulator's mRNA is a weak proxy of its activity, as the factor model predicts.
   - It is not a complex fix: about 0.68 of a3r on complexes, because its final CLR step dilutes as √(G/n). On regulon co-membership it is mixed (×1.12 base, ×0.94 with global programs).

   So C3 is a candidate *edge-level* readout next to a3r, for the paper's stated weak point (attention "is a module ranker, not a better inferrer of direct regulatory edges"). It is not a replacement for a3r. A follow-up plan is in section 6.

The three criteria:

| Criterion | C1 column-CLR | C2 Sinkhorn | C3 latent-column CLR | C4 residual / units |
|---|---|---|---|---|
| Principled operator | yes (Faith 2007 gene × regulator CLR on attention logits) | yes (doubly-stochastic attention; Sinkformer) | yes (factor model: rows observed, columns latent) | yes (order of operations) |
| Mechanism supported for cliques | no (dilution unchanged; ties a3r) | no (dilution unchanged) | no (√(G/n) dilution in the CLR step) | no (removes the clique signal or the programs) |
| Mechanism supported for edges | no | no | **yes** (12/12 seeds; scales with a weak TF proxy) | no |
| Clean test design possible | – | – | **yes**, after the CLR(log Pᵗ) verdict, on edge truths not used to build C3 | – |

## 2. What is already asymmetric, and the four candidate operators

**The current operator is already directed.** The softmax makes row i a budget over keys j, and Pᵗ is row-stochastic; the readout symmetrizes only at the end. The a3r readout z-scores each row of log Pᵗ against that gene's own background (over reachable pairs), clips at 0 and combines the two directions with CLR's Stouffer rule: S = (Z + Zᵀ)/√2, where Z_ij = z_row i(log Pᵗ_ij). In a3r, the reverse direction of pair (i, j) is judged in row j, against j's background.

The PI's proposal treats columns as regulators or latent factors and rows as observed genes, as in CLR's original gene × TF form (the target is z-scored over the regulators, the regulator over all genes). Four operators make this concrete. All use the method's CLR, softmax (keff 10, cap 50) and lazy walk (α = 0.5) unchanged.

| Code | Operator | Directed score D_ij (i observed row, j column/latent) | Symmetric readout | How it differs from a3r |
|---|---|---|---|---|
| **C1** a3c | CLR, gene × regulator form, on attention logits | (z_row i(L_ij) + z_col j(L_ij))/√2, L = log Pᵗ, reachable pairs, negatives clipped | (D + Dᵀ)/2 | Each entry is also judged against the column's background (how much everyone reaches j), not only the row's. Since log Pᵗ_kj → log π_j for distant k, the column z is a standardized log lift log(Pᵗ_ij/π_j): an in-degree (hub) correction. |
| **C2** Sinkhorn | Doubly-stochastic attention | P_ds = D_u · ½(I + A) · D_v, with rows and columns summing to 1 (the positive diagonal guarantees convergence; error 4e-16) | (P_dsᵗ + P_dsᵗᵀ)/2 and a3r on P_dsᵗ | Every key receives the same total attention, so π is uniform and over-attended columns are scaled down before the walk. |
| **C3** latent-column CLR | Factor model: rows observed, columns latent | f_j = Σ_k W_kj x_k with W = Pᵗᵃ column-normalized (j's attenders, plus j itself through the lazy diagonal). R_ij = cor(x_i, f_j^(−i)), leaving gene i out in closed form. D = Stouffer of the row z (gene i over latents) and the column z (latent j over genes) of Gaussian MI(R). t_a ∈ {1, 4} | (D + Dᵀ)/2 | Not a readout of attention mass at all. Attention only defines the latent of each column; the score is a dependence on that latent, so it has no budget. Its ablation without aggregation (W = I) is plain CLR. |
| **C4** residual | Order of operations | Tight units (average linkage on 1 − r, r ≥ 0.6, ≥ 3 genes); unit factors (means of members) regressed out of every gene; MI, CLR, softmax and walk on residuals. C4b ranks within-unit pairs first (by CLR) | C4a: attention on residuals; C4b: units first, then C4a | Hard latent removal before attending. C4b matches the K562 memo's variants (d) and (b0). |

## 3. Analysis: clique, regulator, hub, satellite

The notation follows a factor model: x_i = Σ_f b_if F_f + e_i. A **clique** is n genes with loading λ on one factor (pairwise r = λ²). A **regulon** is a regulator with mRNA loading λ_R on its activity a_R, plus m targets with loose loadings b (target–target r = b², regulator–target r = λ_R·b). A **hub** loads weakly on several modules. A **satellite** loads on a clique's factor without being a member (INTS8 or EIF4B on the 60S ribosome).

### 3.1 Mass conservation: column normalization cannot remove dilution

- **Lemma.** Let C be n genes exchangeable under the operator, meaning a permutation of C leaves it invariant. Then for any row-stochastic Pᵗ, the within-C block is constant off the diagonal, and each member's pair weight is m_i/(n − 1), where m_i ≤ 1 is the mass retained in C.
- **Scalings keep the bound.** A column scaling D_v, or the Sinkhorn scaling D_u K D_v, is constant on C by the same symmetry, so the block stays constant and the bound still holds. Normalizing rows, columns or both changes how much mass enters or leaves C, never how it is shared inside it.
- **Simulation (base, mean of 4 seeds).**
  - Softmax pair weight: 0.241, 0.161, 0.103, 0.052 and 0.0255 for n = 4, 6, 10, 20 and 40 (λ = 0.6), and 0.0126 for the 80-gene ribosome-like clique. This is exactly m/(n − 1), with in-clique mass m = 0.72, 0.80, 0.92, 0.98 and 0.99 (0.999 for the ribosome), so 1/(n − 1) for n ≥ 20. It is the same at λ = 0.9 (0.252 … 0.0256).
  - The softmax therefore also erases coupling strength. Pairs at r = 0.36 and r = 0.81 receive the same share.
  - At t = 9 the plain walk gives 0.0122 per ribosome pair and the doubly-stochastic walk 0.0115.
- **Small cliques leak.** At keff = 10, a row with fewer than 10 partners must place mass outside its clique, so the perplexity floor forces a leak.
  - In-clique mass at t = 1 is 0.72 (n = 4), 0.80 to 0.86 (n = 6) and ≥ 0.92 (n ≥ 10).
  - At t = 9 the walk has spread most of a small clique's mass outside it: 0.23 to 0.26 is retained for n = 4, against ≥ 0.93 for n ≥ 20.
  - So small complexes lose to leakage and large ones to dilution. Per-pair mass at t = 9 never exceeds about 0.09, whatever the coupling.
- **Consequence.** The fix must make the readout's scale independent of n. A different normalization of the operator cannot do it.

### 3.2 Linear-domain z-scores dilute as √(G/n); log-domain z-scores as log n

- **Linear domain.** Suppose row i holds n − 1 entries at a and G − n at b, with noise small. Then the member z-score is z = (1 − p)/√(p(1 − p)) = √((1 − p)/p), with p = (n − 1)/(G − 1). This is independent of the signal strength a − b.
  - For CLR (G = 1,000), the bilateral pair score is √2·√((G − n)/(n − 1)): 25.8 (n = 4), 7.0 (n = 40) and 4.8 (n = 80).
  - Simulated: 25.5, 7.0 and 4.7. CLR's own complex score is therefore set by clique size, not by coupling: 21.7 at λ = 0.6 vs 25.5 at λ = 0.9 for n = 4.
  - This is why |r| beats CLR on complexes here (72.3× vs 43.7×), as on K562 (7.00× vs 6.39×).
- **Log domain (a3r).** The clique sits at log(m/n) against a background of reached-but-unrelated genes that spans orders of magnitude. The 1/n becomes an additive −log n against a wide background SD.
  - The a3r pair score is 5.4 to 8.3 for every n from 4 to 40 and every λ, with no trend in n below 40 (ribosome 4.4). This is the "inherit, don't divide" behaviour of the K562 memo.
- **C3 and C1.** C3's final CLR step is a linear-domain z over latents, so it inherits CLR's √(G/n) dilution. C1's row half is a3r's, so it keeps a3r's behaviour; its column half adds the lift (section 3.3).

### 3.3 Column normalization penalizes regulators, because regulators are the true hubs

- **C1's column z is the lift.** Column j of log Pᵗ is shifted up by log π_j for every row the walk has mixed: log Pᵗ_kj ≈ log π_j. Column standardization therefore scores (log Pᵗ_ij − log π_j)/σ_j, a standardized log lift.
- **The lift has already failed.** On K562 it failed as a readout (memo table: CLR on log lift 4.26×; diffusion density 1.35×). On M3D, hubness carries real regulon signal: π_i + π_j alone gives AUROC 0.546 (paper).
- **Regulators are high in-degree columns.** Their m targets attend to them. The in-degree correction removes exactly that, so C1 costs regulons: ×0.94 base, ×0.86 with global programs (0/4 seeds better in each).
- **Sinkhorn (C2) does the same, more strongly.** It forces every column sum to 1, so a regulator attended by 100 targets is scaled to the same total as a gene nobody attends: regulons ×0.85 and ×0.74 of a3r.
- **What C2 does fix.** It fixes the one-sided satellite problem for the *probability* readout: satellite AUROC 0.996 vs 0.798 for attention. An over-attended clique column is scaled down and the satellite's budget moves elsewhere. But a3r already fixes that (1.000) through the reverse row.

### 3.4 Satellites and hubs

- **Satellites (one-sided attention).** A satellite's row puts most of its budget on the clique; the members' rows put nothing on it.
  - The probability readout averages the two directions, so satellite–member pairs compete with member–member pairs: attention's satellite AUROC is 0.80 (base) and 0.65 (global).
  - The bilateral log readout requires both rows to find the pair surprising: a3r scores 1.00 (base) and 0.75 (global). With strong global programs no attention variant fully separates satellites. C3 and CLR do somewhat better (0.80 and 0.78), and units-first (C4b) best (0.998), because it hard-codes the clique.
- **Hubs.** A hub that loads 0.3 on six modules has r ≈ 0.3λ with module members. Its attention is spread, it is rarely in members' top-10, and the walk does not favour it over true partners.
  - a3r separates within-module pairs from hub–member pairs best (AUROC 0.781 base, 0.749 global).
  - The column-normalized operators do worse (C1 0.770 and 0.714; C2 0.740 and 0.695), and CLR worse still (0.672 and 0.643).
  - The hypothesis that column normalization separates hubs from cliques is not supported. In-degree in attention is not where module hubs show up. They have modest in-degree; regulators and large-module members have high in-degree.

### 3.5 Latent columns (C3): inherit, don't divide, and denoise the regulator

- **The latent.** For a clique member j, W's column j holds about 1/n weight on each attending member, so f_j ≈ the clique mean. Its correlation with the factor is ρ_f with ρ_f² = nλ²/(nλ² + 1 − λ²) → 1. Hence R_ij ≈ λρ_f → λ, whereas gene–gene r = λ².
- **Simulated gene–latent correlation:** 0.42 to 0.47 at λ = 0.6, 0.63 to 0.67 at λ = 0.75 and 0.85 to 0.87 at λ = 0.9 (ribosome 0.80), flat in n. This is the soft version of the PI's latent-variable picture (PI decision 3: invariance to the number of noisy copies), obtained without collapsing anything.
- **Regulator–target pairs.** f_R averages the regulator with the targets that attend to it, so R_target,R ≈ b·ρ_f(a_R) instead of b·λ_R.
- **Prediction.** The gain over gene–gene scoring should grow as the regulator's mRNA becomes a weaker proxy of its activity (smaller λ_R). The simulation agrees:
  - The edge gain over CLR is ×1.06 when λ_R ∈ [0.5, 0.8] and ×1.32 when λ_R is drawn like the targets' [0.25, 0.6].
  - This is the same logic as regulator-activity inference from targets (network component analysis; VIPER).
- **Costs.**
  1. The CLR step's √(G/n) dilution (section 3.2): complexes ×0.68 of a3r.
  2. Under strong global programs, the attention-weighted latent also averages the members' global loadings: regulons ×0.94 of a3r.
- **Leave-one-out is necessary.** Without it, a gene that puts a large share of its attention on j would partly define f_j and correlate with itself.

### 3.6 Residualizing a clique factor (C4)

- **C4a** removes what defines the clique, so within-clique pairs lose their signal: complexes 9.9× vs a3r 74.5×.
- **Reinstating units first (C4b)** gives the hard-collapse behaviour: complexes 66.2× (×0.89) with regulons ×0.30.
- **Why regulons collapse.** Data-driven units are not pure cliques. With N = 400, targets of two regulators at loadings near 0.6 reach r ≥ 0.6 and form units (24 targets in units in seed 1). Their factor is a mixture of regulator activities, and removing it from every gene removes those programs from all their targets: within-regulon |r| 0.20 → 0.07, and the |r| regulon AUPR 7.8× → 2.3× (seed 1).
- **Relation to K562.** This is the "leak into the wider program" in its hard form, and the synthetic counterpart of the K562 results: residual attention 2.25×; units first 6.30×; coarse-graining costs DoRothEA −0.013.

## 4. Synthetic test

**Generator** (`sim_asymmetry.R`): G = 1,000 genes, N = 400 samples, Gaussian factor model, unit-variance genes.

- 15 cliques: sizes 4, 6, 10, 20 and 40, each at loadings 0.60, 0.75 and 0.90.
- One ribosome-like clique of 80 genes at 0.85.
- 25 satellites (0.6 on the ribosome factor, not members).
- 10 regulators: mRNA loading 0.5–0.8 on their own activity, 10–150 targets each, drawn from a pool of 380 (so regulons overlap), target loadings U(0.25, 0.60).
- 20 module hubs, each loading 0.30 on six random module factors.
- A weak global factor (0.20 on 30% of genes); the rest is background.

**Scenarios:**

- `base`;
- `global`: two strong global programs added, with gene-specific loadings N(0, 0.25²). This is the per-gene background structure CLR is designed to remove.
- `tf_matched`: base, with the regulator mRNA loading drawn like the targets', U(0.25, 0.60). This is the control for the direction metric, and the same seeds as base.

Four seeds per scenario; every operator is compared on the same simulated data (paired).

**Dependence and depth:**

- Dependence is Gaussian MI, −½ log(1 − r²), which is monotone in |r|; B-spline MI is unnecessary for Gaussian data.
- CLR is the package's `clr_calibrate` (Stouffer); the softmax is the package's `.softmax_rows` (keff 10, cap 50).
- t = 9 is primary (the method's transferred depth); t = 4 is also written out. No depth was selected on truth.

**Metrics:** pooled AUPR over all G(G − 1)/2 pairs, divided by the base rate (×):

- **complexes:** pairs in the same clique;
- **no-ribosome:** the same, with pairs touching the 80-gene clique removed;
- **regulons:** regulator-agnostic co-membership, regulator included;
- **TF–target edges:** universe = pairs with exactly one regulator gene;
- **satellite AUROC:** ribosome–ribosome pairs against satellite–ribosome pairs;
- **hub AUROC:** within-module pairs against (hub, member) pairs, for the modules the hub loads on;
- **target→TF:** the share of true edges with D[target, TF] > D[TF, target].

Per-clique AP saturates at 1.00 for nearly every clique and operator (N = 400 separates every planted clique within its own rows), so the differences below are about calibration across rows, as on K562.

### 4.1 Results at t = 9 (means of 4 seeds)

**Base**

| operator | complexes × | no-ribosome × | regulons × | TF–target edges × | satellite AUROC | hub AUROC | target→TF |
|---|---|---|---|---|---|---|---|
| \|r\| | 72.3 | 108.3 | 8.00 | 15.31 | 1.000 | 0.657 | – |
| CLR | 43.7 | 81.3 | 5.54 | 13.42 | 0.998 | 0.672 | 0.45 |
| attention | 42.9 | 80.5 | 5.80 | 11.78 | 0.798 | 0.768 | 0.65 |
| **a3r** | **74.5** | **134.2** | 6.23 | 12.40 | 1.000 | **0.781** | 0.61 |
| C1 column-CLR on log Pᵗ | 75.1 | 134.3 | 5.86 | 12.50 | 1.000 | 0.770 | 0.59 |
| C2 Sinkhorn, attention readout | 49.5 | 93.6 | 5.05 | 12.12 | 0.996 | 0.758 | 0.59 |
| C2 Sinkhorn, a3r readout | 77.5 | 134.4 | 5.29 | 12.68 | 1.000 | 0.740 | 0.52 |
| C3 latent-column CLR, t_a = 1 | 46.3 | 86.5 | 6.07 | 14.07 | 1.000 | 0.689 | 0.80 |
| **C3 latent-column CLR, t_a = 4** | 50.8 | 95.1 | **6.98** | **14.24** | 1.000 | 0.722 | 0.65 |
| C4a residual attention | 9.9 | 30.1 | 2.25 | 5.72 | 0.506 | 0.681 | 0.64 |
| C4b units first | 66.2 | 83.2 | 1.87 | 5.72 | 1.000 | 0.723 | – |

**Global programs**

| operator | complexes × | no-ribosome × | regulons × | TF–target edges × | satellite AUROC | hub AUROC | target→TF |
|---|---|---|---|---|---|---|---|
| \|r\| | 67.9 | 96.8 | 4.17 | 11.17 | 0.939 | 0.646 | – |
| CLR | 39.3 | 70.9 | 4.11 | 11.14 | 0.782 | 0.643 | 0.42 |
| attention | 32.0 | 63.9 | 5.79 | 10.47 | 0.645 | 0.735 | 0.60 |
| **a3r** | **66.1** | **127.8** | **5.87** | 11.15 | 0.753 | **0.749** | 0.58 |
| C1 column-CLR on log Pᵗ | 63.8 | 124.7 | 5.07 | 10.86 | 0.672 | 0.714 | 0.57 |
| C2 Sinkhorn, attention readout | 37.5 | 74.4 | 5.15 | 11.34 | 0.691 | 0.709 | 0.60 |
| C2 Sinkhorn, a3r readout | 51.8 | 99.2 | 4.36 | 10.81 | 0.648 | 0.695 | 0.53 |
| C3 latent-column CLR, t_a = 1 | 41.9 | 77.0 | 4.62 | 11.82 | 0.809 | 0.658 | 0.64 |
| **C3 latent-column CLR, t_a = 4** | 45.2 | 85.2 | 5.52 | **12.43** | 0.804 | 0.691 | 0.59 |
| C4a residual attention | 6.5 | 16.9 | 2.07 | 5.42 | 0.595 | 0.678 | 0.60 |
| C4b units first | 62.5 | 74.0 | 1.69 | 5.42 | 0.998 | 0.715 | – |

**Regulator loading matched to targets (direction control):**

- target→TF: CLR 0.14, attention 0.21, a3r 0.23, C1 0.21, C2 0.52 to 0.58, C3 0.70 and 0.79.
- TF–target edges: CLR 9.91, attention 10.29, a3r 10.53, C3 (t_a = 4) 13.13, \|r\| 13.50.
- The complex, regulon, satellite and hub columns are within 1% of base, as expected (only the regulator loading changed).

**Paired ratio to a3r** (mean over seeds; seeds better than a3r; base / global / tf_matched):

| operator | complexes | regulons | TF–target edges |
|---|---|---|---|
| CLR | 0.59 (0/4) / 0.59 (0/4) / 0.59 (0/4) | 0.89 (0/4) / 0.70 (0/4) / 0.88 (0/4) | 1.08 (4/4) / 1.00 (2/4) / 0.94 (1/4) |
| attention | 0.58 (0/4) / 0.48 (0/4) / 0.58 (0/4) | 0.93 (0/4) / 0.99 (1/4) / 0.93 (0/4) | 0.95 (0/4) / 0.94 (0/4) / 0.98 (1/4) |
| C1 | 1.01 (4/4) / 0.97 (0/4) / 1.01 (4/4) | 0.94 (0/4) / 0.86 (0/4) / 0.94 (0/4) | 1.01 (3/4) / 0.98 (1/4) / 1.02 (4/4) |
| C2, a3r readout | 1.04 (4/4) / 0.78 (0/4) / 1.04 (4/4) | 0.85 (0/4) / 0.74 (0/4) / 0.84 (0/4) | 1.02 (3/4) / 0.97 (1/4) / 1.03 (3/4) |
| C3, t_a = 4 | 0.68 (0/4) / 0.68 (0/4) / 0.69 (0/4) | 1.12 (4/4) / 0.94 (1/4) / 1.12 (4/4) | **1.15 (4/4) / 1.12 (4/4) / 1.25 (4/4)** |
| C4b units first | 0.89 (0/4) / 0.95 (1/4) / 0.89 (0/4) | 0.30 (0/4) / 0.29 (0/4) / 0.30 (0/4) | 0.46 (0/4) / 0.49 (0/4) / 0.41 (0/4) |

Full tables: `out/summary_<scenario>.csv`, `out/paired_vs_a3r.csv`, `out/<scenario>/metrics_by_seed.csv`. Mechanism numbers by clique: `out/<scenario>/mechanism_by_clique.csv`.

### 4.2 What the simulation does and does not show

- **It reproduces the K562 failure and its rescue.**
  - Attention ties CLR on complexes (base 42.9× vs 43.7×) or loses to it (global 32.0× vs 39.3×). It also ranks satellites into the clique (satellite AUROC 0.80).
  - a3r roughly doubles the complex score over attention (×1.74 base, ×2.07 global) and slightly raises regulons and edges.
  - The generator is therefore a fair model of the mechanism the asymmetry was meant to address.
- **A linear Gaussian world favours |r|.** It wins on complexes, and on regulons in base, because nothing in the generator makes CLR's per-gene background normalization pay as it does on real compendia (heterogeneous noise, platform effects, non-linear dependences). Read the tables as *relative* evidence among the attention operators and their readouts, not as absolute rankings against |r| or CLR.
- **Small scale.** Four seeds per scenario, one generator family, one N, and loadings chosen by the author. The paired ratios are tight (seed ranges in `paired_vs_a3r.csv`: most within ±0.05 of the mean, edges up to ±0.15), but they are synthetic.

## 5. Bearing on the CORUM order-of-operations question

The PI's earlier question was whether to collapse cliques into latents before attending, or to attend on residuals after removing a clique's shared factor. The asymmetry proposal is the same move expressed through columns. The synthetic evidence and the K562 memo now agree:

1. **Hard removal or collapse costs regulons.** Data-driven units are mixtures of programs (C4, K562 (b) and (d)).
2. **Soft column latents (C3) inherit without dividing,** but a linear-domain CLR readout re-introduces the √(G/n) dilution, so they do not fix complexes.
3. **What fixes complexes is a scale-free readout of the attention walk:** the log-domain bilateral z (a3r), which is already preregistered. No asymmetric operator beats it on complexes without a regulon cost: C1 ties it and costs regulons ×0.94; C2 gains ×1.04 in base, loses ×0.78 under global programs and costs regulons ×0.74 to 0.85.

For the expert item, the open order-of-operations question reduces to the edge-level role of soft latents (C3). No complex-level operator is left to test beyond the CLR(log Pᵗ) preregistration.

## 6. Follow-up plan: latent-column CLR as an edge-level readout (the part that passes)

*Not to be started until the CLR(log Pᵗ) verdict is written (`verdict.csv`, truth phase 0910–0915). No C3 score is computed on any real dataset before its preregistration is committed.*

**Operator (to fix in the preregistration):**

- Upstream: the method's operator, unchanged and fixed at the K10 operator: B-spline MI → Stouffer CLR → softmax (keff 10, cap 50) → P = ½(I + A).
- Latents: W = P^{t_a}, column-normalized; f_j = Σ_k W_kj x̃_k, where x̃ are normal scores (Blom) of each gene.
- Dependence: R_ij = cor(x̃_i, f_j^(−i)) in closed form (leave i out), with MI = −½ log(1 − R²) (Gaussian-copula MI).
  - This keeps the closed form, and CLR on normal-score |r| already behaves like CLR on MI on M3D (paper, Table 1).
  - B-spline MI on latents is a descriptive secondary, because leave-one-out makes each pair's latent different, so it needs G² separate estimates.
- Directed score: D_ij = (z_row i + z_col j)/√2 with negatives clipped.
- Readouts: (D + Dᵀ)/2 for pairs; D[target, TF] restricted to TF columns as a descriptive directed form. No direction claims (section 3, finding 4).
- Aggregation depth: t_a chosen label-free on the existing lab splits with criterion A (AP for the other half's top-1% CLR pairs), grid {1, 2, 4, 9}; committed before truths, flagged if t_a = 1.

**Endpoints and truths:**

- **Primary (confirmatory): direct regulator→target edges**, E. coli M3D and EcoMAC, RegulonDB SC.
  - Universe: (TF, gene) pairs for RegulonDB TFs with ≥ 1 SC target in the gene universe; autoregulation excluded.
  - Comparators: CLR (the Faith 2007 edge method; superiority) and attention at t = 9 (superiority). Also a3r if the CLR(log Pᵗ) verdict adopts it (descriptive, or a third primary if adopted).
  - The intersection–union rule follows the CLR(log Pᵗ) preregistration.
- **No-harm:** RegulonDB SC co-membership against the main method, margin −0.005. Expect a pass in base-like data and a small loss where global programs dominate (synthetic ×0.94 of a3r). Preregister the margin, not the expectation.
- **Complexes: descriptive only**, with the expected direction stated in advance (C3 below a3r, synthetic ×0.68). C3 is not a complex method and must not be adopted for complexes.
- **Key secondary:** RPE1 CollecTRI edges, only after the CLR(log Pᵗ) truth phase for RPE1 is complete. BEELINE edge truths (cell-type ChIP-seq) as descriptive.
- K562 is in-sample for about 20 operators (descriptive only, never confirmatory). C3 has never been computed on it; keep it so until the E. coli scoring is done.

**Truth hygiene and disclosure:**

- RegulonDB edges underlie the co-membership truths already scored many times for CLR and attention; the edge form has not been scored for any method as far as this memo knows.
  - Before preregistering, check `results/` *by file name only* for edge-level metrics (the planned DREAM5-style directed evaluation has not been run) and record the finding in the disclosure table.
  - Merge with that planned directed evaluation (with ARACNE) into one preregistration, so that the edge truth is touched once for all methods.
- **Order of discovery, to disclose:**
  1. K562 CORUM loss;
  2. exploratory a3r on K562 (22 Oct);
  3. CLR(log Pᵗ) preregistered (22 Oct);
  4. PI's asymmetry proposal (26 Oct), evaluated on synthetic data only (this memo; operators, t_a ∈ {1, 4} and leave-one-out chosen on synthetic data);
  5. C3 found to be an edge-level, not a complex-level, readout;
  6. preregistration before any real-data computation of C3.

**Interactions:**

- **CLR(log Pᵗ) preregistration:** untouched. C3 runs only after its verdict; the comparator set includes whichever readout the verdict makes the method. No CLR(log Pᵗ) endpoint is re-scored.
- **keff preregistration:** untouched. C3 is defined on the keff 10 operator, and no width or multi-head variant is combined with C3. If the keff test later changes the method's width, C3 keeps keff 10 in its confirmatory run, and any wider variant is descriptive.

**Power.**

- The synthetic edge effect is ×1.06 to ×1.32 over CLR (median about ×1.12) and ×1.12 to ×1.25 over a3r.
- Edge AUPRs are small (base rates under 1%), so an absolute smallest effect of interest of 0.005 is too coarse. Preregister a relative one instead: 10% of CLR's edge AUPR, tested on the log ratio with the paired gene jackknife and the M3D cluster bootstrap.
- Estimate the paired SE beforehand, without real scores: plant RegulonDB's out-degree distribution (truth structure only, as the CLR(log Pᵗ) evaluability counts did) into this generator at M3D dimensions (G 4,297, N 466), and run the jackknife on the synthetic truth.
- The risk is that the real regulator mRNA is a strong proxy for most RegulonDB TFs. Then the factor model predicts a small gain (synthetic ×1.06 over CLR at λ_R 0.5–0.8). Say in advance that a null result is informative.

**Compute:** one extra G × G product and elementwise operations per dataset (seconds to minutes); the jackknife as in the CLR(log Pᵗ) script. Jobs go through `jobs_pc/sched` after the truth phase finishes.

## 7. Files

- `analysis/synthetic_asymmetry/sim_asymmetry.R`: generator, operators C1–C4, metrics. Usage: `Rscript sim_asymmetry.R 1:4 <base|global|tf_matched>`, run with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19`.
- `analysis/synthetic_asymmetry/summarize.R`: summaries and paired ratios.
- `analysis/synthetic_asymmetry/out/`: summaries (≈150 KB in total) and run logs (`run_<scenario>.log`).

The simulation seeds are 20261026 + s, for s = 1–4. The base and tf_matched scenarios share seeds and differ only in the regulator loading.
