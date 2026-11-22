# Spectral-embedding exploratory pass: pre-specified plan (21 Nov 2026)

*Written before any real score of any variant below exists. Exploratory; pre-specified, not registered externally. Code, smoke tests on synthetic truths and the controls of section 7 exist; no real truth has been scored by this code, and the sealed dataset (Jurkat) has not been touched.*

## 1. Question
The scGPT token-embedding cosine (a data-free similarity of pretrained gene embeddings) has a higher complex AUPR/base (about 3-6) than the CLR family (about 2) in all four datasets scored so far (`docs/plans/scgpt_extension_plan.md`, Revisions 1-3; `results/scgpt_ext/`). The PI asks, as preparation for adding pretraining later, whether structure in the spectral (principal-component) space of tcaCLR2 recovers cliques and complexes, with BIC choosing how many dimensions and how many components. Pretraining may be unreliable in biology, so nothing here presupposes that the spectral variants fall short of, or match, the token cosine. Whatever comes out is reported.

## 2. Data, units, truths (no new datasets, nothing changed)
- **Development set:** hESC, hHep, K562, HepG2, RPE1; scGPT gene set; units exactly as in `analysis/scgpt_ext/score_ext.R` (hESC: hESC_500, hESC_1000; hHep: hHep_500, hHep_1000; K562, HepG2, RPE1: subsamples 1-5). The expression matrix of each unit is rebuilt by the same construction; the pass refuses to score a unit whose permuted-data edge and t2 differ from `results/scgpt_ext/<dataset>/tcaclr2_edge_t2.csv` by more than 1e-9 (gate).
- **Sealed dataset:** Jurkat is the held-out test for the FROZEN variants only (section 8). It is not in the development set, no development code path opens it, and its job is clearly labelled and must not be submitted before the lead freezes the variants.
- **Truths (as `score_ext.R`, same code):** DoRothEA A-B and CollecTRI (directed TF-node `score_truth`; symmetric score used in both directions), Complex Portal complexes (pooled co-membership AUPR, `ap_fast`, pair universe of variant `main`, undirected), and, where a bulk file exists, KD response (secondary). Metric: AUPR / base rate ("AUPR ratio", column `auprc_ratio`). All of these truths have already been scored in this project for the comparison columns and for tcaCLR2 (a "used" truth); the spectral variants are new methods scored on used truths.
- **Comparison columns (read from `results/scgpt_ext/<dataset>/metrics_ext.csv`, not recomputed):** tcaCLR2, scGPT token cosine (`scgpt_tok`), scGPT attention (`scgpt_attn`); for context also aCLR, tcaCLR, CLR-NS, scGPT hidden cosine.

## 3. Embedding and read-outs (all settings fixed here)
**Operator and walk.** Per unit, on the CLR-NS matrix `ns_clr(X)`: the aCLR-style operator P = `gl_make_P` (keff 10, cap 50, lazy 0.5; the operator of tcaCLR and tcaCLR2), its degree-symmetrised form P_sym = `ad_Psym(P)` (symmetric, spectrum in [0, 1]), the top 61 eigenpairs by `ad_eigs`. The trivial eigenpair (value 1) is dropped. Walk depth t2 and the edge are those of tcaCLR2 (`ext_edge`: mean over 10 per-gene permutations of the second eigenvalue of the permuted-data P_sym; t2 = max(1, round(ln 0.05 / ln edge))).
**Embedding.** Y_k[, j] = lambda_j^(t2/2) v_j for the first k non-trivial eigenpairs (so Y_k Y_k^T is the rank-k part of P_sym^t2).
**k rules.** (a) edge rule: k = number of non-trivial eigenvalues, among the first 60, that exceed the permuted-data edge, capped at 50 (if fewer than 2, the edge-rule variants are not scored and the unit is reported as having no structure above the edge); (b) fixed grid k = 5, 10, 20, 50.
**Read-outs** (each gives a gene-by-gene pair score with zero diagonal; every seed is fixed in `sp_lib.R`, `SP$seed = 20261121`):
1. `cos`: cosine similarity of the rows of Y_k.
2. `cos_clr2`: `sk_c2z` (the existing second-CLR read-out: bilateral row standardisation, negatives to 0, Stouffer combination) applied to the cosine matrix.
3. `mclust`: `mclust::mclustBIC` on Y_k, G = 1, 2..30, all 14 covariance models. The scored model is the best BIC with G >= 2 (as specified); whether G = 1 had the better BIC is recorded (`g1_wins`), as is whether G hit 30 (`at_max`). Pair score = sum_k P(k|i) P(k|j) from the posterior of the scored model.
4. `hddc`: `HDclassif::hddc` on Y_k, K = 2..20, `model = "ALL"` (14 models), `d_select = "BIC"` (BIC chooses the number of components and, for every component, its intrinsic dimension), k-means start, 60 EM iterations; pair score as above from the posterior.
5. `pgmm`: `pgmm::pgmmEM` (mixtures of parsimonious Gaussian mixture models = factor analysers), G = 2..4, q = 1..3 (q is common to the components, limited by the Ledermann bound for the number of columns), all 8 families, k-means start, BIC; pair score as above from `zhat`. Narrower than the other mixtures because it is the slowest (3-5 s per fit at k = 10 on the smoke inputs, 72 fits per call); `at_max` is recorded. **pgmm is run only for k <= 30:** on HepG2 unit 1 the k = 50 call ran 43 min and returned no model, so `pgmm_k50` is listed in the ledger but not run (status `not_run`), and an edge-rule k above 30 is not run either. The family is 35 listed variants, 34 scored (fewer on a unit whose edge-rule k exceeds 30).
6. `ring` (family A, polar transform; in each plane the coordinates are centred eigenvector coordinates, centre = origin). **Planes:** all adjacent planes (j, j+1) of the weighted embedding, j = 1..k-1 (an unsupervised criterion; eigenvalue pairs are adjacent in this order; coordinate weighting equalises near-degenerate pairs). In a plane, r = sqrt(x^2 + y^2), theta = atan2(y, x). *Ring model* (H = 0..3 rings): one isotropic Gaussian bulk at the origin plus H rings, ring h with radial Gaussian N(r; mu_h, s_h), uniform angle, planar density N(r; mu_h, s_h) / (2 pi r), constraint s_h <= mu_h / 3 (thin ring) and effective size >= 8 genes; fitted by EM (radius-quantile start plus 3 seeded random starts); parameters 1 + 3H. *Cluster model:* mclust BIC on (x, y), G = 1..6, 14 models. BIC = 2 log L - npar log n, both on the planar likelihood. **A plane is accepted if the best ring model with H >= 1 exceeds the best cluster model by more than 10 BIC units.** Genes get ring responsibilities gamma_ih; pair score = sum over accepted planes (j <= k-1) and rings of gamma_ih gamma_jh (1 + cos(theta_i - theta_j)) / 2. No accepted plane: the score is identically 0 (AUPR ratio 1 under the tie convention) and the number of rings found is 0.
7. `rmix` (family B, mixture of radial Gaussians with a fitted centre; added on the PI's instruction). Same planes, same acceptance rule and pair score. Component m has planar density f_m(x) = exp(-(r - mu_m)^2 / (2 s_m^2)) / Z(mu_m, s_m), r the distance to the common centre c, Z = 2 pi [s^2 exp(-mu^2 / 2s^2) + mu s sqrt(2 pi) Phi(mu / s)], mu >= 0. This is a proper density on the plane for every (mu, s) (mu = 0 is an isotropic Gaussian blob, mu > 0 and small s a thin ring, mu > 0 and large s a fat doughnut, two components concentric rings) and equals the PI's N(r; mu, s) / (2 pi r) to first order for thin rings; **the literal form N(r; mu, s) / (2 pi r) is not integrable at r = 0 when the radial Gaussian does not vanish there, so the normalised form is used** (recorded for the PI). M = 1 or 2 components with positive weights; parameters counted honestly: 2 (centre) + 2M (mu, s) + (M - 1) (weights) = 3M + 1; fitted by L-BFGS-B on the exact likelihood from two starts (centroid; best of 6,000 random three-point circles). A component is a *ring component* if mu >= s and its effective size is >= 8 genes. The plane is accepted if the best model with at least one ring component exceeds the plane's mclust BIC (the same one as family A) by more than 10 units. The angle theta is taken about the fitted centre.
**Validator (not a read-out, never enters a pair score).** In every plane the largest H1 bar (death - birth) of the Vietoris-Rips filtration (`TDAstats::calculate_homology`, dim 1; points scaled by the plane's RMS radius) is reported next to the per-gene permuted-data null of the same unit (same pipeline on the planes of the permuted embedding, seed 90001) and the flag "exceeds the maximum over all null planes". Disagreement between BIC acceptance and the validator is reported, not resolved here.

## 4. Variants and ledger (M = 35 listed, 34 scored)
The family is 5 embeddings x 7 read-outs = 35 listed variants (`pgmm_k50` is not run, see section 3), each scored on every unit and every truth of the development set. **No variant is added or dropped after scoring starts; no tuning on any score.** The table below is the ledger: the definition is fixed now; the result columns are appended by the lead after the run (development ratios against the comparators, section 5; frozen yes/no; Jurkat ratio only for frozen variants).

| ID | Name | Embedding | Read-out | Development result (to be appended) | Frozen? | Jurkat (once, frozen only) |
|---|---|---|---|---|---|---|
| V01 | cos_kE | k = number of non-trivial eigenvalues above the permuted-data edge (cap 50) | cosine of the rows of Y_k (diagonal 0) | | | |
| V02 | cos_clr2_kE | k = number of non-trivial eigenvalues above the permuted-data edge (cap 50) | sk_c2z (the second CLR, bilateral row standardisation of stack_lib.R) applied to the cosine matrix of the same Y_k | | | |
| V03 | mclust_kE | k = number of non-trivial eigenvalues above the permuted-data edge (cap 50) | mclust, BIC over G = 2..30 and all 14 covariance models on Y_k (G = 1 also fitted and recorded); pair score = sum_k P(k|i) P(k|j) | | | |
| V04 | hddc_kE | k = number of non-trivial eigenvalues above the permuted-data edge (cap 50) | HDDC (HDclassif::hddc), K = 2..20, all 14 models, BIC for K and for each component's intrinsic dimension; pair score = sum_k P(k|i) P(k|j) | | | |
| V05 | pgmm_kE | k = number of non-trivial eigenvalues above the permuted-data edge (cap 50) | mixture of parsimonious factor analysers (pgmm::pgmmEM), G = 2..4, q = 1..3 (common q), all 8 families, BIC; pair score = sum_k P(k|i) P(k|j) | | | |
| V06 | ring_kE | k = number of non-trivial eigenvalues above the permuted-data edge (cap 50) | ring family A (polar transform, bulk + thin rings at the origin, BIC against mclust) over planes (j, j+1), j <= k-1; accepted planes summed; pair score = sum gamma_ih gamma_jh (1 + cos(theta_i - theta_j)) / 2 | | | |
| V07 | rmix_kE | k = number of non-trivial eigenvalues above the permuted-data edge (cap 50) | ring family B (mixture of 1-2 radial Gaussians with fitted centre, BIC against mclust) over planes (j, j+1), j <= k-1; same pair score, angle about the fitted centre | | | |
| V08 | cos_k5 | k = 5 | cosine of the rows of Y_k (diagonal 0) | | | |
| V09 | cos_clr2_k5 | k = 5 | sk_c2z (the second CLR, bilateral row standardisation of stack_lib.R) applied to the cosine matrix of the same Y_k | | | |
| V10 | mclust_k5 | k = 5 | mclust, BIC over G = 2..30 and all 14 covariance models on Y_k (G = 1 also fitted and recorded); pair score = sum_k P(k|i) P(k|j) | | | |
| V11 | hddc_k5 | k = 5 | HDDC (HDclassif::hddc), K = 2..20, all 14 models, BIC for K and for each component's intrinsic dimension; pair score = sum_k P(k|i) P(k|j) | | | |
| V12 | pgmm_k5 | k = 5 | mixture of parsimonious factor analysers (pgmm::pgmmEM), G = 2..4, q = 1..3 (common q), all 8 families, BIC; pair score = sum_k P(k|i) P(k|j) | | | |
| V13 | ring_k5 | k = 5 | ring family A (polar transform, bulk + thin rings at the origin, BIC against mclust) over planes (j, j+1), j <= k-1; accepted planes summed; pair score = sum gamma_ih gamma_jh (1 + cos(theta_i - theta_j)) / 2 | | | |
| V14 | rmix_k5 | k = 5 | ring family B (mixture of 1-2 radial Gaussians with fitted centre, BIC against mclust) over planes (j, j+1), j <= k-1; same pair score, angle about the fitted centre | | | |
| V15 | cos_k10 | k = 10 | cosine of the rows of Y_k (diagonal 0) | | | |
| V16 | cos_clr2_k10 | k = 10 | sk_c2z (the second CLR, bilateral row standardisation of stack_lib.R) applied to the cosine matrix of the same Y_k | | | |
| V17 | mclust_k10 | k = 10 | mclust, BIC over G = 2..30 and all 14 covariance models on Y_k (G = 1 also fitted and recorded); pair score = sum_k P(k|i) P(k|j) | | | |
| V18 | hddc_k10 | k = 10 | HDDC (HDclassif::hddc), K = 2..20, all 14 models, BIC for K and for each component's intrinsic dimension; pair score = sum_k P(k|i) P(k|j) | | | |
| V19 | pgmm_k10 | k = 10 | mixture of parsimonious factor analysers (pgmm::pgmmEM), G = 2..4, q = 1..3 (common q), all 8 families, BIC; pair score = sum_k P(k|i) P(k|j) | | | |
| V20 | ring_k10 | k = 10 | ring family A (polar transform, bulk + thin rings at the origin, BIC against mclust) over planes (j, j+1), j <= k-1; accepted planes summed; pair score = sum gamma_ih gamma_jh (1 + cos(theta_i - theta_j)) / 2 | | | |
| V21 | rmix_k10 | k = 10 | ring family B (mixture of 1-2 radial Gaussians with fitted centre, BIC against mclust) over planes (j, j+1), j <= k-1; same pair score, angle about the fitted centre | | | |
| V22 | cos_k20 | k = 20 | cosine of the rows of Y_k (diagonal 0) | | | |
| V23 | cos_clr2_k20 | k = 20 | sk_c2z (the second CLR, bilateral row standardisation of stack_lib.R) applied to the cosine matrix of the same Y_k | | | |
| V24 | mclust_k20 | k = 20 | mclust, BIC over G = 2..30 and all 14 covariance models on Y_k (G = 1 also fitted and recorded); pair score = sum_k P(k|i) P(k|j) | | | |
| V25 | hddc_k20 | k = 20 | HDDC (HDclassif::hddc), K = 2..20, all 14 models, BIC for K and for each component's intrinsic dimension; pair score = sum_k P(k|i) P(k|j) | | | |
| V26 | pgmm_k20 | k = 20 | mixture of parsimonious factor analysers (pgmm::pgmmEM), G = 2..4, q = 1..3 (common q), all 8 families, BIC; pair score = sum_k P(k|i) P(k|j) | | | |
| V27 | ring_k20 | k = 20 | ring family A (polar transform, bulk + thin rings at the origin, BIC against mclust) over planes (j, j+1), j <= k-1; accepted planes summed; pair score = sum gamma_ih gamma_jh (1 + cos(theta_i - theta_j)) / 2 | | | |
| V28 | rmix_k20 | k = 20 | ring family B (mixture of 1-2 radial Gaussians with fitted centre, BIC against mclust) over planes (j, j+1), j <= k-1; same pair score, angle about the fitted centre | | | |
| V29 | cos_k50 | k = 50 | cosine of the rows of Y_k (diagonal 0) | | | |
| V30 | cos_clr2_k50 | k = 50 | sk_c2z (the second CLR, bilateral row standardisation of stack_lib.R) applied to the cosine matrix of the same Y_k | | | |
| V31 | mclust_k50 | k = 50 | mclust, BIC over G = 2..30 and all 14 covariance models on Y_k (G = 1 also fitted and recorded); pair score = sum_k P(k|i) P(k|j) | | | |
| V32 | hddc_k50 | k = 50 | HDDC (HDclassif::hddc), K = 2..20, all 14 models, BIC for K and for each component's intrinsic dimension; pair score = sum_k P(k|i) P(k|j) | | | |
| V33 | pgmm_k50 | k = 50 | NOT RUN (pgmm limited to k <= 30; see section 3) | n/a | n/a | n/a |
| V34 | ring_k50 | k = 50 | ring family A (polar transform, bulk + thin rings at the origin, BIC against mclust) over planes (j, j+1), j <= k-1; accepted planes summed; pair score = sum gamma_ih gamma_jh (1 + cos(theta_i - theta_j)) / 2 | | | |
| V35 | rmix_k50 | k = 50 | ring family B (mixture of 1-2 radial Gaussians with fitted centre, BIC against mclust) over planes (j, j+1), j <= k-1; same pair score, angle about the fitted centre | | | |

## 5. Success criterion (descriptive; no threshold is a pass/fail claim)
For a variant v, a development dataset d and a truth T, with units u of d as in `score_ext.R`:
ratio_tok(v; d, T) = mean_u AUPR/base(v; d, T, u) / mean_u AUPR/base(scGPT token cosine; d, T, u),
ratio_tcaCLR2(v; d, T) = mean_u AUPR/base(v; d, T, u) / mean_u AUPR/base(tcaCLR2; d, T, u),
per truth (DoRothEA A-B, CollecTRI, Complexes; KD where it exists) and per dataset, printed for all scored variants (`results/spectral_pass/summary/ratio_by_dataset_truth.csv`); the same ratio against scGPT attention is printed as a third column. A convenience column, the mean over units and the three primary truths of log2 of the per-unit ratio (the d_i convention of `stats_ext.R`), is in `ratio_primary_mean_log2.csv`. A ratio of 1 means equal, above 1 the spectral variant is higher; no value is a threshold, an equivalence margin or a decision. The number of units in each cell (2 for hESC and hHep, 5 for K562, HepG2, RPE1) is printed and a cell with 2 units is not an interval.

## 6. Multiplicity and claims
This is a family of M = 34 scored variants (35 listed) scored on used truths. Descriptive only: no test, no interval across datasets, no Holm family, no claim that any variant works or that anything transfers. Large ratios among 34 variants and 4 truths are expected by selection alone; the table is shown whole, ordered by name, not by result. Per-variant ratios on the development set are used by the lead to freeze a subset of variants (the freezing choice is a data-dependent selection and is stated as such). The sealed dataset Jurkat is the held-out test for the FROZEN variants only: scored once, after the freeze, with the frozen list written into this ledger first. A statement about transfer is allowed only for the frozen variants on Jurkat, and only descriptively (Jurkat is one dataset with 5 pseudo-replicate units; no interval).

## 7. Controls (run before the real scores; outputs in `results/spectral_pass/controls/`)
- **Plane tests** (`sp_controls.R --mode planes`): thin ring, fat doughnut, two concentric rings with and without a Gaussian bulk, a ring off-centre from a bulk, three Gaussian blobs, one blob, standard normal noise, a uniform disc. Reading in advance: both families should accept thin rings, family B should accept fat doughnuts, both should reject blobs and noise; the uniform disc and the off-centre ring beside a bulk are known misspecifications (a common centre; the uniform disc is not Gaussian-radial) and are reported as found. The validator is shown for each.
- **Positive control** (`--mode positive`): planted 3 Gaussian clusters (40 genes, 4 SD from the origin), 2 subspace clusters (40 genes, SD 3 in a random 2-D subspace, 0.3 elsewhere, centred at the origin) and one thin ring (50 genes, radius 5 SD, radial SD 0.25 SD, coordinates 7-8) in (a) a synthetic 12-D Gaussian embedding and (b) the real spectral embedding of a unit (HepG2 unit 1), by replacing the rows of 250 randomly chosen genes, planted coordinates scaled by the real column SDs. Every variant is scored by pooled pair AUPR/base on the planted co-membership (all pairs), overall and per structure type; recovery diagnostics (G, K, rings, plane, validator) are reported.
- **Negative control** (`--mode negative`): per-gene permuted expression of a real unit (seed 90001), the same pipeline with its own edge; reported: eigenvalues above the edge, k by the edge rule, components found by mclust (G, whether G = 1 had the better BIC), HDDC, pgmm, and the number of rings found by each family, at the edge-rule k and at forced k = 5, 10, 20. The expectation stated in advance is that no structure is reported at the edge rule; at forced k the mixtures may still return components (G >= 2 is imposed on mclust), which is why `g1_wins` and the ring counts are reported and the result is read as found.
Controls are for checking the code and for calibrating how to read the ring and mixture read-outs; they do not tune a variant.
**Control outcomes already seen when this plan was finalised (disclosed; they were produced by the smoke runs of the code, with no real truth):** (i) plane tests: thin ring, fat doughnut and concentric rings without a bulk are accepted by both families (concentric rings plus a bulk need three components and family B, limited to two, does not accept them); blobs and noise are rejected by both; the off-centre ring beside a bulk is not accepted by either family (the validator flags it); the uniform disc is accepted by both families (validator at the null level). (ii) negative control (HepG2 unit 1, hESC_500): 0 eigenvalues above the edge, so no edge-rule variant is scored there; at forced k = 5, 10, 20 mclust, HDDC and pgmm return components (G = 2-17, K = 3-19, pgmm at its limit G = 4; mclust's G = 1 never had the better BIC), family A accepts 0-5 of 19 planes, family B 0-1. (iii) positive control: in the synthetic 12-D embedding both ring families accept only the plane of the planted ring (plane 7; validator 0.61 against a null maximum 0.31), plus one family-A plane at the cluster block; in the real HepG2 embedding with the same structure planted, neither ring family accepts any plane although the validator does flag plane 7 (0.94 against 0.41), while mclust, HDDC and pgmm recover the planted pairs at the edge rule and at k = 10. These facts are used only to read the real scores (a ring read-out reporting no ring in a real embedding is weak evidence of absence); no setting was changed after seeing them except those listed in section 9.

## 8. Order of work and freeze
Commit this plan and the code; run the five development jobs (`analysis/spectral_pass/jobs/7900-7904`), the controls job and the summary job; the lead appends the ledger and decides which variants to freeze (listing them in this document); only then the Jurkat job (7990) may be submitted. Outputs of the development run go to `results/spectral_pass/<dataset>/`; the Jurkat run to `results/spectral_pass/Jurkat/` and is invoked with an explicit flag (`--heldout YES`) so that it cannot run by accident.

## 9. Implementation notes and limitations stated in advance
- Settings fixed or changed during the build, before any real score: mclust G upper limit 20 -> 30 and HDDC K 10 -> 20 (BIC reached the old limits on a real embedding), pgmm G 2-6 and q 1-5 -> G 2-4 and q 1-3 (run time), ring family B added and the H1 validator added on the PI's instruction (below).
- **PI instructions recorded before any real score (21 Nov 2026):** (1) R packages princurve, circular, movMF, dbscan, flexmix, TDAstats, ripserr are installed; TDAstats (H1 persistent homology) is the ring validator in the planted-ring control and in the real planes, against the bars of per-gene permuted data; princurve and dbscan are optional comparators and are not used. (2) A second ring family, a mixture of radial Gaussians with a fitted centre (read-out 7, `rmix`), with its BIC comparable to mclust on the same plane, parameters counted honestly, tested on thin rings, fat doughnuts, Gaussian blobs and permuted data (section 7).
- No second operator or scorer: operator `gl_make_P`, symmetriser `ad_Psym`, eigen-solver `ad_eigs`, edge `ext_edge`, second CLR `sk_c2z`, scorers `score_truth` and `ext_score_complex` (`ap_fast`). New code: `analysis/spectral_pass/` only.
- The pair scores from mixtures are posterior overlaps; many pairs score near 0 and ties follow the registered `ap_fast` convention.
- A mixture with G at the upper limit (`at_max`) means BIC wanted more components than allowed; it is reported, not extended after the fact.
- The ring families assume a centred plane and a common centre per plane; a ring displaced from a bulk, or a non-Gaussian radial profile (a uniform disc), is not modelled (section 7); the validator is the check.
- The edge-rule count is capped at 50 and taken among the first 60 eigenvalues.
- Packages (user library, offline): mclust 6.1.4, HDclassif 2.2.2, pgmm 1.2.8 (installed from CRAN by this pass); the PI's additional princurve, circular, movMF, dbscan, flexmix, TDAstats, ripserr are installed (TDAstats is used by the validator; the others are not used).


## 10. Development results and freeze (22 Nov 2026; written before the sealed dataset is scored)

All 34 scored variants were run on the five development datasets (hESC, hHep, K562, HepG2, RPE1; jobs 7900-7904, controls 7905-7906, summary 7910, all rc = 0). Descriptive only. Mean over the five datasets of the per-dataset mean log2 ratio (variant AUPR/base over comparator AUPR/base, over units and the three primary truths); negative means the variant scored lower than the comparator:

| Variant | vs scGPT token cosine | vs tcaCLR2 | vs scGPT attention |
|---|---|---|---|
| cos_kE | -0.69 | -0.18 | +0.82 |
| cos_k20 | -0.70 | -0.18 | +0.82 |
| cos_k50 | -0.70 | -0.19 | +0.81 |
| cos_k10 | -0.71 | -0.19 | +0.81 |
| cos_clr2_k10 | -0.72 | -0.20 | +0.80 |
| cos_clr2_kE | -0.77 | -0.25 | +0.75 |
| cos_k5 | -0.78 | -0.27 | +0.73 |
| cos_clr2_k5 | -0.79 | -0.28 | +0.72 |
| cos_clr2_k20 | -0.80 | -0.28 | +0.71 |
| mclust_k5 | -0.80 | -0.28 | +0.71 |
| mclust_k10 | -0.84 | -0.32 | +0.68 |
| hddc_k10 | -0.85 | -0.33 | +0.66 |
| cos_clr2_k50 | -0.85 | -0.33 | +0.67 |
| pgmm_k5 | -0.86 | -0.34 | +0.66 |
| hddc_k5 | -0.86 | -0.34 | +0.66 |
| hddc_kE | -0.92 | -0.40 | +0.60 |
| hddc_k20 | -0.92 | -0.40 | +0.60 |
| mclust_kE | -0.94 | -0.42 | +0.58 |
| pgmm_k10 | -0.94 | -0.43 | +0.57 |
| mclust_k20 | -0.95 | -0.43 | +0.57 |
| pgmm_kE | -0.98 | -0.47 | +0.53 |
| pgmm_k20 | -1.01 | -0.50 | +0.50 |
| hddc_k50 | -1.04 | -0.52 | +0.48 |
| mclust_k50 | -1.08 | -0.57 | +0.43 |
| ring_k50 | -1.28 | -0.76 | +0.24 |
| rmix_k10 | -1.30 | -0.79 | +0.21 |
| rmix_k50 | -1.30 | -0.79 | +0.21 |
| rmix_k5 | -1.30 | -0.79 | +0.21 |
| rmix_k20 | -1.30 | -0.79 | +0.21 |
| ring_k10 | -1.30 | -0.79 | +0.21 |
| ring_kE | -1.30 | -0.79 | +0.21 |
| ring_k5 | -1.30 | -0.79 | +0.21 |
| rmix_kE | -1.30 | -0.79 | +0.21 |
| ring_k20 | -1.31 | -0.79 | +0.21 |

Per truth (mean over the five datasets, AUPR/base): Complexes: best variant (plain cosine) 2.07, tcaCLR2 1.68, tcaCLR 1.86, CLR-NS 2.07, scGPT token cosine 4.95. DoRothEA A-B: 1.42, 1.75, 1.55, 1.32, 2.01. CollecTRI: 1.30, 1.86, 1.66, 1.38, 1.61. The ring families (polar-transform and radial-mixture) accepted no ring in any real plane; their scores equal the no-ring floor. The mixture families (mclust, HDDC, pgmm) scored below plain cosine, more so with more eigenmodes. Controls: planted structure recovered by the mixtures and not by the ring families on the real HepG2 embedding; permuted data give 0 eigenvalues above the permutation edge, so edge-rule variants are not scored there; forced-k mixtures still return components.

**Freeze.** The frozen variants for the sealed dataset are cos_kE and cos_k10 (plain cosine similarity in the top-k eigenmodes of the tcaCLR2 walk, k by the permutation-edge rule and k = 10), the two best development variants; tcaCLR2 and the scGPT readouts are the comparators and are scored by the scGPT extension plan. Jurkat is scored once for the frozen variants, descriptively; no variant is added or dropped after this entry. Decided by the PI on 22 Nov 2026 ('freeze and score Jurkat'; the case for pretraining is left open in the paper).
