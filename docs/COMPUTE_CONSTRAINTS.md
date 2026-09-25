# Compute-constrained design decisions

*Started 2026-09-24.  Every run so far has been done on a MacBook Pro (M2 Max, 12 cores, 32 GB).  Several choices below were made to fit that machine or to avoid long waits, not because they are the best design.  Once the dual-socket Xeon workstation (28 physical cores, 128 GB, Ubuntu) is available around the clock, the plan is to revert to the full versions listed here, all at once.  Each entry says what we did, why, and what the full version is.*

Setup on the workstation:  `tools/setup_linux.sh`, then `tools/bench_ht.sh` to decide whether hyper-threading helps (see D28 in `CITATION_LOG.md`).

## A. Resampling depth and replicate counts

1. **Chips bootstrap draws.**
   - **Now:**  40 cluster and 20 replicate draws, being extended on the Mac to 200 and 100.
   - **Why:**  each draw recomputes MI, CLR and attention for 4,297 genes and takes minutes.
   - **Full version:**  at least 1,000 cluster draws, so that the 2.5% and 97.5% percentiles rest on 25 draws each, plus 200 replicate draws.

2. **Depth is fixed inside every bootstrap and perturbation draw.**
   - **Now:**  the chips bootstrap scores attention at fixed depths (t = 4, 8 and 16 for `soft10`; 6 for `pearson_top10`).  The perturbation influence uses a fixed t = 8.  The held-out depth search runs once, on the full data.
   - **Why:**  re-running the held-out lab-split depth search inside every draw multiplies the cost by roughly 20–40 fits per draw.
   - **Consequence:**  the intervals do not include the uncertainty of choosing t, so they are somewhat optimistic for "attention at the selected depth".
   - **Full version:**  a nested bootstrap, in which each draw re-selects t on its own lab split and is scored at its own t*.

3. **One held-out lab split.**
   - **Now:**  one balanced split of the 39 labs into two halves, used in both directions, so t* is a geometric mean of two values.  BEELINE uses one random halving of cells per case.
   - **Why:**  each split costs two full fits plus a depth scan.
   - **Full version:**  20–50 random lab splits (and cell splits for BEELINE), reporting the distribution of t*, choosing t by the median, and checking how flat the criterion is around the optimum.

4. **Depth search by geometric scan and golden section.**
   - **Now:**  t = 1, 2, 4, …, 128, then golden-section refinement, chosen to minimise the number of fits.
   - **Full version:**  evaluate every integer t up to 64 and show the whole curve, which costs little once fits are cheap.

5. **Permutation null size.**
   - **Now:**  B = 100 permutations for the main fit, and B = 50 (`--Btrain`) for the operators built on training halves.
   - **Why:**  each permutation is a full MI + CLR computation.
   - **Full version:**  B = 500–1,000 for the full-data null.  Also, BH q-values near the threshold are resolved by roughly 1/(B × number of pairs), so re-check the edge sets at q = 0.05 and 0.20 with the larger null.

6. **The FDR-thresholded operator is left out of the bootstrap.**
   - **Now:**  `fdr05_top10` and the `fdr20` variants appear only in the point estimates and in the gene jackknife.
   - **Why:**  they need a permutation null inside every draw.
   - **Full version:**  include them in the 907-chip bootstrap, with B = 100 per draw.

7. **Perturbation null.**
   - **Now:**  100 random delete-k subsets per k, shared by every perturbation with the same k (the first run used 30).  Subsets are drawn uniformly from all experiments.
   - **Full version:**  500 or more per k.  Also a *matched* null, in which the random subsets are drawn from the same lab or the same condition type as the perturbation's experiments.  That separates "this perturbation matters" from "removing any k experiments from this lab matters".  It costs about one network fit per null draw, and the matched version needs its own null for each perturbation.

8. **Uncertainty for BEELINE and the TF-node benchmark.**
   - **Now:**  point estimates only.
   - **Full version:**  bootstrap cells (BEELINE) or experiments (M3D) with MI recomputed in every draw, as in the chips bootstrap, to give intervals on the differences between attention and CLR.

## B. Problem size

9. **Gene context for BEELINE.**
   - **Now:**  CLR and attention are built only on BEELINE's selected genes (TFs + 500 or TFs + 1,000, that is 514–1,528 genes, 30–60% of them TFs).
   - **Why:**  that is BEELINE's protocol, and it is cheap.  Widening it was deferred because a dense 15,000-gene matrix is 1.8 GB, and a single fit holds about ten of them.
   - **Full version:**  build CLR and the operator on all expressed genes, and score only the BEELINE gene set.  This gives each gene's CLR background a representative context instead of a thin, TF-heavy one.  The intermediate step is the 5,000 most variable genes.

10. **Softmax candidate cap.**
    - **Now:**  each row's softmax is over its top 50 CLR scores (`softmax_cap = 50`), with the temperature set so that exp(entropy) = 10.
    - **Why:**  partly for sparsity (a sparse P keeps the attention products cheap), partly as a modelling choice.
    - **Full version:**  test caps of 100, 200 and no cap (full softmax), to confirm that the cap does not drive the results.

11. **Worker count limited by RAM.**
    - **Now:**  on the Mac, 7 workers × 1 OpenMP thread (3.3 GB per chips worker, 6 GB left to the user), so about 3 of the 10 usable cores sit idle.
    - **Full version:**  on the Xeon, about 26–32 workers (see `tools/bench_ht.sh`), with no change to the code.  Only `CLR_HT` and `CLR_RESERVE_RAM_GB` are set in the environment.

## C. Analyses deferred for time

12. **Graded ChIP-seq truth.**  Score BEELINE ChIP benchmarks by the rank correlation between each TF's CLR row and its ChIP-Atlas binding scores (MACS2, per experiment, cell-type matched), instead of 0/1 edges.  This is deferred for the download size (per-TF tables for hundreds of TFs), not the compute.
13. **Evidence-filtered STRING.**  Rebuild the STRING truth from the channel-level files, keeping curated-database and experimental evidence and dropping edges supported only by co-expression, text mining or genome context.  If the current STRING version publishes directed transcriptional-regulation actions, also build a TF-to-target subset from them.
14. **SSEM-Lasso re-run.**  Re-run the Cosgrove et al. (2008) network filter on the same perturbations, with the same truth and gene universe, for a like-for-like comparison of top-100 sensitivity.
15. **Depth diagnostics on BEELINE.**  Score attention at t = 1, 2, 4 and 8 against each truth, as a diagnostic only, never for choosing t.

## D. What is *not* a compute compromise

These choices were made on scientific grounds and should stay unless the evidence changes:

- The 466 replicate-averaged experiments as the primary set, to avoid pseudo-replication.
- HG common bin counts on raw values.
- The CLR-score permutation null (rather than an MI null).
- Stouffer combination.
- The regulator-agnostic co-membership benchmark as the headline.
- Evaluation-only variance stratification in the perturbation test.  Its planned refinement is to compute each gene's variance without the perturbation's own experiments.
