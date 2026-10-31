# Trace: self-loop for rows without positive CLR candidates

Date 2026-10-31, branch accuracy-review. Read-only trace; no code, paper text or results changed.
Question: Methods ("Attention operator", paper.md ~l.184) says rows without candidates keep a self-loop;
`R/packed.R:91` says rows without candidates "stay empty". Which is true for the numbers in the paper?

## Verdict

**No effect on any reported number.** The packed path is never used by a results-producing script; the dense
path that is used gives empty rows a self-loop (or the script adds one); and no dataset or split checked
contains a row without a positive candidate. The Methods sentence matches the code that produced the results.
The only discrepancy is a stale/misleading docstring in `R/packed.R` (and an inconsistency in five scripts that omit
the self-loop, which is moot because no row is empty). Nothing needs a fix for the paper; optional hygiene below.

## (1) Code path per script

No file under `analysis/`, `tools/`, `jobs_pc/` or the figure scripts calls `clr_packed()`, `softmax_operator()` or
any `src/*.cpp` operator code. `grep -rl "clr_packed\|softmax_operator"` outside R/, src/ and tests finds only
PORT_NOTES.md (a memory benchmark on K562, 8,563 genes, 16 Oct) and docs/plans/fig1_fresh_eyes_review.md.
The packed path is exercised only by `tests/testthat/test-packed.R`. The C++ core computes MI and CLR calibration
(`cpp_clr_calibrate`); the operator is built in R.

Every results-producing script builds the operator through one of two R routes, both in `R/clr_attention.R`:

| Route | Self-loop on empty row | Scripts |
|---|---|---|
| `ClrAttention$build_operator(softmax_keff=, softmax_cap=)` -> `.softmax_rows()` -> `.finish_operator()` (l.568-575 sets `A[i,i]=1` for rows with sum 0) | built in | baselines.R, tfnode_attention.R, clr_logpt_test.R (build_operator calls l.579, 789, 932), chips_bootstrap.R, nested_bootstrap.R, lab_splits.R (soft10 ops), perturbation.R, followup.R, run_m3d_regulondb.R, fig_method.R, beeline.R and beeline_boot.R (`soft_op`, the MI-based soft10 used for the reported BEELINE CLR attention) |
| `clr:::.softmax_rows()` called directly, then a script-local self-loop | added by script | keff_test.R (`selfloop`, l.313-314, softmax_custom), clr_logpt_test.R (l.238, 247, 316), shapley.R (l.71-74), granularity.R (l.31), perturbseq.R (l.102), perturbseq_bulk.R (l.97), multiscale_pilot.R (l.38), lab_splits.R ablation ops (l.148-149), synthetic_asymmetry/sim_asymmetry.R (l.95), synthetic_width/sim_width.R, sim_gate.R, sim_round2.R, sim_rclr.R (`selfloop`), rclr_k562/rclr_k562.R (l.91, 96) |
| `clr:::.softmax_rows()` called directly, **no** self-loop | no (row stays empty) | ecomac_kernels.R (l.48), pc_removal.R (l.78), beeline_blom_test.R (l.142), and the CLR-NS operator in beeline.R (l.325, 353) and beeline_boot.R (l.188) |
| Reads results only | n/a | macro_aupr.R, sota_table.R, rerun_gate.R, fig_graphical_abstract.R, fig_shapley.R |

`analysis/synthetic_width/realdata_gates.R` and `m3d_label_free.R` do not build a softmax operator.
So the dense route is used everywhere; the packed route nowhere.

## (2) Behaviour of a row with no positive candidate

- `.softmax_rows()` (`R/clr_attention.R` l.690): `if (!length(top)) next` leaves the row all zero.
- `build_operator()` -> `.finish_operator()` (l.568-577): rows with sum 0 get `A[i,i] = 1`. Row-stochastic.
- script-local `selfloop()`: same effect, `A[i,i] = 1`.
- `softmax_operator()` (packed, R/packed.R l.91-100): `if (!any(ok)) next`, row stays empty, no self-loop. This is the
  documented behaviour of the packed path, unused.
- No path produces NaN: the C++ calibration sets z = 0 when a row's sd is 0 (clr_core.cpp l.285), so a degenerate
  gene gives an all-zero CLR row, not NaN. `clr_packed()` refuses constant genes outright.

Walk consequences, P = (1-a)I + aA, S = (P^t + P^t')/2 with zero diagonal:
- With self-loop: row i of P is e_i, the gene is absorbing and keeps its mass.
- Without self-loop (empty row): row i of P is (1-a)e_i; its mass leaks (1-a)^t of it remains after t steps, the
  rest is lost (substochastic). The gene does not distribute mass to anyone (isolated source).
- CLR is symmetric (Stouffer combination), so C_ij <= 0 for all j implies C_ji <= 0 for all j: an empty row is a
  fully isolated gene (no gene lists it as a candidate either). Then off-diagonal entries of P^t in row and
  column i are 0 in both variants, so **S is identical off the diagonal** and the diagonal is zeroed. Only the
  stationary distribution pi / hub diagnostics (baselines.R l.196-200) and row-sum checks could differ.
- A row with all off-diagonal scores <= 0 requires every z_ij <= 0 (and z_ji <= 0); a row with non-constant MI
  always has an entry above its own row mean, hence a positive CLR. Empty rows therefore need a gene with identical
  dependence to all others, i.e. a degenerate (constant) gene. All analyses drop constant genes first
  (e.g. beeline.R l.~287, ecomac_kernels.R l.65).

## (3) Do real datasets have such rows? (count of rows with zero positive off-diagonal CLR)

Cached CLR matrices (`C`), counted directly from the .rds files (min = smallest number of positive candidates in any row):

| Dataset | Cached file | Genes | Empty rows | Min positives per row |
|---|---|---|---|---|
| M3D | results/clr_logpt_test/m3d/cache/clr_full.rds (also keff_test/m3d/cache/fulldata.rds, clr_logpt_test/truth/m3d/cache/fit.rds) | 4,297 | 0 | 1,276 |
| EcoMAC | results/clr_logpt_test/ecomac/cache/clr_full.rds (also keff_test/ecomac) | 4,189 | 0 | 1,097 |
| K562, 5,671 genes | results/clr_logpt_test/truth/k562ess/cache/fit.rds | 5,671 | 0 | 1,698 |
| K562, 8,559 genes | results/keff_test/k562ess/cache/fulldata.rds | 8,559 | 0 | 2,396 |
| RPE1 (extra) | results/clr_logpt_test/rpe1/cache/clr_full.rds | 8,747 | 0 | 2,454 |

No NaN in any of them. Not cached: the CLR matrices of the 14 BEELINE cases, and those of the bootstrap draws and
held-out halves of the bulk data. For BEELINE I recomputed (cheap: 514-1,528 genes), with the same gene
selection, constant-gene filter, MI (hg bins, Stouffer) and CLR-NS construction as `analysis/beeline.R`, on the full
data and on one random cell split (both halves, each half's constant genes dropped as in `heldout_depth`):

| Case | Genes | Cells | Empty rows (MI full / NS full / MI halves / NS halves) | Min positives (MI full) |
|---|---|---|---|---|
| hESC 500 / 1000 | 873 / 1,323 | 758 | 0 / 0 / 0 / 0 | 283 / 421 |
| hHep 500 / 1000 | 911 / 1,367 | 425 | 0 / 0 / 0 / 0 | 300 / 426 |
| mDC 500 / 1000 | 807 / 1,270 | 383 | 0 / 0 / 0 / 0 | 321 / 517 |
| mESC 500 / 1000 | 1,071 / 1,528 | 421 | 0 / 0 / 0 / 0 | 323 / 495 |
| mHSC-E 500 / 1000 | 666 / 1,135 | 1,071 | 0 / 0 / 0 / 0 | 215 / 347 |
| mHSC-GM 500 / 1000 | 593 / 1,055 | 889 | 0 / 0 / 0 / 0 | 180 / 255 |
| mHSC-L 500 / 1000 | 514 / 692 | 847 | 0 / 0 / 0 / 0 | 184 / 215 |

(Smallest min over any half-split matrix was 187 for MI and 193 for NS.) Not checked: the 5,000-gene and
all-expressed-gene BEELINE contexts (heavier, secondary, and they use build_operator for the MI soft10), the ~1,000
bootstrap-draw CLR matrices, and the platform subsets in ecomac_kernels.R. For those the structural argument in
section 2 applies (constant genes are removed in each); I did not count them.

## (4) Do the paths agree where both apply?

Smoke test (`systemd-run --user --scope -p MemoryMax=8G`, nice 19): M3D, random 800 of 4,297 genes, 466 experiments.
Scratch scripts in the session scratchpad, not committed.

- `build_operator(softmax_keff=10, softmax_cap=50)` vs `.softmax_rows()` + `selfloop()`: max abs difference 1.1e-16.
- Packed `softmax_operator()` (float32 CLR) vs dense: identical sparsity pattern, max abs difference 3.2e-7 in A;
  max abs difference 2.9e-8 in S at t = 9 (lazy walk, a = 0.5, symmetrized, zero diagonal).
- There were no empty rows in this subset, so the two paths could not differ on the point in question. A search
  for a naturally empty row on 2,000 tiny random datasets found none (consistent with section 2).
- Forced test: set 8 genes' rows and columns of the dense C to zero (isolated genes, as a symmetric CLR would give).
  `.softmax_rows()` leaves 8 zero rows (sum 0); with `selfloop` all rows sum to 1. Over the 800-gene walk at t = 9,
  S is identical (max abs difference 0, Spearman 1) for all pairs, including the 6,364 pairs involving the isolated genes
  (all 0 in both variants). P^9 row mass of an isolated gene: 1 with self-loop, 0.00195 = 0.5^9 without; total mass
  per row averaged over all genes 1 versus 0.990. So the two conventions differ only in the retained mass of isolated
  genes, which is not read out by S (diagonal zeroed) but would alter pi if pi were computed with such rows.

## (5) Does Methods match the results-producing path?

Methods: "Rows sum to one, rows without candidates keep a self-loop, and P = (1-a)I + aA with a = 0.5 is read out
as 1/2(P^t + P^t')." This matches `build_operator`/`.finish_operator` and the script-local `selfloop` wrappers, which
produced the headline results. Five scripts (ecomac_kernels.R, pc_removal.R, beeline_blom_test.R, the CLR-NS
operator in beeline.R and beeline_boot.R) skip the self-loop; that does not match the Methods wording literally,
but it changes nothing because those operators have no empty row (and for a symmetric CLR it would not change S
even if one existed). The packed docstring (R/packed.R l.91) contradicts Methods but the packed path produced no
result.

## Optional hygiene (not done; no change needed for any number)

1. Either give `softmax_operator()` the same self-loop as `.finish_operator()` or state in its docstring that the caller
   must add one (the K562 packed benchmark in PORT_NOTES.md predates any use).
2. Add the `selfloop` wrapper to the five scripts above if they are ever rerun with a different gene filter.
3. If the hub-stratification pi is ever computed on data with constant genes, fix the self-loop convention first.
