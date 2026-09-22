# Diffusion-depth design review: fixed-MI (A) vs re-estimated-MI (B)

Independent technical review of two experimental designs for the diffusion-depth
study in the CLR-as-attention project. Written 2026-09-22. Citation hold is ON:
citation needs are flagged as [CITATION NEEDED]; none are invented here.

## 1. Shared setup and notation

- Data: `X`, genes x samples (toy: 24 x 100 synthetic; real: M3D E. coli v4 build 6,
  4,297 genes x 907 arrays).
- `estimate_mi()`: B-spline MI matrix `M(X)`, per-gene adaptive bins ("fd").
- `calibrate()`: CLR z-scores `S` from `M` (normal/euclidean default).
- `select_threshold(B=100)`: independent per-gene shuffles, rebuild MI+CLR per
  bootstrap, pooled empirical null via streaming histogram, cutoff `tau` by
  Tukey's higher criticism (default) or BH-FDR. **Cost driver of the pipeline.**
- `build_operator()`: `A_hat` = row-stochastic sparsified `S` (tau or top-k).
- `diffuse(steps)`: `E^{(t+1)} = P E^{(t)}`, `P = (1-a)I + a*A_hat`, trajectory
  `E^{(0..T)}` cached; `E^{(t)} = P^t X`.
- `find_modules()`: connected components of `((A + A')/2 > 0)`; `min_size=1`
  keeps singletons.
- `ModuleTester`: Hungarian F1 headline per reference module vs planted/RegulonDB
  modules; soft Sinkhorn layer for uncertainty only. Headline aggregates: weighted
  mean F1, recovery rate (F1 >= 0.5), all per-module rows retained.

Standing decisions this review takes as fixed:

- (S1) The module finder must run **independently on every E^{(t)}**; components of
  one fixed operator cannot demonstrate improvement across depth.
- (S2) Diffusion is defined with a **fixed** `A_hat`: `E^{(t+1)} = ((1-a)I +
  a*A_hat) E^{(t)}`.
- (S3) Depth selection is **unresolved**; the plan is a convolution-based statistic
  with a test/validate log ratio (not train/validate), and validation must not
  choose the stopping depth.
- (S4) Toy results need not be exact; the toy is a machinery shakedown.

## 2. The two designs, precisely

**Design A ("fixed MI", user's model: representation = f(MI, t)).**
MI estimated once from `X`. `A_hat`, `P` fixed. Trajectory `E^{(t)} = P^t X`.
At each `t`, modules are found independently **from `E^{(t)}` itself**: build a fresh
gene-gene similarity graph over the diffused profiles (rows of `E^{(t)}`), e.g.
absolute correlation + a fresh HC threshold from a correlation null + connected
components. MI is never re-estimated. The only scientific knob varying with `t`
is the iteration count.

**Design B ("re-estimated MI").** At each `t`, treat `E^{(t)}` as fresh input and
rerun the full pipeline independently: B-spline MI on `E^{(t)}` -> CLR ->
permutation threshold (`B=100` MI bootstraps) -> new operator -> components ->
ModuleTester. The diffusion trajectory itself still uses the `t=0` operator
(per S2), so the smoother is fixed while the evaluated graph is rebuilt per `t`.

**Design C (iterative, discussed in section 8).** Rebuild the operator from
`E^{(t)}`, then diffuse one step with the *new* operator: `E^{(t+1)} = P_t E^{(t)}`.

## 3. Causal attribution: what does F1-vs-t actually measure?

This is the decisive dimension, and it favors A.

- **Design A.** The data-generating process has exactly one moving part: `t`.
  `M(X)`, `S`, `A_hat`, `P` are frozen. The per-t module extraction (similarity
  graph + fresh threshold + components) is a *fixed measurement procedure*
  applied identically at each `t`; the threshold re-selection adapts to the
  smoothed null, but that is part of the instrument, not a second experimental
  knob. An F1-vs-t curve under A therefore measures the effect of diffusion depth
  on module recoverability, full stop. If the curve rises then falls, the rise is
  smoothing helping and the fall is over-smoothing hurting.
- **Design B.** Three things move with `t`: (i) the diffused representation,
  (ii) the re-estimated MI/CLR graph, (iii) the re-selected sparsity threshold
  (each `t` gets its own permutation null, hence its own `tau_t`). An F1
  improvement from `t=0` to `t=3` cannot be attributed to depth: it may be
  smoothing, or better MI estimation on smoothed profiles, or a luckier
  threshold draw. The curve is a joint function F1(t, graph_t, tau_t), and the
  depth-selection question "how far to diffuse" is not identified from it.

A secondary confound specific to B, worth naming explicitly: **self-fulfilling
recovery**. Diffusion with the module structure bakes within-module similarity
into `E^{(t)}` by construction; re-estimating MI on `E^{(t)}` then "discovers"
inflated within-module MI and finds the modules more easily. Part of any F1(t)
rise under B is the pipeline recovering structure it injected. Under A this
cannot happen to the diffusion operator (it is fixed); the per-t similarity graph
does see smoothed profiles, but the smoothing kernel is not re-fit to flatter
itself.

## 4. Consistency with the standing decisions

- **S1 (independent finder per E^{(t)}).** Both satisfy it, differently. A runs a
  fresh similarity-graph + threshold + components procedure on each `E^{(t)}`.
  B reruns the entire pipeline per `t`. Note the S1 rationale ("components of one
  fixed operator cannot demonstrate improvement") is respected by A: A never
  takes components of the fixed `A_hat`; the fixed operator is only the diffusion
  kernel, while modules come from per-t graphs over `E^{(t)}`.
- **S2 (fixed-`A_hat` diffusion).** A is exactly this. B is consistent with the
  letter (its trajectory also uses the `t=0` operator) but conceptually muddy:
  the kernel that produced `E^{(t)}` is frozen at `t=0`'s graph while the
  evaluated graph is rebuilt per `t` — you end up judging graphs that played no
  role in producing the representations they are judged on.
- **S3 (depth selection unresolved).** Neither design settles it, but they demand
  different things from the eventual statistic (section 6).

Implementation note for A: the package's `select_threshold()` is MI-specific
(it rebuilds `bspline_mi` + `clr_calibrate` internally per bootstrap). A's per-t
correlation graph needs a parallel path: correlation null + HC cutoff. That is a
small, cheap-to-run code addition, not a research question. (A hybrid A':
fixed diffusion operator, per-t module graph built with MI on `E^{(t)}` — keeps
the diffusion clean per the user's model but reintroduces the weaker form of B's
self-fulfilling concern at the extraction step. Prefer pure A with correlation
unless there is a principled reason MI is required for extraction.)

## 5. The degenerate limit t -> infinity

`P^t -> 1 pi'` (rank-1; `pi` the stationary distribution), so every row of
`E^{(t)}` converges to the same `pi`-weighted mean profile. Modular signal is
provably gone in the limit under *both* designs — the fall of the F1(t) curve is
expected, which is good (it is what makes depth selection non-trivial). The
designs differ in *how* the fall happens and whether it is clean or erratic:

- **Design A.** Pairwise correlations of rows -> 1 for all pairs (up to floating
  point noise). Two sub-cases: (i) if the threshold selector still separates
  signal from null, the graph goes complete -> one hairball component ->
  precision collapses; (ii) more likely, the permutation null collapses too
  (shuffled near-identical profiles are still near-identical), observed scores
  stop being extreme vs the null, HC selects ~nothing -> all singletons ->
  recall collapses. Either way F1 falls, but the mechanism (hairball vs
  fragmentation) depends on fragile threshold behavior in the degenerate regime.
  **Recommendation:** put a collapse guard on the experiment — stop (or flag)
  when the median pairwise row-correlation of `E^{(t)}` exceeds ~0.99 or the
  effective rank of `E^{(t)}` drops below a floor. Record the guard value; do not
  let `t` run into pure numerical noise.
- **Design B.** Strictly worse. B-spline MI with per-gene FD-adaptive bins on
  near-constant vectors is numerically unstable (bin widths degenerate as
  marginal variance -> 0); CLR z-scores of garbage MI are garbage; the per-t
  permutation null is degenerate as in A. The F1(t) tail under B is not just
  falling but potentially **erratic and non-monotonic** — poison for any
  depth-selection statistic that assumes a smooth rise-then-fall shape (e.g. the
  planned convolution-based one). Under B the usable `t` range must be cut off
  well before collapse, and "well before" is itself a tuning choice that the
  experiment was supposed to inform — circular.

## 6. Statistical issues

**Per-t threshold re-selection.** Legitimate under both designs as finder
machinery (S1 requires the finder to be independent per `t`; a stale `tau_0`
applied to smoothed data would be the real error, since smoothing shifts the
null). It is not a confound under A (fixed procedure, one knob) but it *is* a
confound under B (each `t` gets its own MI-estimation noise *and* its own null
*and* its own threshold — three noise sources per point on the curve).

**Multiple comparisons across t.** Scanning `t = 0..T` and reporting
`max_t F1(t)` is selection-biased under either design; the headline "best depth"
F1 is optimistic. This is precisely what S3's machinery (convolution-based
statistic, test/validate log ratio, validation not choosing depth) is meant to
fix, and it remains open. What the designs demand of that future statistic:

- Under A, `F1(t)` is a noisy measurement of a univariate function of one knob;
  smoothing the curve (convolution) and doing test/validate comparisons on it
  is well-posed.
- Under B, `F1(t)` mixes depth effects with per-t estimation noise; any
  selection statistic must disentangle them, which requires modeling the
  estimator noise per `t` — a strictly harder problem, arguably intractable
  without strong assumptions.

**What "validation must not choose depth" means here.** The depth must be picked
by a statistic computed on test data (or a test/validate log ratio), not by
tuning on the validation modules. Under A this is a clean prescription: pick `t`
by the statistic, report F1(t) once. Under B, even the *statistic itself* would
be computed from per-t re-estimated pipelines, reintroducing the tuning the rule
was meant to prevent.

## 7. Computational cost

Let `C_MI` = one B-spline MI matrix, `C_boot` = B=100 MI rebuilds for the null
(`C_boot` >> `C_MI` dominates), `C_corr` = one correlation matrix (cheap).

- **Design A:** `1 x (C_MI + C_boot)` for the `t=0` operator, then per `t`:
  `C_corr` + a correlation-null HC selection (100 *cheap* correlation
  bootstraps, or a parametric/approximate null — open implementation choice).
  Total ~ `C_boot + T * small`.
- **Design B:** per `t`: `C_MI + C_boot` (full MI bootstraps each time). Total ~
  `T x C_boot`.

On the toy (24 genes) both are seconds; the toy cannot adjudicate cost. On M3D
(4,297 genes -> ~9.2M pairs; 907 samples), `C_boot` is 100 full MI matrices at
that scale — Design B multiplies it by `T` (e.g. 10), which is likely prohibitive
outside the user's 28-core desktop and infeasible in this workspace. Design A
pays `C_boot` once. This alone nearly decides the real-data run: **B at full
`T` is not runnable where the real experiment needs to run.**

## 8. Failure modes, one per design (steelman each side's worst case)

**A's worst case: garbage in, smoothed garbage out.** A inherits *all* errors of
the `t=0` graph. False edges in `A_hat` create spurious neighborhoods that
diffusion then reinforces; missing true edges can never be recovered by
smoothing, because `P` has no mass where `A_hat` has none (up to the
`(1-a)I` self-term). If the `t=0` threshold was too stringent, diffusion
propagates within fragments and the F1(t) curve is flat — correctly measuring
"depth doesn't help *this* graph," but saying nothing about whether a better
graph would have benefited. Mitigation: sensitivity runs over `tau_0`
(a small grid, still under A) and over `a`; cheap under A, prohibitive under B.

**B's worst case: the experiment measures itself.** Beyond the self-fulfilling
concern (section 3) and the erratic degenerate tail (section 5), B has a
quiet failure: because each `t` re-optimizes sparsity via its own null, a depth
can "win" purely by landing a favorable threshold draw. With per-t pipeline
noise, `argmax_t F1(t)` is as much a lottery over bootstrap seeds as a statement
about depth. Mitigation (fixed seeds, averaging over seeds per `t`) multiplies
the already-prohibitive cost.

**Design C (iterative rebuild + diffuse): verdict — more confounded, not more
principled *for a depth study*.** `E^{(t+1)} = P_t E^{(t)}` with `P_t` rebuilt
from `E^{(t)}` is a genuinely interesting *method* (adaptive diffusion; echoes
the "iterative CLR aggregates concepts" line), and it can repair early graph
errors, which is its real attraction. But as an experiment about *depth*, `t`
no longer indexes one thing — it indexes smoothing *and* graph evolution
jointly, with a feedback loop (smoothing -> denser re-estimated graph -> faster
collapse). It also breaks S2. Recommendation: keep C as a labeled future
direction (it may be the right *method*), but do not use it to answer "how deep."

## 9. Recommendation

**Toy shakedown: Design A.** It is the only design whose F1-vs-t curve answers
the question being asked, it is cheap, and it matches S2 exactly. Rewrite the
depth script as: MI once -> `A_hat` -> trajectory -> per `t`, independent
correlation-similarity graph + fresh HC threshold + components + ModuleTester
(Hungarian headline; soft layer recorded but not headline).

**Real M3D x RegulonDB run: Design A as primary.** Cost makes this nearly
forced (section 7). Optionally, run Design B at **2-3 depths only**
(e.g. `t in {0, t*_A, t_large}`) as a robustness check on whether the useful
depth range survives graph re-estimation — a diagnostic, not the curve.

**Record now to keep both options open** (cheap, one-time):

1. The full trajectory matrices `E^{(0..T)}` (or code + seeds to regenerate
   bit-identically).
2. `t=0` artifacts: MI matrix, CLR scores, `tau_0`, `A_hat`, `alpha`, bin rule,
   `B`, RNG seeds, combine method.
3. Per `t`: `tau_t`, mean degree, # modules, full ModuleTester outputs
   (per-module Hungarian rows + soft uncertainty), all as CSVs.
4. Collapse diagnostics per `t`: median pairwise row-correlation of `E^{(t)}`,
   effective rank, row-variance floor — whichever guard is adopted.
5. With (1)-(2) saved, Design B can be run later at any subset of depths without
   redoing diffusion; nothing is foreclosed.

## 10. Open questions (not settled by this review)

- The per-t similarity for A: correlation is recommended (cheap, keeps MI fixed
  as the user's model requires); needs the small correlation-null + HC code path.
- The collapse guard threshold (section 5): pick a number, record it, sensitivity-check.
- S3's depth-selection statistic itself: still open; this review only establishes
  that A gives it a well-posed input and B does not.
- Whether Design C (iterative rebuild) should become a separate method track:
  yes, later, as its own study — not as the depth experiment.

## 11. Citation needs (flagged, not filled; hold is ON)

- [CITATION NEEDED] Tukey's higher criticism (Donoho & Jin) — threshold rule.
- [CITATION NEEDED] Lazy random walks / diffusion maps (Coifman & Lafon;
  Zhou & Scholkopf 2004 already logged per project notes — confirm scope).
- [CITATION NEEDED] B-spline MI estimator used in `bspline_mi`.
- [CITATION NEEDED] Hungarian algorithm (Kuhn / Munkres) and Sinkhorn
  (Cuturi) for the tester — likely already covered by tester notes; confirm.
- [CITATION NEEDED] M3D database and RegulonDB for the real-data run.
- [CITATION NEEDED] Freedman-Diaconis bin rule (already logged per notes —
  confirm it covers the per-gene adaptive use).
