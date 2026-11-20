# tcaCLR2 against tcaCLR on the five used truths (exploratory)

Rule and reading as written in `docs/plans/tcaclr2_definition.md` before any score existed: depth t2 = round(ln 0.05 / ln edge), the edge being the mean over 10 per-gene permutations of the largest nontrivial eigenvalue modulus of the symmetrised lazy aCLR operator (keff 10, cap 50); "not worse" if pooled AUPR is at least 0.95 times the registered tcaCLR's on that truth. Outputs: `results/tcaclr2/` (summary.csv, edge_*.tsv, score_*/). Reference check: tcaCLR at t = 9 reproduces the stored Arm B scores exactly (maximum difference 0). Descriptive only: the depth-24 and depth-32 curves of the grid had been seen on these truths.

| Truth | t2 | AUPR tcaCLR (t = 9) | AUPR tcaCLR2 | Ratio (gene interval) | Reading |
|---|---|---|---|---|---|
| K562 CORUM | 29 | 0.1054 | 0.1038 | 0.985 (0.945 to 1.024) | not worse |
| Yeast net4 Complex Portal | 21 | 0.2278 | 0.2244 | 0.985 (0.965 to 1.005) | not worse |
| HepG2 Arm A | 32 | 0.0669 | 0.0615 | 0.919 (0.875 to 0.963) | **worse** (below 0.95; the interval excludes 1) |
| Kemmeren Arm B | 23 | 0.2076 | 0.2087 | 1.005 (0.977 to 1.033) | not worse |
| M3D TF regulons (SC) | 25 | 0.1622 | 0.1681 | 1.036 (1.018 to 1.055) | not worse (higher) |
| M3D TF regulons (all) | 25 | 0.1874 | 0.1937 | 1.034 (1.020 to 1.047) | not worse (higher) |

Edge and depth by dataset: M3D 0.886 (t2 25), K562 0.903 (29), net4 0.868 (21), Kemmeren 0.876 (23), HepG2 0.910 (32). The rule gives the deepest walk to the largest matrix (HepG2, 9,622 genes), and HepG2 is where the grid curves had deeper walks scoring lower (depth 4 best at 0.0714 against 0.0656 at depth 9).

**Reading.** tcaCLR2 meets the pre-written "not worse" reading on four of five truths (both M3D scopes are higher, with intervals above 1), and fails it on HepG2 by about 8% relative. The HepG2 interval is gene-resampling only. So the data-derived depth is not shown to be safer than the constant on the one truth that carries the registered H1; "more likely to work on a new dataset" is not supported, and the wording rule of the definition note stands (no claim that it will work). The edge rule's depth grows with matrix size and with how sparse the operator is; whether it should is the question a fresh truth would have to answer.
