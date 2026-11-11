# DREAM5 truth stage, part 1: extraction of the scoring archive and facts from the official scorer

Written by Claude (Sonnet 5.5), 2026-11-06, for the lead. Facts only: no score of any of our methods was computed, no
score file, verdict or refit file was opened, nothing in the design, the methods, the seeds or the endpoints was changed.
Lead decides the switches and writes the dated change-log entry.

## 1. Extraction (PI authorised opening the scoring archive)

- Archive: `data/dream5/synapse/scoring/Evaluation scripts.zip` (Synapse syn2787219), mode 000. Temporarily `chmod 600`,
  SHA-256 recomputed = `937b039cb50f42c086464605cb317a82a70ad15e02bfd4a4b1b40533ed78d60a` = the value in
  `download_log.tsv` (row of 2026-11-03T19:19:26Z); member list read with `unzip -Z1`; extracted with
  `unzip -d data/dream5/truth_extracted -x '__MACOSX/*' '*.DS_Store'`; then `chmod 000` again (verified).
- Output directory: `data/dream5/truth_extracted/` (new, empty before; `git check-ignore` confirms it is ignored via `data`).
  Nothing from it was executed (no MATLAB or Octave exists on this machine; see section 4).
- 90 archive members; 30 extracted files (the rest are `__MACOSX` resource forks and `.DS_Store`, skipped). SHA-256 of the
  extracted files:

```
2b5b808846e553e229723cb50f0be324f0250e722965665978aeda6c7fc12703  INPUT/gold_standard_complete/DREAM5_NetworkInference_GoldStandard_Network1.tsv
657c9d413e24279466df6ed3219c12c42b59068b0d4d92b28ed4f1e078ba78d3  INPUT/gold_standard_complete/DREAM5_NetworkInference_GoldStandard_Network3.tsv
8a2dfa1278c062e73b6b4511d0d84b055f3982e467e92a48467b05574ea9a0a0  INPUT/gold_standard_complete/DREAM5_NetworkInference_GoldStandard_Network4.tsv
4201d8e7fa475f45516eb454b6e5b8b9cd15504c4ab58263850f9d44fbf9b574  INPUT/gold_standard_complete/README.txt
5c6867c91e3b2584e0f9ae82aa9b1f67ca5cc52b122e0707917e52027e09ad82  INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network1.tsv
f6947dab958acb40858479bef3b0535fbeebdd5bca26edc1433cdcf4b72a2f59  INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network3.tsv
a8e1eec57707d129cf5354f732d92ce90540b954b0fb82c3d809a9151f786295  INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network4.tsv
d2cb89f304355733a17b6818fdf8dd65c14614fd6b618075d088ab171fca46f8  INPUT/gold_standard_edges_only/README.txt
04e7bcb4e1712c31aff225aeb114af99de6b607c78d2ddec8d6d6331ee1109fc  INPUT/predictions/myteam/DREAM5_NetworkInference_myteam_Network1.txt
3e0affca2e3b204f0c77d16d2e58028c0a00b6d1e1a1c34c90ddd7b2db6d0b77  INPUT/predictions/myteam/DREAM5_NetworkInference_myteam_Network3.txt
e536ca28a62a9c576f0e7563854c1ba0af83c4d1c32ff34a5aa8634fbc71112d  INPUT/predictions/myteam/DREAM5_NetworkInference_myteam_Network4.txt
a20ae86a92dcbd34f1a99e37793da9795c7e5d6321303e876245da40470fc917  INPUT/probability_densities/AUPR.pdf
53f0f7f35b3cceb20f1d3e090eae4548f38f58b2ff9d81721f6cc9230d932f1b  INPUT/probability_densities/AUROC.pdf
0d3c12b7542508b8c7ea51bac1c45fc923e6527fcc0b1bf5a235b33d0a13d64a  INPUT/probability_densities/Network1_AUPR.mat
9d2108dfe4e452d271673c59338aaa38c5e44bdb9a30c76b1454ca3ef5d96df4  INPUT/probability_densities/Network1_AUROC.mat
2e930a9d0e0fe531149bb6761b85d47332a40691fd91cfdfad48ea50efae9e6e  INPUT/probability_densities/Network3_AUPR.mat
9b7119e0c6341128a7e3ab3ae07e1d8ef6fb9d140f750c27d461c312e0d405a7  INPUT/probability_densities/Network3_AUROC.mat
8494de776737d42ba3bf2e4c2668ae255155c1086ef992a903d9c8a52b8e7732  INPUT/probability_densities/Network4_AUPR.mat
462f37fde5ed54ae0b118e833adb07d7171a31efff041524b15010bb520c7cd9  INPUT/probability_densities/Network4_AUROC.mat
1a9e0ccf6a465061af966a3c7a127d6a5ea8b55dba20288908429ccd245aa916  INSTRUCTIONS.txt   (identical to the separately downloaded INSTRUCTIONS.txt)
995d9f9cf497b313011bebadb6c009c6e79207d471535834132a9f360d987b45  matlab/directory_list.m
3122783accce0e2c05a52717a83eb94e8d1842e7d6acef4885497a1167c92dd4  matlab/DREAM5_Challenge4_Evaluation.m
57a04469e1772746ff43e0b4c25f2266e6fab67c2379fe323f5613f1ee4e007d  matlab/edgelist_2_D5C4_gs.m
ec6f7945cdbbd9bad7c6435a6818047b07d3979213e1014bdd6a7e410890bb7a  matlab/edgelist2sparse.m
7cb3925d0916b47495bee8cc7fd0a7cc35cc68c3c98fe5e81ed55b75d90b84d5  matlab/go_all.m
1889b5aaaf51fba420c81e75e50f6b65091974da8922aa2b73f14745ccf88558  matlab/go.m
6c1e8673129e7b20713779fe0b47e02ba7e5e30826ff8e58dc9674e7141ac8d6  matlab/load_dream_network.m
1c45c8e1acee8724f5b511962beac350c5e91f91946f7f2250939e29b343fe70  matlab/probability.m
829a596e25039c70c7c628f4227e9a0b2b24ddf91df117102a2d3c3f298934cb  matlab/remove_edges_not_in_gs.m
27920db03f1668921a1b002156fedab37473aab92c3cd4e5b84b9f6df9245de6  OUTPUT/myteam.txt
```

(`INPUT/predictions/myteam/*` and `OUTPUT/myteam.txt` are the organisers' SAMPLE files shipped for demonstration; they are
not the public challenge predictions, which are Synapse syn2787251 and are NOT in the archive.)

## 2. Facts from the official scripts (file:line in `Evaluation scripts/matlab/`)

- **Gold file used by the scorer:** `gold_standard_edges_only` (`go_all.m:24`, `go.m:25`; README of the `complete`
  directory: "the evaluation scripts actually use the files in the directory called gold_standard_edges_only").
- **Edge counts (data rows of the edges-only files; every row has label 1):** net1 4,012; net3 2,066; net4 3,940. (The
  design's recalled "3,940 for the in silico network" is the net4 count; net1 is 4,012.)
- **Cap:** no cap, no sort, no tie rule anywhere in the scripts (grep for 100000, sort, rand: no hit). The scorer uses the
  prediction file in FILE ORDER as the ranking; column 3 (confidence) is loaded (`load_dream_network.m:9,23`) and never
  used (`DREAM5_Challenge4_Evaluation.m:15-26` read only columns 1 and 2). The cap of 100,000 is therefore a property of
  the submitted file (each sample file has exactly 100,000 lines, sorted by confidence descending), applied BEFORE the
  scorer drops pairs outside the gold universe (sample file lengths after dropping, L: 86,225 / 12,032 / 12,305).
  Ties at rank 100,000 do not exist for the scorer (positions are ranks); a tie is broken by whatever order the submitter
  wrote. Duplicated pairs are not removed (`remove_edges_not_in_gs.m:10-20` copies every row whose cell is non-zero).
- **TF-TF pairs:** scored in both directions, neither collapsed nor dropped: the universe is regulators x targets with
  regulators included among the targets (`edgelist_2_D5C4_gs.m:3-4,13-19`); a cell (i, j) is positive iff the directed
  edge i to j is in the gold list, else negative unless i == j (`:15-19`).
- **Pair universe:** `regulators = unique(gold(:,1))` (`edgelist_2_D5C4_gs.m:3`) = TFs with at least one gold edge;
  `targets = unique(gold(:,1:2))` (`:4`) = every gene appearing in a gold edge in either column (so regulators are targets
  too). Universe = regulators x targets minus self (`:11-19`); positives P = `sum(sum(G>0))`, negatives N = `sum(sum(G<0))`,
  T = P + N (`DREAM5_Challenge4_Evaluation.m:10-12`). Predictions outside that box are discarded
  (`remove_edges_not_in_gs.m:13-19`, G(i,j) must be non-zero). Decoy genes and TFs without gold edges never enter.
  Sizes: net1 T = 278,392 (178 regulators, 1,565 genes); net3 T = 152,280 (141, 1,081); net4 T = 227,202 (114, 1,994);
  these T equal the row counts of the `complete` files exactly (278,392 / 152,280 / 227,202).
- **Unlisted universe pairs:** completed with the expected fraction of positives, p = (P - TPL)/(T - L), one fractional
  increment per remaining position (`DREAM5_Challenge4_Evaluation.m:33-45`): the tail is a uniform random order in
  expectation (the paper's "randomly ordered at the end").
- **AUROC:** `trapz(FPR, TPR)` over the completed list of T points (`:76`; the point (0,0) is not added; a loop version
  is commented out `:69-73`).
- **AUPR:** `trapz(REC, PREC) / (1 - 1/P)` over the same T points (`:104`), normalised "by max possible value". This is
  the trapezoid rule on the (recall, precision) points, NOT step-wise average precision and not interpolated precision. An
  alternative integration after Stolovitzky et al. 2009 is commented out (`:79-101`).
- **P-values:** `probability.m:1-5`: P(X >= x) = `sum((X>=x) .* Y * dx)` with dx = X(2) - X(1), from the precomputed null
  density `INPUT/probability_densities/Network{n}_{AUPR,AUROC}.mat` (variables X, Y, C; 100,000 grid points; dx = 1e-5;
  loaded in `go_all.m:49-52`; applied `DREAM5_Challenge4_Evaluation.m:110-116`). The code contains no random network
  generation or stretched-exponential fit: the null is a lookup table. Y is exactly 0 above a support top, so p = 0 exactly
  for a score above it: AUPR support top net1 0.4134, net3 0.3033, net4 0.1217; AUROC net1 0.9082, net3 0.8199, net4 0.6174.
  (Total mass of the density is 0.9994 to 0.9999, not 1, so p at the bottom of the support is 0.9994 to 0.9999.)
- **Overall score** (`go_all.m:23,66-68`): `aupr_score = -mean(log10(P_AUPR))`, `auroc_score = -mean(log10(P_AUROC))` over
  networks [1,3,4], `overall_score = mean([aupr_score auroc_score])`. Equal to minus the mean over the 3 networks of
  log10 of the geometric mean of the two p-values (the design's working definition). An exact p of 0 gives Inf.
- **Inputs expected** (scripts run from `matlab/`, relative paths): `../INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network{1,3,4}.tsv`;
  `../INPUT/predictions/<team>/DREAM5_NetworkInference_<team>_Network{1,3,4}.txt` (tab-separated `G<i> G<j> confidence`,
  already sorted best first); `../INPUT/probability_densities/Network{1,3,4}_{AUPR,AUROC}.mat`; `go_all.m` requires all
  three networks and takes the first directory of `INPUT/predictions` as the team.

## 3. Structure-only counts from the gold standards (no score; training-file ID tables used only for names and TF lists)

| | net1 | net3 | net4 |
|---|---|---|---|
| edges (distinct; self-edges 0) | 4,012 | 2,066 | 3,940 |
| distinct TFs with edges (all in the TF list) | 178 of 195 | 141 of 334 | 114 of 333 |
| targets per TF min / median / mean / max | 1 / 8 / 22.5 / 432 | 1 / 5 / 14.7 / 275 | 1 / 18.5 / 34.6 / 217 |
| genes in gold edges (either end) / distinct targets | 1,565 / 1,498 | 1,081 / 999 | 1,994 / 1,934 |
| TF-TF edges (target in TF list); reciprocal pairs | 246; 16 | 143; 11 | 324; 5 |
| edges with a decoy TF / decoy target | 0 / 0 | 0 / 0 | 0 / 0 |
| TFs with 5-500 targets (= with >= 5; none above 500) | 130 | 82 | 90 |
| TFs with 10-500 targets | 76 | 43 | 77 |
| genes in >= 1 regulon (5-500) | 1,475 | 948 | 1,918 |
| positive co-membership pairs, unordered distinct (5-500) | 208,688 | 67,929 | 141,943 |
| same, regulons with >= 10 targets | 207,968 | 67,544 | 141,774 |

Y-M evaluability gate (>= 20 regulons and >= 1,000 positive pairs on the 5-500 set), net4: 90 regulons, 141,943 pairs:
passes on both counts (computed directly from the gold standard, not with the project's gate function; the pair count here is
unordered distinct pairs among genes sharing a regulon, so the project's count may differ by definition).

## 4. Scorer gate inputs

- `Network_predictions.zip` (syn2787251, dataFileHandleId 759759, modified 2014-10-29) and
  `DREAM5_NetworkInference_MethodScoringResults.xls` (syn2787254, handle 759763, modified 2014-10-29) are NOT in the archive
  and NOT on disk under `data/dream5`. They exist only on Synapse (the `truth` stage of `analysis/dream5_download.py`); sizes
  are not visible without a login. They need the PI's Synapse login; nothing was attempted.
- No MATLAB or Octave on this machine (`which matlab octave octave-cli`: nothing; no octave package). The official scripts
  could not be run. A faithful Python port (`data/dream5/truth_extracted_tools/d5_official_port.py`, run with `python3 -I`)
  reproduces every number of the shipped `OUTPUT/myteam.txt` from the shipped sample predictions: AUPR 0.027128 / 0.027424 /
  0.019338 (file: 0.027 / 0.027 / 0.019), AUROC 0.681817 / 0.599660 / 0.515037 (0.682 / 0.600 / 0.515), P_AUPR 9.999e-01 /
  9.992e-01 / 9.129e-01, P_AUROC 6.268e-01 / 2.533e-04 / 1.197e-01, overall 0.794, AUPR score 0.013, AUROC score 1.574, all
  equal at the printed precision. This validates the port, not the published-method numbers (the gate still needs the
  public predictions and the xls).
- Convenience files (gitignored): `data/dream5/truth_extracted_universe/net{1,3,4}_universe_{tf,gene}.txt` (one ID per
  line: gold regulators / genes in gold edges), in the format `--universe_tf` / `--universe_gene` expect.

## 5. Where the official scorer and the implementation differ (see the hand-back report for the switch mapping)

Facts only: (a) AUPR integration: trapezoid on PR points / (1 - 1/P), implementation uses step-wise average precision;
(b) the official list has no tie rule (file order is the rank); the implementation offers end_of_block / expected only;
(c) the official universe is the gold box (universe = gs, cap before universe); (d) the p-values come from a lookup of the
precomputed density, with p = 0 above the support top.

## 6. Scorer gate (2026-11-10)

Written by Claude (Sonnet 5.5) for the lead. This changes the VALIDATION harness only, before any score of ours exists: no method,
endpoint, switch, seed or the official scorer itself was changed. No score file, verdict or refit file was opened and
`analysis/dream5_truth_run.R` was not run.

**Files received** (the PI downloaded them with their own login; rows of 2026-11-11T00:31Z UTC in `docs/dream5_selections/download_log.tsv`,
`truth` stage; stored under `data/dream5/synapse/published/`):

| file | Synapse id | bytes | SHA-256 |
|---|---|---|---|
| `DREAM5_NetworkInference_MethodScoringResults.xls` | syn2787254 | 58,368 | `244714913510e9bfc429aebe7efb2b3aa12a1c49afebb77eddcce21f16adc6ce` |
| `Network_predictions.zip` | syn2787251 | 79,484,816 | `1d6a2d7bba57a55f1dbff8fc2bf8c704466266c6f8766947af430c85d4c1ef2e` |

Both were extracted by the lead into the gitignored `data/dream5/published_extracted/` (`predictions/Network_predictions/`, `results/`;
the extracted `r.xls` has the logged SHA-256). The extracted prediction files are untrusted data: only read.

**Table structure.** The xls has one sheet, `survey4.xls` (a survey of the teams, exported to `r-survey4.xls.csv` with LibreOffice). Row 1-2 are
titles; then, for each category, a header row (`Regression`, `Mutual information`, `Correlation`, `Bayesian networks`, `Other`, `Meta`,
`Community`) followed by rows whose column 2 is the within-category number (1..8) or `RankSum` (Community). Columns (0-based): 2 overall score,
3 AUPR score, 4 AUROC score (the organisers' -log10 p scores averaged over networks 1, 3, 4; overall = mean of the two), 5-7 AUPR networks 1, 3, 4,
8-10 AUROC networks 1, 3, 4, 11-12 network 2 (S. aureus, not scored; ignored); columns 14 on are the survey answers. There are NO per-network
p-values in the table, so p-values enter the gate only through the three aggregate scores. 36 rows: Regression 8, MI 5, Correlation 3, Bayesian 6,
Other 8, Meta 5, Community 1. A re-export with "save cell contents as shown" off gave identical values: the sheet itself stores 3 decimals (a few
cells 4), so the printed rounding is the stored precision. The archive holds 29 team folders in `Challenge participants`, 6 in `Off-the-shelf
methods` (Correlation2-3, MI1-3, Regression8) and the four `Community` files: 36 names, the same 36 names as the table (matched by name only).
Each team has `DREAM5_NetworkInference_<team>_Network{1,2,3,4}.txt`; every file has at most 100,000 lines (net1 files 99,999).

**Converter and flat directory.** `analysis/dream5_published_table.py` (run with `python3 -I`): `tidy` writes
`data/dream5/published_extracted/published_table_tidy.csv` (team = category prefix + number, `Regression3`, `MI4`, ..., `Community`; score_overall,
score_aupr, score_auroc, aupr1/3/4, auroc1/3/4), asserting the row counts above, that every number parses, and that overall is within 2e-3 of the
mean of the two scores. `flat-pred` builds `data/dream5/published_extracted/flat_pred/<team>/DREAM5_NetworkInference_<team>_Network<n>.txt` as symlinks (36
teams x 4); the extracted files are not modified. Both are tested (synthetic sheet and tree, and the real sheet when present).

**What was adapted in the harness** (`analysis/dream5_official.R`, `analysis/dream5_scorer_gate.R`): the results reader is strict (the tidy
columns, duplicates and non-numbers are errors); the comparison covers per-network AUPR and AUROC (6 values) and the three aggregate scores, which
our side computes as go_all.m does (`d5_official_overall`: aupr_score = -mean(log10 P_AUPR) over networks 1, 3, 4, auroc_score likewise, overall =
their mean); R vs port vs matrix agreement (1e-9) now also covers -log10 of both p-values.
Confirmation of the aggregate definition on the shipped demo (`OUTPUT/myteam.txt`, scores 0.794 / 0.013 / 1.574): from the printed p-values
(P_AUPR 9.999e-01, 9.992e-01, 9.129e-01; P_AUROC 6.268e-01, 2.533e-04, 1.197e-01) the aggregates are AUPR score 0.013323, AUROC score 1.573714,
overall 0.793518 = their mean; the file shows 0.013 / 1.574 / 0.794 (go_all.m:66-68 and :79).

**Pass rule, fixed before the full run.** For a team, every one of the six per-network AUPR/AUROC values and the three aggregate scores must satisfy
|ours - published| <= 1e-3, for R, for the Python port and, where it ran, for the matrix path. The OR "within the printed rounding" is not used as a
second alternative: the table stores at most 3 decimals (rounding error <= 5e-4), and trailing zeros cannot be recovered from a stored number
(8.8 may be 8.800), so the printed precision is never coarser than 3 decimals and 1e-3 already covers it. Gate PASS = at least 3 teams match fully
(R and port), R vs port vs matrix agree to 1e-9 in every team. Teams are matched by name only. (A 3-team smoke run, Regression1, Community and Other7,
was done to time the pipeline before the full run; all three matched; nothing was changed afterwards except the explanatory comment of the rule
and an empty-result return value, neither of which affects any number.)

**p = 0 finding.** The official code has no floor and no cap anywhere on p or on -log10 p: `probability.m:2-5` is `dx = X(2)-X(1); P = sum(double(X>=x) .* Y * dx)`;
`DREAM5_Challenge4_Evaluation.m:112` and `:116` assign it to P_AUPR and P_AUROC unchanged; `go_all.m:66-67` take `-mean(log10(P))` directly, so an
exact p of 0 gives `Inf` (MATLAB `log10(0) = -Inf`) in the aggregate score. p is exactly 0 only above the support top of the null density
(AUPR 0.4134 / 0.3033 / 0.1217, AUROC 0.9082 / 0.8199 / 0.6174 for nets 1 / 3 / 4); below it the density is strictly positive (down to denormals,
about 1e-322), so p can be tiny (smallest positive tail 5e-324) but is not 0. In the published table no value is Inf, the largest score is 58.6
(Community AUPR score), and none of the 36 teams reached p = 0 in any network (all our aggregates are finite); so the gate never needed the p = 0 case. The gate mimics the official
code (no floor, no cap; Inf compared as Inf, a finite published value would then fail and be reported). The design pipeline does not cap either
(change-log decision 8).

**Gate result: PASS** (`results/dream5/scorer_gate/gate_report.md`, `gate_comparison.csv`; force-added to git as other result files are). 36 teams, all
fully matching for R, the Python port and the matrix path (required: 3). Command:
`nice -n 19 Rscript analysis/dream5_scorer_gate.R --pred_dir data/dream5/published_extracted/flat_pred --results_csv data/dream5/published_extracted/published_table_tidy.csv --out_dir results/dream5/scorer_gate`
(about 2 minutes wall, one process). Max |R - port| = 8.8e-12, max |R - matrix path| = 3.9e-15 over AUPR, AUROC, both p-values and their -log10.
No file repeats a pair (matrix path skipped 0 of 108 team-networks), so the repeated-pair handling (R and port do not de-duplicate, as the
official scorer; the matrix path would skip the file) was never exercised on real data. Largest differences to the published table: per-network
metric 5.8e-4 (Correlation3 net3 AUROC: ours 0.575579, published 0.575, i.e. the sheet appears truncated rather than rounded in this cell; MI1 net1 AUROC
0.773513 vs 0.773 likewise); 2 of 324 metric values exceed 5e-4, all within 1e-3; aggregate scores at most 4.9e-4. Differences are consistent with the sheet's 3-decimal
storage, not with a systematic offset. The table below shows the three aggregate scores (published / ours) and the maximum absolute differences;
port and matrix path equal R to 1e-11.

| team | published overall / AUPR score / AUROC score | ours (R) | max abs diff, 6 per-net metrics | max abs diff, 3 scores | all 3 implementations pass |
|---|---|---|---|---|---|
| Bayesian1 | 4.973 / 9.947 / 0 | 4.9734 / 9.9466 / 0.0002 | 2.95e-04 | 4.19e-04 | yes |
| Bayesian2 | 2.266 / 4.531 / 0 | 2.2656 / 4.5311 / 0.0002 | 3.31e-04 | 3.70e-04 | yes |
| Bayesian3 | 0.331 / 0.173 / 0.489 | 0.3307 / 0.1726 / 0.4888 | 4.61e-04 | 3.68e-04 | yes |
| Bayesian4 | 0.257 / 0.04 / 0.473 | 0.2566 / 0.0403 / 0.4729 | 3.54e-04 | 4.01e-04 | yes |
| Bayesian5 | 0.256 / 0.04 / 0.472 | 0.2563 / 0.0403 / 0.4724 | 3.38e-04 | 4.07e-04 | yes |
| Bayesian6 | 0.035 / 0.003 / 0.068 | 0.0352 / 0.0026 / 0.0679 | 4.32e-04 | 4.28e-04 | yes |
| Community | 52.4271 / 58.5993 / 46.2549 | 52.4271 / 58.5993 / 46.2549 | 4.57e-04 | 2.69e-05 | yes |
| Correlation1 | 8.126 / 5.632 / 10.619 | 8.1260 / 5.6324 / 10.6195 | 4.08e-04 | 4.87e-04 | yes |
| Correlation2 | 1.191 / 1.717 / 0.666 | 1.1914 / 1.7171 / 0.6656 | 3.20e-04 | 4.27e-04 | yes |
| Correlation3 | 0.834 / 1.422 / 0.246 | 0.8337 / 1.4217 / 0.2456 | 5.79e-04 | 3.63e-04 | yes |
| MI1 | 20.462 / 23.906 / 17.018 | 20.4623 / 23.9062 / 17.0184 | 5.13e-04 | 3.67e-04 | yes |
| MI2 | 12.013 / 8.477 / 15.549 | 12.0134 / 8.4775 / 15.5494 | 4.96e-04 | 4.93e-04 | yes |
| MI3 | 9.2351 / 6.0926 / 12.3776 | 9.2351 / 6.0926 / 12.3776 | 4.90e-04 | 3.89e-05 | yes |
| MI4 | 8.8 / 3.984 / 13.616 | 8.7998 / 3.9836 / 13.6160 | 3.39e-04 | 3.90e-04 | yes |
| MI5 | 1.475 / 2.372 / 0.579 | 1.4754 / 2.3722 / 0.5787 | 4.39e-04 | 4.04e-04 | yes |
| Meta1 | 22.711 / 36.795 / 8.627 | 22.7112 / 36.7953 / 8.6271 | 4.76e-04 | 3.01e-04 | yes |
| Meta2 | 16.686 / 14.538 / 18.833 | 16.6856 / 14.5379 / 18.8333 | 4.83e-04 | 3.92e-04 | yes |
| Meta3 | 10.586 / 10.491 / 10.681 | 10.5856 / 10.4905 / 10.6807 | 3.76e-04 | 4.87e-04 | yes |
| Meta4 | 6.42 / 8.183 / 4.658 | 6.4204 / 8.1830 / 4.6578 | 4.66e-04 | 4.15e-04 | yes |
| Meta5 | 6.332 / 8.337 / 4.327 | 6.3322 / 8.3375 / 4.3269 | 4.67e-04 | 4.79e-04 | yes |
| Other1 | 40.279 / 41.295 / 39.263 | 40.2790 / 41.2953 / 39.2626 | 3.59e-04 | 3.53e-04 | yes |
| Other2 | 34.023 / 30.908 / 37.137 | 34.0225 / 30.9082 / 37.1369 | 4.61e-04 | 4.65e-04 | yes |
| Other3 | 8.887 / 3.076 / 14.699 | 8.8874 / 3.0757 / 14.6992 | 4.69e-04 | 4.22e-04 | yes |
| Other4 | 2.484 / 2.133 / 2.835 | 2.4839 / 2.1327 / 2.8352 | 4.93e-04 | 2.88e-04 | yes |
| Other5 | 1.897 / 0.066 / 3.727 | 1.8966 / 0.0665 / 3.7266 | 4.31e-04 | 4.78e-04 | yes |
| Other6 | 0.144 / 0.25 / 0.038 | 0.1441 / 0.2502 / 0.0380 | 4.11e-04 | 2.02e-04 | yes |
| Other7 | 0 / 0 / 0 | 0.0002 / 0.0002 / 0.0002 | 3.90e-04 | 1.65e-04 | yes |
| Other8 | 0 / 0 / 0 | 0.0002 / 0.0002 / 0.0002 | 4.54e-04 | 1.65e-04 | yes |
| Regression1 | 31.099 / 41.368 / 20.831 | 31.0990 / 41.3676 / 20.8305 | 4.11e-04 | 4.67e-04 | yes |
| Regression2 | 28.747 / 44.689 / 12.806 | 28.7475 / 44.6890 / 12.8060 | 4.70e-04 | 4.87e-04 | yes |
| Regression3 | 21.398 / 32.838 / 9.957 | 21.3979 / 32.8383 / 9.9575 | 4.09e-04 | 4.80e-04 | yes |
| Regression4 | 20.694 / 31.965 / 9.422 | 20.6937 / 31.9652 / 9.4222 | 3.99e-04 | 3.26e-04 | yes |
| Regression5 | 13.397 / 13.781 / 13.013 | 13.3969 / 13.7809 / 13.0128 | 3.94e-04 | 1.97e-04 | yes |
| Regression6 | 6.026 / 8.924 / 3.128 | 6.0265 / 8.9244 / 3.1285 | 4.02e-04 | 4.91e-04 | yes |
| Regression7 | 0.997 / 0 / 1.995 | 0.9974 / 0.0002 / 1.9945 | 4.11e-04 | 4.74e-04 | yes |
| Regression8 | 0.045 / 0.058 / 0.032 | 0.0446 / 0.0577 / 0.0316 | 2.77e-04 | 4.37e-04 | yes |
