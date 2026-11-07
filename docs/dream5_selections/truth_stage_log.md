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
