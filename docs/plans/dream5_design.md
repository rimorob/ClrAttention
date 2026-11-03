# Preregistration (DRAFT): CLR, CLR-NS, CLR attention and baselines on the DREAM5 network-inference challenge data

Status: DRAFT rev 1, 2026-11-03. Not fixed. (The truth stage of `analysis/dream5_download.py` refuses to run until a line starting `Status: FIXED <date>` replaces this one, the file is committed and pushed, and the data log is committed and pushed; see "Protocol".)

*Drafted by Claude (Sonnet 5.5) for the PI; not yet reviewed.*

- Items marked **[V]** are recalled from Marbach et al. 2012 (*Nat. Methods* 9:796) or from file listings and **must be verified** against the downloaded training files and the challenge's `INSTRUCTIONS.txt` and evaluation scripts (download stage `scoring`, no truth) **before this file is fixed**. Any difference goes in the change log.
- No DREAM5 gold standard has been downloaded, opened or scored. The only DREAM5 content seen so far is public metadata (the Synapse file tree and names), the organisers' login-free expression archive (dimensions and replicate structure only), and the GitHub file lists of two third-party mirrors (names and sizes only).
- No score of any method has been computed on any DREAM5 data.
- After this file is fixed, any change is a dated change-log entry. A change after the first truth-phase run is a deviation and is reported as one.

## Question

The headline result (attention beats CLR on *E. coli* RegulonDB co-membership, t = 9, keff 10) was found on M3D and EcoMAC, which are one compendium in effect (`docs/plans/ecomac_m3d_overlap.md`). The yeast held-out test ended at preregistered verdict 1 (not run; `docs/STATUS.md`). DREAM5 is a public, community-standard benchmark with an *S. cerevisiae* network whose data and truth are independent of everything in this project. Two questions, both on the same fixed methods:

1. **Replication.** Does CLR attention beat CLR on an organism and a truth the project has never seen, for the headline readout (regulon co-membership) and for the benchmark's own readout (directed TF to target edges)?
2. **Standing.** Where do CLR, CLR-NS and CLR attention sit next to GENIE3, GRNBoost2, PIDC and |r| on the three DREAM5 networks, in DREAM5's own scoring? (Descriptive.)

## Data: what DREAM5 provides, and how it is obtained

Challenge project: Synapse `syn2787209`. Training data are one expression matrix per network (rows = chips, columns = genes with anonymised IDs G1…), a chip-features table, a gene-ID table and a TF list. The challenge expression data are already normalised by the organisers; **no further normalisation, log transform or batch correction is applied here.**

| Network | Organism | Genes x chips (replicate-level) | Experiments (replicates averaged) | TFs **[V]** | Gold-standard edges **[V]** | Truth source **[V]** |
|---|---|---|---|---|---|---|
| net1 | in silico (GeneNetWeaver) | 1,643 x 805 | 487 | 195 | about 4,000 | the simulated network (known) |
| net3 | *E. coli* | 4,511 x 805 | 487 | 334 | about 2,100 | RegulonDB, experimentally supported interactions |
| net4 | *S. cerevisiae* | 5,950 x 536 | 321 | 333 | about 3,900 | ChIP-chip binding plus conserved motifs (MacIsaac et al. 2006 lineage) |

Gene and chip counts are read from the organisers' login-free expression archive (3 Nov 2026). Net2 (*S. aureus*) is not used: its challenge truth is not a usable gold standard, and it is outside the overall score. Chips of one experiment are replicates; the experiment index is in the chip labels and in `net*_chip_features`.

**Availability verdict.**

- **Login-free (done):** the organisers' alternative-formats archive on the GeneNetWeaver site (`gnw.sourceforge.net/resources/DREAM5_NetworkInferenceChallenge_AlternativeDataFormats.zip`, 59 MB, SHA-256 `217e8869b379bbeeb3b1ed2fc7951021083a2c0711047fe6ff13cf2f850fd634`, 9 Jun 2010) holds the expression matrices of nets 1–4 transposed and replicate-averaged, with chip features. It has no TF lists, no gene-ID maps and no truth.
- **Synapse (PI's own login):** everything else. The tree is public metadata (no login needed to list it; `analysis/dream5_download.py listing` checks all 27 manifest IDs against it). Gold standards are Synapse files `syn2787240`–`syn2787244`; scoring scripts `syn2787219` and `syn2787218`.
- **Third-party GitHub mirrors** (`JBris/dream5_grn_data`, `k3jph/coms4761`) hold the data **together with the gold standards**, in renamed formats and with no provenance guarantee. They are not used.
- Download is staged: `data` (expression, chip features, gene IDs, TF lists, real names for *E. coli* and yeast; no truth), `scoring` (scoring code only), `truth` (gold standards and the challenge's published results; gated). Every file's size, SHA-256 and Synapse MD5 are logged to `docs/dream5_selections/download_log.tsv`. No credential is handled, printed or written by the script.

## Independence of each network from existing evidence (disclosure)

| Network | Data | Truth | Status and label |
|---|---|---|---|
| net3 *E. coli* | **M3D lineage.** The compendium derives from M3D (Faith et al. 2007); our 466 replicate-averaged M3D experiments and DREAM5's 487 experiments are expected to overlap heavily. EcoMAC already overlaps M3D in 37.7% of arrays. The depth t = 9, keff 10 and alpha were chosen on M3D. | **RegulonDB**, the source of our headline truths. Our M3D edge-level (TF-node) results exist: TF-node AUPR over random M3D CLR 6.75, GENIE3 6.54, attention 5.17 (`docs/whitepapers/paper.md`, "Direct edges"). | **NOT INDEPENDENT.** Labelled "same lineage as the headline" wherever reported. Never pooled with net4, never cited as a replication. It is a re-scoring of known material under another scorer. |
| net4 *S. cerevisiae* | GEO-compiled microarray chips (challenge data released June 2010). Predates the Kemmeren deleteome (2014) in `data/yeast`, so no overlap with it. Never used in this project. | Binding-based (ChIP-chip plus conservation), not expression-derived. A different file from the YEASTRACT truth of the failed yeast preregistration (that truth is unconsumed). | **The independent test.** One caveat: binding truths are weak proxies for regulation (known low scores for all methods in the challenge **[V]**). A later yeast test restricted to DNA-binding evidence would share provenance with this truth; this must be disclosed there. |
| net1 in silico | Simulated by GeneNetWeaver, with network structure from biological networks **[V]**; replicate and perturbation design mirrors the *E. coli* compendium (487 experiments, 805 chips). | Known exactly. | **Synthetic benchmark.** Descriptive mechanism check only (not organism-independent evidence); it can show whether a pattern is a property of the method or of real-data noise. |

**Method leakage and priors (disclosure).**

- GENIE3 and CLR (challenge-era implementations) took part in the original challenge; GENIE3 was the best single method, and the challenge's full results and predictions are public. This run is not blind, and published numbers are known to the team in outline. We reimplement all methods and make no claim to compare with the challenge ranking beyond a descriptive overall score.
- Our implementations differ from the challenge-era ones: the R `GENIE3` package, `tools/grnboost2.py` (arboreto routines, legacy early stopping), NetworkInference.jl PIDC, and the project's CLR (`ClrAttention`, B-spline MI, `bins = "hg"`, calibrated with `method = "normal"`, `combine = "stouffer"`). No statement is made about the origin of the combination step.
- Methods use expression only. Some challenge participants used the perturbation annotations in the chip features (e.g. which gene was deleted); we do not.
- **Prior from M3D:** attention was below CLR on direct TF to gene edges on M3D and EcoMAC (ratio about 0.77 vs CLR on M3D) and above on regulon co-membership. The design therefore states, before any DREAM5 truth, the expected pattern: module readout positive, edge readout not positive. It is a prediction to be tested, not an assumption.
- Nothing was tuned on DREAM5. No DREAM5 truth, score or chip-level quantity has entered any choice below.

## Hypotheses

All on net4 (S. cerevisiae), CLR attention vs CLR, fixed settings. Effect size is the log2 ratio of AUPRs, d = log2(AUPR_attention / AUPR_CLR). Smallest effect of interest: 10% relative AUPR, delta = log2(1.10) = 0.1375.

- **Y-M (module readout; the headline construct).** Regulon co-membership AUPR on regulons built from the yeast gold standard: the target set of each TF with 5 to 500 targets in the expression gene set (the project's `make_regulons` limits; targets only; no operon exclusion in yeast); pairs among genes in at least one regulon, positive if they share a regulon. Prediction: d > 0.
- **Y-E (edge readout; DREAM5's own).** Directed TF to target AUPR under the official scorer. Prediction: d <= 0 (the M3D pattern).
- **Context, not tested:** net3 and net1 with the same two readouts; all other methods.

## Methods run (fixed settings; the same expression matrix for every method)

Input per network: the official replicate-level `net*_expression_data.tsv` (all chips), genes with zero variance dropped (count reported). Candidate edges: TF list x all genes, no self-edges, TF to TF pairs scored in both directions. For gene-by-gene methods the scores of the TF rows are the predictions; undirected methods give equal scores to (a, b) and (b, a), as in the challenge.

| Method | Setting | Directed? | Gene-gene scores for Y-M? |
|---|---|---|---|
| \|r\| | absolute Pearson | no | yes |
| CLR | `estimate_mi(bins = "hg", transform = "none")`, `calibrate("normal", "stouffer")` as `analysis/baselines.R` | no | yes |
| CLR-NS | Blom normal scores with average ranks per gene (`blom_avg` in `analysis/beeline.R`), \|r\| of the scores, CLR (Stouffer) | no | yes |
| **CLR attention** | `build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50)` on CLR, lazy walk, **t = 9**, score (P^t + t(P^t))/2 on all genes, as `analysis/baselines.R` | no | yes |
| GENIE3 | R `GENIE3`, RF, K = "sqrt", 1,000 trees, seed 20261016; **regulators = the TF list** (DREAM5-native, directed, no symmetrisation) | yes | no |
| GRNBoost2 | `tools/grnboost2.py`, legacy early stopping, seed 20261016; **regulators = the TF list** (needs a `--regulators` option, see "Implementation") | yes | no |
| PIDC | NetworkInference.jl, Bayesian blocks, all genes (as `tools/pidc.jl`) | no | yes |

- Tree methods also get one run with **all genes as regulators**, then symmetrised by max (exactly as `analysis/baselines.R`), to provide a module readout for them. It is descriptive.
- t = 9 is a transfer value from M3D, not re-selected. No label-free depth rule is applied here.
- Nothing is chosen on DREAM5 truth, at any step.

## Scoring

All scorers are written, tested and committed before the truth stage, and used once.

- **Official edge scorer (Y-E, E-arm, in silico) [V]:** the challenge's own definitions of AUPR and AUROC over the ranked TF to target list, including its pair universe, its treatment of unlabelled pairs and its cap on the number of predictions (recalled: 100,000 ranked edges). Network score = minus log10 of the geometric mean of the AUPR p-value and the AUROC p-value, each against the organisers' random-network null; overall score = mean of the network scores over nets 1, 3 and 4. **The primary effect uses raw AUPR, not the p-value score** (the p-value score saturates at the null's resolution and is not a difference measure). The overall score, AUROC and the p-values are co-reported.
- **Scorer gate (before any score of ours):** the scorer reproduces the published AUPR and AUROC of at least three challenge methods, computed from the organisers' public predictions (`Network_predictions.zip`) against the published values (`MethodScoringResults.xls`), to 1e-3. If the official script can run, it is used and our reimplementation must match it. If neither matches, no scoring happens; the cause is logged.
- **Module scorer (Y-M, and the E-arm and in silico analogues):** pooled co-membership AUPR with the project's `comembership_pairs` and `ap_fast` (ties counted at the end of their block), regulons as defined above. Sensitivities (descriptive): expected-tie AUPR; regulons with at least 10 targets.
- **Replicates for intervals (shared by all methods):** see "Statistics".

## Endpoints and decision rule

Two co-primary endpoints, both on net4, each with its own fixed rule. **Any text that cites either must cite both** (the paper, abstract notes, any page for a reader). They answer different questions; there is no combined claim, and the familywise error across the two is not controlled (stated).

| Code | Readout | Contrast | Prediction |
|---|---|---|---|
| Y-M | co-membership AUPR, yeast regulons | attention vs CLR, d | d > 0 |
| Y-E | official directed edge AUPR, yeast | attention vs CLR, d | d <= 0 |

For each endpoint, with SE from "Statistics", one-sided tests at 0.025 per side (in effect a two-sided 0.05 test with directional conclusions, as in the yeast preregistration):

| Verdict | Condition |
|---|---|
| Superior | p_sup < 0.025 **and** d >= delta |
| Inferior (loss) | p_harm < 0.025 **and** d <= -delta |
| Equivalent within 10% | both non-inferiority tests (d > -delta, d < +delta) rejected at 0.025 (the 95% interval lies inside +-delta) |
| Inconclusive | none of the above; interval reported, no claim |

where p_sup = 1 - Phi(d / SE), p_harm = Phi(d / SE), and the two equivalence tests use (d + delta)/SE and (delta - d)/SE. Verdicts are written by the script to `verdict.csv` before anyone interprets them. Wording follows the earlier preregistrations: "not confirmed" with estimate and interval, never "no effect".

- **Evaluability gate for Y-M (structure only, before any score):** at least 20 regulons with >= 5 targets and at least 1,000 positive co-membership pairs in the expression gene set. Otherwise Y-M is "not evaluable" and Y-E is unaffected.
- **Interpretation lines fixed in advance.** Y-M superior and Y-E loss or inconclusive: "module gain with a direct-edge cost replicates in yeast". Y-M superior and Y-E superior: "attention improved both". Y-M not superior: "the co-membership headline did not replicate in yeast", with the interval. Y-M inferior: said plainly.
- **What does not count:** a win at another t, on another readout, on another network, or with a method variant; net3 or net1 results as support; pooling networks; results within +-10%.

## Statistics

- **Truth jackknife (SE_J).** Delete-half over TFs and, jointly, over target genes: each replicate keeps a random half of the TFs and a random half of the genes; the metric is recomputed on kept pairs (for Y-M, on the kept genes' pairs). 200 replicates, the same halves for every method; seeds fixed here: **net4 Y-E 20261203, Y-M 20261204; net3 20261213 and 20261214; net1 20261223 and 20261224.** SE_J = SD of the 200 paired differences d_b.
- **Experiment jackknife (SE_B).** Delete-half over **experiments** (chips of an experiment stay together; 321 experiments in net4). 100 replicates, seed 20263000 + b. CLR and attention (and |r|, CLR-NS) are refitted in every replicate; the truth and gene universe are fixed. Tree methods and PIDC are not refitted. Delete-half makes the delete-d factor 1, so SE_B = SD of the paired differences; it is conservative (half-size networks are noisier).
- **Deciding statistic:** SE = sqrt(SE_J^2 + SE_B^2) (a crossed genes x experiments structure; the sum counts the interaction twice and is conservative, Owen 2007). The interval for d is d +- 1.96 SE.
- **Other methods:** AUPR with the truth jackknife only, descriptive intervals, no p-values.
- The AUPR ratio and difference are both reported for every pair.
- Seeds are never changed; a rerun is allowed only after a crash or for a bug shown to violate a named line of this file, and both runs are reported.

## What is exploratory (descriptive; intervals only, never p-values, outside the endpoints)

- net3 (*E. coli*) and net1 (in silico), both readouts, all methods, with intervals, labelled as above.
- All other contrasts on net4 (CLR vs |r|, CLR-NS vs CLR, GENIE3, GRNBoost2 and PIDC against CLR and attention), AUROC, the p-value network score and the overall score over nets 1, 3 and 4.
- CLR-NS attention (soft10 on CLR-NS); directed attention (unsymmetrised P^t[TF, target]); attention at t = 3 and 27; hub-stratified scores (the project's `hub_stratify`).
- Tree methods' module readout (all-genes-regulator runs); precision at the top 100 and 1,000 edges; per-TF AUROC.
- Input sensitivity: replicate-averaged chips (the organisers' `_avg` file; our M3D analysis used replicate-averaged data).
- Label-free: overlap of net3 chips with our M3D chips (profile correlation after gene matching; descriptive, no gate); chip counts per experiment.
- Comparison with the published challenge table, as context only (not blind, not like-for-like).

## Protocol and gates (in order)

1. **Download, data stage (PI runs the command; no truth).** Commit the log. Verify dimensions against the table above, TF counts **[V]**, no duplicate gene IDs, matrices finite, and that the login-free archive's replicate-level matrices equal the Synapse ones after transposition (integrity cross-check).
2. **Download, scoring stage (scoring code only).** Read `INSTRUCTIONS.txt` and the scripts; record the pair universe, the prediction cap, the null construction and the p-value formula. Resolve every **[V]** in this file by a change-log entry.
3. **Fix this file** (replace the status line by `Status: FIXED <date>`), commit and push. Truth stage checks (in `dream5_download.py`): the design is tracked, identical to HEAD, carries the FIXED line, and its last commit and the data-log's last commit are on the upstream branch.
4. **Implementation, label-free.** `analysis/dream5_run.R` (one pipeline for all three networks; caches each method's TF-row and gene-gene scores with SHA-256), `tools/grnboost2.py --regulators`, the scorer and the replicate drivers. Smoke tests on **decoy truths only** (random TF to target edges with the real degree structure of a made-up network, never the real truth), and an equivalence test that the new GRNBoost2 option equals the old output when regulators = all genes. Commit.
5. **Compute all method scores** (nets 1, 3, 4; nets 1 and 4 first) and the refit replicates. Commit the score-file checksums. **No truth is opened during this step.** Truth-phase code refuses to run unless these are committed and unchanged.
6. **Code review** of the committed scoring and test scripts against this file, by a reviewer who has seen no DREAM5 score. Required changes are logged deviations.
7. **Truth stage** (`truth`, gated). Structure-only count (regulon counts, positive pairs, edges per TF) for the evaluability gate; then the **scorer gate** on the published predictions; then **score ours once**. Commit `verdict.csv` and all tables.
8. Report: both endpoints together, the labelling of net3 and net1, and the disclosure table.

## Implementation, compute and storage

- **New code:** `analysis/dream5_run.R` and scorer; `--regulators <file>` for `tools/grnboost2.py` (candidate regulators from the TF list; every target's own removal as before). Existing pieces reused unchanged: `ClrAttention`, `blom_avg`, `hub_stratify`, `make_regulons`, `ap_fast`.
- **Compute (estimates, from M3D scaling; to be replaced by measured values in the first smoke run).** CLR and attention per network: minutes (net4 about 15–20 min at 16 threads); refit replicates for net4: 100 x that, about 25–35 h; GENIE3 TF-regulators: hours; GENIE3 all regulators on net4: about a day; GRNBoost2: similar; PIDC (cubic in genes; the K562 subset run took about 40 h for 5,672 genes): about 2 days single-threaded on net4 and about 1 day on net3. All through `jobs_pc/sched` by the lead, with `# mem` and `# cpus`, caches resuming; PIDC on its own lane.
- **Storage:** raw data about 0.5 GB; PIDC edge lists about 1.5 GB; cached scores (TF-row matrices small, gene-gene for symmetric methods 0.3 GB each on net4); replicate tables small. Expected total well under 5 GB (`docs/disk_budget.md`).

## Limitations to state in any report

- net4's truth is binding evidence, incomplete and a weak proxy for regulation; every method scored low in the challenge **[V]**, so differences may be undetectable. Inconclusive is a likely outcome and a legitimate one.
- Power was not computed (no yeast truth structure is available before the truth stage); the realised SE is reported, and nothing is claimed from intervals that include both +-delta.
- net3 shares data and truth lineage with the headline; net1 is simulated. Neither supports a replication claim.
- Replicate-level chips are used (as in the challenge), so replicate chips are not independent samples; the experiment jackknife respects the grouping.
- Not blind: published results and the challenge ranking are known in outline.

## Change log

- 2026-11-03: DRAFT rev 1. Not fixed. Decisions needed from the PI are listed in the hand-off report of this draft (readout choice, margin, replicate-level versus averaged input, tree module readout, scope of comparators, PIDC compute).
- (skeleton) YYYY-MM-DD: [V] items resolved after the `data` and `scoring` stages: network sizes, TF counts, gold-standard counts, pair universe, prediction cap, null construction, p-value formula. Which items changed: ...
- 2026-11-03: PI decisions received (all recommendations accepted, still DRAFT): (1) two co-primary yeast endpoints Y-M and Y-E, always reported together; (2) margin 10% relative AUPR, with the stated prediction that Y-M is positive and Y-E is not; (3) replicate-level chips as the primary input, replicate-averaged file as a sensitivity run; (4) extra all-genes-regulator GENIE3 and GRNBoost2 module-readout runs, descriptive; (5) no C3 latent-column edge readout in this file; (6) PIDC and the net4 experiment-jackknife refits approved to run via `jobs_pc/sched`; (7) proceed with DREAM5. Next: the PI runs the two download commands; the [V] items are then verified against INSTRUCTIONS.txt and the scoring scripts; then FIXED.
- 2026-11-03: Synapse download stages `data` and `scoring` run by the main session with the PI's own temporary token (supplied by the PI for this purpose; never written to a file; the PI revokes it afterwards). Logged in `download_log.tsv`. **Disclosure:** the Synapse file `Evaluation scripts.zip` (syn2787219, SHA-256 937b039c...) contains the three gold-standard files (`INPUT/gold_standard_complete/` and `gold_standard_edges_only/`, networks 1, 3, 4) alongside the scoring code. The script refused to extract it (gold-like member names), so the gold standards are on disk inside the zip but **unextracted and unopened**; the zip was then sealed with `chmod 000`. Nothing was read from it beyond the member list printed by the script. The scoring code will be extracted only at the truth stage (after FIXED), together with the gold standards, so the `scoring` stage cannot be used to read scoring scripts earlier; the [V] items that depend on `INSTRUCTIONS.txt` (377 bytes; read, it only describes running the MATLAB scripts) are checked against the published paper instead, and against the scripts at the truth stage. Not a change to any endpoint.
- (skeleton) YYYY-MM-DD: FIXED. PI decisions recorded: ... No other change from DRAFT rev N.
- (skeleton) YYYY-MM-DD: DEVIATION n (before / after the first truth-phase run): where, what was found, what changed, who decided.
- (skeleton) YYYY-MM-DD: scorer gate result (methods, published vs reproduced values); truth stage run: commit, design hash, data-log hash, score checksums.
