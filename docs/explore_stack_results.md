# Exploratory: stacked a3r layers (does a3r lend itself to multiple CLR passes?): results

**Status.** EXPLORATORY. K562 is spent and M3D is the tuning set. a3r itself was chosen from about 20 readouts on K562 truths, so K562 is optimistic, and the stacked variants add to that selection. No claim, no p-values, no winner selected, nothing for the paper or the pre-read. Intervals are descriptive paired delete-half gene jackknife intervals (d +/- 1.96 SD of the paired replicate differences, networks held fixed). Plan: `docs/plans/explore_stack_plan.md` (commit 98f3b0d, written before any compute; deviations in the last section). Code: `analysis/explore_stack/` (`stack_lib.R`, `stack_k562.R`, `stack_m3d.R`, `stack_summarize.R`). Outputs: `results/explore_stack/{k562,m3d}/` (full tables with every truth, scope and interval in `tables.md`; `pass_fail.csv`, `contrasts.csv`, `pooled_levels.csv`, `diagnostics_*.csv`, `equivalence_gate_stack.csv`; rds caches not committed).

## 1. Literal pass / fail (plan section 6)

**All four stacked readouts (a3r2, a3r3 on K562 and on M3D) FAIL the literal rule, with both hub stratifiers, in the primary and in the strict reading.** Layer 2 and layer 3 both fail.

### K562 (rule: CORUM-all gain over a3r >= +0.005 with interval excluding 0; DoRothEA A-B and CollecTRI each >= -0.005; gain over the matched-depth control >= +0.005 with interval excluding 0; same sign under the hub stratifier). AUPR differences, raw unless stated.

| readout | (1) CORUM_all minus a3r | (2) DoRothEA A-B / CollecTRI minus a3r | (3) CORUM_all minus matched c1 | (4) hub MI / hub pi, CORUM_all minus a3r | pass (MI) | pass (pi) | strict MI | strict pi |
|---|---|---|---|---|---|---|---|---|
| a3r2 (2 layers), c1 = t 18 | +0.0031 [-0.0027, +0.0089] no | -0.0023 ok / -0.0061 no | +0.0025 [-0.0036, +0.0087] no | -0.0009 no / -0.0005 no | **FAIL** | **FAIL** | FAIL | FAIL |
| a3r3 (3 layers), c1 = t 27 | +0.0027 [-0.0062, +0.0115] no | -0.0077 no / -0.0108 no | +0.0039 [-0.0063, +0.0140] no | +0.0020 ok / +0.0017 ok | **FAIL** | **FAIL** | FAIL | FAIL |

Reasons: neither CORUM gain reaches +0.005 or has an interval excluding 0; both stacked readouts lose on CollecTRI by more than 0.005 (a3r3 also on DoRothEA A-B); neither beats the matched-depth control by +0.005; for two layers the hub-stratified sign is negative. "Strict" = conditions 1-3 on the hub-stratified scores (X, a3r and c1 all stratified); every strict reading also fails (CORUM_all minus a3r hub MI: a3r2 -0.0009, a3r3 +0.0020; hub pi: -0.0005, +0.0017).

### M3D (rule: RegulonDB SC and all each non-inferior to a3r at margin -0.005, jackknife lower bound; complexes not decreased)

| readout | scope | SC d [lower bound] | all d [lower bound] | complexes d | pass |
|---|---|---|---|---|---|
| a3r2 | raw | +0.0006 [-0.0043] ok | -0.0055 [-0.0092] **no** | +0.0362 ok | **FAIL** |
| a3r2 | hub MI | +0.0003 [-0.0043] ok | -0.0062 [-0.0098] no | +0.0184 ok | FAIL |
| a3r2 | hub pi | +0.0015 [-0.0030] ok | -0.0034 [-0.0070] no | +0.0353 ok | FAIL |
| a3r3 | raw | -0.0056 [-0.0111] no | -0.0120 [-0.0163] no | +0.0468 ok | **FAIL** |
| a3r3 | hub MI | -0.0060 [-0.0111] no | -0.0131 [-0.0174] no | +0.0222 ok | FAIL |
| a3r3 | hub pi | -0.0045 [-0.0095] no | -0.0088 [-0.0129] no | +0.0518 ok | FAIL |

a3r2 is non-inferior on SC but not on "all"; a3r3 is non-inferior on neither. The complex truth rises for both (point estimates), which satisfies that clause.

## 2. Top-line numbers (pooled AUPR, enrichment in brackets)

K562 raw (base rates: CORUM_all 0.0133, DoRothEA A-B 0.1830, CollecTRI 0.1278, KD 0.0787):

| readout | CORUM_all | CORUM merged | CORUM noRP | DoRothEA A-B | CollecTRI | KD response |
|---|---|---|---|---|---|---|
| K10 (t 9) | 0.0465 (3.51x) | 0.0464 (3.39x) | 0.0243 (2.05x) | 0.2003 | 0.1402 | 0.1870 |
| CLR | 0.0847 (6.39x) | 0.0843 (6.16x) | 0.0259 (2.18x) | 0.1939 | 0.1315 | 0.1908 |
| **a3r (1 layer)** | **0.1054 (7.95x)** | 0.1050 (7.66x) | 0.0300 (2.53x) | 0.2036 | 0.1450 | 0.2460 |
| a3r2 | 0.1086 (8.19x) | 0.1078 (7.87x) | 0.0270 (2.28x) | 0.2013 | 0.1388 | 0.2025 |
| a3r3 | 0.1081 (8.15x) | 0.1074 (7.84x) | 0.0285 (2.40x) | 0.1959 | 0.1342 | 0.1675 |
| c1, t 18 (matched to a3r2) | 0.1060 (8.00x) | 0.1056 (7.71x) | 0.0283 (2.38x) | 0.2009 | 0.1486 | 0.2704 |
| c1, t 27 (matched to a3r3) | 0.1042 (7.86x) | 0.1039 (7.59x) | 0.0264 (2.23x) | 0.1986 | 0.1502 | 0.2783 |
| c2_z (repeat CLR on S1) | 0.0776 (5.85x) | 0.0773 (5.65x) | 0.0250 (2.11x) | 0.2020 | 0.1391 | 0.1867 |
| c2_log (a3r(S1)) | 0.0906 (6.83x) | 0.0903 (6.59x) | 0.0265 (2.24x) | 0.2031 | 0.1421 | 0.2248 |
| a3r2, t2 = 2 | 0.0216 (1.63x) | 0.0219 | 0.0167 | 0.1881 | 0.1300 | 0.1117 |
| a3r2, t2 = 4 | 0.0411 (3.10x) | 0.0410 | 0.0262 | 0.2014 | 0.1372 | 0.1696 |

K562 differences to a3r (raw, AUPR, interval): a3r2 CORUM_all +0.0031 [-0.0027, +0.0089], merged +0.0028 [-0.0022, +0.0078], noRP -0.0030 [-0.0047, -0.0013], DoRothEA A-B -0.0023 [-0.0053, +0.0006], CollecTRI -0.0061 [-0.0087, -0.0036], KD -0.0436 [-0.0481, -0.0390]. a3r3: CORUM_all +0.0027 [-0.0062, +0.0115], noRP -0.0015 [-0.0050, +0.0020], DoRothEA A-B -0.0077 [-0.0125, -0.0030], CollecTRI -0.0108 [-0.0148, -0.0068], KD -0.0785 [-0.0860, -0.0711]. Controls against a3r: c1 t18 CORUM_all +0.0006 [-0.0021, +0.0032], CollecTRI +0.0036, KD +0.0243; c1 t27 CORUM_all -0.0012 [-0.0052, +0.0028], CollecTRI +0.0053, KD +0.0323 (KD response is not in the rule; near-duplicate caveat). c2_z minus a3r on CORUM_all -0.0279 [-0.0369, -0.0188]; c2_log -0.0149 [-0.0236, -0.0061] (the pure repeat is not idempotent here: it loses on CORUM; on DoRothEA it is within about 0.002 of a3r).

K562 hub-stratified (MI row mean), CORUM_all AUPR: a3r 0.0854 (6.44x), a3r2 0.0844 (6.37x), a3r3 0.0874 (6.60x), c1 t18 0.0838, c1 t27 0.0810; pi-stratified values and all other truths in `results/explore_stack/k562/tables.md`.

M3D raw (base rates: complexes 0.00736, SC 0.0943, all 0.1205):

| readout | complexes | RegulonDB SC | RegulonDB all |
|---|---|---|---|
| K10 | 0.1664 (22.6x) | 0.1579 | 0.1733 |
| CLR | 0.1605 (21.8x) | 0.1359 | 0.1648 |
| **a3r** | **0.1921 (26.1x)** | **0.1622** | **0.1874** |
| a3r2 | 0.2283 (31.0x) | 0.1629 | 0.1819 |
| a3r3 | 0.2389 (32.5x) | 0.1566 | 0.1753 |
| c1, t 18 | 0.1931 | 0.1671 | 0.1920 |
| c1, t 27 | 0.1926 | 0.1682 | 0.1940 |
| c2_z | 0.1878 | 0.1554 | 0.1798 |
| c2_log | 0.1347 | 0.1576 | 0.1804 |
| a3r2, t2 = 2 / t2 = 4 | 0.0607 / 0.1386 | 0.1199 / 0.1625 | 0.1390 / 0.1737 |

M3D differences to a3r (raw): a3r2 complexes +0.0362 [-0.0028, +0.0752], SC +0.0006 [-0.0043, +0.0055], all -0.0055 [-0.0092, -0.0018]; a3r3 complexes +0.0468 [+0.0035, +0.0902], SC -0.0056 [-0.0111, -0.0001], all -0.0120 [-0.0163, -0.0077]. Matched-depth controls against a3r: c1 t18 complexes +0.0010 [-0.0062, +0.0081], SC +0.0049 [+0.0026, +0.0071], all +0.0046 [+0.0029, +0.0064]; c1 t27 complexes +0.0005, SC +0.0059, all +0.0066. Stacked minus matched-depth control: a3r2 - c1_t18 complexes +0.0352 [-0.0032, +0.0737], SC -0.0042, all -0.0101 [-0.0140, -0.0062]; a3r3 - c1_t27 complexes +0.0463 [+0.0033, +0.0893], SC -0.0115 [-0.0178, -0.0053], all -0.0186 [-0.0239, -0.0133]. Hub-stratified complexes gains are smaller (MI: a3r2 +0.0184 [-0.0123, +0.0490], a3r3 +0.0222 [-0.0098, +0.0542]; pi: +0.0353, +0.0518, intervals include 0).

## 3. Reading of the tables (limited to what they show)

- **K562: the stacked CORUM gain is a point estimate of +0.003 (8.19x against 7.95x) whose interval includes 0, and one layer at matched depth gives +0.0006 (t 18) and -0.0012 (t 27).** So two or three layers are not distinguishable from one layer on CORUM_all, and not distinguishable from the same depth in one layer (a3r2 minus c1_t18 +0.0025 [-0.0036, +0.0087]). On everything else stacking is at or below a3r: CollecTRI -0.006 and -0.011, KD response -0.044 and -0.079, noRP -0.003 and -0.0015 (a3r2 interval excludes 0), DoRothEA A-B -0.002 and -0.008.
- **K562 depth alone is roughly neutral on CORUM_all, slightly negative on noRP (-0.0017 at t 18, -0.0036 at t 27) and positive on CollecTRI and KD**; the stacked readouts lose on KD where depth alone gains (c1_t27 minus a3r +0.032, a3r3 minus a3r -0.079). Stacking is not the same operation as depth.
- **The pure repeat is not idempotent.** On K562 CORUM_all c2_z falls to 5.85x and c2_log to 6.83x (both intervals against a3r exclude 0); on M3D complexes c2_z is within noise (-0.004) and c2_log falls (-0.057).
- **M3D: complexes rise with stacking and not with depth** (a3r2 +0.036, a3r3 +0.047 against c1 +0.001), but the interval for two layers includes 0 and the complex truth has 1,150 positive pairs in 54 complexes (the 26 Oct verdict gave SE 0.018 for this truth; the interval here is from the same style of gene jackknife). At the same time RegulonDB "all" falls for stacking (-0.006, -0.012) and rises with depth (+0.005, +0.007); SC is flat for two layers (+0.0006) and falls for three (-0.0056). Hub stratification shrinks the M3D complex gain by about half.
- **Secondary t2 grid (descriptive):** a shorter second layer is worse everywhere: t2 = 2 gives CORUM 1.63x on K562 (the P^2 walk leaves 90% of pairs unreachable, so the log readout has few rows of support) and complexes 0.061 on M3D; t2 = 4 gives 3.10x and 0.139. The t = 9 transfer is the best of the three layer-2 depths on both datasets, but t2 > 9 was not run.

## 4. Where things degrade (label-free diagnostics; layers 1-5, layers 4 and 5 never scored against a truth)

Per layer (operator built from S_{k-1}; walk at t = 9). Full table: `diagnostics_layers.csv` in each results folder.

| quantity | K562 layers 1 / 2 / 3 / 4 / 5 | M3D layers 1 / 2 / 3 / 4 / 5 |
|---|---|---|
| rows with no positive candidate (self-loop rows) | 0 / 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 |
| rows with fewer than 10 candidates; median candidates | 0 ; 50 at every layer | 0 ; 50 at every layer |
| rows of all zero in S_k ("collapsed rows") | 0 at every layer | 0 at every layer |
| unreachable share of P^9 pairs | 0.28% / 5.7% / 12.3% / 17.2% / 18.1% | 0.07% / 3.0% / 2.9% / 5.5% / 11.1% |
| rows with < 2 reachable entries | 0 at every layer | 0 at every layer |
| Spearman of stationary pi_k with pi_1 | 1 / 0.74 / 0.52 / 0.40 / 0.29 | 1 / 0.84 / 0.64 / 0.53 / 0.40 |
| Spearman of pi_k with the MI row mean | 0.20 / -0.02 / -0.04 / -0.03 / -0.04 | 0.11 / 0.03 / 0.04 / 0.05 / 0.07 |
| pi_k max / median | 18,460 / 24,090 / 13,210 / 4,763 / 1,484 | 311 / 2,887 / 4,856 / 2,309 / 680 |
| share of in-weight held by the top 1% of genes | 0.136 / 0.189 / 0.179 / 0.170 / 0.156 | 0.045 / 0.076 / 0.095 / 0.102 / 0.102 |
| Spearman of the row mean of S_k with the MI row mean | 0.30 / 0.01 / 0.00 / -0.03 / -0.03 | 0.25 / 0.04 / 0.04 / 0.12 / 0.06 |
| top-1% pair overlap (Jaccard) of S_k with S_{k-1} | 0.32 / 0.53 / 0.58 / 0.57 / 0.62 | 0.46 / 0.64 / 0.69 / 0.68 / 0.60 |
| top-1% pair overlap of S_k with the CLR score | 0.32 / 0.21 / 0.16 / 0.12 / 0.10 | 0.46 / 0.33 / 0.28 / 0.24 / 0.17 |
| median effective support 1 / sum A_ij^2 | 4.8 / 5.3 / 6.0 / 6.6 / 7.1 | 4.5 / 5.4 / 6.2 / 6.9 / 7.4 |

What the numbers show: the candidate lists never empty and no row collapses at any of five layers, so nothing breaks outright. What changes progressively is (i) reachability, the unreachable share of P^9 grows from 0.3% to 17-18% by layer 4-5 on K562 (to 11% at layer 5 on M3D), so the log readout rests on fewer pairs per row; (ii) the stationary distribution drifts away from the layer-1 one (Spearman 0.74 at layer 2, 0.29 at layer 5 on K562), and its concentration first rises (K562 pi max/median 18,460 to 24,090, the top-1% in-weight share 0.136 to 0.189 at layer 2) and then falls; on M3D the concentration rises through layer 3; (iii) the MI-hub dependence of the score disappears at layer 2 (Spearman of the row mean of S_k with the MI row mean 0.30 to 0.01 on K562), so the stacked score is less a function of the MI hub degree, while the hub-stratified differences to a3r stay about equal to the raw ones; (iv) the top-1% pairs of S_k keep drifting from the CLR evidence (overlap 0.32, 0.21, 0.16, 0.12, 0.10). The scored degradation is already visible at layer 2 (CollecTRI, KD, noRP on K562; RegulonDB all on M3D) and larger at layer 3 (DoRothEA, SC on M3D). These are descriptive and are not tied to the truth scores beyond coincidence of ordering.

## 5. Plain-language summary

Feeding a3r's output back through the same attention step does not give a measurable gain on the K562 complexes (about +0.003 AUPR, 8.2x against 8.0x, interval including 0), the same as just walking one layer deeper (c1), and it costs the regulon truths (CollecTRI, DoRothEA A-B) and the knockdown-response truth. A third layer is no better than the second. On E. coli, stacking raises the RegulonDB-style complex truth (+0.04 to +0.05 AUPR, one interval excluding 0, one including it) and that rise is not reproduced by walking deeper, but it costs RegulonDB all-evidence regulons and, at three layers, SC regulons too, so the non-inferiority rule on regulons fails. Neither two nor three layers meet the literal rule, on either dataset, with either hub stratifier. The pure repeat of the CLR on a3r's output is not a no-op: it loses on K562 complexes. Layers do not break mechanically (no empty rows, no collapse), but reachability falls and the hub structure of the walk is rewritten from layer 2 onward.

## 6. Caveats

- a3r was chosen from about 20 readouts on K562 truths (K562 is optimistic) and was not adopted by the PI as the default readout; K562 CORUM and DoRothEA were used in earlier tests; stacked variants add to that selection. M3D is a tuning set. Nothing here is confirmatory and nothing was selected.
- The M3D complex truth is small (1,150 positive pairs, 54 complexes) and its gene-jackknife intervals are wide; the E. coli complex source was the Complex Portal fallback of the 26 Oct preregistration and cannot stand in for the preregistered EcoCyc truth.
- Jackknife intervals hold the networks fixed and are descriptive. The 100-replicate stream is the one of the source pipelines; there are about 70 contrasts per table, so intervals that exclude 0 should not be counted as findings.
- t = 9 was transferred to every layer without tuning, and layer-2 depths longer than 9 and a layer-3 depth grid were not run, so "stacking fails" means "stacking at t = 9 per layer (t2 in {2, 4, 9}) fails the literal rule".
- The matched-depth control matches total reach (t_1 + t_2 and so on), not the number of softmax operators; a3r2 uses two top-50 sparsifications, c1 uses one.
- c2 has two forms fixed in the plan (c2_z primary, c2_log); the "nearly idempotent" expectation was not met by either on K562 CORUM.

## 7. Equivalence gates (all PASS, tolerance 1e-9 unless stated)

- **K562** (`equivalence_gate.csv`, `equivalence_gate_stack.csv`): (a) CLR and K10 against the `robust_z` V0 caches, 7 truths x 3 scopes x (AUPR + 100 replicates), max |diff| 0; (b) CLR and attention against `baselines_k562ess_truthgenes_genie3/metrics.csv`, 3.9e-16; (c) CLR, K10, a3r against the `corum_expert2` jackknife record, 4.2e-16, a3r CORUM_all enrichment 7.9533 (prints 7.95); (d) a3r minus K10 point difference and jackknife interval against `metrics_robust.csv` for DoRothEA A-B, CollecTRI, CORUM, KD response, 4.2e-17, the CORUM line printing +0.059 [+0.033, +0.084] as in `docs/plans/corum_analysis.md`; for CORUM noRP the point difference only, 8.7e-19 (see deviation 1); (e) the generic stacking code with N = 1 equals the package-path a3r and operator exactly (max |diff| 0).
- **M3D** (`equivalence_gate_stack.csv`): K10 AUPR against `baselines_m3d/metrics.csv` 2.5e-16 and K10 enrichments (raw, hub MI, hub pi) against `robust_z` V0 att9 3.6e-15; CLR, K10, a3r pooled AUPR (complexes, SC, all) against `results/clr_logpt_test/truth/m3d/metrics.csv` 5.0e-16; all 100 stored jackknife replicates of those three methods x three truths 5.0e-16; hub-pi AUPR and replicates against the `_hubstrat` rows 4.2e-16 / 5.0e-16; a3r AUPR, d and SE against `docs/clr_logpt_verdict/endpoints_all.csv` (P1, P3, S1, S3) 5.0e-16; (e) N = 1 identity, max |diff| 0. a3r complexes 0.19209 and SC 0.16222 as in the verdict.

## 8. Jobs

| job | run | declared | peak | note |
|---|---|---|---|---|
| `5490_stk_smoke_k562`, `5491_stk_smoke_m3d` | 20:35-20:37 | 8 GB each | 0.8, 0.4 GB | code smoke test on a 600-gene subset (truths read on the subset, in-sample; nothing reported; output deleted) |
| `5500_stk_k562` | 20:38-20:47, rc 1 | 34 GB | 13.3 GB | aborted at gate (d): see deviation 1; no stacked readout was scored |
| `5500_stk_k562_retry` | 20:48-21:14 (26 min) | 34 GB | 25.7 GB | the K562 run reported here; CLR, K10, a3r caches reused from the first run |
| `5510_stk_m3d` | 20:38-20:49 (11 min) | 17 GB | 12.4 GB | the M3D run reported here |
| `5520_stk_summarize` | 21:14-21:15 | 6 GB | < 0.1 GB | contrasts, pass/fail, tables |

No OOM kill and no scheduler retry occurred (`_retry` of 5500 was a manual resubmission after the gate abort).

## 9. Deviations from the plan

1. **Gate (d), CORUM noRP.** The plan required the a3r minus K10 point difference and interval for CORUM noRP to match `metrics_robust.csv` to 1e-9. The point difference matches (8.7e-19) but the interval differs by 2.1e-4: that record drew the noRP halves right after the five original truths, while the `robust_z` / ideas 4-7 stream (used here, and gated at 0 against the stored `robust_z` replicates in gate (a)) has 100 delete-20% draws in between. This is a different jackknife stream, not an implementation difference. The job aborted at that gate (as the plan requires) before any stacked readout was scored; the check was changed to compare the point difference only for noRP (commit 1ee8ff8) and the job resubmitted. Timing relative to results: before any stacked score existed.
2. Layers 4 and 5 (label-free diagnostics only) and the `diagnostics_c1_depth` / `diagnostics_t2_grid` tables are as the plan section 7 describes; the plan's "13 readouts" wording was corrected to 11 before the commit.
3. The K562 truth list scored is the ten truths of `ideas47_lib.R` (including CORUM_RP, merged and merged noRP), a superset of the pass-rule truths. M3D scores complexes, SC and all.
4. Nothing else changed after the first stacked score: no readout, depth, truth, stratifier or pass rule was altered.

## 10. Disk

`results/explore_stack/k562` 0.6 MB and `results/explore_stack/m3d` 0.2 MB (CSV, md and logs force-added to git; rds readout caches under `cache/`, about 0.3 MB, not committed; no matrices kept). Row in `docs/disk_budget.md`.
