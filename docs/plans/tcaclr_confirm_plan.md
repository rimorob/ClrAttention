# Preregistration DRAFT (confirmatory): does tcaCLR beat aCLR on protein-complex co-membership in data and truths it was not selected or estimated on, and how does it compare with CLR, |r| and the established methods?

**Status: draft for review. Nothing in this plan has been scored. No tcaCLR, aCLR, CLR, |r|, GENIE3, PIDC or GRNBoost2 value exists on any dataset named here.** Items marked OPEN are decided in the review round and fixed before any run. Names are temporary: aCLR = CLR with attention (internal: CLR attention), tcaCLR = twice cooked aCLR (internal: a3r).

## 0. Why this test, and what is already used

tcaCLR has no fitted parameters, but it was selected among about 20 read-outs on K562 truths, so every K562 number is optimistic. Its evidence elsewhere:

| Evidence | Result | Status |
|---|---|---|
| *E. coli* preregistered test (Complex Portal 83333), tcaCLR minus aCLR | M3D +0.0257 (p 0.079); EcoMAC +0.0285 (p 0.051) | missed one-sided p < 0.025 |
| DREAM5 net4 yeast with Complex Portal yeast, tcaCLR minus aCLR | +0.096, total 95% interval +0.067 to +0.125; complex-level interval −0.023 to +0.215 | exploratory estimation |
| net4 by complex size | 20 or more genes +0.137; 7 to 19 genes −0.003; 3 to 6 genes −0.004 | exploratory, post hoc strata |
| RPE1 merged CORUM, tcaCLR against CLR | 6.14 against 6.77 | unexplained |

**Consumed (never to be used as untouched):** K562 CORUM, RPE1 CORUM, M3D and EcoMAC with Complex Portal 83333, the DREAM5 net4 compendium with Complex Portal yeast (559292), the DREAM5 network 3 and 4 gold standards.

## 1. Hypotheses (fixed now)

- **H1 (primary, one per arm).** On pooled complex co-membership AUPR, tcaCLR exceeds aCLR.
- **H2 (size-stratified prediction, from the theory memo; registered now but derived after the net4 look, so it is a prediction made with that knowledge).** The tcaCLR minus aCLR difference is positive in complexes of 20 or more genes, near zero in complexes of 7 to 19, and not positive in complexes of 3 to 6. A tie is |difference| < 0.01 AUPR.
- **H3 (comparators, descriptive plus one gated claim).** tcaCLR against each of CLR, |r|, GENIE3, PIDC, GRNBoost2 and aCLR, per arm. The statement "tcaCLR is the best of the seven methods in this arm" is made only if tcaCLR exceeds each of the other six with Holm correction over six one-sided tests at overall 0.025. |r| was first on *E. coli* complexes (29.2) and ahead of CLR and aCLR on yeast net4, so this statement may well fail, and the report will say so.

## 2. Arms

**Arm A (primary): human HepG2 Perturb-seq pseudobulk with EBI Complex Portal human (taxon 9606).**
- Data: GEO GSE264667, HepG2 h5ad (5.6 GB). Not yet downloaded; needs a pseudobulk build with the same recipe as the K562 profiles (recipe and counts frozen and committed before any score; OPEN: exact recipe and gene universe).
- Truth: Complex Portal 9606 (downloaded, 5.4 MB), never scored. On a proxy universe (K562 gene names) it gives 329 complexes and n_eff 31; 60% of its positive pairs are CORUM positives, 40% are not. Both numbers are recomputed on the real HepG2 universe before any score.
- Independence: new cell type and new truth set, but the same laboratory, screen design and pipeline as K562 and RPE1, and shared positives with CORUM. So this arm tests transfer across cell types and truths, not across technologies.
- Power (worker's calculation; n_eff 31): detects +0.03 (needs n_eff 15 gene-only, 27 with the net4-sized experiment component) and +0.05 (5.5 or 10). Adequate for both unless the experiment component is as large as 1.78 times the gene component (n_eff 63 needed for +0.03).
- Fallback if HepG2 fails the data gates below: Jurkat (9.4 GB), same truth.

**Arm B (secondary): yeast Kemmeren deletion compendium with CYC2008.**
- Data: GEO superseries GSE42528 (GSE42527, 700 responsive mutants, and GSE42526, 784 non-responsive; GPL11232), already in `data/yeast`: 6,169 ORFs by 1,283 profiles (the full matrix; the 652-profile responsive-only file is not used, chosen now so the choice cannot depend on a result). Independent of net4.
- Truth: CYC2008 (408 complexes, 1,627 ORFs; taken from an Internet Archive capture of 16 Feb 2020 because the authors' server did not answer; licence not verified, OPEN). On the Kemmeren universe: 213 merged complexes, 10,308 positive pairs, n_eff 10.0.
- **Independence of the truth is poor.** 89.5% of CYC2008 positive pairs are also Complex Portal yeast positives; only 41 complexes (19%) have no match; the 1,079 novel pairs have n_eff about 7. So this arm is a test on independent *data* with a largely reused *truth*, and will be described that way. The novel-pairs-only AUPR difference is reported as descriptive and is expected to be underpowered.
- Power: n_eff 10 gives power 0.63 (gene-only) or 0.40 (with the net4 experiment component) for +0.03, and 0.97 or 0.81 for +0.05. Adequate for +0.05 and for the +0.096 seen on net4; marginal for +0.03; the 20-or-more stratum has 4 to 8 effective units and cannot be tested alone.
- Known data caveats: gate F of an earlier reciprocity preregistration failed on this matrix (median chance-corrected reproducibility 0.122 for mutants against 0.296 for wild types); 376 of the 1,600 truth genes are themselves deleted in the panel; the Kemmeren expression and gene names were opened before label-free (counts, reproducibility), and no method was scored on it.

**Claim rule.** The confirmatory claim "tcaCLR exceeds aCLR on complexes in untouched data" rests on Arm A, H1. Arm B is reported with the same estimand and the limits above and is not needed for the claim. If both arms are positive the report says so; if A is positive and B is not, the report says the evidence is from one arm and gives B's interval. OPEN: whether to make both arms co-primary with a Holm split instead.

## 3. Methods (all fixed; nothing tuned)

- aCLR, tcaCLR: the exact R implementation and settings used on K562 and net4 (B-spline MI, row z-score clipped at 0, Stouffer, softmax over the top 50 with effective partners 10, lazy walk P = 0.5 I + 0.5 A, t = 9, symmetric read-out; tcaCLR adds log P^9 with unreachable pairs NA, row statistics over reachable pairs, second CLR clipped at 0, Stouffer). The commit id of the code and the operator hash are recorded in the freeze file.
- CLR: Stouffer version, as used throughout. |r|: absolute Pearson correlation.
- GENIE3: all genes as regulators, random forest, K = sqrt, 1,000 trees, `set.seed(20261016)`, `doRNG`, so results do not depend on the worker count. PIDC: `tools/pidc_chunked.jl` (bit-identical to the package). GRNBoost2: `tools/grnboost2.py`, legacy early stopping, as in the scoreboard.
- Equivalence gates, the same as in the scoreboard competitor study: on a 300-to-600-gene subset each re-implementation must match its reference exactly (or to the declared tolerance) before the full run; any failure stops the job.
- No method parameter is changed for any arm. Depth t = 9 is the label-free choice already made.

## 4. Truth build (frozen before any score)

The label-free builder conventions of the yeast net4 study apply: merged complexes by overlapping membership, positives are within-complex pairs among universe genes, complexes with fewer than three universe genes are dropped, and the paralog and cross-subunit exclusions of that study. For human, OPEN: whether a paralog exclusion applies and how it is defined (it must be label-free). The members file, excluded-pairs file and their SHA-256 hashes are committed with the freeze; the score stage refuses to run on a mismatch.

## 5. Estimand, inference and reporting (fixed before any score)

- Estimand: pooled AUPR over all in-universe pairs, with positives from the frozen truth; difference tcaCLR minus aCLR.
- Primary inference follows the *E. coli* preregistered test: resampling of genes in blocks that respect the truth's n_eff clusters, combined with resampling of experiments (profiles). One-sided p < 0.025 for H1.
- Complex-level resampling is reported as a descriptive and conservative interval. On net4 it was 5.6 times the gene-level standard error and its interval included 0 for the +0.096 estimate, so it will probably include 0 here too; this is stated in advance. OPEN: whether this convention is acceptable as primary (the review should challenge it).
- Experiment resampling unit: for HepG2, perturbation; for Kemmeren, profile (there is no series structure, so the experiment-to-gene standard error ratio is unknown; the plan uses three planning values, 0, 0.89 and 1.78).
- H2: the same difference in three size classes (3 to 6, 7 to 19, 20 or more genes), reported with intervals; no p-values for H2 in Arm B.
- All pairwise differences among the seven methods are reported as estimates with total intervals; no unadjusted p-values other than H1 and the Holm family in H3.

## 6. Gates (in order; each writes a CSV; any failure stops the job before any truth is read)

1. Data gates, label-free: gene universe and profile counts match the frozen build; expression-matrix hash matches; no duplicated gene names (as for RPE1, where HSPA14 was duplicated and dropped by a rule fixed beforehand).
2. Truth gates: n_eff and complex counts match the freeze; Arm A needs n_eff of at least 15 on the real universe; below that, the arm is reported as underpowered and not as a test.
3. Method gates: equivalence checks of Section 3; deterministic seeds; peak memory logged.
4. Score gates: AUPR reproduces on a fixed smoke subset to 1e-12 against the registered reference.

## 7. Compute and scheduling

Everything runs through the resource-aware scheduler at 70 to 80% utilisation, one run, no retries beyond the automatic out-of-memory retry. Estimates at 32 threads for Kemmeren (6,169 genes by 1,283 profiles): |r|, CLR, aCLR and tcaCLR under 3 minutes combined (peak 6 to 8 GB); PIDC 10 to 15 minutes (5 to 6 GB); **GENIE3 about 5 hours and GRNBoost2 about 2.5 to 4.2 hours (range 3 to 10), both over the 2-hour limit and requiring the PI's approval before they run.** Arm A has similar costs (GENIE3 took an estimated 6.6 hours at 8,745 genes on RPE1) plus the pseudobulk build; re-estimated once the universe is fixed and again submitted for approval.

## 8. Disclosure and what is not claimed

- Drafter's exposure: the drafter knows the net4 and *E. coli* estimates and the net4 size strata; H2 is therefore a prediction made with that knowledge. The drafter has seen structure counts only for CYC2008, Complex Portal human and Kemmeren; no method has been scored on any of them.
- Not claimed, whatever the result: that tcaCLR is a theoretically justified or optimal estimator (the theory review calls it a heuristic with derivable properties); a replication of the *E. coli* test; any statement about direct regulatory edges; any statement about a method on a dataset where it was not run.
- Limitations stated in any report: shared laboratory and pipeline in Arm A; reused truth in Arm B; the CYC2008 licence; the Kemmeren reproducibility caveat; a K562-selected method.

## 9. Open items for the review round

1. Is Arm A the right primary, and is a separate secondary Arm B worth running given its power and truth overlap?
2. The primary inference convention (gene blocks plus experiments) against the conservative complex-level interval.
3. Whether H2 is worth registering as stated (post hoc strata), and the tie definition of 0.01.
4. Multiplicity: Holm over six comparators for the "best method" statement, and co-primary against primary-plus-secondary.
5. The HepG2 pseudobulk recipe and paralog rule for human.
6. Whether the simulation in the theory memo (R1 to R12) should be run before the real data, and whether any of its eight cheap checks should be a gate.
