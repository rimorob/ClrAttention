# Learned-attention comparator, step 3: production extraction (2026-11-02) -- RUNNING LOG

**No truth file was opened and nothing was scored in this step.** Preregistration change log: step-2 decisions E1-E12 committed 2026-11-01 (17ca611). The E13 "step-3 freeze" entry is NOT yet committed; it follows the final QC.

## Progress log (updated as steps finish)

| Step | Status | Commit / job |
|---|---|---|
| 3.0 change-log entry E1-E12 | done | 17ca611 |
| 3.1 production scripts (`analysis/scfm_extract_scprint.py`, `scfm_extract_scgpt.py`, `scfm_controls.py`) | done; both reproduce the step-1 hESC_1000 64-cell arrays bitwise (max abs diff 0) | 84586ef |
| 3.2 gates | (a) shim 1.0.13 vs scprint 2.2.0 on mHSC-E_500 arm A, 1071 cells: max abs diff 0, bitwise; (b) determinism hESC_500 x2: bitwise; (c) large-v1/small-v1 load strictly; (d) GNInfer returns a G x G matrix. All pass | d6a3b32, jobs 1001, 1002 |
| 3.3 scPRINT medium-v1.5 arm A, 14 cases | done, QC pass (7 min) | job 1003 |
| 3.4 arm B, 8 cases | done, QC pass (5 min) | job 1004 |
| 3.5 R1 (GNInfer) A 14 + B 8 | done, QC pass (8 min) | job 1005 |
| 3.6 controls | done, QC pass (63 min) | job 1006 |
| 3.7 scGPT + seed-1 rerun | done, QC pass (8 min) | job 1007 |
| 3.8 K562 | **BLOCKED, not run** (see below) | |
| 3.9 large-v1 / small-v1 | queued/running | job 1008 |
| 3.10 fine-tuning | not done (report only) | |
| 3.11 QC document and E13 | pending | |

## 3.8 K562 is blocked by the no-truth rule

The preregistered K562 gene set is "CollecTRI TFs with >= 10 targets detected in >= 5% of control cells, plus the 1,000 most variable detected genes". The TF list needs the CollecTRI regulon sizes (`data/human_truth/collectri.tsv`, a truth file under this task's rules). `results/baselines_k562ess/expression.tsv` is not a substitute: its columns are perturbation pseudobulks, not control cells. Also the raw single-cell Replogle file is a truth source for the knockdown-response modules, so it was not opened. Not run; needs a PI decision (for example: allow reading CollecTRI only to compute regulon sizes for the gene-set definition, before any scoring, or give a gene list).
