. jobs_pc/common.sh
# mem: 4
# cpus: 1
# Spectral-embedding exploratory pass (docs/plans/spectral_embedding_pass_plan.md), development dataset RPE1: the 35 variants (5 embeddings x 7 read-outs; ledger in the plan) on every unit with the real truths (DoRothEA A-B, CollecTRI, Complex Portal complexes, KD response where it exists) by the metrics code of score_ext.R, after the gate that the permuted-data edge and t2 equal results/scgpt_ext/RPE1/tcaclr2_edge_t2.csv. NOT SUBMITTED by the builder: submit only when the lead has committed the plan. One R process, no forked workers (tree = main process); measured peak RSS about 1.0 GB. Output results/spectral_pass/RPE1/.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/spectral_pass/sp_score.R --datasets RPE1 --out results/spectral_pass --truth real --variants all --threads 1
