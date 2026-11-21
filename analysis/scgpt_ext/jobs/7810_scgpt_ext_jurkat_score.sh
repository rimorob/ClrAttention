. jobs_pc/common.sh
# mem: 8
# cpus: 4
# after: 7802_scgpt_ext_jurkat_freeze_export
# HELD FOR THE PI 'GO': real truths (DoRothEA A-B, CollecTRI, Complex Portal complexes; knockdown response secondary). Submit only after the freeze manifest is committed. Run once. Scores Jurkat: |r|, CLR, CLR-NS, tcaCLR (t = 9), tcaCLR2 (t2 from the permuted-data edge rule recomputed on each unit's matrix) and the four scGPT readouts on the scGPT gene set, all primary truths, with the gates of score_ext.R (stops before writing on a failure). No forked workers (main-process peak equals the process tree). Output results/scgpt_ext/<dataset>/. No knockdown-response truth (no bulk file).
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/scgpt_ext/score_ext.R --datasets Jurkat --out results/scgpt_ext --truth real --methods all --threads 4
