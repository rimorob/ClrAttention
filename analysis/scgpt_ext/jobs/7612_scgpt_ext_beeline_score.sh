. jobs_pc/common.sh
# mem: 8
# cpus: 4
# HELD FOR THE PI 'GO': real truths (DoRothEA A-B, CollecTRI, Complex Portal complexes; knockdown response secondary). Submit only after the freeze manifest is committed. Run once. Scores hESC and hHep (BEELINE cases, scGPT extraction of step 3, step4_scoring exports; all cells): same methods and truths as the single-cell datasets. Adds aCLR (descriptive column). Gates first: aCLR (at the stored held-out depth t* and at t = 9), tcaCLR and tcaCLR2 on the full hESC_500 and hHep_500 BEELINE matrices reproduce the stored scoreboard values, edge files and t2 to 1e-9, and the label-free depth function reproduces the stored t*. Output results/scgpt_ext/hESC, hHep.
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/scgpt_ext/score_ext.R --datasets hESC,hHep --out results/scgpt_ext --truth real --methods all --threads 4
