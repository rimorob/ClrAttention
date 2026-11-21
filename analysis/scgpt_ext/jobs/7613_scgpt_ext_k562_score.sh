. jobs_pc/common.sh
# mem: 8
# cpus: 4
# HELD FOR THE PI 'GO': real truths (DoRothEA A-B, CollecTRI, Complex Portal complexes; knockdown response secondary). Submit only after the freeze manifest is committed. Run once. Scores K562 (5 subsamples of 1,024 controls, existing scGPT extraction, step4_k562 exports): adds aCLR (descriptive), tcaCLR and tcaCLR2 and the Complex Portal truth to the earlier DoRothEA A-B / CollecTRI / knockdown-response scores. Output results/scgpt_ext/K562. (The legacy-row check against results/scfm_score_k562/metrics_long.csv is a separate --methods legacy call, see the README.)
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/scgpt_ext/score_ext.R --datasets K562 --out results/scgpt_ext --truth real --methods all --threads 4
