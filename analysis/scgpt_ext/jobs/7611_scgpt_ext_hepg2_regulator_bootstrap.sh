. jobs_pc/common.sh
# mem: 24
# cpus: 6
# after: 7602_scgpt_ext_hepg2_freeze_export
# HELD FOR THE PI 'GO': real truths (DoRothEA A-B, CollecTRI, Complex Portal complexes; knockdown response secondary). Submit only after the freeze manifest is committed. Run once. Regulator bootstrap within HepG2 (B = 1,000, descriptive; scGPT-only contrasts) by analysis/scfm_score_k562.R --tag HepG2, as the K562 job 1023 (same resources: 6 forked workers; 24 GB declared as that job, measured tree peak PSS 9.2 GB at 4 workers). Output results/scfm_score_hepg2/.
set -e
export OMP_NUM_THREADS=6 MKL_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6
mkdir -p results/scfm_score_hepg2
nice -n 19 Rscript analysis/scfm_score_k562.R --tag HepG2 --out results/scfm_score_hepg2 --truth real --boot 1000 --workers 6
