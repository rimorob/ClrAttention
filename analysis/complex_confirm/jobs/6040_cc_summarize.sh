. jobs_pc/common.sh
# mem: 4
# cpus: 2
# after: 6030_cc_net4_expjack
# complex_confirm summariser (docs/plans/complex_confirm_design.md, job 6040): writes estimate.csv ONCE (refuses if it exists) and estimate_final.csv, tables.md.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
Rscript analysis/complex_confirm/cc_run.R --mode summarize
