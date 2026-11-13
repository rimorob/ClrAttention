. jobs_pc/common.sh
# mem: 16
# cpus: 8
# after: 6010_cc_net4_labelfree
# complex_confirm TRUTH STAGE (docs/plans/complex_confirm_design.md, job 6020): primary and all secondary scores, 500 gene and 500 complex delete-half jackknife replicates, seeds 20261114 / 20261115. ONE RUN: the driver refuses if any truth-stage output exists.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript analysis/complex_confirm/cc_run.R --mode truth --workers 8 --jack 500 --cjack 500 --t 9 --seed 20261114 --cseed 20261115
