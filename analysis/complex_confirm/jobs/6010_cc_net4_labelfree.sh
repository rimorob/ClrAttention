. jobs_pc/common.sh
# mem: 10
# cpus: 8
# after: 6000_cc_gate_a
# complex_confirm label-free stage on net4 (docs/plans/complex_confirm_design.md, Gates Ga parts 4/4b/5 and Gb, job 6010): guard (no prior a3r score on net4), operator, P^9, CLR/attention/a3r caches with SHA-256, shared-operator check, feasibility gate, decoy diagnostic, one timed label-free refit. Nothing is scored against any truth.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript analysis/complex_confirm/cc_run.R --mode labelfree --workers 8
Rscript analysis/complex_confirm/cc_run.R --mode refit_timing --workers 8
