. jobs_pc/common.sh
# mem: 4
# cpus: 1
# Spectral pass controls 1: plane-level tests of the two ring families and the planted-structure positive control (synthetic 12-D embedding and the real HepG2 unit-1 embedding with planted rows; no truth file). Output results/spectral_pass/controls/. Run before or with the development jobs.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/spectral_pass/sp_controls.R --mode planes --out results/spectral_pass/controls --threads 1
nice -n 19 Rscript analysis/spectral_pass/sp_controls.R --mode positive --dataset HepG2 --unit 1 --expr files --out results/spectral_pass/controls --krules E,5,10 --threads 1
