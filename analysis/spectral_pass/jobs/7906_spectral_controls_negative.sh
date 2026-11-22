. jobs_pc/common.sh
# mem: 4
# cpus: 1
# Spectral pass controls 2: per-gene permuted-expression negative control on unit 1 of each development dataset (eigenvalues above the edge, components and rings found at the edge-rule k and forced k = 5, 10, 20). No truth file. Output results/spectral_pass/controls/.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for d in hESC:hESC_500 hHep:hHep_500 K562:1 HepG2:1 RPE1:1; do
  nice -n 19 Rscript analysis/spectral_pass/sp_controls.R --mode negative --dataset ${d%%:*} --unit ${d##*:} --out results/spectral_pass/controls --krules E,5,10,20 --threads 1
done
