. jobs_pc/common.sh
# mem: 26
# cpus: 8
# complex_confirm gate Ga parts 1-3 (docs/plans/complex_confirm_design.md): library identity vs d8fdb58; K562 and M3D rebuild and scoring of spent truths against stored values; spent data only, no net4, no yeast truth.
rc=0
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1
Rscript analysis/complex_confirm/gate_a_k562.R --out results/complex_confirm/gate_a/k562 --threads 8 --workers 8 || rc=1
Rscript analysis/complex_confirm/gate_a_m3d.R --out results/complex_confirm/gate_a/m3d --threads 8 --workers 8 || rc=1
Rscript analysis/complex_confirm/cc_run.R --mode gate_a || rc=1
exit $rc
