. jobs_pc/common.sh
# mem: 12
# cpus: 8
# after: 8104_gt_edge_kemmeren
# Global-temperature ladder (docs/plans/global_temperature_ladder_plan.md), score stage: gate 1 (registered operator) and gate 2 (V2, M3D and EcoMAC only) inline, stop on failure; then the 13 temperatures + V2 reference under R9 (t = 9) and RE (t from the edge table), multi-scale scores MS1-3 for both rules, comparators CLR / tcaCLR / tcaCLR2, every metric of the plan. No forked workers (BLAS threads only). Exploratory, descriptive; spent truths. Main-process peak only (no workers).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/global_temp_ladder/gt_run.R --mode score --ds kemmeren --threads 8
