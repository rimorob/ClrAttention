. jobs_pc/common.sh
# mem: 24
# cpus: 8
# Global-temperature ladder (docs/plans/global_temperature_ladder_plan.md), edge stage: real-data temperature table (V2's tau and the effective partners per tau; no truth), permutation edge of the registered reference operator, V2's tau and the 13 ladder temperatures (one fit per permuted matrix, 10 permutations, registered seeds), t table for the edge rule RE, reference edge checked against the stored edge file (stops on a difference). Exploratory, descriptive; spent truths. Forked permutation workers: declared memory at least 2.5 x the measured whole-tree peak (README).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/global_temp_ladder/gt_run.R --mode edge --ds m3d --workers 4 --threads 2
