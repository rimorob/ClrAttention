. jobs_pc/common.sh
# mem: 30
# cpus: 6
# Width ladder (docs/plans/width_ladder_plan.md), edge stage: permutation edge of the reference operator and of the 13 ladder widths (one fit per permuted matrix, 10 permutations, registered seeds), t table for the edge rule RE, reference edge checked against the stored edge file (stops on a difference). No truth is read. Exploratory. Estimated about 13 min wall, tree peak interpolated about 12 GB (3 workers x 2 threads).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/width_ladder/wl_run.R --mode edge --ds yeast --workers 3 --threads 2
