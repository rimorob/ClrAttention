. jobs_pc/common.sh
# mem: 30
# cpus: 8
# Width ladder (docs/plans/width_ladder_plan.md), edge stage: permutation edge of the reference operator and of the 13 ladder widths (one fit per permuted matrix, 10 permutations, registered seeds), t table for the edge rule RE, reference edge checked against the stored edge file (stops on a difference). No truth is read. Exploratory. Estimated about 50 min wall, tree peak measured 11.3 GB (2 workers x 4 threads, 2 permutations x 3 operators). OPTIONAL (separate job).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/width_ladder/wl_run.R --mode edge --ds hepg2 --workers 2 --threads 4
