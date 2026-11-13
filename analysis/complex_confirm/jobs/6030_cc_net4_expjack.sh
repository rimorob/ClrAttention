. jobs_pc/common.sh
# mem: 30
# cpus: 40
# after: 6020_cc_net4_truth
# complex_confirm experiment-level jackknives (docs/plans/complex_confirm_design.md, job 6030): 100 series-level and 100 group-level refits of CLR, attention and a3r, each scored on all truths; CLR and attention hashes checked against both DREAM5 manifests.
# Implementation: the 200 independent, deterministic refits (seeds fixed by the design) are split over 5 lanes of 8 threads each (40 of 56 threads, about 71%); each lane measured about 3.3 GB at the refit peak in 6010 (1.4 x 5 lanes x ~4.3 GB incl. scoring = 30 GB declared; 1.4 x one refit peak = 4.7 GB, below the 16 GB trigger).
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
rc=0; pids=()
for l in 0 1 2 3 4; do
  Rscript analysis/complex_confirm/cc_run.R --mode expjack --lane $l --nlanes 5 --workers 8 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || rc=1; done
Rscript analysis/complex_confirm/cc_run.R --mode expjack --collect 1 || rc=1
exit $rc
