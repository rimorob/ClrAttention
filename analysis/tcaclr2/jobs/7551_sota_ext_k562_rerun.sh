. jobs_pc/common.sh
# mem: 36
# cpus: 4
# Extended SOTA comparison, k562: tcaCLR (t = 9) and tcaCLR2 (t2 = 29) columns next to aCLR, CLR, GENIE3, PIDC, GRNBoost2, |r| (and Hub only) on the 18-cell rows of k562 (5 truths: DoRothEA AB, ABC, CollecTRI, CORUM, KD response), with the paired delete-half gene jackknife (100 replicates, the baselines.R halves). Gate (stops before writing): tcaCLR and tcaCLR2 reproduce results/tcaclr2_scoreboard/cells_reference_tcaclr_all.csv and cells_tcaclr2_all.csv, the old methods reproduce results/baselines_*/metrics.csv and jack_replicates.csv, all to 1e-9. tcaCLR2 is exploratory; intervals resample genes only. Output results/sota_extended/{rows,jack,gate,run}_k562.csv (never overwrites). Rerun of the OOM-killed 7550_sota_ext_k562 with 4 threads and 36 GB (the first runs used 8 forked workers, whose memory the smoke peak did not count). Smoke: 742 s at 100 replicates, main-process peak RSS 9.3 GB, 43 gate checks, max difference 5.1e-15.
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/sota_extended.R --matrix k562 --threads 4
