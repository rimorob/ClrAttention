. jobs_pc/common.sh
# mem: 6
# cpus: 8
# Extended SOTA comparison, m3d: tcaCLR (t = 9) and tcaCLR2 (t2 = 25) columns next to aCLR, CLR, GENIE3, PIDC, GRNBoost2, |r| (and Hub only) on the 18-cell rows of m3d (2 truths: SC, all), with the paired delete-half gene jackknife (100 replicates, the baselines.R halves). Gate (stops before writing): tcaCLR and tcaCLR2 reproduce results/tcaclr2_scoreboard/cells_reference_tcaclr_all.csv and cells_tcaclr2_all.csv, the old methods reproduce results/baselines_*/metrics.csv and jack_replicates.csv, all to 1e-9. tcaCLR2 is exploratory; intervals resample genes only. Output results/sota_extended/{rows,jack,gate,run}_m3d.csv (never overwrites). Measured: 229 s at 100 replicates, peak RSS 3.84 GB, 43 gate checks, max difference 4.9e-15.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/sota_extended.R --matrix m3d --threads 8
