. jobs_pc/common.sh
# mem: 3
# cpus: 2
# after: 7522_tcaclr2_sb_edge_beeline
# tcaCLR2 on the scoreboard rows, step 2e: BEELINE median AUPRC ratios (44 primary case-truth cells and the 14 cell-type ChIP cases) for tcaCLR2 (t2 per case from 7522) and the tcaCLR t = 9 gate. Exploratory. Output results/tcaclr2_scoreboard/cells_beeline.csv, gate_beeline.csv, beeline_metrics.csv, beeline_t2.csv. Measured: 197 s, peak RSS 1.2 GB, both cells gated (max difference 3e-15).
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/score_scoreboard.R --family beeline --threads 2
