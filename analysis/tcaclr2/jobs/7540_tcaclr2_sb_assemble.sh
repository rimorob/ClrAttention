. jobs_pc/common.sh
# mem: 2
# cpus: 1
# after: 7530_tcaclr2_sb_score_m3d 7531_tcaclr2_sb_score_k562 7532_tcaclr2_sb_score_ecomac 7533_tcaclr2_sb_score_rpe1 7534_tcaclr2_sb_score_beeline
# tcaCLR2 on the scoreboard rows, step 3: assemble every tcaCLR2 cell (family cells plus the DREAM5 yeast, HepG2 and Kemmeren cells already scored in results/tcaclr2), cross-check against results/tcaclr2, write t2_by_dataset.csv. Output results/tcaclr2_scoreboard/cells_tcaclr2_all.csv, cells_reference_tcaclr_all.csv, crosscheck.csv, t2_by_dataset.csv. Seconds.
set -e
nice -n 19 Rscript analysis/tcaclr2/assemble_tcaclr2_scoreboard.R
