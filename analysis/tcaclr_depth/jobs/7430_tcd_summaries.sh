. jobs_pc/common.sh
# mem: 8
# cpus: 1
# after: 7410_tcd_stage1_m3d,7411_tcd_stage1_k562,7420_tcd_stage2_k562_corum,7421_tcd_stage2_net4_cp_yeast,7422_tcd_stage2_armA_cp_human,7423_tcd_stage2_armB_cyc2008,7424_tcd_stage2_m3d_tf_regulons
# tcaCLR depth/width/clip grid: stage 1 summary (argmax, plateau, registered comparison, descriptive agreement) and stage 2 heatmaps with the registered and best cells. Seconds.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/stage1_summary.R --files results/tcaclr_depth/stage1_m3d/stage1_M3D.csv,results/tcaclr_depth/stage1_k562/stage1_K562.csv --out results/tcaclr_depth/stage1_summary
Rscript analysis/tcaclr_depth/stage2_summary.R --dirs results/tcaclr_depth/stage2_k562_corum,results/tcaclr_depth/stage2_net4_cp_yeast,results/tcaclr_depth/stage2_armA_cp_human,results/tcaclr_depth/stage2_armB_cyc2008,results/tcaclr_depth/stage2_m3d_tf_regulons --out results/tcaclr_depth/stage2_summary
