. jobs_pc/common.sh
# mem: 16
# cpus: 8
# tcaCLR2 step 1 (docs/plans/tcaclr2_definition.md): the permuted-data noise edge of the registered operator (keff 10, cap 50) and the walk depth t2 = max(1, round(ln 0.05 / ln edge)) for the five datasets (K562, net4 yeast, HepG2 Arm A, Kemmeren Arm B, E. coli M3D); 10 per-gene permutations, fixed seeds, no truth read. Exploratory. Estimated about 35 min wall, 4 workers x about 2 GB (Arm A: 2 workers x about 5.2 GB).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
analysis/tcaclr2/run_tcaclr2.sh edge k562_corum net4_cp_yeast armA_cp_human armB_cyc2008 m3d_tf_regulons
