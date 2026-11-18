. jobs_pc/common.sh
# mem: 8
# cpus: 8
# after: 7351_tc_gates_interim_armB
# Arm B (rerun of 7361, which stopped because the stage wanted all seven methods present; the four refitted methods are all it uses): experiment-component jackknife, 100 delete-half replicates over profiles (seed 20261118 + b - 1), CLR, aCLR, tcaCLR and |r| refitted per replicate (estimate: 2.3 h at 16 threads, about 3.5 h at 8).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript analysis/tcaclr_confirm/score_confirm.R --stage expjack --arm B --scores results/tcaclr_confirm/armB/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armB_cyc2008_members.tsv --excluded docs/tcaclr_confirm_selections/armB_cyc2008_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armB.sha256 --counts docs/tcaclr_confirm_selections/armB_cyc2008_counts.csv --perturbed results/tcaclr_confirm/armB/perturbed_deleted_orfs.txt --out results/tcaclr_confirm/result_B --upstream_gates results/tcaclr_confirm/result_B/gates.csv --expr results/tcaclr_confirm/armB/expression.tsv --workers 8 --threads 8 --allow_missing 1
Rscript analysis/tcaclr_confirm/score_confirm.R --stage collect --arm B --scores results/tcaclr_confirm/armB/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armB_cyc2008_members.tsv --excluded docs/tcaclr_confirm_selections/armB_cyc2008_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armB.sha256 --counts docs/tcaclr_confirm_selections/armB_cyc2008_counts.csv --perturbed results/tcaclr_confirm/armB/perturbed_deleted_orfs.txt --out results/tcaclr_confirm/result_B
