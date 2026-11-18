. jobs_pc/common.sh
# mem: 12
# cpus: 8
# after: 7350_tc_gates_interim_armA
# Arm A: experiment-component jackknife, 100 delete-half replicates over profiles (seed 20261118 + b - 1), CLR, aCLR, tcaCLR and |r| refitted per replicate (estimate: 5.4 h at 16 threads, about 8 h at 8).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript analysis/tcaclr_confirm/score_confirm.R --stage expjack --arm A --scores results/tcaclr_confirm/armA/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armA_complex_portal_members.tsv --excluded docs/tcaclr_confirm_selections/armA_complex_portal_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armA.sha256 --counts docs/tcaclr_confirm_selections/armA_complex_portal_counts.csv --perturbed docs/tcaclr_confirm_selections/armA_perturbed_genes_networkhalf.txt --out results/tcaclr_confirm/result_A --upstream_gates results/tcaclr_confirm/result_A/gates.csv --expr results/tcaclr_confirm/armA/expression.tsv --workers 8 --threads 8
Rscript analysis/tcaclr_confirm/score_confirm.R --stage collect --arm A --scores results/tcaclr_confirm/armA/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armA_complex_portal_members.tsv --excluded docs/tcaclr_confirm_selections/armA_complex_portal_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armA.sha256 --counts docs/tcaclr_confirm_selections/armA_complex_portal_counts.csv --perturbed docs/tcaclr_confirm_selections/armA_perturbed_genes_networkhalf.txt --out results/tcaclr_confirm/result_A
