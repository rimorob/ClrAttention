. jobs_pc/common.sh
# mem: 6
# cpus: 2
# after: 7310_tc_fit_armA
# Arm A: interim gates (hash guard, universe, truth counts and n_eff, power class, AP-implementation identity) on the score files present (CLR, aCLR, tcaCLR, |r|, hub variable); label-free except for the AP identity check, which prints only a difference.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript analysis/tcaclr_confirm/score_confirm.R --stage gates --arm A --scores results/tcaclr_confirm/armA/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armA_complex_portal_members.tsv --excluded docs/tcaclr_confirm_selections/armA_complex_portal_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armA.sha256 --counts docs/tcaclr_confirm_selections/armA_complex_portal_counts.csv --perturbed docs/tcaclr_confirm_selections/armA_perturbed_genes_networkhalf.txt --out results/tcaclr_confirm/result_A --allow_missing 1
