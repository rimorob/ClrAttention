. jobs_pc/common.sh
# mem: 4
# cpus: 1
# Arm A: collect stage alone (waits for the Arm A jackknife job to end; the collect line of the expjack job lacked --allow_missing 1 and stopped before writing anything).
set -e
until grep -q 'END 7360_tc_expjack_armA ' jobs_pc/runner.log; do sleep 60; done
Rscript analysis/tcaclr_confirm/score_confirm.R --stage collect --arm A --scores results/tcaclr_confirm/armA/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armA_complex_portal_members.tsv --excluded docs/tcaclr_confirm_selections/armA_complex_portal_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armA.sha256 --counts docs/tcaclr_confirm_selections/armA_complex_portal_counts.csv --perturbed docs/tcaclr_confirm_selections/armA_perturbed_genes_networkhalf.txt --out results/tcaclr_confirm/result_A --allow_missing 1
