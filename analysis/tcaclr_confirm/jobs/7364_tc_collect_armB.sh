. jobs_pc/common.sh
# mem: 4
# cpus: 1
# Arm B: collect stage alone (the 100 experiment-jackknife replicates are complete; the collect line of the expjack job lacked --allow_missing 1 and stopped before writing anything).
set -e
Rscript analysis/tcaclr_confirm/score_confirm.R --stage collect --arm B --scores results/tcaclr_confirm/armB/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armB_cyc2008_members.tsv --excluded docs/tcaclr_confirm_selections/armB_cyc2008_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armB.sha256 --counts docs/tcaclr_confirm_selections/armB_cyc2008_counts.csv --perturbed results/tcaclr_confirm/armB/perturbed_deleted_orfs.txt --out results/tcaclr_confirm/result_B --allow_missing 1
