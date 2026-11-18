. jobs_pc/common.sh
# mem: 6
# cpus: 2
# after: 7311_tc_fit_armB
# Arm B: interim gates (hash guard, universe, truth counts and n_eff, power class, AP-implementation identity) on the score files present (CLR, aCLR, tcaCLR, |r|, hub variable); label-free except for the AP identity check, which prints only a difference.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -f results/tcaclr_confirm/armB/perturbed_deleted_orfs.txt ] || tail -n +2 docs/tcaclr_confirm_selections/armB_panel_deleted_orfs.tsv | cut -f2 | sort -u > results/tcaclr_confirm/armB/perturbed_deleted_orfs.txt
Rscript analysis/tcaclr_confirm/score_confirm.R --stage gates --arm B --scores results/tcaclr_confirm/armB/scores --truth_dir docs/tcaclr_confirm_selections --members docs/tcaclr_confirm_selections/armB_cyc2008_members.tsv --excluded docs/tcaclr_confirm_selections/armB_cyc2008_excluded_pairs.tsv --sha docs/tcaclr_confirm_selections/armB.sha256 --counts docs/tcaclr_confirm_selections/armB_cyc2008_counts.csv --perturbed results/tcaclr_confirm/armB/perturbed_deleted_orfs.txt --out results/tcaclr_confirm/result_B --allow_missing 1
