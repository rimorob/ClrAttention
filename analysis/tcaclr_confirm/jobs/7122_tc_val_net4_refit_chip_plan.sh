. jobs_pc/common.sh
# mem: 8
# cpus: 8
# after: 7132_tc_scale_refit_A
# tcaCLR code validation, consumed net4: 100 profile-level refits (each of the 536 chips its own unit, replicate b with seed 20261118 + b - 1: the plan, draft 3b) through the new expjack stage. Gives the plan-convention SE_exp on net4 next to the series- and group-level values of the net4 study.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
IN=results/tcaclr_confirm/_inputs_net4; OUT=results/tcaclr_confirm/validation_net4; SEL=docs/complex_confirm_selections
COMMON="--scores $IN/scores --truth_dir $SEL --members $SEL/yeast_complexes_members.tsv --excluded $SEL/yeast_complexes_excluded_pairs.tsv --sha $SEL/yeast_complexes.sha256 --out $OUT --validation 1 --order net4 --jack_seed 20261114 --cjack_seed 20261115 --arm A --workers 8 --threads 8"
Rscript analysis/tcaclr_confirm/score_confirm.R --stage expjack $COMMON --expr $IN/expr.rds --exp_name chip_plan --exp_scheme plan --exp_seed 20261118
Rscript analysis/tcaclr_confirm/score_confirm.R --stage collect $COMMON --exp_name chip_plan
