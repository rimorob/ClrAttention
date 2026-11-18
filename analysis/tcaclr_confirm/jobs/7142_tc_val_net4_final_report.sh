. jobs_pc/common.sh
# mem: 6
# cpus: 2
# after: 7122_tc_val_net4_refit_chip_plan,7102_tc_val_net4_score_final
# tcaCLR code validation, consumed net4: final report stage with every experiment source side by side (stored series and group replicates, fresh series refits, fresh profile-level refits under the plan's convention), then the comparison with the stored study including the fresh refits and their hashes.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
IN=results/tcaclr_confirm/_inputs_net4; OUT=results/tcaclr_confirm/validation_net4; SEL=docs/complex_confirm_selections
COMMON="--scores $IN/scores --truth_dir $SEL --members $SEL/yeast_complexes_members.tsv --excluded $SEL/yeast_complexes_excluded_pairs.tsv --sha $SEL/yeast_complexes.sha256 --counts $SEL/yeast_complexes_counts.csv --out $OUT --validation 1 --order net4 --jack_seed 20261114 --cjack_seed 20261115 --arm A --workers 2"
FIXED=$(Rscript -e 'e <- read.csv("results/complex_confirm/estimate.csv", stringsAsFactors = FALSE); cat(e$S_fixed[e$row == "P"])')
Rscript analysis/tcaclr_confirm/score_confirm.R --stage report $COMMON --exp_import "series=results/complex_confirm/net4/expjack/replicates_series.csv,group=results/complex_confirm/net4/expjack/replicates_group.csv" --exp_primary series --names_new a3r --names_ref attention --context "yeast complexes" --kind "preregistered exploratory estimation; not a confirmatory test" --fixed_sentence "$FIXED"
Rscript analysis/tcaclr_confirm/validate_net4_compare.R $OUT series_refit
