. jobs_pc/common.sh
# mem: 14
# cpus: 8
# tcaCLR scoring/reporting code validation on the CONSUMED DREAM5 net4 yeast study (Complex Portal yeast): inputs from the stored net4 caches, gates, score stage (net4 truth order, net4 seeds 20261114 / 20261115), report with the stored series- and group-level replicates imported, comparison with results/complex_confirm. Consumed data only (no HepG2, Kemmeren, CYC2008, Complex Portal human).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
IN=results/tcaclr_confirm/_inputs_net4; OUT=results/tcaclr_confirm/validation_net4; SEL=docs/complex_confirm_selections
mkdir -p $OUT
[ -f $IN/scores/GRNBoost2.scores ] || Rscript analysis/tcaclr_confirm/validate_net4_prep.R $IN
COMMON="--scores $IN/scores --truth_dir $SEL --members $SEL/yeast_complexes_members.tsv --excluded $SEL/yeast_complexes_excluded_pairs.tsv --sha $SEL/yeast_complexes.sha256 --counts $SEL/yeast_complexes_counts.csv --out $OUT --validation 1 --order net4 --jack_seed 20261114 --cjack_seed 20261115 --arm A --workers 8"
Rscript analysis/tcaclr_confirm/score_confirm.R --stage gates $COMMON
Rscript analysis/tcaclr_confirm/score_confirm.R --stage score $COMMON --overwrite 1
FIXED=$(Rscript -e 'e <- read.csv("results/complex_confirm/estimate.csv", stringsAsFactors = FALSE); cat(e$S_fixed[e$row == "P"])')
Rscript analysis/tcaclr_confirm/score_confirm.R --stage report $COMMON --exp_import "series=results/complex_confirm/net4/expjack/replicates_series.csv,group=results/complex_confirm/net4/expjack/replicates_group.csv" --exp_primary series --names_new a3r --names_ref attention --context "yeast complexes" --kind "preregistered exploratory estimation; not a confirmatory test" --fixed_sentence "$FIXED"
Rscript analysis/tcaclr_confirm/validate_net4_compare.R $OUT
