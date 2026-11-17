. jobs_pc/common.sh
# mem: 8
# cpus: 8
# after: 7101_tc_val_net4_score
# tcaCLR code validation, consumed net4: 100 series-level refits (59 series, seeds 20264000 + b, the net4 convention) through the new expjack stage; CLR, aCLR, tcaCLR, |r| refitted and scored on every truth; packed-score hashes compared with the DREAM5 manifest by the comparer. One lane, 8 threads (about 57 s per refit in the net4 study).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
IN=results/tcaclr_confirm/_inputs_net4; OUT=results/tcaclr_confirm/validation_net4; SEL=docs/complex_confirm_selections
COMMON="--scores $IN/scores --truth_dir $SEL --members $SEL/yeast_complexes_members.tsv --excluded $SEL/yeast_complexes_excluded_pairs.tsv --sha $SEL/yeast_complexes.sha256 --out $OUT --validation 1 --order net4 --jack_seed 20261114 --cjack_seed 20261115 --arm A --workers 8 --threads 8"
Rscript analysis/tcaclr_confirm/score_confirm.R --stage expjack $COMMON --expr $IN/expr.rds --exp_units $IN/units_series.tsv --exp_scheme perrep --exp_seed 20264000 --exp_name series_refit --hash_refit 1
Rscript analysis/tcaclr_confirm/score_confirm.R --stage collect $COMMON --exp_name series_refit
