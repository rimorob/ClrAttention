#!/bin/bash
# Equivalence gate (27 Oct 2026): stock arboreto 0.1.6 + scikit-learn 1.2.2 (Python 3.10 (3.11 fails to import dask.dataframe with old dask), ~/venvs/gb2-legacy-check)
# vs tools/grnboost2.py legacy mode under scikit-learn 1.9.1 (~/venvs/clr-baselines), on one small BEELINE case, 3 seeds.
# Arms: stock_dask (arboreto.algo.grnboost2 itself), stock_pool (stock per-target routine, gives trees per target), legacy (new default).
set -e; cd "$(dirname "$0")/.."
IN=${1:-results/beeline_baselines/mDC_500/expression.tsv}; OUT=${2:-results/grnboost2_equivalence}; W=${3:-6}
OLD=~/venvs/gb2-legacy-check/bin/python; NEW=~/venvs/clr-baselines/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p $OUT
for s in 1 2 3; do
  [ -s $OUT/stock_pool_s$s.tsv ] || nice -n 19 $OLD tools/grnboost2.py $IN $OUT/stock_pool_s$s.tsv --workers $W --seed $s --early_stop sklearn | tee -a $OUT/log.txt
  [ -s $OUT/legacy_s$s.tsv ]     || nice -n 19 $NEW tools/grnboost2.py $IN $OUT/legacy_s$s.tsv --workers $W --seed $s | tee -a $OUT/log.txt
  [ -s $OUT/stock_dask_s$s.tsv ] || nice -n 19 $OLD tools/grnboost2_variants.py $IN $OUT/stock_dask_s$s.tsv --mode dask --workers $W --seed $s | tee -a $OUT/log.txt
done
