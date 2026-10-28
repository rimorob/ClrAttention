#!/bin/bash
# Move aside every cached GRNBoost2 output (nothing is deleted) into results/grnboost2_superseded/<run>/, and
# copy the result files that the rerun will overwrite, so the superseded numbers stay available for the disclosure.
# baselines.R and beeline_baselines.R reuse grnboost2_scores.rds / grnboost2_edges.tsv when present, hence the move.
set -e; cd "$(dirname "$0")/.."
S=results/grnboost2_superseded; mkdir -p $S
for d in results/baselines_m3d results/baselines_ecomac results/baselines_k562ess results/baselines_k562ess_genie3 \
         results/baselines_k562ess_truthgenes results/baselines_k562ess_truthgenes_genie3; do
  [ -d $d ] || continue; t=$S/$(basename $d); mkdir -p $t
  for f in grnboost2_edges.tsv grnboost2_scores.rds; do [ -e $d/$f ] && mv -n $d/$f $t/; done
  for f in metrics.csv jack_replicates.csv per_regulon.csv per_regulon_jack.csv regulon_members.csv baselines.log; do
    [ -e $d/$f ] && cp -an $d/$f $t/; done
  echo "$d -> $t"
done
b=results/beeline_baselines; t=$S/beeline_baselines; mkdir -p $t
for c in $b/*_[0-9]*/; do n=$(basename $c); mkdir -p $t/$n
  for f in grnboost2_edges.tsv grnboost2.log; do [ -e $c$f ] && mv -n $c$f $t/$n/; done; done
for f in metrics.csv comembership.csv summary.csv beeline_baselines.log; do [ -e $b/$f ] && cp -an $b/$f $t/; done
echo "$b -> $t"
cat > $S/README.txt <<'TXT'
GRNBoost2 outputs superseded on 27 Oct 2026 (early stopping under scikit-learn 1.9.1 broke arboreto's stopping rule;
see docs/plans/grnboost2_diagnostic.md). Moved here (edges, scores) or copied (metrics, jackknife, per-regulon tables)
before the legacy-early-stopping rerun. Nothing here is used by any script; kept for the supplement's disclosure.
TXT
