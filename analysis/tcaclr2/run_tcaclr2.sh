#!/bin/bash
# tcaCLR2 runner (docs/plans/tcaclr2_definition.md). Run from the repo root.
#   analysis/tcaclr2/run_tcaclr2.sh edge  [dataset ...]   edge + t2 of each dataset  -> results/tcaclr2/edge_<dataset>.tsv
#   analysis/tcaclr2/run_tcaclr2.sh score [dataset ...]   registered tcaCLR (t = 9) and tcaCLR2 side by side, ratio, jackknife -> results/tcaclr2/score_<dataset>/ and results/tcaclr2/summary.csv
#   analysis/tcaclr2/run_tcaclr2.sh files [dataset ...]   tcaCLR2.scores (+ genes.txt) in the make_scores.R format -> results/tcaclr2/scores_<dataset>/ (not part of the two job files)
# Datasets (default all five, in this order): k562_corum net4_cp_yeast armA_cp_human armB_cyc2008 m3d_tf_regulons. Output only under results/tcaclr2/.
# Env: TC2_WORKERS (edge permutation workers, default 4; armA_cp_human uses 2 for memory), TC2_THREADS (clr threads per edge worker, default 2; armA 4), TC2_SCORE_THREADS (default 8).
set -e
stage=$1; shift || true
ds=${*:-k562_corum net4_cp_yeast armA_cp_human armB_cyc2008 m3d_tf_regulons}
out=results/tcaclr2; mkdir -p $out
for d in $ds; do
  case $stage in
    edge)
      w=${TC2_WORKERS:-4}; th=${TC2_THREADS:-2}; if [ "$d" = armA_cp_human ]; then w=2; th=4; fi
      Rscript analysis/tcaclr2/edge_depth.R --dataset $d --perms 10 --workers $w --threads $th --out $out ;;
    score)
      Rscript analysis/tcaclr2/score_tcaclr2.R --truth $d --edge $out/edge_$d.tsv --out $out/score_$d --threads ${TC2_SCORE_THREADS:-8} ;;
    files)
      Rscript analysis/tcaclr2/make_scores_tcaclr2.R --dataset $d --edge $out/edge_$d.tsv --out $out/scores_$d --threads ${TC2_SCORE_THREADS:-8} ;;
    *) echo "stage must be edge, score or files"; exit 2 ;;
  esac
done
if [ "$stage" = score ]; then Rscript analysis/tcaclr2/collect_tcaclr2.R --out $out; fi
