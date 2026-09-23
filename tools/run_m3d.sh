#!/bin/bash
# Full M3D x RegulonDB run: 907-chip compendium (primary) and the 466
# replicate-averaged experiments (sensitivity). Usage:
#   tools/run_m3d.sh [RegulonDB extract dir] [B]
# The extract dir must contain RISet.tsv, TUSet.tsv, OperonSet.tsv and
# NetworkSigmaGene.tsv (RegulonDB "Datasets" downloads).
set -euo pipefail
cd "$(dirname "$0")/.."
RDB="${1:-data/RegulonDBExtract}"
B="${2:-100}"
for f in RISet.tsv TUSet.tsv OperonSet.tsv NetworkSigmaGene.tsv; do
  [ -f "$RDB/$f" ] || { echo "missing $RDB/$f"; exit 1; }
done
for SET in chips avg; do
  Rscript analysis/run_m3d_regulondb.R --m3d data/E_coli_v4_Build_6 \
    --rdb "$RDB" --set "$SET" --B "$B" --mi-null --out "results/$SET"
done
