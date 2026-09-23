#!/bin/bash
# Full M3D x RegulonDB run: 907-chip compendium (primary) and the 466
# replicate-averaged experiments (sensitivity). Usage:
#   tools/run_m3d.sh data/NetworkRegulatorGene.tsv [B]
set -euo pipefail
cd "$(dirname "$0")/.."
REG="${1:?path to RegulonDB TF-gene file}"
B="${2:-100}"
for SET in chips avg; do
  Rscript analysis/run_m3d_regulondb.R --m3d data/E_coli_v4_Build_6 \
    --regulondb "$REG" --set "$SET" --B "$B" --mi-null --out "results/$SET"
done
