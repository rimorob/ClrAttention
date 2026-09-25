#!/bin/bash
# STRING files for evidence-filtered BEELINE truths (analysis/beeline_truths.R):
#   v12.0 channel-level links (neighborhood, fusion, cooccurence, coexpression,
#   experimental, database, textmining) and protein names, and v11.0 "actions"
#   (mode = expression is transcriptional regulation, with direction), the
#   last STRING release that published actions. Human (9606) and mouse (10090).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/string && cd data/string
B=https://stringdb-downloads.org/download
for sp in 9606 10090; do
  for f in "protein.links.detailed.v12.0/$sp.protein.links.detailed.v12.0.txt.gz" \
           "protein.info.v12.0/$sp.protein.info.v12.0.txt.gz" \
           "protein.actions.v11.0/$sp.protein.actions.v11.0.txt.gz" \
           "protein.info.v11.0/$sp.protein.info.v11.0.txt.gz"; do
    o=$(basename "$f"); [ -s "$o" ] || curl -L --fail -sS -o "$o" "$B/$f"
  done
done
ls -la
