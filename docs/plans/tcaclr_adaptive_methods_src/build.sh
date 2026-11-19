#!/bin/bash
# Rebuild the explainer: figures (python3 + matplotlib, seeded) then HTML -> PDF with headless Chrome.
# Usage: nice -n 19 bash docs/plans/tcaclr_adaptive_methods_src/build.sh        (add "pdf" to skip the figures)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
if [ "$1" != "pdf" ]; then
  for f in fig1_filter fig2_pa fig3_graphsignal fig4_islands fig5_closed fig6_flow; do nice -n 19 python3 $f.py; done
fi
nice -n 19 google-chrome --headless=new --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf="$HERE/../tcaclr_adaptive_methods_explainer.pdf" "file://$HERE/explainer.html"
