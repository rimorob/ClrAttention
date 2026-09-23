#!/bin/bash
# Full M3D x RegulonDB run on the 466 replicate-averaged experiments (the
# 907-chip set counts technical replicates as samples, inflating N; run it
# explicitly with SET=chips if wanted). Usage:
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
# Always run the code in this checkout: reinstall clr (OpenMP build on macOS)
# so an older installed version can never be used by mistake.
if [ "$(uname)" = "Darwin" ]; then
  LIBOMP=""
  for p in /opt/homebrew/opt/libomp /usr/local/opt/libomp /usr/local; do
    if [ -f "$p/lib/libomp.dylib" ] && [ -f "$p/include/omp.h" ]; then LIBOMP="$p"; break; fi
  done
  if [ -n "$LIBOMP" ]; then
    sed "s#@LIBOMP@#$LIBOMP#g" tools/Makevars.macos-openmp > tools/.Makevars.generated
    R_MAKEVARS_USER="$PWD/tools/.Makevars.generated" R CMD INSTALL --preclean . > results_install.log 2>&1 \
      || { echo "package install failed; see results_install.log"; exit 1; }
  else
    R CMD INSTALL --preclean . > results_install.log 2>&1 || { echo "install failed"; exit 1; }
  fi
else
  R CMD INSTALL --preclean . > results_install.log 2>&1 || { echo "install failed"; exit 1; }
fi
Rscript -e 'library(clr); i <- clr_openmp_info(); cat("clr", as.character(packageVersion("clr")), "OpenMP threads:", i[["openmp_max_threads"]], "\n")'

SET="${SET:-avg}"
Rscript analysis/run_m3d_regulondb.R --m3d data/E_coli_v4_Build_6 \
  --rdb "$RDB" --set "$SET" --B "$B" --out "results/$SET"
