#!/bin/bash
# One-time setup on macOS (Apple Silicon): dependencies, OpenMP build,
# M3D download, test suite. Run from anywhere; operates on this checkout.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"

command -v Rscript >/dev/null || { echo "R not found: install from https://cran.r-project.org/bin/macosx/"; exit 1; }
echo "== R: $(Rscript -e 'cat(R.version.string)')"

if [ ! -d /opt/homebrew/opt/libomp ]; then
  command -v brew >/dev/null || { echo "Homebrew not found (needed for libomp): https://brew.sh"; exit 1; }
  brew install libomp
fi

Rscript -e 'for (p in c("R6","Rcpp","Matrix","testthat","pkgload")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")'

echo "== building clr with OpenMP"
R_MAKEVARS_USER="$ROOT/tools/Makevars.macos-openmp" R CMD INSTALL --preclean "$ROOT"
Rscript -e 'library(clr); i <- clr_openmp_info(); print(i); if (i[["openmp_max_threads"]] == 0) stop("clr was built WITHOUT OpenMP; check libomp / tools/Makevars.macos-openmp")'

mkdir -p data
if [ ! -f data/E_coli_v4_Build_6/E_coli_v4_Build_6_chips907probes4297.tab ]; then
  echo "== downloading M3D E. coli v4 build 6 (117 MB)"
  curl -L --fail -o data/E_coli_v4_Build_6.tar.gz http://m3d.mssm.edu/norm/E_coli_v4_Build_6.tar.gz
  tar -xzf data/E_coli_v4_Build_6.tar.gz -C data
fi

echo "== test suite"
Rscript -e 'library(clr); testthat::test_dir("tests/testthat", reporter = "summary")'
echo "== setup OK"
