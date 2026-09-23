#!/bin/bash
# One-time setup on macOS (Apple Silicon): dependencies, OpenMP build,
# M3D download, test suite. Run from anywhere; operates on this checkout.
# Everything is also logged to tools/setup_mac.log (useful for debugging).
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
LOG="$ROOT/tools/setup_mac.log"
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
die() { echo "!! $*"; exit 1; }

command -v Rscript >/dev/null || die "R not found: install from https://cran.r-project.org/bin/macosx/"
echo "== R: $(Rscript -e 'cat(R.version.string, R.version$arch)')"
echo "== macOS $(sw_vers -productVersion), $(uname -m); clang: $(clang --version | head -1)"

Rscript -e 'for (p in c("R6","Rcpp","Matrix","testthat","pkgload")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")' \
  || die "installing R dependencies failed"

# --- locate libomp -----------------------------------------------------------
LIBOMP=""
for p in /opt/homebrew/opt/libomp /usr/local/opt/libomp /usr/local; do
  if [ -f "$p/lib/libomp.dylib" ] && [ -f "$p/include/omp.h" ]; then LIBOMP="$p"; break; fi
done
if [ -z "$LIBOMP" ] && command -v brew >/dev/null; then
  brew install libomp && LIBOMP="$(brew --prefix libomp)"
fi

try_openmp_build() {
  local mv="$ROOT/tools/.Makevars.generated"
  sed "s#@LIBOMP@#$LIBOMP#g" "$ROOT/tools/Makevars.macos-openmp" > "$mv"
  echo "== building clr with OpenMP (libomp at $LIBOMP)"
  echo "   libomp: $(file -b "$LIBOMP/lib/libomp.dylib")"
  echo "   install name: $(otool -D "$LIBOMP/lib/libomp.dylib" | tail -1)"
  R_MAKEVARS_USER="$mv" R CMD INSTALL --preclean --no-test-load "$ROOT" || return 1
  local so
  so="$(Rscript -e 'cat(system.file("libs", "clr.so", package = "clr"))')"
  echo "   built: $so"; otool -L "$so" | sed 's/^/     /'
  Rscript -e 'library(clr); i <- clr_openmp_info(); print(i); if (i[["openmp_max_threads"]] == 0) quit(status = 2)' || return 1
}

serial_build() {
  echo "== building clr WITHOUT OpenMP (serial fallback)"
  R CMD INSTALL --preclean "$ROOT" || die "even the serial build failed; send tools/setup_mac.log"
  Rscript -e 'library(clr); print(clr_openmp_info())'
  echo "!! clr is SERIAL: MI builds will use one core (the full run will be ~8-10x slower)."
  echo "!! Please send tools/setup_mac.log so the OpenMP build can be fixed."
}

if [ -n "$LIBOMP" ]; then
  if ! try_openmp_build; then
    echo "!! OpenMP build or load FAILED. Diagnostics:"
    Rscript -e 'tryCatch(library(clr), error = function(e) cat("   load error:", conditionMessage(e), "\n"))'
    ls -l /opt/R/arm64/lib/libomp* /usr/local/lib/libomp* 2>/dev/null | sed 's/^/     /'
    serial_build
  fi
else
  echo "!! no libomp found and no Homebrew to install it"
  serial_build
fi

# --- data --------------------------------------------------------------------
mkdir -p data
if [ ! -f data/E_coli_v4_Build_6/E_coli_v4_Build_6_chips907probes4297.tab ]; then
  echo "== downloading M3D E. coli v4 build 6 (117 MB)"
  curl -L --fail -o data/E_coli_v4_Build_6.tar.gz http://m3d.mssm.edu/norm/E_coli_v4_Build_6.tar.gz \
    || die "M3D download failed"
  tar -xzf data/E_coli_v4_Build_6.tar.gz -C data || die "M3D unpack failed"
fi

echo "== test suite"
Rscript -e 'library(clr); r <- as.data.frame(testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)); if (sum(r$failed)) quit(status = 1)' \
  || die "tests failed"
echo "== setup OK (log: $LOG)"
