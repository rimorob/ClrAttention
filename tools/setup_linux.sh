#!/bin/bash
# One-time setup on Linux (tested target: Ubuntu, dual Xeon E5-2690 v4):
# R dependencies, package build (OpenMP via R's own flags), data downloads,
# test suite. Run from anywhere; operates on this checkout.
# Log: tools/setup_linux.log.
#
# Prerequisites (once, needs sudo; not done by this script):
#   sudo apt-get install -y r-base r-base-dev build-essential git curl unzip \
#        libcurl4-openssl-dev libssl-dev
# Optional, faster dense algebra (tcrossprod, solve) than the reference BLAS;
# the scripts pin it to 1 thread per worker:
#   sudo apt-get install -y libopenblas0-pthread
# The RegulonDB extract is not downloadable by script: copy the folder
# data/RegulonDBExtract from the Mac checkout (4.6 MB), e.g.
#   scp -r <mac>:"<mac checkout>/data/RegulonDBExtract" data/
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
LOG="$ROOT/tools/setup_linux.log"
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
die() { echo "!! $*"; exit 1; }

command -v Rscript >/dev/null || die "R not found: sudo apt-get install -y r-base r-base-dev build-essential"
echo "== R: $(Rscript -e 'cat(R.version.string, R.version$arch)')"
echo "== $(. /etc/os-release 2>/dev/null; echo "$PRETTY_NAME"), $(uname -m); $(g++ --version | head -1)"
echo "== CPU: $(lscpu | grep -E '^Model name' | sed 's/.*: *//'); sockets $(lscpu | awk -F: '/^Socket/ {gsub(/ /,"",$2); print $2}'), cores/socket $(lscpu | awk -F: '/^Core\(s\) per socket/ {gsub(/ /,"",$2); print $2}'), threads/core $(lscpu | awk -F: '/^Thread\(s\) per core/ {gsub(/ /,"",$2); print $2}')"
echo "== RAM: $(free -g | awk '/^Mem/ {print $2}') GB"

# R's default per-user library (R puts it on the search path by itself once
# it exists), so no sudo and no shell configuration are needed.
export R_LIBS_USER="$(Rscript -e 'cat(path.expand(Sys.getenv("R_LIBS_USER")))')"
mkdir -p "$R_LIBS_USER"
Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths())); for (p in c("R6","Rcpp","Matrix","foreach","doParallel","testthat","pkgload")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib, repos = "https://cloud.r-project.org", Ncpus = 8)' \
  || die "installing R dependencies failed"
R CMD INSTALL --preclean -l "$R_LIBS_USER" "$ROOT" || die "package build failed"
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); library(clr); print(clr_openmp_info()); cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")'

mkdir -p data
if [ ! -f data/E_coli_v4_Build_6/E_coli_v4_Build_6_chips907probes4297.tab ]; then
  echo "== downloading M3D E. coli v4 build 6 (117 MB)"
  curl -L --fail -o data/E_coli_v4_Build_6.tar.gz http://m3d.mssm.edu/norm/E_coli_v4_Build_6.tar.gz \
    || die "M3D download failed"
  tar -xzf data/E_coli_v4_Build_6.tar.gz -C data || die "M3D unpack failed"
fi
[ -d data/beeline/Networks ] || tools/get_beeline.sh || die "BEELINE download failed"
[ -f data/RegulonDBExtract/RISet.tsv ] || echo "!! data/RegulonDBExtract missing: copy it from the Mac (see header)"

echo "== test suite"
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); library(clr); r <- as.data.frame(testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)); if (sum(r$failed)) quit(status = 1)' \
  || die "tests failed"
echo "== setup OK (log: $LOG); packages in $R_LIBS_USER"
