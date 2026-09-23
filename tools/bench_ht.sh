#!/bin/bash
# Does hyper-threading help? Runs the same chips-bootstrap draws twice, with
# the MI kernel's OpenMP threads on physical cores only (CLR_HT unset) and on
# both hyper-threads of each core (CLR_HT=1). One draw per worker, so each
# run is one full round of the real workload (G = 4,297 genes). Prints the
# wall time of each; the numbers themselves are checked to be identical.
# Usage (repo root, on the workstation):  tools/bench_ht.sh
cd "$(dirname "$0")/.." || exit 1
W=$(Rscript -e 'source("analysis/parallel.R"); u <- max(1L, .par_cores() - 2L); r <- .par_total_ram_gb(); cat(if (is.na(r)) u else min(u, max(1L, floor((r * 0.8 - 1.5) / 3.5))))' 2>/dev/null)
echo "workers: $W  (physical cores: $(Rscript -e 'source("analysis/parallel.R"); cat(.par_cores(), "x", .par_threads_per_core(), "threads")' 2>/dev/null))"
for ht in 0 1; do
  out=results/bench_ht_$ht; rm -rf "$out"
  t0=$(date +%s)
  CLR_HT=$ht Rscript analysis/chips_bootstrap.R --Bcluster "$W" --Brep 0 --out "$out" > "$out.log" 2>&1 || { echo "run CLR_HT=$ht failed; see $out.log"; exit 1; }
  t1=$(date +%s)
  echo "CLR_HT=$ht: $((t1 - t0)) s for $W draws + 2 point estimates"
  grep "compute:" "$out.log"
done
Rscript -e 'a <- read.csv("results/bench_ht_0/draws.csv"); b <- read.csv("results/bench_ht_1/draws.csv"); m <- merge(a, b, by = c("draw","evidence","method")); cat("identical results:", isTRUE(all.equal(m$aupr.x, m$aupr.y)), "\n")'
