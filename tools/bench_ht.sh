#!/bin/bash
# Does hyper-threading raise throughput? Runs the real chips-bootstrap
# workload (G = 4,297 genes) under two plans and reports draws per hour:
#   CLR_HT=0  one worker per physical core (cores - 2), RAM permitting;
#   CLR_HT=1  workers may also use hyper-threads, up to the RAM limit, and
#             spare hyper-threads go to the MI kernel's OpenMP threads.
# Each plan runs exactly two full rounds (2 x its worker count draws), so
# neither gets a partial last round. The draws are the same bootstrap
# draws in both plans; their overlap is checked to give identical numbers.
# RAM left to the user: CLR_RESERVE_RAM_GB (default max(6 GB, 15% of RAM)).
# Usage (repo root, on the workstation):  tools/bench_ht.sh
cd "$(dirname "$0")/.." || exit 1
for ht in 0 1; do
  W=$(CLR_HT=$ht Rscript -e 'source("analysis/parallel.R"); cat(par_plan(3.3)$W)' 2>/dev/null)
  N=$((2 * W))
  out=results/bench_ht_$ht; rm -rf "$out"
  t0=$(date +%s)
  CLR_HT=$ht Rscript analysis/chips_bootstrap.R --Bcluster "$N" --Brep 0 --out "$out" > "$out.log" 2>&1 \
    || { echo "run CLR_HT=$ht failed; see $out.log"; exit 1; }
  t1=$(date +%s)
  grep "compute:" "$out.log"
  # the two point estimates also run; count them as draws
  echo "CLR_HT=$ht: $((N + 2)) draws in $((t1 - t0)) s = $(( (N + 2) * 3600 / (t1 - t0) )) draws/hour"
done
Rscript -e 'a <- read.csv("results/bench_ht_0/draws.csv"); b <- read.csv("results/bench_ht_1/draws.csv"); m <- merge(a, b, by = c("draw","evidence","method")); cat("shared draws:", length(unique(m$draw)), " identical numbers:", isTRUE(all.equal(m$aupr.x, m$aupr.y)), "\n")'
