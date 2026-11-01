#!/usr/bin/env Rscript
# SYNTHETIC ONLY. Summary of sim_p1_rescore.R: per cell (gen, f, b), over 4 paired seeds:
# geometric-mean ratios and mean absolute differences of R vs K10, R vs CLR, K10 vs CLR on P1 / S1 / REG,
# at t = 9 and at each operator's own criterion-A t*; seeds better; t* values.
OUT <- "analysis/synthetic_width/out/p1_rescore"
fs <- list.files(file.path(OUT, "runs"), pattern = "\\.csv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, utils::read.csv))
stopifnot(all(x$sanity_maxreldiff < 1e-8))
x$cell <- sprintf("%s f%02d b%d", x$gen, round(100 * x$f), x$b)
get <- function(op, depth) { y <- x[x$op == op & (depth == "na" | x$depth == depth), ]; y[order(y$cell, y$seed), ] }
cmp <- function(a, b, lab) {
  do.call(rbind, lapply(split(seq_len(nrow(a)), a$cell), function(ix) {
    A <- a[ix, ]; B <- b[match(paste(A$cell, A$seed), paste(b$cell, b$seed)), ]
    do.call(rbind, lapply(c("P1", "S1", "REG"), function(e) data.frame(cell = A$cell[1], contrast = lab, endpoint = e,
      ratio_gm = exp(mean(log(A[[e]] / B[[e]]))), d_mean = mean(A[[e]] - B[[e]]), d_min = min(A[[e]] - B[[e]]),
      seeds_better = sum(A[[e]] > B[[e]]), level_ref = mean(B[[e]]),
      t_a = paste(A$t, collapse = "/"), t_b = paste(B$t, collapse = "/"))))
  }))
}
res <- rbind(cmp(get("R", "t9"), get("K10", "t9"), "R vs K10, t = 9"),
             cmp(get("R", "tstar"), get("K10", "tstar"), "R vs K10, own t*"),
             cmp(get("R", "tstar"), get("CLR", "na"), "R (own t*) vs CLR"),
             cmp(get("K10", "tstar"), get("CLR", "na"), "K10 (own t*) vs CLR"),
             cmp(get("R", "t9"), get("CLR", "na"), "R (t = 9) vs CLR"),
             cmp(get("K10", "t9"), get("CLR", "na"), "K10 (t = 9) vs CLR"))
utils::write.csv(res, file.path(OUT, "summary_p1.csv"), row.names = FALSE)
tst <- aggregate(t ~ cell + op, data = x[x$depth == "tstar", ], FUN = function(v) paste(v, collapse = "/"))
utils::write.csv(tst, file.path(OUT, "tstar_by_cell.csv"), row.names = FALSE)
cat(sprintf("%d runs; sanity max rel diff %.1e\n", nrow(x) / 5, max(x$sanity_maxreldiff)))
