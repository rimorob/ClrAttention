#!/usr/bin/env Rscript
# EXPLORATORY (docs/plans/explore_idea3_plan.md, sections 3, 5, 7). Step 0: hub deciles (MI row mean and attention pi),
# per-gene knee of log k_i(t) on the real full-data operator and on 50 degree-preserving rewirings of it, per-decile knees,
# and the between-decile statistics of the null. No truth and no held-out split is used.
# Usage (repo root): Rscript analysis/explore_idea3/idea3_knee_null.R --dataset m3d|k562 [--nrew 50] [--workers 6] [--quick 0]
args <- commandArgs(TRUE)
opt <- list(dataset = "m3d", out = NULL, nrew = 50L, workers = 6L, quick = 0L, threads = 8L, alpha = 0.5, seed = 20261112L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); if (!key %in% names(opt)) stop("unknown argument: ", args[k])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[k + 1]) else args[k + 1] }
for (k in c("nrew", "workers", "quick", "threads", "seed")) opt[[k]] <- as.integer(opt[[k]])
if (is.null(opt$out)) opt$out <- file.path("results/explore_idea3", opt$dataset)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
script_dir <- (function() { f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis/explore_idea3" })()
source(file.path(script_dir, "idea3_lib.R")); i3_source_all()
logf <- file.path(opt$out, "knee_null.log"); say <- function(...) i3_say(logf, ...)
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

D <- i3_load(opt$dataset, opt$quick); X <- D$X; genes <- D$genes; G <- length(genes); rm(D)
say(sprintf("%s: %d genes x %d samples", opt$dataset, G, ncol(X)))
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
f$calibrate(method = "normal", combine = "stouffer")
M <- f$mi; f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L); A <- f$operator
f$release(); rm(f)
A0 <- rz_operator(rz_scores(M, "V0")); stopifnot(max(abs(A0 - A)) < 1e-8); rm(A0)
H_mi <- (rowSums(M) - diag(M)) / (G - 1); rm(M)
P <- rz_lazy_P(A, opt$alpha); pi_att <- rz_stationary(P); rm(P)
dec <- data.frame(gene = genes, H_mi = H_mi, pi_att = pi_att, dec_mi = i3_decile(H_mi), dec_pi = i3_decile(pi_att))
utils::write.csv(dec, file.path(opt$out, "deciles.csv"), row.names = FALSE)
say(sprintf("deciles written; Spearman(MI row mean, pi) = %.3f; decile agreement (same decile) = %.3f",
            stats::cor(H_mi, pi_att, method = "spearman"), mean(dec$dec_mi == dec$dec_pi)))

knee_of <- function(Aop) {
  P <- rz_lazy_P(Aop, opt$alpha); st <- i3_knee_new(G)
  i3_walk(P, I3_TMAX, function(Pt, t) st <<- i3_knee_add(st, Pt, t))
  i3_knee_idx(st$logk)
}
dec_knee <- function(kidx, d) vapply(1:10, function(k) stats::median(log2(I3_GRID[kidx[d == k]])), 0)
stats_of <- function(l2) c(sd = stats::sd(l2), range = diff(range(l2)), spearman = suppressWarnings(stats::cor(1:10, l2, method = "spearman")))

t0 <- Sys.time(); kr <- knee_of(A); say(sprintf("real operator knee walk: %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
real_gene <- data.frame(gene = genes, knee_t = I3_GRID[kr], dec_mi = dec$dec_mi, dec_pi = dec$dec_pi)
utils::write.csv(real_gene, file.path(opt$out, "knee_gene_real.csv"), row.names = FALSE)

invisible(i3_rewire_cpp()); invisible(gc(verbose = FALSE))
nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))
say(sprintf("%d rewirings, %d workers", opt$nrew, nc))
res <- parallel::mclapply(seq_len(opt$nrew), function(r) {
  Ar <- i3_rewire(A, opt$seed + r); k <- knee_of(Ar); rm(Ar); invisible(gc(verbose = FALSE)); k }, mc.cores = nc, mc.preschedule = FALSE)
bad <- vapply(res, function(x) !is.integer(x) || length(x) != G, NA)
if (any(bad)) stop("rewiring workers failed: ", sum(bad), " ", paste(unique(vapply(res[bad], function(x) paste(class(x), collapse = "/"), "")), collapse = ","))

out <- list(); nullst <- list()
for (sname in c("mi", "pi")) {
  d <- dec[[paste0("dec_", sname)]]
  real <- dec_knee(kr, d); rw <- do.call(rbind, lapply(res, function(k) dec_knee(k, d)))
  out[[sname]] <- data.frame(stratifier = sname, decile = 1:10, real_median_log2_knee = real, real_knee_snapped = vapply(real, i3_snap, 0),
                             rewired_mean_log2_knee = colMeans(rw), rewired_q025 = apply(rw, 2, stats::quantile, 0.025), rewired_q975 = apply(rw, 2, stats::quantile, 0.975))
  rs <- t(apply(rw, 1, stats_of)); r0 <- stats_of(real)
  nullst[[sname]] <- rbind(data.frame(stratifier = sname, rep = 0L, t(r0)), data.frame(stratifier = sname, rep = seq_len(nrow(rs)), rs))
}
# ADDED after the 500-gene smoke test showed a degenerate rewired null (every gene's knee identical on a rewired operator);
# logged as a deviation in the results document: label-permutation null on the REAL per-gene knees (2000 random
# gene-to-decile labellings, decile sizes preserved), the sampling noise of a decile median of per-gene knees.
set.seed(opt$seed)
for (sname in c("mi", "pi")) {
  d <- dec[[paste0("dec_", sname)]]
  pm <- t(vapply(1:2000, function(r) stats_of(dec_knee(kr, sample(d))), numeric(3)))
  nullst[[paste0(sname, "_perm")]] <- data.frame(stratifier = paste0(sname, "_labelperm"), rep = 0:2000, rbind(stats_of(dec_knee(kr, d)), pm))
}
out <- do.call(rbind, out); nullst <- do.call(rbind, nullst)
utils::write.csv(out, file.path(opt$out, "knee_decile.csv"), row.names = FALSE)
utils::write.csv(nullst, file.path(opt$out, "knee_null_stats.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(c("mi", "pi", "mi_labelperm", "pi_labelperm"), function(s) { x <- nullst[nullst$stratifier == s, ]; r <- x[x$rep == 0, ]; n <- x[x$rep > 0, ]
  data.frame(stratifier = s, real_sd = r$sd, null_sd_mean = mean(n$sd), null_sd_q95 = stats::quantile(n$sd, 0.95), frac_null_sd_ge_real = mean(n$sd >= r$sd),
             real_range = r$range, null_range_mean = mean(n$range), frac_null_range_ge_real = mean(n$range >= r$range),
             real_spearman = r$spearman, null_spearman_mean = mean(n$spearman), n_rewirings = nrow(n),
             real_sd_above_q95 = r$sd > stats::quantile(n$sd, 0.95)) }))
utils::write.csv(summ, file.path(opt$out, "knee_null_summary.csv"), row.names = FALSE)
print(out, row.names = FALSE); print(summ, row.names = FALSE)
say("knee null done")
