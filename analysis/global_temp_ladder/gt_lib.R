# Global-temperature ladder (docs/plans/global_temperature_ladder_plan.md): shared helpers. EXPLORATORY, DESCRIPTIVE; spent truths. Source from the repo root.
# The operator with ONE global temperature is variant V2 of the recipe ladder (analysis/recipe_ladder/01_ladder.R): the candidate set (each gene's top-50 POSITIVE CLR scores,
# self excluded, order(decreasing = TRUE) ties as there), the weights exp((s - max s)/tau) normalised per row, gl_finish_A, lazy weight 0.5. gt_topmat / gt_effv / gt_tau_for_eff /
# gt_weights are the recipe-ladder lines moved into functions UNCHANGED (checked by gate 2). Everything else (walk, read-out, metrics, truths, edge machinery) is the width-ladder
# code (analysis/width_ladder/wl_lib.R and the files it sources); the new code only loops over temperatures. No second operator or scorer.
suppressPackageStartupMessages(source("analysis/width_ladder/wl_lib.R"))

# half-octave ladder of the global temperature tau (units of the clipped CLR score); endpoints and the number of points are fixed in the plan from the probe (gt_probe.R)
GTL <- list(log2_lo = -4, log2_hi = 2, step = 0.5, cap = 50L, eff_target = 10, ms3_log2 = seq(-2.5, 1, by = 0.5), tmax = WL$tmax, t_fixed = WL$t_fixed, ref = WL$ref, out = "results/global_temp_ladder")
GTL$taus <- 2^seq(GTL$log2_lo, GTL$log2_hi, by = GTL$step)
gt_tau_id <- function(tau) sprintf("T%+.2f", log2(tau))   # id = log2 tau

## ---- candidate set and weights (recipe ladder, verbatim logic) --------------------------------------------------------------------
gt_topmat <- function(C, cap = GTL$cap) {
  G <- nrow(C); idx <- matrix(NA_integer_, G, cap); sc <- matrix(NA_real_, G, cap)
  for (g in seq_len(G)) { r <- C[g, ]; r[g] <- -Inf; o <- order(r, decreasing = TRUE)[seq_len(min(cap, G - 1L))]; o <- o[r[o] > 0]; if (length(o)) { idx[g, seq_along(o)] <- o; sc[g, seq_along(o)] <- r[o] } }
  list(idx = idx, sc = sc, rowmax = apply(sc, 1, max, na.rm = TRUE)) }
# unnormalised weights (NA -> 0), the recipe-ladder expression
gt_w <- function(Tm, tau) { w <- exp((Tm$sc - apply(Tm$sc, 1, max, na.rm = TRUE)) / tau); w[is.na(w)] <- 0; w }
# effective partners exp(entropy) of every row that has a weight (the recipe-ladder effv)
gt_effv <- function(Tm, tau) { w <- gt_w(Tm, tau); has <- rowSums(w) > 0; w <- w / pmax(rowSums(w), 1e-300); exp(-rowSums(ifelse(w > 0, w * log(w), 0)))[has] }
# V2's temperature: bisection on log tau over [log 1e-4, log 1e4], 80 iterations, mean effective partners = target (the recipe-ladder loop)
gt_tau_for_eff <- function(Tm, target = GTL$eff_target) { lo <- log(1e-4); hi <- log(1e4)
  for (it in 1:80) { mid <- (lo + hi) / 2; if (mean(gt_effv(Tm, exp(mid))) < target) lo <- mid else hi <- mid }; exp((lo + hi) / 2) }
# lazy operator P of the global-temperature softmax (recipe-ladder fromW + gl_finish_A + rz_lazy_P / gl_make_P final step)
gt_make_P <- function(Tm, tau, alpha = GL$alpha) {
  w2 <- gt_w(Tm, tau); W <- w2 / pmax(rowSums(w2), 1e-300); G <- nrow(W)
  A2 <- matrix(0, G, G); ok <- !is.na(Tm$idx); A2[cbind(row(Tm$idx)[ok], Tm$idx[ok])] <- W[ok]; A <- gl_finish_A(A2); rm(A2, W, w2)
  P <- Matrix::Diagonal(G, 1 - alpha) + alpha * Matrix::Matrix(A, sparse = TRUE); rm(A); P }
gt_eff_row <- function(Tm, tau, id, dataset) { e <- gt_effv(Tm, tau)
  data.frame(dataset = dataset, op = id, tau = tau, log2_tau = log2(tau), n_rows = length(e), eff_mean = mean(e), eff_median = stats::median(e), eff_q10 = unname(stats::quantile(e, 0.1)), eff_q90 = unname(stats::quantile(e, 0.9)),
             eff_min = min(e), eff_max = max(e), stringsAsFactors = FALSE) }

## ---- operator list (ids, taus). 'ref' = registered per-gene operator (keff 10, cap 50); 'V2' = V2's own temperature of this dataset (reference point); T<log2 tau> = ladder.
gt_ops <- function(tau_v2) {
  o <- list(list(id = "ref", kind = "ref", tau = NA_real_), list(id = "V2", kind = "tau", tau = tau_v2))
  c(o, lapply(GTL$taus, function(tau) list(id = gt_tau_id(tau), kind = "tau", tau = tau))) }
gt_P <- function(o, C, Tm) if (o$kind == "ref") gl_make_P(C, GTL$ref$keff, GTL$ref$cap) else gt_make_P(Tm, o$tau)

## ---- edge (permutation) stage: wl_edges with the operator builder swapped -----------------------------------------------------------
# One fit per permuted matrix (seeds ED seed + 7000 + p). The permuted-data operator is built with the SAME global tau (not re-solved) and the same top-50 rule.
gt_edges <- function(X, ops, perms, workers, threads, say = function(...) invisible(NULL)) {
  t0 <- Sys.time()
  one <- function(p) {
    f <- ad_fit(ad_permute_genes(X, ED$seed + 7000L + p), threads); on.exit(f$release()); C <- f$clr_scores; Tm <- gt_topmat(C)
    m <- t(vapply(ops, function(o) { v <- ad_eigs(ad_Psym(gt_P(o, C, Tm)), ED$kvals)$values; c(top = v[1L], edge = v[2L]) }, c(top = 0, edge = 0)))
    list(m = m, rss = gl_vmhwm_gb()) }
  all <- ad_lapply(seq_len(perms), one, workers)
  bad <- vapply(all, function(x) !is.list(x), NA); if (any(bad)) stop("permutation task failed: ", all[bad][[1]])
  say(sprintf("%d permutations x %d operators: %.0f s (worker peak RSS up to %.1f GB)", perms, length(ops), as.numeric(difftime(Sys.time(), t0, units = "secs")), max(vapply(all, `[[`, 0, "rss"))))
  per <- do.call(rbind, lapply(seq_along(all), function(p) data.frame(op = vapply(ops, `[[`, "", "id"), perm = p, top = all[[p]]$m[, "top"], edge = all[[p]]$m[, "edge"])))
  tab <- do.call(rbind, lapply(seq_along(ops), function(k) { o <- ops[[k]]; e <- per$edge[per$op == o$id]; tt <- wl_t_from_edge(mean(e))
    data.frame(op = o$id, tau = o$tau, edge = mean(e), edge_min = min(e), edge_max = max(e), tstar = tt[["tstar"]], t = as.integer(tt[["t"]]), capped = as.integer(tt[["capped"]]),
               top_eig_min = min(per$top[per$op == o$id]), perms = perms, stringsAsFactors = FALSE) }))
  list(tab = tab, per = per)
}
