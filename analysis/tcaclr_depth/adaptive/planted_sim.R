#!/usr/bin/env Rscript
# Semi-synthetic planted-module simulation (plan Revision 2, items 9, 11, 12; review section 4). EXPLORATORY; synthetic modules in real null data, no truth is read.
#   Rscript analysis/tcaclr_depth/adaptive/planted_sim.R --matrix results/baselines_m3d/expression.tsv --dataset M3D --genes 3000 --out results/tcaclr_depth/adaptive_planted_m3d
#
# Background: a real matrix (--genes random genes) with every gene's values permuted across experiments independently (real marginals, sparsity and zeros kept,
# all dependence destroyed). Planting (copula plant): for a module of m genes and dose rho, y_g = sqrt(1 - rho) z_g + sqrt(rho) f, z_g the gene's rank-normal
# background score, f one shared factor over the experiments (linear link: f ~ N(0, 1); "nonlin" link: f = (f0^2 - 1)/sqrt(2), a U-shaped link that is
# invisible to Pearson correlation but not to MI); the gene's own background values are then re-assigned in the rank order of y_g, so every planted gene keeps
# EXACTLY its real marginal. Gaussian-space pair correlation inside a module = rho. One dataset = one (module size, dose, link, replicate); --counts modules of
# that size are planted (non-overlapping), so oracle cells are per module size.
# Methods per cell (keff, t): CLR (no cell), aCLR = (P^t + P^t')/2, tcaCLR (clipped, registered read-out; grid_lib gl_variants), and diffPMI (EXPLORATORY ONLY, plan item 12):
#   W_ij = log max(P^t_ij / pi_j, 0.01) (pi = stationary mass of the lazy walk; unreachable entries floored at log 0.01), each row centred on its median over all
#   off-diagonal entries and divided by 1.4826 MAD (the row standard deviation if the MAD is 0), clipped at 0, Stouffer (Z + Z')/sqrt(2).
# Scoring: pooled co-membership AUPR of every cell, positives = the pairs inside the planted modules, negatives = every other pair, the registered tie convention
# (gl_ap_T = ap_fast, checked in the smoke test). Rules (label-free; computed from the planted matrix alone): the edge rule t* = ln(eps)/ln(lambda_edge) (eps 0.05, 0.01;
# spectral_edge.R) and the graph-signal likelihood (graph_signal_lik.R: best t, 1-SE paired, 1-SE marginal), each snapped to the nearest grid depth (log scale);
# reference "fixed_t9" = the registered depth. Regret = 1 - AUPR(chosen cell) / AUPR(oracle), against the best t at the same keff (same_keff) and against the best cell over
# the whole (keff, t) grid (all_cells), per method.
# Output (<out>/): ps_cells.tsv, ps_rules.tsv, ps_oracle.tsv, ps_regret.tsv, ps_summary.tsv (per rule and method: share of configurations with regret <= 5%, Spearman of the
# chosen depth against module size), ps_datasets.tsv (timings); datasets/<id>.rds restartable.
.ps_main <- sys.nframe() == 0L
if (!exists("ad_Psym")) source("analysis/tcaclr_depth/adaptive/adaptive_lib.R")
if (!exists("se_run")) source("analysis/tcaclr_depth/adaptive/spectral_edge.R")
if (!exists("gs_split")) source("analysis/tcaclr_depth/adaptive/graph_signal_lik.R")

## ---- planting -----------------------------------------------------------------------------------------------------------------
ps_pair_pos <- function(a, b) { lo <- pmin(a, b); hi <- pmax(a, b); as.integer((hi - 1) * (hi - 2) / 2 + lo) }   # position in the upper-triangle (gl_ut) vector
ps_within <- function(m) { m <- sort(m); ij <- utils::combn(length(m), 2L); ps_pair_pos(m[ij[1L, ]], m[ij[2L, ]]) }
ps_plant <- function(X0, size, count, dose, link, seed) {
  G <- nrow(X0); N <- ncol(X0); Xp <- ad_permute_genes(X0, seed); Zb <- ad_rank_normal(Xp)
  set.seed(seed + 1L); genes <- sample.int(G, size * count); mods <- split(genes, rep(seq_len(count), each = size))
  X <- Xp
  for (m in mods) {
    f0 <- stats::rnorm(N); fU <- (f0^2 - 1) / sqrt(2)      # linear loading / U-shaped loading (unit variance each)
    for (k in seq_along(m)) { g <- m[k]; f <- if (link == "nonlin" && k %% 2L == 0L) fU else f0
      y <- sqrt(1 - dose) * Zb[g, ] + sqrt(dose) * f; X[g, order(y)] <- sort(Xp[g, ]) }
  }
  list(X = X, modules = mods, truth = sort(unlist(lapply(mods, ps_within))))
}

## ---- diffusion PMI read-out (exploratory) ------------------------------------------------------------------------------------------------
ps_pmi_S <- function(Pt, pi, eps = 0.01) {
  G <- nrow(Pt); W <- log(pmax(Pt * rep(1 / pi, each = G), eps)); diag(W) <- NA_real_
  med <- apply(W, 1L, stats::median, na.rm = TRUE); W <- W - med
  sc <- 1.4826 * apply(abs(W), 1L, stats::median, na.rm = TRUE); sdv <- sqrt(rowSums(W^2, na.rm = TRUE) / (G - 2)); sc[!(sc > 0)] <- sdv[!(sc > 0)]
  Z <- W / sc; Z[is.na(Z) | Z < 0] <- 0; Z
}

## ---- one dataset -----------------------------------------------------------------------------------------------------------------
ps_snap <- function(t, depths) depths[which.min(abs(log(depths) - log(pmax(t, 1))))]
ps_dataset <- function(X0, size, count, dose, link, rep, keffs, depths, perms, r, threads, seed, id, say) {
  t0 <- Sys.time(); pl <- ps_plant(X0, size, count, dose, link, seed); X <- pl$X; G <- nrow(X); U <- gl_ut(G)
  lab <- logical(length(U$ut)); lab[pl$truth] <- TRUE; np <- length(lab)
  f <- ad_fit(X, threads); cells <- list()
  ap <- function(v) { v[!is.finite(v)] <- 0; gl_ap_T(v, gl_sorted(v), pl$truth) }
  add <- function(method, keff, t, v) cells[[length(cells) + 1L]] <<- data.frame(id = id, size = size, dose = dose, link = link, rep = rep, method = method, keff = keff, t = t, aupr = ap(v), base = length(pl$truth) / np, stringsAsFactors = FALSE)
  add("CLR", NA_real_, NA_integer_, f$clr_scores[U$ut])
  for (keff in keffs) {
    P <- ad_P_from_fit(f, keff); pi <- gl_stationary(P)
    gl_walk(P, depths, function(t, Pt) {
      v <- gl_variants(Pt, U); add("aCLR", keff, t, v$aCLR); add("tcaCLR", keff, t, v$tca_clip); rm(v)
      Z <- ps_pmi_S(Pt, pi); add("diffPMI", keff, t, (Z[U$ut] + Z[U$utT]) / sqrt(2)); rm(Z) })
    rm(P)
  }
  cells <- do.call(rbind, cells); f$release(); rm(f); invisible(gc(FALSE))
  ## label-free rules from the planted matrix alone
  se <- se_run(X, keffs, perms, kobs = 40L, kperm = 30L, threads = threads, workers = 1L, seed = seed + 5L, dataset = id)$summary
  un <- ad_units(X, "random", 1L, seed + 6L); gs <- gs_split(X, un$splits[[1]], un$unit, keffs, depths, r, threads, 1L)
  rules <- do.call(rbind, lapply(keffs, function(k) {
    s <- se[se$keff == k, ]; g <- gs$selected[gs$selected$keff == k, ]
    mk <- function(rule, tc) data.frame(id = id, size = size, dose = dose, link = link, rep = rep, keff = k, rule = rule, t_cont = tc, t_grid = ps_snap(tc, depths), stringsAsFactors = FALSE)
    rbind(mk("edge_eps05", s$tstar_eps05), mk("edge_eps01", s$tstar_eps01), mk("gs_best", g$t_best), mk("gs_1se", g$t_selected_1se), mk("gs_1se_marginal", g$t_selected_1se_marginal),
      mk("fixed_t9", 9)) }))
  rules$n_above_max <- se$n_above_max[match(rules$keff, se$keff)]; rules$edge_mean <- se$edge_mean[match(rules$keff, se$keff)]
  rules$gs_gain_best <- gs$selected$gain_best[match(rules$keff, gs$selected$keff)]
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")); say(sprintf("%s: %.0f s, peak RSS %.2f GB", id, secs, ad_mem()))
  list(cells = cells, rules = rules, secs = secs, peak_gb = ad_mem(), n_planted_pairs = length(pl$truth))
}

## ---- oracle, regret, summary ---------------------------------------------------------------------------------------------------------------
ps_regret <- function(cells, rules) {
  cc <- cells[cells$method != "CLR", ]; key <- paste(cc$id, cc$method)
  orc_all <- tapply(cc$aupr, key, max); orc_k <- tapply(cc$aupr, paste(key, cc$keff), max)
  oracle <- do.call(rbind, lapply(split(cc, key), function(d) { b <- d[which.max(d$aupr), ]; data.frame(id = b$id, size = b$size, dose = b$dose, link = b$link, rep = b$rep, method = b$method, best_keff = b$keff, best_t = b$t, best_aupr = b$aupr) }))
  rg <- do.call(rbind, lapply(unique(cc$method), function(m) {
    r <- merge(rules, cc[cc$method == m, c("id", "keff", "t", "aupr")], by.x = c("id", "keff", "t_grid"), by.y = c("id", "keff", "t"), all.x = TRUE)
    r$method <- m; r$aupr_chosen <- r$aupr; r$aupr <- NULL
    r$aupr_oracle_same_keff <- as.numeric(orc_k[paste(r$id, r$method, r$keff)]); r$aupr_oracle_all <- as.numeric(orc_all[paste(r$id, r$method)])
    r$regret_same_keff <- 1 - r$aupr_chosen / r$aupr_oracle_same_keff; r$regret_all_cells <- 1 - r$aupr_chosen / r$aupr_oracle_all; r }))
  clr <- cells[cells$method == "CLR", c("id", "aupr", "base")]; rg$aupr_clr <- clr$aupr[match(rg$id, clr$id)]
  rg$base <- clr$base[match(rg$id, clr$id)]; rg$oracle_over_base <- rg$aupr_oracle_all / rg$base      # near 1: nothing detectable, regret is noise
  list(oracle = oracle, regret = rg[order(rg$id, rg$method, rg$keff, rg$rule), ])
}
ps_summary <- function(rg, oracle) {
  keys <- unique(rg[, c("method", "rule", "keff")]); rows <- list()
  for (i in seq_len(nrow(keys))) {
    d <- rg[rg$method == keys$method[i] & rg$rule == keys$rule[i] & rg$keff == keys$keff[i], ]
    sp <- function(x, y) if (length(unique(x)) > 2L && length(unique(y)) > 2L) suppressWarnings(stats::cor(x, y, method = "spearman")) else NA_real_
    rows[[i]] <- data.frame(keys[i, ], configs = nrow(d), frac_regret_same_keff_le5pct = mean(d$regret_same_keff <= 0.05, na.rm = TRUE), frac_regret_all_cells_le5pct = mean(d$regret_all_cells <= 0.05, na.rm = TRUE),
      median_regret_same_keff = stats::median(d$regret_same_keff, na.rm = TRUE), median_regret_all_cells = stats::median(d$regret_all_cells, na.rm = TRUE),
      spearman_tgrid_vs_size = sp(d$t_grid, d$size))
  }
  sm <- do.call(rbind, rows)
  o <- oracle; sm_o <- do.call(rbind, lapply(split(o, o$method), function(d) data.frame(method = d$method[1], rule = "ORACLE_best_t", keff = NA_real_, configs = nrow(d), spearman_tgrid_vs_size = suppressWarnings(stats::cor(d$best_t, d$size, method = "spearman")))))
  plyr_fill <- function(a, b) { for (n in setdiff(names(a), names(b))) b[[n]] <- NA_real_; b[names(a)] }
  rbind(sm, plyr_fill(sm, sm_o))
}

if (.ps_main) {
  opt <- gl_options(commandArgs(TRUE), list(matrix = "", dataset = "X", genes = 3000L, sizes = "5,10,20,50,100", counts = "20,10,5,2,2", doses = "0.1,0.2,0.3,0.4,0.5,0.6", links = "linear,nonlin", reps = 1L,
    keffs = "10,28", depths = "1,2,4,9,18,32,64", perms = 2L, r = 300L, threads = 4L, workers = 1L, seed = AD$seed, out = "", tag = ""))
  if (!nzchar(opt$matrix) || !nzchar(opt$out)) stop("--matrix and --out are required"); if (!nzchar(opt$tag)) opt$tag <- opt$dataset
  dir.create(file.path(opt$out, "datasets"), recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, sprintf("ps_%s.log", opt$tag)))
  say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
  X0 <- ad_load(opt$matrix, opt$genes, opt$seed); say(sprintf("%s: %d genes x %d experiments (background)", opt$dataset, nrow(X0), ncol(X0)))
  num <- function(s) as.numeric(strsplit(s, ",")[[1]]); sizes <- as.integer(num(opt$sizes)); counts <- as.integer(num(opt$counts)); stopifnot(length(sizes) == length(counts))
  keffs <- num(opt$keffs); depths <- as.integer(num(opt$depths)); doses <- num(opt$doses); links <- strsplit(opt$links, ",")[[1]]
  grid <- expand.grid(size = sizes, dose = doses, link = links, rep = seq_len(opt$reps), stringsAsFactors = FALSE); grid$count <- counts[match(grid$size, sizes)]
  grid$id <- sprintf("%s_m%d_d%.2f_%s_r%d", opt$tag, grid$size, grid$dose, grid$link, grid$rep); say(sprintf("%d datasets, keffs %s, depths %s", nrow(grid), opt$keffs, opt$depths))
  one <- function(i) { fo <- file.path(opt$out, "datasets", paste0(grid$id[i], ".rds")); if (file.exists(fo)) return(readRDS(fo))
    res <- ps_dataset(X0, grid$size[i], grid$count[i], grid$dose[i], grid$link[i], grid$rep[i], keffs, depths, opt$perms, opt$r, opt$threads, opt$seed + 100L * i, grid$id[i], say)
    saveRDS(res, paste0(fo, ".tmp")); file.rename(paste0(fo, ".tmp"), fo); res }
  t0 <- Sys.time(); res <- ad_lapply(seq_len(nrow(grid)), one, opt$workers)
  bad <- vapply(res, is.character, NA); if (any(bad)) { say("ERROR: ", res[bad][[1]]); quit(status = 1) }
  cells <- do.call(rbind, lapply(res, `[[`, "cells")); rules <- do.call(rbind, lapply(res, `[[`, "rules")); rr <- ps_regret(cells, rules); sm <- ps_summary(rr$regret, rr$oracle)
  w <- function(x, nm) utils::write.table(x, file.path(opt$out, sprintf("ps_%s_%s.tsv", nm, opt$tag)), sep = "\t", quote = FALSE, row.names = FALSE)
  w(cells, "cells"); w(rules, "rules"); w(rr$oracle, "oracle"); w(rr$regret, "regret"); w(sm, "summary")
  w(data.frame(id = grid$id, secs = vapply(res, `[[`, 0, "secs"), peak_gb = vapply(res, `[[`, 0, "peak_gb"), planted_pairs = vapply(res, `[[`, 0, "n_planted_pairs")), "datasets")
  print(sm[sm$method == "tcaCLR", ], digits = 3, row.names = FALSE)
  say(sprintf("done: %.1f s total, peak RSS (this process) %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), ad_mem()))
}
