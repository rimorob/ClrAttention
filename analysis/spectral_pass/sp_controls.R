#!/usr/bin/env Rscript
# Spectral-embedding exploratory pass: controls (docs/plans/spectral_embedding_pass_plan.md section 6). No truth file is read.
#  POSITIVE: planted structure (3 Gaussian clusters, 2 subspace clusters, 1 ring) in (a) a purely synthetic embedding and (b) the real spectral embedding of a unit (the planted genes REPLACE the rows of
#            randomly chosen genes; planted coordinates scaled by the real column SDs). Each variant is scored by pooled AUPR / base rate (ap_fast) on the co-membership of the planted groups (all genes pairs),
#            for all planted groups together and for each structure type; recovery diagnostics (mclust G, HDDC K, pgmm G, rings found and their plane and radius).
#  NEGATIVE: per-gene permuted expression of a real unit (each gene's values permuted across cells, seed 90001), the same pipeline (own edge from 10 further per-gene permutations): number of eigenvalues above
#            the edge, components found (mclust G and whether G = 1 had the better BIC; HDDC K; pgmm G) and rings found, at the edge-rule k (if >= 2) and at the forced k = 5, 10, 20.
# Usage (repo root): Rscript analysis/spectral_pass/sp_controls.R --mode positive|negative|both --dataset HepG2 --unit 1 --out results/spectral_pass/controls [--only ...] [--krules ...] [--threads 4]
a0 <- commandArgs(trailingOnly = TRUE)
so <- list(mode = "both", dataset = "HepG2", unit = "1", out = "results/spectral_pass/controls", only = "cos,cos_clr2,mclust,hddc,pgmm,ring,rmix", krules = "E,5,10,20", threads = 4L, expr = "auto", gres = 12L, seed = 20261121L)
i <- 1L; while (i <= length(a0)) { k <- sub("^--", "", a0[i]); stopifnot(k %in% names(so)); so[[k]] <- if (k %in% c("threads", "gres", "seed")) as.integer(a0[i + 1L]) else a0[i + 1L]; i <- i + 2L }
SP_OUT <- so$out; SP_THREADS <- so$threads; source("analysis/spectral_pass/sp_lib.R")
ONLY <- strsplit(so$only, ",")[[1]]; KR <- strsplit(so$krules, ",")[[1]]; dir.create(so$out, recursive = TRUE, showWarnings = FALSE)
say <- tsb_say(file.path(so$out, sprintf("sp_controls_%s_%s_%s.log", so$mode, so$dataset, so$unit))); say("options: ", paste(names(so), unlist(so), sep = "=", collapse = " "))
tag <- sprintf("%s_%s", so$dataset, so$unit)

## ---- planted structure in a G x d coordinate matrix (rows replaced) -----------------------------------------------------------------------
# groups: cl1..cl3 (40 genes each, Gaussian clusters at distance 4 sd from the origin, sd 0.5), sub1, sub2 (40 genes each, centred at the origin, SD 3 in a random 2-D subspace and 0.3 elsewhere),
# ring (50 genes in coordinates (7, 8): radius 5, radial SD 0.25, uniform angle; SD 0.5 elsewhere). Units: multiples of the background SD `sdv` (per column).
sp_plant <- function(Y, sdv, seed) {
  set.seed(seed); G <- nrow(Y); d <- ncol(Y); stopifnot(d >= 8L); idx <- sample.int(G, 250L); grp <- character(G); grp[] <- "none"
  take <- function(n) { r <- idx[seq_len(n)]; idx <<- idx[-seq_len(n)]; r }
  sdm <- function(r) matrix(rep(sdv, each = length(r)), length(r), d)
  for (c in 1:3) { r <- take(40L); mu <- stats::rnorm(6); mu <- 4 * mu / sqrt(sum(mu^2)); M <- matrix(0, length(r), d); M[, 1:6] <- rep(mu, each = length(r))
    Y[r, ] <- (M + matrix(stats::rnorm(length(r) * d, sd = 0.5), length(r), d)) * sdm(r); grp[r] <- paste0("cl", c) }
  for (c in 1:2) { r <- take(40L); Q <- qr.Q(qr(matrix(stats::rnorm(d * 2), d, 2))); inside <- matrix(stats::rnorm(length(r) * 2, sd = 3), length(r), 2) %*% t(Q)
    Y[r, ] <- (inside + matrix(stats::rnorm(length(r) * d, sd = 0.3), length(r), d)) * sdm(r); grp[r] <- paste0("sub", c) }
  r <- take(50L); th <- stats::runif(length(r), 0, 2 * pi); rad <- 5 + stats::rnorm(length(r), sd = 0.25); M <- matrix(stats::rnorm(length(r) * d, sd = 0.5), length(r), d)
  M[, 7] <- rad * cos(th); M[, 8] <- rad * sin(th); sc8 <- mean(sdv[7:8]); sdm_r <- sdm(r); sdm_r[, 7:8] <- sc8; Y[r, ] <- M * sdm_r; grp[r] <- "ring"
  list(Y = Y, grp = grp)
}
sp_fake_emb <- function(Y) list(V = Y, lam = rep(1, ncol(Y)), edge = NA_real_, t2 = 0L, n_above = ncol(Y), k_edge = ncol(Y), Y = function(k) Y[, seq_len(min(k, ncol(Y))), drop = FALSE])
sp_pair_ap <- function(S, lab) { ut <- upper.tri(S); ap <- ap_fast(S[ut], lab[ut]); c(ratio = ap / mean(lab[ut]), base = mean(lab[ut])) }
sp_lab <- function(grp, which) { g <- ifelse(grp %in% which, grp, NA); L <- outer(g, g, function(a, b) !is.na(a) & !is.na(b) & a == b); L * 1 }
run_positive <- function(Y, sdv, label, seed, null_fun) {
  pl <- sp_plant(Y, sdv, seed); emb <- sp_fake_emb(pl$Y); grp <- pl$grp; G <- nrow(Y)
  say(sprintf("positive control %s: G = %d, d = %d, planted: %s", label, G, ncol(Y), paste(names(table(grp)), table(grp), collapse = ", ")))
  krp <- KR[KR == "E" | suppressWarnings(as.integer(KR)) <= ncol(Y)]
  sc <- sp_unit_scores(emb, say, only = ONLY, krules = krp)
  types <- list(all = c("cl1", "cl2", "cl3", "sub1", "sub2", "ring"), clusters = c("cl1", "cl2", "cl3"), subspace = c("sub1", "sub2"), ring = "ring")
  res <- do.call(rbind, lapply(names(sc$S), function(vn) do.call(rbind, lapply(names(types), function(ty) { x <- sp_pair_ap(ext_clean(sc$S[[vn]]), sp_lab(grp, types[[ty]])); data.frame(control = label, variant = vn, structure = ty, aupr_ratio = x[["ratio"]], base = x[["base"]]) }))))
  rp <- NULL; if (!is.null(sc$ring_planes)) { nullh1 <- null_fun(length(sc$ring_planes)); ringg <- which(grp == "ring")
    rp <- do.call(rbind, lapply(seq_along(sc$ring_planes), function(j) { p <- sc$ring_planes[[j]]; ing <- function(f) if (f$accepted && f$H >= 1) max(colSums(f$gamma[ringg, , drop = FALSE] > 0.5)) else NA
      data.frame(control = label, plane = j, A_accepted = p$ring$accepted, A_H = p$ring$H, A_dbic = p$ring$bic_ring - p$ring$bic_cluster, A_mu = paste(round(p$ring$mu, 3), collapse = ";"), A_s = paste(round(p$ring$s, 3), collapse = ";"), A_ring_genes_in_comp = ing(p$ring),
                 B_accepted = p$rmix$accepted, B_H = p$rmix$H, B_dbic = p$rmix$bic_ring - p$rmix$bic_cluster, B_mu = paste(round(p$rmix$mu, 3), collapse = ";"), B_s = paste(round(p$rmix$s, 3), collapse = ";"), B_ring_genes_in_comp = ing(p$rmix),
                 h1_max_persistence = p$h1, h1_null_plane = nullh1[j], h1_null_max_all_planes = max(nullh1), h1_exceeds_null_max = p$h1 > max(nullh1))})) }
  list(res = res, info = data.frame(control = label, sc$info), rings = rp)
}
if (so$mode == "planes") {     # plane-level tests of the two ring families on planted 2-D structures (no embedding, no data)
  set.seed(so$seed); n <- 1250L
  ring <- function(m, rad, sdr, cx = 0, cy = 0) { th <- stats::runif(m, 0, 2 * pi); r <- rad + stats::rnorm(m, 0, sdr); cbind(cx + r * cos(th), cy + r * sin(th)) }
  bulk <- function(m, sd = 0.4) cbind(stats::rnorm(m, 0, sd), stats::rnorm(m, 0, sd))
  cases <- list(thin_ring = rbind(bulk(1150), ring(100, 3, 0.1)), fat_doughnut = rbind(bulk(1000), ring(250, 3, 1.0)), concentric_rings_with_bulk = rbind(bulk(1050), ring(100, 2, 0.1), ring(100, 4, 0.1)),
    concentric_rings_no_bulk = rbind(ring(600, 2, 0.15), ring(650, 4, 0.15)), ring_off_centre_from_bulk = rbind(bulk(1150), ring(100, 3, 0.1, 5, 2)),
    three_gaussian_blobs = rbind(bulk(400), sweep(bulk(400), 2, c(3, 0), "+"), sweep(bulk(450), 2, c(-2, 3), "+")), one_gaussian_blob = bulk(n, 1), standard_normal_noise = cbind(stats::rnorm(n), stats::rnorm(n)),
    uniform_disc_not_a_ring = { r <- sqrt(stats::runif(n)) * 3; th <- stats::runif(n, 0, 2 * pi); cbind(r * cos(th), r * sin(th)) })
  res <- do.call(rbind, lapply(names(cases), function(nm) { xy <- cases[[nm]]; bc <- sp_plane_cluster_bic(xy[, 1], xy[, 2]); A <- sp_ring_plane(xy[, 1], xy[, 2], bc = bc); B <- sp_rmix_plane(xy[, 1], xy[, 2], bc = bc)
    say(sprintf("%-27s A accepted %d (rings %d, dBIC %.1f) | B accepted %d (rings %d, dBIC %.1f)", nm, A$accepted, A$H, A$bic_ring - bc, B$accepted, B$H, B$bic_ring - bc))
    data.frame(case = nm, mclust_bic = bc, A_accepted = A$accepted, A_rings = A$H, A_dbic = A$bic_ring - bc, A_mu = paste(round(A$mu, 2), collapse = ";"), B_accepted = B$accepted, B_rings = B$H, B_dbic = B$bic_ring - bc,
               B_mu = paste(round(B$mu, 2), collapse = ";"), B_s = paste(round(B$s, 2), collapse = ";"), B_centre = paste(round(B$centre, 2), collapse = ";"), h1_max_persistence = sp_tda_h1(xy[, 1], xy[, 2])) }))
  utils::write.csv(res, file.path(so$out, "planes_ring_tests.csv"), row.names = FALSE)
}
if (so$mode %in% c("positive", "both")) {
  set.seed(so$seed); Gs <- 1250L; ds <- so$gres
  null_syn <- function(np) { Z <- matrix(stats::rnorm(Gs * ds), Gs, ds); vapply(seq_len(np), function(j) sp_tda_h1(Z[, j], Z[, j + 1L]), 0) }
  r1 <- run_positive(matrix(stats::rnorm(Gs * ds), Gs, ds), rep(1, ds), "synthetic", so$seed, null_syn)
  D <- sp_dataset(so$dataset, so$expr); u <- D$unit_X(D$units[match(so$unit, as.character(D$units))]); e <- ext_edge(u$X); emb <- sp_embed(ns_clr(u$X), e$edge, e$t2); Yr <- emb$Y(ds)
  say(sprintf("real embedding %s: k_edge %d, column SD of the first %d weighted coordinates: %s", tag, emb$k_edge, ds, paste(round(apply(Yr, 2, stats::sd), 4), collapse = " ")))
  r2 <- run_positive(Yr, apply(Yr, 2, stats::sd), paste0("real_", tag), so$seed + 1L, function(np) sp_tda_null(u$X, e$t2, e$edge, np))
  for (nm in c("res", "info", "rings")) utils::write.csv(dplyr_bind(Filter(Negate(is.null), list(r1[[nm]], r2[[nm]]))), file.path(so$out, sprintf("positive_%s_%s.csv", nm, tag)), row.names = FALSE)
  say("positive control written")
}
if (so$mode %in% c("negative", "both")) {
  D <- sp_dataset(so$dataset, so$expr); u <- D$unit_X(D$units[match(so$unit, as.character(D$units))]); Xp <- ad_permute_genes(u$X, 90001L); e <- ext_edge(Xp)
  emb <- sp_embed(ns_clr(Xp), e$edge, e$t2); say(sprintf("negative control %s (per-gene permuted): permuted-data edge %.5f, t2 %d, eigenvalues above the edge (first %d) = %d, k_edge = %d", tag, e$edge, e$t2, SP$neig, emb$n_above, emb$k_edge))
  sc <- sp_unit_scores(emb, say, only = ONLY, krules = KR)
  nc <- sc$info; nc$dataset <- so$dataset; nc$unit <- so$unit; nc$n_above <- emb$n_above; nc$k_edge <- emb$k_edge; nc$top_lambda <- emb$lam[1]; nc$edge <- e$edge
  nc$rings_found <- ifelse(nc$readout %in% c("ring", "rmix"), nc$n_comp, NA)
  utils::write.csv(nc, file.path(so$out, sprintf("negative_info_%s.csv", tag)), row.names = FALSE)
  if (!is.null(sc$ring_planes)) utils::write.csv(do.call(rbind, lapply(seq_along(sc$ring_planes), function(j) { p <- sc$ring_planes[[j]]; data.frame(plane = j, A_accepted = p$ring$accepted, A_H = p$ring$H, A_dbic = p$ring$bic_ring - p$ring$bic_cluster, B_accepted = p$rmix$accepted, B_H = p$rmix$H, B_dbic = p$rmix$bic_ring - p$rmix$bic_cluster, h1_max_persistence = p$h1) })), file.path(so$out, sprintf("negative_ring_planes_%s.csv", tag)), row.names = FALSE)
  utils::write.csv(data.frame(lam = emb$lam, edge = e$edge), file.path(so$out, sprintf("negative_eigenvalues_%s.csv", tag)), row.names = FALSE)
  say("negative control written")
}
