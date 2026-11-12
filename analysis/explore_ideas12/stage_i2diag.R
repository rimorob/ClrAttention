#!/usr/bin/env Rscript
# EXPLORATORY. Step 2 of docs/plans/explore_ideas12_plan.md: idea-2 detectors (cross-compendium top-20 hypergeometric overlap) and the
# label-free diagnostic (fraction of genes with rho < 0.5). No truth is used. Writes i2_flag.txt: stop | proceed | proceed_flag.
args <- commandArgs(TRUE); threads <- if (length(args)) as.integer(args[1]) else 8L
source("analysis/explore_ideas12/ideas12_lib.R")
out <- file.path(EI_OUT, "k562"); say <- ei_say(file.path(out, "stage_i2diag.log"))
core <- readRDS(file.path(out, "cache", "k562_core.rds")); genes <- core$genes; G <- length(genes); C <- core$C
## compendium Y exactly as analysis/perturbseq_bulk.R (half A of the significant knockdowns)
load_Y <- function(h5) {
  set.seed(20261015L); rd <- function(p) rhdf5::h5read(h5, p)
  cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field)); if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL); if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
  Z <- rd("X"); g <- cat_of("/var", "gene_name"); rownames(Z) <- g; pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
  pval <- as.numeric(rd("/obs/energy_test_p_value")); ok <- pval < 1e-3 & !grepl("non-targeting", pert)
  K <- pert[ok]; K <- K[sample(length(K))]; Ah <- sort(K[seq_len(floor(length(K) / 2))])
  Y <- Z[, Ah]; rm(Z); invisible(gc(verbose = FALSE))
  bad <- rowMeans(!is.finite(Y)) > 0.01; Y <- Y[!bad, , drop = FALSE]; n_imp <- sum(!is.finite(Y)); Y[!is.finite(Y)] <- 0
  Y <- Y[apply(Y, 1, function(x) diff(range(x))) > 0, ]
  list(Y = Y, n_profiles = ncol(Y), n_genes = nrow(Y), n_dropped = sum(bad), n_imputed = n_imp) }
res <- list(); diag_rows <- list()
for (nm in c(gwps = "data/perturbseq/K562_gwps_normalized_bulk_01.h5ad", rpe1 = "data/perturbseq/rpe1_normalized_bulk_01.h5ad")) {
  key <- names(nm)[1]; yl <- load_Y(nm); Y <- yl$Y
  uidx <- which(genes %in% rownames(Y)); U <- genes[uidx]; Yu <- Y[U, , drop = FALSE]; rm(Y); invisible(gc(verbose = FALSE))
  say(sprintf("%s: %d profiles, %d genes (%d dropped >1%% missing, %d values set to 0); universe U = %d of %d essential-compendium genes", key, yl$n_profiles, yl$n_genes, yl$n_dropped, yl$n_imputed, length(U), G))
  f <- ClrAttention$new(Yu)$estimate_mi(bins = "hg", transform = "none", threads = threads); f$calibrate(method = "normal", combine = "stouffer")
  Cy <- f$clr_scores; f$release(); rm(f, Yu); invisible(gc(verbose = FALSE))
  Tx <- ei_topk(C, uidx, 20L)                                         # essential half A, columns restricted to U (essential-gene index space)
  TyU <- ei_topk(Cy, seq_along(uidx), 20L); Ty <- matrix(0L, G, 20L)   # compendium Y in its own U space -> map to essential index space
  mapped <- matrix(0L, length(uidx), 20L); nz <- TyU > 0; mapped[nz] <- uidx[TyU[nz]]; Ty[uidx, ] <- mapped
  r <- ei_rho(Tx, Ty, length(uidx)); rho <- r$rho; inU <- seq_len(G) %in% uidx; rho[!inU] <- 1; r$p[!inU] <- NA; r$overlap[!inU] <- NA_integer_
  res[[key]] <- list(rho = rho, p = r$p, overlap = r$overlap, in_universe = inU, n_universe = length(uidx), n_profiles = yl$n_profiles)
  diag_rows[[key]] <- data.frame(detector = key, n_profiles_Y = yl$n_profiles, n_universe = length(uidx), n_outside = sum(!inU), frac_rho_lt_0.5 = mean(rho[inU] < 0.5),
    frac_rho_lt_0.5_all_genes = mean(rho < 0.5), frac_overlap0 = mean(r$overlap[inU] == 0), frac_p_le_0.05 = mean(r$p[inU] <= 0.05), median_overlap = stats::median(r$overlap[inU]),
    mean_overlap = mean(r$overlap[inU]), chance_overlap = 20 * 20 / (length(uidx) - 1))
  say(sprintf("%s: fraction rho < 0.5 = %.4f (U genes), overlap 0: %.4f, p <= 0.05: %.4f, mean overlap %.2f (chance %.2f)", key, diag_rows[[key]]$frac_rho_lt_0.5, diag_rows[[key]]$frac_overlap0,
              diag_rows[[key]]$frac_p_le_0.05, diag_rows[[key]]$mean_overlap, diag_rows[[key]]$chance_overlap))
  rm(Cy, Tx, Ty, TyU, mapped); invisible(gc(verbose = FALSE)) }
saveRDS(res, file.path(out, "cache", "i2_rho.rds"))
dg <- do.call(rbind, diag_rows); utils::write.csv(dg, file.path(out, "i2_diagnostic.csv"), row.names = FALSE)
fr <- dg$frac_rho_lt_0.5[dg$detector == "gwps"]
flag <- if (fr >= 0.30) "stop" else if (fr <= 0.20) "proceed" else "proceed_flag"
writeLines(flag, file.path(out, "i2_flag.txt")); say(sprintf("gwps fraction rho < 0.5 = %.4f -> %s (rule: <= 0.20 proceed, >= 0.30 stop, between: proceed flagged)", fr, flag)); say("done")
