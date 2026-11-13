# Scoreboard (docs/plans/scoreboard_plan.md): shared helpers for the four scoreboard jobs. DESCRIPTIVE: point estimates only, no
# intervals, no tests. Sources the existing libraries read-only (rz_*, sk_*, ap_fast, regulon_ap_table, truth builders).
sb_say <- function(logf, ...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

sb_source <- function() {
  suppressPackageStartupMessages({ library(clr); library(Matrix) })
  source("analysis/robust_z/robust_z_lib.R")
  source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
  source("analysis/explore_stack/stack_lib.R")
}

# ---- gate bookkeeping ------------------------------------------------------------------------------------------------------
# rows: list of data.frame(check, max_abs_diff, tol, pass); finish writes the csv and stops on failure
sb_gate_add <- function(rows, check, value, tol = 1e-9) { rows[[length(rows) + 1L]] <- data.frame(check = check, max_abs_diff = value, tol = tol, pass = is.finite(value) && value <= tol, stringsAsFactors = FALSE); rows }
sb_gate_finish <- function(rows, file, say) {
  g <- do.call(rbind, rows); utils::write.csv(g, file, row.names = FALSE); print(g, row.names = FALSE)
  if (!all(g$pass)) stop("EQUIVALENCE GATE FAILED (", file, "): nothing new is reported")
  say("EQUIVALENCE GATE PASSED: ", basename(file)); invisible(g)
}

# ---- E. coli data and truths: analysis/baselines.R, verbatim logic for m3d and ecomac -----------------------------------
sb_load_ecoli <- function(dataset) {
  if (dataset == "m3d") {
    m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
    raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
    X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
  } else {
    X <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
    rownames(X) <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
    plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
    if (length(unique(plat)) > 1L) X <- suppressMessages(sva::ComBat(X, batch = factor(plat)))
  }
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  genes <- rownames(X); G <- length(genes)
  to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
  rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
  operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
    ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck)
    list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes), co = comember_matrix(reg, G), operon = operon_of) })
  tfr <- rtab[rtab$reg_class == "TF" & !is.na(rtab$bnumber), ]
  tf_parts <- strsplit(.classic_tf_to_genes(tfr$regulator), ";", fixed = TRUE)
  tfnet <- data.frame(tf_gene = to_bn(unlist(tf_parts)), target = rep(tfr$bnumber, lengths(tf_parts)), confidence = rep(tfr$confidence, lengths(tf_parts)), stringsAsFactors = FALSE)
  tfnet <- tfnet[!is.na(tfnet$tf_gene), ]
  U <- list(SC = edge_universe(tfnet, genes, conf_keep = c("S", "C")), all = edge_universe(tfnet, genes))
  list(X = X, genes = genes, bench = bench, U = U)
}

# per-TF median AUROC (analysis/baselines.R::per_tf_median_auroc, verbatim; G = number of genes)
sb_per_tf_median_auroc <- function(S, u, G) stats::median(vapply(names(u$targets), function(tn) {
  t <- as.integer(tn); y <- integer(G); y[unique(u$targets[[tn]])] <- 1L; y[t] <- NA
  if (sum(y, na.rm = TRUE) < 3) return(NA_real_)
  s <- S[t, ]; ok <- !is.na(y); r <- rank(s[ok]); n1 <- sum(y[ok]); n0 <- sum(ok) - n1
  (sum(r[y[ok] == 1]) - n1 * (n1 + 1) / 2) / (as.numeric(n1) * n0) }, 0), na.rm = TRUE)

# ---- operator, layers, stratifiers (package code path; stack_lib.R) -------------------------------------------------------
# returns the CLR score C, MI M, operator A, readout matrices S9 (attention, t = 9), S1, S2, S3 (a3r layers, t = 9 each), H_att, H_mi
sb_prep <- function(X, threads, say) {
  G <- nrow(X); t0 <- Sys.time()
  fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = threads)
  fit$calibrate(method = "normal", combine = "stouffer")
  M <- fit$mi; C <- fit$clr_scores
  fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A_pkg <- fit$operator; fit$release()
  A <- rz_operator(C); eq_op <- max(abs(A - A_pkg)); rm(A_pkg)
  say(sprintf("MI + CLR + operator %.0f s; max|A(rz_operator) - package operator| = %.3g", as.numeric(difftime(Sys.time(), t0, units = "secs")), eq_op))
  stopifnot(eq_op < 1e-12)
  P <- rz_lazy_P(A); q_att <- rz_stationary(P); H_att <- outer(q_att, q_att, "+")
  mi_rm <- (rowSums(M) - diag(M)) / (G - 1); H_mi <- outer(mi_rm, mi_rm, "+"); rm(M)
  Pt9 <- as.matrix(P); for (s in 2:9) Pt9 <- as.matrix(Pt9 %*% P)
  S9 <- (Pt9 + t(Pt9)) / 2; diag(S9) <- 0; S1 <- sk_a3r(Pt9); rm(Pt9, P); invisible(gc(verbose = FALSE))
  L1 <- sk_layer(C, ts_out = 9L); d_e <- max(abs(L1$outs[["9"]] - S1)); rm(L1)
  say(sprintf("generic stack N = 1 vs a3r of the package path: max|diff| = %.3g", d_e)); stopifnot(d_e < 1e-12)
  L2 <- sk_layer(S1, ts_out = 9L); S2 <- L2$outs[["9"]]; rm(L2); invisible(gc(verbose = FALSE))
  L3 <- sk_layer(S2, ts_out = 9L); S3 <- L3$outs[["9"]]; rm(L3); invisible(gc(verbose = FALSE))
  list(C = C, S = list(CLR = C, attention = S9, a3r = S1, a3r2 = S2, a3r3 = S3), H_att = H_att, H_mi = H_mi)
}

# pooled AUPR on every truth of `bench` (list with pairs$i, pairs$j, pairs$label), raw / hub_att / hub_mi.
# returns data.frame(readout, scope, truth, aupr, base_rate, enrichment)
sb_pooled <- function(name, S, bench, H_att, H_mi, scopes = c("raw", "hub_att", "hub_mi")) {
  S[!is.finite(S)] <- 0; diag(S) <- 0; out <- list()
  for (sc in scopes) {
    Ss <- switch(sc, raw = S, hub_att = rz_hub_stratify(S, H_att), hub_mi = rz_hub_stratify(S, H_mi))
    for (e in names(bench)) { b <- bench[[e]]; br <- mean(b$pairs$label); a <- ap_fast(Ss[cbind(b$pairs$i, b$pairs$j)], b$pairs$label)
      out[[length(out) + 1L]] <- data.frame(readout = name, scope = sc, truth = e, aupr = a, base_rate = br, enrichment = a / br, stringsAsFactors = FALSE) }
    if (sc != "raw") rm(Ss); invisible(gc(verbose = FALSE))
  }
  do.call(rbind, out)
}
