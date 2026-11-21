#!/usr/bin/env Rscript
# Extended method comparison ("SOTA" table with tcaCLR and tcaCLR2 columns). EXPLORATORY for tcaCLR2 (docs/plans/tcaclr2_definition.md; outside the preregistration).
# Same statistics as analysis/sota_table.R: pooled co-membership AUPR over base rate, paired delete-half gene jackknife (100 replicates; the halves of analysis/baselines.R:
# set.seed(20261016), drawn per truth in the baselines.R order), SE = sd of the replicate differences (delete-d factor 1), a method is bold when the Bonferroni-corrected
# 95% interval of (best - method) includes 0, dagger = tied only after the correction. Intervals resample GENES only (fitted networks held fixed; no experiments, cells or seeds).
#
#   Rscript analysis/sota_extended.R --matrix m3d|ecomac|k562 [--out results/sota_extended] [--threads 8] [--jack 100] [--tag ""]
#       scores the nine columns on the 18-cell rows of one matrix, gates, writes rows_<m>.csv, jack_<m>.csv, gate_<m>.csv, run_<m>.csv
#   Rscript analysis/sota_extended.R --matrix assemble [--matrices m3d,ecomac,k562] [--out results/sota_extended] [--tag ""]
#       marks over eight / six / legacy-seven methods, paired differences, the old-table reproduction gate; writes table_extended.md/.csv, marks_by_row.csv,
#       paired_differences.csv, gate_assemble.csv
#
# Columns: aCLR (= "Attention" of sota_table.R, symmetrised P^9), tcaCLR (t = 9, registered: second CLR on log P^9, sk_a3r), tcaCLR2 (same at t2 = depth of the permuted-data edge rule),
# CLR, GENIE3, PIDC, GRNBoost2, |r|, Hub only (raw rows only; a degree baseline, not one of the eight). Everything is reused: operator rz_operator / rz_lazy_P (checked equal to the package
# operator), registered sequential chain (tsb_walk / the K562 chain of score_scoreboard.R), read-out sk_a3r, stratifier rz_stationary + rz_hub_stratify, truths of sb_load_ecoli /
# ideas47_lib.R, ap_fast. GENIE3 / PIDC / GRNBoost2 matrices are the stored results/baselines_*/{genie3,pidc,grnboost2}_scores.rds. No operator or scorer is written here.
# Gate (a matrix job stops before writing anything on failure): tcaCLR (t = 9) and tcaCLR2 reproduce results/tcaclr2_scoreboard/cells_reference_tcaclr_all.csv and cells_tcaclr2_all.csv
# (raw and hub_att pooled enrichment) to 1e-9; the six old methods (and Hub only) reproduce results/baselines_*/metrics.csv AUPR and jack_replicates.csv replicates to 1e-9.
# The assemble step additionally reproduces docs/whitepapers/sota_table.md (printed values and bold / dagger marks, legacy convention) and stops on any difference.
# Never overwrites an existing output file (a failed gate writes gate_<m>_FAILED.csv only).
suppressPackageStartupMessages(source("analysis/tcaclr2/tsb_lib.R"))
SE <- tsb_args(list(matrix = "", matrices = "m3d,ecomac,k562", out = "results/sota_extended", threads = 8L, jack = 100L, tag = "", md = "docs/whitepapers/sota_table.md"))
MX <- SE$matrix; if (!MX %in% c("m3d", "ecomac", "k562", "assemble")) stop("--matrix must be m3d, ecomac, k562 or assemble")
tsb_refuse_frozen(SE$out); dir.create(SE$out, recursive = TRUE, showWarnings = FALSE)
say <- tsb_say(file.path(SE$out, sprintf("sota_extended_%s%s.log", MX, SE$tag))); say("options: ", paste(names(SE), unlist(SE), sep = "=", collapse = " "))
METH <- c("aCLR", "tcaCLR", "tcaCLR2", "CLR", "GENIE3", "PIDC", "GRNBoost2", "abs_r", "hub_only")
SHOW <- c(aCLR = "aCLR", tcaCLR = "tcaCLR", tcaCLR2 = "tcaCLR2", CLR = "CLR", GENIE3 = "GENIE3", PIDC = "PIDC", GRNBoost2 = "GRNBoost2", abs_r = "\\|r\\|", hub_only = "Hub only")
OLDKEY <- c(aCLR = "attention", CLR = "clr", GENIE3 = "genie3", PIDC = "pidc", GRNBoost2 = "grnboost2", abs_r = "abs_pearson", hub_only = "hub")
SIX <- c("aCLR", "CLR", "GENIE3", "PIDC", "GRNBoost2", "abs_r"); EIGHT <- c("aCLR", "tcaCLR", "tcaCLR2", SIX[-1])
RUN <- list(m3d = "results/baselines_m3d", ecomac = "results/baselines_ecomac", k562 = "results/baselines_k562ess_truthgenes_genie3")
EDGE <- list(m3d = "results/tcaclr2/edge_m3d_tf_regulons.tsv", ecomac = file.path(TSB$out, "edge_ecomac.tsv"), k562 = "results/tcaclr2/edge_k562_corum.tsv")
TRU <- list(m3d = c("SC", "all"), ecomac = c("SC", "all"), k562 = c("DoRothEA_AB", "DoRothEA_ABC", "CollecTRI", "CORUM", "KD_response"))
out_file <- function(stem, m) file.path(SE$out, sprintf("%s_%s%s.csv", stem, m, SE$tag))
write_new <- function(x, f) { if (file.exists(f)) stop("refusing to overwrite an existing result file: ", f); tmp <- paste0(f, ".tmp", Sys.getpid()); utils::write.csv(x, tmp, row.names = FALSE); invisible(file.rename(tmp, f)) }
write_new_lines <- function(x, f) { if (file.exists(f)) stop("refusing to overwrite an existing result file: ", f); tmp <- paste0(f, ".tmp", Sys.getpid()); writeLines(x, tmp); invisible(file.rename(tmp, f)) }

# ============================================================================ one matrix ==================================================================================
if (MX != "assemble") {
  T0 <- Sys.time(); ds <- MX; NJ <- SE$jack; stopifnot(NJ >= 2L, NJ <= 100L)
  for (f in c(out_file("rows", ds), out_file("jack", ds), out_file("gate", ds))) if (file.exists(f)) stop("output exists (use another --tag): ", f)
  et <- tsb_edge_read(EDGE[[ds]]); t2 <- et$t2; say(sprintf("t2 = %d from %s (edge %.5f)", t2, EDGE[[ds]], et$edge))
  nc <- max(1L, min(8L, SE$threads))
  if (ds %in% c("m3d", "ecomac")) {
    source("analysis/scoreboard/sb_lib.R"); sb_source()
    D <- sb_load_ecoli(ds); X <- D$X; genes <- D$genes; G <- length(genes); bench <- D$bench[TRU[[ds]]]; rm(D)
    fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = SE$threads); fit$calibrate(method = "normal", combine = "stouffer")
    C <- fit$clr_scores; fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A_pkg <- fit$operator; fit$release(); rm(fit)
    A <- rz_operator(C); eq_op <- max(abs(A - A_pkg)); rm(A_pkg); say(sprintf("max|A(rz_operator) - package operator| = %.3g", eq_op)); stopifnot(eq_op < 1e-12)
    P <- rz_lazy_P(A); rm(A); q <- rz_stationary(P); H_att <- outer(q, q, "+")
    W <- tsb_walk(P, c(9L, t2), function(t, Pt) if (t == 9L) list(att = (Pt + t(Pt)) / 2, a3r = sk_a3r(Pt)) else list(a3r = sk_a3r(Pt)))
    rm(P); invisible(gc(verbose = FALSE))
    SC <- list(aCLR = W[["9"]]$att, tcaCLR = W[["9"]]$a3r, tcaCLR2 = W[[as.character(t2)]]$a3r); rm(W)
  } else {
    opt <- list(out = file.path(SE$out, "work_k562"), threads = SE$threads, workers = nc, jack = 1L, quick = 0L)               # ideas47_lib.R set-up (as score_scoreboard.R --family k562)
    source("analysis/explore_ideas47/ideas47_lib.R"); source("analysis/explore_stack/stack_lib.R"); source("analysis/scoreboard/sb_lib.R"); unlink(file.path(opt$out, "cache"), recursive = TRUE)
    say <- tsb_say(file.path(SE$out, sprintf("sota_extended_%s%s.log", MX, SE$tag))); stopifnot(is.list(RUN), is.list(EDGE), is.list(TRU), is.list(SE), identical(MX, "k562"), is.character(METH))
    stopifnot(G == 5671L); abc <- dorothea_sets(genes, c("A", "B", "C"))
    bench <- list(DoRothEA_AB = bench$DoRothEA_AB, DoRothEA_ABC = list(reg = abc, pairs = comembership_pairs(abc, none, genes)), CollecTRI = bench$CollecTRI, CORUM = bench$CORUM_all, KD_response = bench$KD_response)
    Pt <- Pt9; rm(Pt9, halves); invisible(gc(verbose = FALSE)); S9a <- S9; S9b <- a3r(Pt)                                    # t = 9 chain of ideas47_lib.R, continued sequentially to t2
    if (t2 < 9L) stop("t2 < 9 not supported for K562 (chain starts at t = 9)")
    if (t2 > 9L) for (s in 10:t2) { Pt <- as.matrix(Pt %*% P); if (s %% 8L == 0L) invisible(gc(verbose = FALSE)) }
    SC <- list(aCLR = S9a, tcaCLR = S9b, tcaCLR2 = if (t2 > 9L) a3r(Pt) else S9b); rm(Pt, P, S9a, S9b, S9, M); invisible(gc(verbose = FALSE))
  }
  say(sprintf("%s: %d genes, t2 = %d; chain done (%.0f s, peak RSS %.1f GB); truths: %s", ds, G, t2, as.numeric(difftime(Sys.time(), T0, units = "secs")), tsb_hwm(),
              paste(sprintf("%s %d pairs (%d pos)", names(bench), vapply(bench, function(b) length(b$pairs$label), 0L), vapply(bench, function(b) sum(b$pairs$label), 0)), collapse = "; ")))
  stopifnot(identical(names(bench), TRU[[ds]]))
  # the other columns: CLR (fit$clr_scores), |r|, stored baselines matrices
  SC$CLR <- C; rm(C); r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; SC$abs_r <- r; rm(r, X)
  for (m in c("GENIE3", "PIDC", "GRNBoost2")) { S <- readRDS(file.path(RUN[[ds]], paste0(tolower(m), "_scores.rds"))); stopifnot(identical(rownames(S), genes), identical(colnames(S), genes)); SC[[m]] <- S }
  SC$hub_only <- H_att; rm(S); invisible(gc(verbose = FALSE))

  # jackknife halves: baselines.R, set.seed(seed) then per truth in bench order, `jack` draws of floor(n/2) genes of genes_in (the first NJ are used)
  set.seed(20261016L)
  KEEP <- lapply(bench, function(b) lapply(seq_len(100L), function(r) sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2))))
  IX <- list()
  ix_of <- function(tr) { if (is.null(IX[[tr]])) { b <- bench[[tr]]; IX[[tr]] <<- lapply(KEEP[[tr]][seq_len(NJ)], function(keep) { inh <- logical(G); inh[keep] <- TRUE; which(inh[b$pairs$i] & inh[b$pairs$j]) }) }; IX[[tr]] }
  rows <- list(); jk <- list(); tt <- c()
  for (m in METH) {
    t1 <- Sys.time(); S <- SC[[m]]; SC[[m]] <- NULL; S[!is.finite(S)] <- 0; diag(S) <- 0
    for (sc in if (m == "hub_only") "raw" else c("raw", "hub_att")) {
      Ss <- if (sc == "raw") S else rz_hub_stratify(S, H_att)
      for (tr in names(bench)) { b <- bench[[tr]]; scv <- Ss[cbind(b$pairs$i, b$pairs$j)]; lab <- b$pairs$label; br <- mean(lab); a <- ap_fast(scv, lab)
        jv <- unlist(parallel::mclapply(ix_of(tr), function(ix) ap_fast(scv[ix], lab[ix]), mc.cores = nc)); stopifnot(length(jv) == NJ)
        rows[[length(rows) + 1L]] <- data.frame(dataset = ds, scope = sc, truth = tr, method = m, aupr = a, base_rate = br, enrichment = a / br, t2 = t2, stringsAsFactors = FALSE)
        jk[[length(jk) + 1L]] <- data.frame(dataset = ds, scope = sc, truth = tr, method = m, rep = seq_len(NJ), aupr = jv, stringsAsFactors = FALSE) }
      if (sc != "raw") rm(Ss); invisible(gc(verbose = FALSE)) }
    rm(S); invisible(gc(verbose = FALSE)); tt[m] <- as.numeric(difftime(Sys.time(), t1, units = "secs")); say(sprintf("scored %-10s %.0f s (peak RSS %.1f GB)", m, tt[m], tsb_hwm()))
  }
  rows <- do.call(rbind, rows); jk <- do.call(rbind, jk)

  # ---- gate -----------------------------------------------------------------------------------------------------------------------------------------------------------
  g <- list(); gadd <- function(check, d, tol = 1e-9) g[[length(g) + 1L]] <<- data.frame(check = check, n = length(d), max_abs_diff = max(abs(d)), tol = tol, pass = all(is.finite(d)) && max(abs(d)) <= tol, stringsAsFactors = FALSE)
  for (spec in list(c("tcaCLR", "cells_reference_tcaclr_all.csv", "tcaCLR_t9"), c("tcaCLR2", "cells_tcaclr2_all.csv", "tcaCLR2"))) {
    cf <- utils::read.csv(file.path(TSB$out, spec[2]), stringsAsFactors = FALSE); cf <- cf[cf$dataset == ds & cf$metric == "pooled_aupr_over_base" & cf$method == spec[3], ]
    for (sc in c("raw", "hub_att")) { x <- rows[rows$method == spec[1] & rows$scope == sc, ]; y <- cf[cf$scope == sc, ]; v <- y$enrichment_over_random[match(x$truth, y$truth)]; stopifnot(!anyNA(v))
      gadd(sprintf("%s %s %s pooled AUPR/base vs %s (%d cells)", ds, spec[1], sc, spec[2], nrow(x)), x$enrichment - v) } }
  mt <- utils::read.csv(file.path(RUN[[ds]], "metrics.csv"), stringsAsFactors = FALSE); jr <- utils::read.csv(file.path(RUN[[ds]], "jack_replicates.csv"), stringsAsFactors = FALSE)
  for (m in names(OLDKEY)) for (sc in if (m == "hub_only") "raw" else c("raw", "hub_att")) {
    key <- if (sc == "raw") OLDKEY[[m]] else paste0(OLDKEY[[m]], "_hubstrat"); x <- rows[rows$method == m & rows$scope == sc, ]
    gadd(sprintf("%s %s %s AUPR vs metrics.csv (%s)", ds, m, sc, key), x$aupr - mt$aupr[match(paste(key, x$truth), paste(mt$method, mt$truth))])
    gadd(sprintf("%s %s %s base rate vs metrics.csv", ds, m, sc), x$base_rate - mt$base_rate[match(paste(key, x$truth), paste(mt$method, mt$truth))])
    j0 <- jk[jk$method == m & jk$scope == sc, ]; v <- jr$aupr[match(paste(key, j0$truth, j0$rep), paste(jr$method, jr$truth, jr$rep))]
    gadd(sprintf("%s %s %s replicates 1..%d vs jack_replicates.csv (%s)", ds, m, sc, NJ, key), j0$aupr - v) }
  g <- do.call(rbind, g); print(g, row.names = FALSE, digits = 4)
  if (!all(g$pass)) { ff <- out_file("gate", paste0(ds, "_FAILED")); write_new(g, ff); say("GATE FAILED; see ", ff); stop("GATE FAILED: nothing else written") }
  say(sprintf("GATE PASSED: %d checks, max abs diff %.3g", nrow(g), max(g$max_abs_diff)))
  write_new(g, out_file("gate", ds)); write_new(rows, out_file("rows", ds)); write_new(jk, out_file("jack", ds))
  write_new(data.frame(dataset = ds, G = G, t2 = t2, edge = et$edge, jack = NJ, n_rows = nrow(rows), wall_s = as.numeric(difftime(Sys.time(), T0, units = "secs")), peak_rss_gb = tsb_hwm(), threads = SE$threads,
                       t(round(tt, 1))), out_file("run", ds))
  say(sprintf("done: %.0f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), tsb_hwm()))
  quit(save = "no")
}

# ============================================================================ assemble ====================================================================================
mats <- strsplit(SE$matrices, ",")[[1]]; stopifnot(all(mats %in% c("m3d", "ecomac", "k562")))
for (f in c("table_extended.md", "table_extended.csv", "marks_by_row.csv", "paired_differences.csv", "gate_assemble.csv")) { ff <- file.path(SE$out, sub("(\\.[a-z]+)$", paste0(SE$tag, "\\1"), f)); if (file.exists(ff)) stop("output exists (use another --tag): ", ff) }
RW <- do.call(rbind, lapply(mats, function(m) utils::read.csv(out_file("rows", m), stringsAsFactors = FALSE)))
JK <- do.call(rbind, lapply(mats, function(m) utils::read.csv(out_file("jack", m), stringsAsFactors = FALSE)))
DSN <- c(m3d = "M3D", ecomac = "EcoMAC", k562 = "K562")
cell_rows <- unique(RW[, c("dataset", "truth")]); cell_rows <- cell_rows[order(match(cell_rows$dataset, mats)), ]
jrep <- function(d, sc, tr, m) { x <- JK[JK$dataset == d & JK$scope == sc & JK$truth == tr & JK$method == m, ]; stats::setNames(x$aupr[order(x$rep)], sort(x$rep)) }
aup <- function(d, sc, tr, m) RW$aupr[RW$dataset == d & RW$scope == sc & RW$truth == tr & RW$method == m]
# marks exactly as sota_table.R: best = max AUPR in the set, Bonferroni over the (set size - 1) comparisons, SE = sd of the paired replicate differences
marks_for <- function(d, sc, tr, ms) {
  a <- stats::setNames(vapply(ms, function(m) aup(d, sc, tr, m), 0), ms); best <- names(which.max(a)); zc <- stats::qnorm(1 - 0.025 / max(1L, length(ms) - 1L))
  z <- vapply(ms, function(m) if (m == best) NA_real_ else (a[[best]] - a[[m]]) / stats::sd(jrep(d, sc, tr, best) - jrep(d, sc, tr, m)), 0)
  vapply(ms, function(m) if (m == best || z[[m]] < 1.96) "bold" else if (z[[m]] < zc) "bold_dagger" else "none", "")
}
mk <- list(); tabrows <- list()
for (sc in c("raw", "hub_att")) for (k in seq_len(nrow(cell_rows))) { d <- cell_rows$dataset[k]; tr <- cell_rows$truth[k]
  leg <- if (sc == "raw") c(SIX, "hub_only") else SIX
  sets <- list(eight = EIGHT, six = SIX, legacy = leg); if (sc == "raw") sets$eight_plus_hub <- c(EIGHT, "hub_only")
  for (vn in names(sets)) { mm <- marks_for(d, sc, tr, sets[[vn]]); mk[[length(mk) + 1L]] <- data.frame(dataset = d, scope = sc, truth = tr, variant = vn, method = names(mm), mark = unname(mm), stringsAsFactors = FALSE) } }
mk <- do.call(rbind, mk)
getmark <- function(d, sc, tr, vn, m) { x <- mk$mark[mk$dataset == d & mk$scope == sc & mk$truth == tr & mk$variant == vn & mk$method == m]; if (length(x)) x else NA_character_ }
cell_txt <- function(d, sc, tr, m, vn) { if (m == "hub_only" && sc != "raw") return(NULL); v <- sprintf("%.2f", RW$enrichment[RW$dataset == d & RW$scope == sc & RW$truth == tr & RW$method == m]); mm <- getmark(d, sc, tr, vn, m)
  if (is.na(mm) || mm == "none") v else if (mm == "bold") paste0("**", v, "**") else paste0("**", v, "**†") }
wide <- do.call(rbind, lapply(c("raw", "hub_att"), function(sc) do.call(rbind, lapply(seq_len(nrow(cell_rows)), function(k) { d <- cell_rows$dataset[k]; tr <- cell_rows$truth[k]
  x <- data.frame(scope = sc, dataset = DSN[[d]], truth = gsub("_", " ", tr), base_rate = RW$base_rate[RW$dataset == d & RW$truth == tr][1], t2 = RW$t2[RW$dataset == d][1], stringsAsFactors = FALSE)
  for (m in METH) { rr <- RW[RW$dataset == d & RW$scope == sc & RW$truth == tr & RW$method == m, ]; x[[paste0("enrich_", m)]] <- if (nrow(rr)) rr$enrichment else NA_real_
    for (vn in c("eight", "six", "legacy")) x[[paste0("mark_", vn, "_", m)]] <- if (nrow(rr)) getmark(d, sc, tr, vn, m) else NA_character_ }
  x }))))

# ---- paired differences with 95% jackknife intervals (normal: diff +/- 1.96 sd of the replicate differences; same SE as the marks) ------------------------------------------------------
CMP <- list(c("tcaCLR", "CLR"), c("tcaCLR", "aCLR"), c("tcaCLR2", "CLR"), c("tcaCLR2", "aCLR"), c("tcaCLR2", "tcaCLR"))
pd <- list()
for (sc in c("raw", "hub_att")) for (k in seq_len(nrow(cell_rows))) { d <- cell_rows$dataset[k]; tr <- cell_rows$truth[k]; br <- RW$base_rate[RW$dataset == d & RW$truth == tr][1]
  for (cp in CMP) { ja <- jrep(d, sc, tr, cp[1]); jb <- jrep(d, sc, tr, cp[2]); stopifnot(identical(names(ja), names(jb))); a <- aup(d, sc, tr, cp[1]); b <- aup(d, sc, tr, cp[2])
    se <- stats::sd(ja - jb); df <- a - b
    pd[[length(pd) + 1L]] <- data.frame(scope = sc, dataset = DSN[[d]], truth = tr, comparison = paste(cp, collapse = " minus "), aupr_a = a, aupr_b = b, base_rate = br, diff_aupr = df, se_aupr = se, lo_aupr = df - 1.96 * se, hi_aupr = df + 1.96 * se,
      diff_enrich = df / br, lo_enrich = (df - 1.96 * se) / br, hi_enrich = (df + 1.96 * se) / br, z = df / se, interval_excludes_0 = (df - 1.96 * se > 0) | (df + 1.96 * se < 0), n_rep = length(ja), stringsAsFactors = FALSE) } }
pd <- do.call(rbind, pd)

# ---- reproduction of the old table (docs/whitepapers/sota_table.md): printed values and marks, legacy convention ---------------------------------------------------------------------------
md <- readLines(SE$md); hdr <- grep("^\\| Dataset", md); stopifnot(length(hdr) == 2L); g <- list()
md <- gsub("\\|", "@@", md, fixed = TRUE)                                                                        # the escaped pipes of the |r| header
parse_tab <- function(h) { cols <- trimws(strsplit(md[h], "|", fixed = TRUE)[[1]]); cols <- cols[nzchar(cols)]; out <- list(); i <- h + 2L
  while (i <= length(md) && grepl("^\\|", md[i])) { p <- trimws(strsplit(md[i], "|", fixed = TRUE)[[1]]); p <- p[seq_along(cols) + 1L]; names(p) <- cols; out[[length(out) + 1L]] <- p; i <- i + 1L }; out }
mdmap <- c(Attention = "aCLR", CLR = "CLR", GENIE3 = "GENIE3", PIDC = "PIDC", GRNBoost2 = "GRNBoost2", "@@r@@" = "abs_r", "Hub only" = "hub_only"); nchk <- 0L; bad <- character(0)
for (ti in 1:2) { sc <- c("raw", "hub_att")[ti]
  for (p in parse_tab(hdr[ti])) { dd <- names(DSN)[match(p[["Dataset"]], DSN)]; tr <- gsub(" ", "_", p[["Truth"]]); if (is.na(dd) || !dd %in% mats) next
    if (!any(RW$dataset == dd & RW$truth == tr)) { bad <- c(bad, paste("row not found:", dd, tr)); next }
    for (cn in names(mdmap)) { if (!cn %in% names(p)) next; m <- mdmap[[cn]]; txt <- p[[cn]]; if (txt == "-") next; nchk <- nchk + 1L
      mine <- cell_txt(dd, sc, tr, m, "legacy"); if (!identical(mine, txt)) bad <- c(bad, sprintf("%s %s %s %s: old table '%s', recomputed '%s'", sc, dd, tr, m, txt, mine)) }
    if (abs(as.numeric(p[["Base rate"]]) - round(RW$base_rate[RW$dataset == dd & RW$truth == tr][1], 3)) > 1e-12) bad <- c(bad, paste("base rate differs:", dd, tr)) } }
g <- data.frame(check = sprintf("old six methods + Hub only: printed value (2 decimals) and bold/dagger mark vs %s, %d cells", basename(SE$md), nchk), n_mismatch = length(bad), pass = length(bad) == 0L && nchk > 0L)
print(g, row.names = FALSE); if (length(bad)) print(bad)
if (!g$pass) { write_new(data.frame(mismatch = bad), file.path(SE$out, sprintf("gate_assemble_FAILED%s.csv", SE$tag))); stop("OLD-TABLE REPRODUCTION GATE FAILED (", length(bad), " cells): see gate_assemble_FAILED*.csv") }
write_new(g, file.path(SE$out, sprintf("gate_assemble%s.csv", SE$tag)))
write_new(wide, file.path(SE$out, sprintf("table_extended%s.csv", SE$tag))); write_new(mk, file.path(SE$out, sprintf("marks_by_row%s.csv", SE$tag)))
write_new(pd, file.path(SE$out, sprintf("paired_differences%s.csv", SE$tag)))

# ---- markdown table ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tbl <- function(sc) { ms <- if (sc == "raw") METH else METH[METH != "hub_only"]
  out <- c(paste0("| Dataset | Truth | Base rate | ", paste(SHOW[ms], collapse = " | "), " |"), paste0("|---|---|---|", paste(rep("---", length(ms)), collapse = "|"), "|"))
  for (k in seq_len(nrow(cell_rows))) { d <- cell_rows$dataset[k]; tr <- cell_rows$truth[k]
    out <- c(out, sprintf("| %s | %s | %.3f | %s |", DSN[[d]], gsub("_", " ", tr), RW$base_rate[RW$dataset == d & RW$truth == tr][1],
      paste(vapply(ms, function(m) { x <- cell_txt(d, sc, tr, m, if (m == "hub_only") "none" else "eight"); if (is.null(x)) "-" else x }, ""), collapse = " | "))) }
  out }
chg <- function(sc) { out <- c("| Scope | Dataset | Truth | Bold among the original six | Bold among all eight | Original-six marks lost |", "|---|---|---|---|---|---|")
  for (k in seq_len(nrow(cell_rows))) { d <- cell_rows$dataset[k]; tr <- cell_rows$truth[k]
    b6 <- SIX[vapply(SIX, function(m) getmark(d, sc, tr, "six", m) != "none", NA)]; b8 <- EIGHT[vapply(EIGHT, function(m) getmark(d, sc, tr, "eight", m) != "none", NA)]
    out <- c(out, sprintf("| %s | %s | %s | %s | %s | %s |", if (sc == "raw") "raw" else "hub-controlled", DSN[[d]], gsub("_", " ", tr), paste(b6, collapse = ", "), paste(b8, collapse = ", "),
      if (length(setdiff(b6, b8))) paste(setdiff(b6, b8), collapse = ", ") else "none")) }
  out }
t2s <- unique(RW[, c("dataset", "t2")]); t2txt <- paste(sprintf("%s t2 = %d", DSN[t2s$dataset], t2s$t2), collapse = "; ")
lines <- c("# Extended method comparison (tcaCLR and tcaCLR2 columns)", "",
  paste0("*Generated by `analysis/sota_extended.R`.* Pooled co-membership AUPR divided by the base rate (1.00 = random). aCLR = CLR attention at t = 9 (the \"Attention\" column of the earlier table); tcaCLR = registered second-CLR read-out at t = 9; tcaCLR2 = the same at the depth set by the permuted-data edge rule (", t2txt, "). **Bold**: best of the eight methods in the row or tied with it (Bonferroni-corrected 95% interval of the paired gene-jackknife difference over the 7 comparisons includes 0); **x**† tied only after the correction. \"Hub only\" is a degree baseline and does not compete for the mark."), "",
  "**tcaCLR2 is exploratory** (rule fixed in `docs/plans/tcaclr2_definition.md`, outside the preregistration); the depth-24 and depth-32 grid curves had been seen on the spent truths. **The intervals resample genes only** (delete-half jackknife, 100 replicates, fitted networks held fixed); they do not resample experiments or cells or vary the seeds of GENIE3 and GRNBoost2.", "",
  "## Raw scores", "", tbl("raw"), "", "## Hub-controlled (within-stratum percentile in 10 strata of pi_i + pi_j)", "", tbl("hub_att"), "",
  "## Bold status with and without the new columns", "", chg("raw"), chg("hub_att")[-(1:2)], "",
  sprintf("Paired differences (tcaCLR minus CLR, tcaCLR minus aCLR, tcaCLR2 minus CLR, tcaCLR2 minus aCLR, tcaCLR2 minus tcaCLR) with 95%% gene-jackknife intervals: `paired_differences.csv` (%d rows). Old-table reproduction gate: %d cells matched.", nrow(pd), nchk))
write_new_lines(lines, file.path(SE$out, sprintf("table_extended%s.md", SE$tag))); cat(lines, sep = "\n")
say("assembled: ", nrow(cell_rows), " rows x ", length(METH), " columns; paired differences ", nrow(pd))
