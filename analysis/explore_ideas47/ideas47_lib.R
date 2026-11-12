# EXPLORATORY, IN-SAMPLE, NOT PREREGISTERED (docs/plans/explore_ideas47_plan.md). Shared set-up for
# ideas 4 and 7: data, K562 truths, jackknife streams, K10/CLR/a3r, stratifiers, scoring, equivalence gate.
# Everything below follows analysis/robust_z/robust_z.R (K562 branch) and
# analysis/corum_expert2/corum_width_readout.R; the caller defines `opt` (out, quick, threads, workers, jack).
suppressPackageStartupMessages({ library(clr); library(Matrix) })
e47_dir <- (function() { f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis/explore_ideas47" })()
source("analysis/robust_z/robust_z_lib.R")
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
dir.create(file.path(opt$out, "cache"), recursive = TRUE, showWarnings = FALSE)
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n")
  cat(m, "\n", file = file.path(opt$out, "run.log"), append = TRUE) }
none <- stats::setNames(character(0), character(0))
QUICK <- opt$quick > 0L

## ---- data ------------------------------------------------------------------------------------
raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (QUICK) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]; rm(v) }
genes <- rownames(X); G <- length(genes); lab <- colnames(X)
UT <- which(upper.tri(matrix(FALSE, G, G)))
say(sprintf("k562: %d genes x %d profiles (quick = %d)", G, ncol(X), opt$quick))

lib_truths <- function() {
  h5 <- "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"; set.seed(20261015L)
  rd <- function(p) rhdf5::h5read(h5, p)
  cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
    if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
    if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
  pert <- cat_of("/obs", "gene_transcript"); pval <- as.numeric(rd("/obs/energy_test_p_value"))
  K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
  A_half <- sort(K[seq_len(floor(length(K) / 2))]); B_half <- setdiff(K, A_half); rm(rd, cat_of, pert, pval, K)
  tr <- list(DoRothEA_AB = dorothea_sets(genes), CORUM = corum_sets(genes), KD_response = response_sets(genes, path = h5, use = B_half))
  cr <- tr$CORUM; cx <- utils::read.csv("analysis/corum_exploratory/corum_complexes.csv", stringsAsFactors = FALSE)
  sz <- lengths(cr$members); nm <- names(cr$members); rp <- nm %in% cx$complex[cx$group == "RP-dominated"]
  if (!QUICK) stopifnot(sum(rp) == 8L, all(sz == cx$size[match(nm, cx$complex)]))
  sub_of <- function(keep) list(members = cr$members[keep], class = cr$class[keep])
  mk <- function(reg) list(reg = reg, pairs = comembership_pairs(reg, none, genes))
  bench <- list(DoRothEA_AB = mk(tr$DoRothEA_AB), CORUM_all = mk(cr), KD_response = mk(tr$KD_response),
                CollecTRI = mk(collectri_sets(genes)),
                CORUM_noRP = mk(sub_of(!rp)), CORUM_RP = mk(sub_of(rp)),
                CORUM_small = mk(sub_of(sz >= 3 & sz <= 6)), CORUM_large = mk(sub_of(sz >= 20)))
  # merged CORUM as analysis/corum_expert2/corum_width_readout.R
  RP_RE <- "^RP[LS][0-9]+[A-Z]?[0-9]*$|^RPLP[0-2]$|^RPSA$"
  merge_sets <- function(s) {                                     # analysis/clr_logpt_test.R, verbatim
    if (length(s) == 0) return(s)
    s <- s[order(names(s))]
    key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]
    drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j)
      length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
    s <- s[!drop]
    if (length(s) < 2) return(s)
    all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
    inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
    cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
    out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
    names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), "")
    out }
  as_reg <- function(sets, cls) list(members = sets, class = stats::setNames(rep(cls, length(sets)), names(sets)))
  raw_m <- cr$members; raw_m <- raw_m[lengths(raw_m) >= 3L & lengths(raw_m) <= 200L]
  mg <- merge_sets(lapply(raw_m, function(ix) genes[ix])); mg <- lapply(mg, function(g) sort(match(g, genes)))
  rpfrac_mg <- vapply(mg, function(ix) mean(grepl(RP_RE, genes[ix])), 0)
  bench$CORUM_merged <- mk(as_reg(mg, "complex")); bench$CORUM_merged_noRP <- mk(as_reg(mg[rpfrac_mg < 0.3], "complex"))
  list(bench = bench, raw_m = raw_m, RP_RE = RP_RE)
}
TT <- lib_truths(); bench <- TT$bench; raw_m <- TT$raw_m; RP_RE <- TT$RP_RE; rm(TT)
for (e in names(bench)) say(sprintf("truth %-18s %5d sets, %9d pairs, base rate %.4f", e, length(bench[[e]]$reg$members),
  length(bench[[e]]$pairs$label), mean(bench[[e]]$pairs$label)))
base_rate <- vapply(bench, function(b) mean(b$pairs$label), 0)
TRUTHS <- names(bench)

## ---- jackknife streams: robust_z.R k562 order, then merged truths appended --------------------------
set.seed(20261016L)
gin_of <- function(reg) sort(unique(unlist(reg$members)))
draw_half <- function(gin) lapply(seq_len(opt$jack), function(r) sort(sample(gin, floor(length(gin) / 2))))
halves <- list()
halves$DoRothEA_AB <- draw_half(bench$DoRothEA_AB$pairs$genes_in)
invisible(draw_half(gin_of(dorothea_sets(genes, c("A", "B", "C")))))
halves$CollecTRI <- draw_half(gin_of(collectri_sets(genes)))           # the draw robust_z makes and discards
halves$CORUM_all <- draw_half(bench$CORUM_all$pairs$genes_in); halves$KD_response <- draw_half(bench$KD_response$pairs$genes_in)
for (r in seq_len(100L)) invisible(sample(G, G - round(0.2 * G)))     # robust_z's jack_reg delete-20% stream (advanced, unused)
for (e in c("CORUM_noRP", "CORUM_RP", "CORUM_small", "CORUM_large", "CORUM_merged", "CORUM_merged_noRP")) halves[[e]] <- draw_half(bench[[e]]$pairs$genes_in)
nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))

## ---- CLR, operator, walk (package code path), stratifiers --------------------------------------------
t0 <- Sys.time()
fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
fit$calibrate(method = "normal", combine = "stouffer")
M <- fit$mi; C <- fit$clr_scores
fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A_pkg <- fit$operator
fit$release(); rm(fit)
A <- rz_operator(C); P <- rz_lazy_P(A)
eq_op <- max(abs(A - A_pkg)); rm(A_pkg)
say(sprintf("MI + CLR + operator %.0f s; max|A(rz_operator) - package operator| = %.3g", as.numeric(difftime(Sys.time(), t0, units = "secs")), eq_op))
stopifnot(eq_op < 1e-12)
walk_to <- function(P, t) { Pt <- as.matrix(P); if (t > 1L) for (s in 2:t) Pt <- as.matrix(Pt %*% P); Pt }
sym0 <- function(Pt) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S }
Pt9 <- walk_to(P, 9L)
S9 <- sym0(Pt9)
a3r <- function(Pt) { L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA       # docs/plans/corum_analysis.md reference code
  mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  Z <- (L - mu) / sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1)); Z[is.na(Z) | Z < 0] <- 0
  S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; S }
minr <- function(Pt) { S <- pmin(Pt, t(Pt)); diag(S) <- 0; S }
q_att <- rz_stationary(P); H_att <- outer(q_att, q_att, "+")
mi_rm <- (rowSums(M) - diag(M)) / (G - 1); H_mi <- outer(mi_rm, mi_rm, "+")
# degree for the degree-scaled readout (corum_expert2): 1 + number of gene i's pairs in the global top 1% of CLR (upper triangle)
obs <- C[UT]; q99 <- unname(stats::quantile(obs, 0.99)); E <- C >= q99 & C > 0; diag(E) <- FALSE; deg <- 1 + rowSums(E); rm(E, obs)
degr <- function(Pt) { D <- deg * Pt; S <- (D + t(D)) / 2; diag(S) <- 0; S }
co <- matrix(FALSE, G, G); for (mm in raw_m) co[mm, mm] <- TRUE; diag(co) <- FALSE
big <- which(vapply(raw_m, function(mm) length(mm) * (length(mm) - 1) / 2 >= 10, NA))
rp_raw <- vapply(raw_m, function(ix) mean(grepl(RP_RE, genes[ix])), 0)
invisible(gc(verbose = FALSE))

## ---- scoring -----------------------------------------------------------------------------------------
jack_one <- function(sc, b, keeps) {
  i <- b$pairs$i; j <- b$pairs$j; lab_ <- b$pairs$label
  unlist(parallel::mclapply(keeps, function(keep) { inh <- logical(G); inh[keep] <- TRUE; ix <- which(inh[i] & inh[j]); ap_fast(sc[ix], lab_[ix]) }, mc.cores = nc)) }
pooled <- function(S) lapply(stats::setNames(TRUTHS, TRUTHS), function(e) { b <- bench[[e]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]
  list(aupr = ap_fast(sc, b$pairs$label), jack = jack_one(sc, b, halves[[e]])) })
macro_of <- function(S) { x <- regulon_ap_table(S, list(members = raw_m[big]), co, none, genes)
  x$size <- lengths(raw_m[big]); x$rp_dominated <- rp_raw[big] >= 0.3; x }
# score one readout (symmetric finite matrix) in the three scopes + per-complex AP (cached)
score_readout <- function(name, S, macro = TRUE) {
  cf <- file.path(opt$out, "cache", paste0("readout__", name, ".rds"))
  if (file.exists(cf)) { say("cached: ", name); return(invisible(readRDS(cf))) }
  t1 <- Sys.time(); S[!is.finite(S)] <- 0; diag(S) <- 0
  out <- list(name = name, raw = pooled(S)); invisible(gc(verbose = FALSE))
  Sh <- rz_hub_stratify(S, H_att); out$hub_att <- pooled(Sh); rm(Sh); invisible(gc(verbose = FALSE))
  Sh <- rz_hub_stratify(S, H_mi); out$hub_mi <- pooled(Sh); rm(Sh); invisible(gc(verbose = FALSE))
  if (macro) out$macro <- macro_of(S)
  out$frac_nonzero <- mean(S[UT] > 0)
  saveRDS(out, cf); say(sprintf("scored %-22s (%.0f s)", name, as.numeric(difftime(Sys.time(), t1, units = "secs")))); invisible(out)
}

## ---- equivalence gate ------------------------------------------------------------------------------
# docs/plans/explore_ideas47_plan.md section 0: (a) robust_z V0 caches (all replicates, 3 scopes), (b) baselines metrics.csv,
# (c) corum_expert2 jackknife csv (CLR, K10, a3r on the seven expert2 truths). Stops the script on failure.
e47_gate <- function(res_clr, res_k10, res_a3r) {
  rows <- list(); add <- function(check, value, ok = value < 1e-9) rows[[length(rows) + 1L]] <<- data.frame(check = check, max_abs_diff = value, pass = ok)
  rz_truths <- c("DoRothEA_AB", "CORUM_all", "KD_response", "CORUM_noRP", "CORUM_RP", "CORUM_small", "CORUM_large")
  for (pr in list(c("V0__pair.rds", "clr"), c("V0__att9.rds", "k10"))) {
    ref <- readRDS(file.path("results/robust_z/k562/cache", pr[1])); mine <- if (pr[2] == "clr") res_clr else res_k10
    for (sc in c("raw", "hub_att", "hub_mi")) { rs <- if (sc == "raw") ref$pooled else ref[[sc]]
      d <- max(vapply(rz_truths, function(e) max(abs(mine[[sc]][[e]]$aupr - rs[[e]]$aupr), abs(mine[[sc]][[e]]$jack - rs[[e]]$jack)), 0))
      add(sprintf("(a) %s vs robust_z cache %s, scope %s, 7 truths x (AUPR + 100 replicates)", pr[2], pr[1], sc), d) } }
  bm <- utils::read.csv("results/baselines_k562ess_truthgenes_genie3/metrics.csv", stringsAsFactors = FALSE); bm <- bm[bm$dataset == "k562ess", ]
  mp <- c(DoRothEA_AB = "DoRothEA_AB", CollecTRI = "CollecTRI", CORUM = "CORUM_all", KD_response = "KD_response")
  for (m in list(c("clr", "clr"), c("attention", "k10"))) { mine <- if (m[2] == "clr") res_clr else res_k10
    d <- max(vapply(names(mp), function(tn) abs(mine$raw[[mp[[tn]]]]$aupr - bm$aupr[bm$method == m[1] & bm$truth == tn]), 0))
    add(sprintf("(b) %s vs baselines metrics.csv, 4 truths", m[2]), d) }
  jk <- utils::read.csv("results/corum_expert2/k562/t9_jackknife_EXPLORATORY.csv", stringsAsFactors = FALSE)
  mp2 <- c(CORUM_raw = "CORUM_all", CORUM_merged = "CORUM_merged", CORUM_merged_noRP = "CORUM_merged_noRP", CORUM_small_3to6 = "CORUM_small",
           CORUM_large_ge20 = "CORUM_large", DoRothEA_AB = "DoRothEA_AB", CollecTRI = "CollecTRI")
  for (m in list(c("CLR", "clr"), c("K10c50", "k10"), c("K10c50_a3r", "a3r"))) { mine <- switch(m[2], clr = res_clr, k10 = res_k10, a3r = res_a3r)
    d <- max(vapply(names(mp2), function(tn) abs(mine$raw[[mp2[[tn]]]]$aupr - jk$aupr[jk$truth == tn & jk$variant == m[1]]), 0))
    add(sprintf("(c) %s vs corum_expert2 t9_jackknife, 7 truths", m[2]), d) }
  enr <- res_a3r$raw$CORUM_all$aupr / base_rate[["CORUM_all"]]
  add(sprintf("(c) a3r CORUM_all enrichment %.4f prints as 7.95", enr), abs(round(enr, 2) - 7.95), abs(round(enr, 2) - 7.95) < 1e-9)
  g <- do.call(rbind, rows); utils::write.csv(g, file.path(opt$out, "equivalence_gate.csv"), row.names = FALSE); print(g, row.names = FALSE)
  if (!all(g$pass)) stop("EQUIVALENCE GATE FAILED: no new readout is computed or reported")
  say("EQUIVALENCE GATE PASSED (a)-(c)"); invisible(g)
}
