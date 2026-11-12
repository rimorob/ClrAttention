#!/usr/bin/env Rscript
# EXPLORATORY, NOT PREREGISTERED, NOTHING CONFIRMATORY (docs/plans/explore_idea3_plan.md, sections 5, 8, 9).
# Scoring stage: reads the FROZEN selection written by idea3_select.R (selected_depths.csv, selected_controls.csv) and the knee tables of
# idea3_knee_null.R, builds the readouts on the full-data operator, scores them against the truths (pooled AUPR, raw and hub-stratified with
# H = MI row mean and H = attention pi, paired delete-half gene jackknife 100 reps, seed 20261016), runs the equivalence gate and writes the
# tables and the literal pass/fail. No p-values.
# Usage (repo root): Rscript analysis/explore_idea3/idea3_score.R --dataset m3d|k562 [--workers 8] [--threads 8]
args <- commandArgs(TRUE)
opt <- list(dataset = "m3d", out = NULL, workers = 8L, threads = 8L, jack = 100L, jack_reg = 100L, alpha = 0.5, seed = 20261016L, quick = 0L, perm = 1L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); if (!key %in% names(opt)) stop("unknown argument: ", args[k])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[k + 1]) else args[k + 1] }
for (k in c("workers", "threads", "jack", "jack_reg", "seed", "quick", "perm")) opt[[k]] <- as.integer(opt[[k]])
if (is.null(opt$out)) opt$out <- file.path("results/explore_idea3", opt$dataset)
dir.create(file.path(opt$out, "cache"), recursive = TRUE, showWarnings = FALSE)
script_dir <- (function() { f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis/explore_idea3" })()
source(file.path(script_dir, "idea3_lib.R")); i3_source_all()
logf <- file.path(opt$out, "score.log"); say <- function(...) i3_say(logf, ...)
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))
none <- stats::setNames(character(0), character(0))

D <- i3_load(opt$dataset, opt$quick); X <- D$X; genes <- D$genes; G <- length(genes); rm(D)
say(sprintf("%s: %d genes x %d samples", opt$dataset, G, ncol(X)))
UT <- which(upper.tri(matrix(FALSE, G, G)))

## ---- the frozen selection (read BEFORE any truth is built) -------------------------------------------------
dec <- utils::read.csv(file.path(opt$out, "deciles.csv"), stringsAsFactors = FALSE); stopifnot(identical(dec$gene, genes))
sel <- utils::read.csv(file.path(opt$out, "selected_depths.csv")); selc <- utils::read.csv(file.path(opt$out, "selected_controls.csv"))
kdec <- utils::read.csv(file.path(opt$out, "knee_decile.csv")); kgene <- utils::read.csv(file.path(opt$out, "knee_gene_real.csv"))
td <- function(ls) { s <- sel[sel$op == "real" & sel$labelset == ls, ]; stats::setNames(s$t_d, s$decile) }
tg <- selc$param[selc$family == "global_grid"]
ppr_a <- selc$param[selc$family == "ppr"]; hk_tau <- selc$param[selc$family == "hk"]
say(sprintf("frozen selection: global grid t_g = %g; PPR alpha = %g; heat kernel tau = %g", tg, ppr_a, hk_tau))
say("t_d (MI deciles 1..10): ", paste(td("mi"), collapse = " "), " | (pi deciles): ", paste(td("pi"), collapse = " "))
kd <- function(s) { x <- kdec[kdec$stratifier == s, ]; stats::setNames(x$real_knee_snapped, x$decile) }
tvecs <- list(t9 = rep(9L, G), global_grid = rep(as.integer(tg), G),
              form1_mi = as.integer(td("mi")[dec$dec_mi]), form1_pi = as.integer(td("pi")[dec$dec_pi]),
              form2_mi = as.integer(kd("mi")[dec$dec_mi]), form2_pi = as.integer(kd("pi")[dec$dec_pi]), form2raw = as.integer(kgene$knee_t))
if (opt$perm > 0L && opt$dataset == "m3d") { set.seed(20261113L); pl <- list(); for (k in 1:4) pl[[k]] <- sample(dec$dec_mi)   # perm_seed stream of idea3_select.R
  tvecs$form1_perm1 <- as.integer(td("perm1")[pl[[1]]]) }
stopifnot(all(unlist(lapply(tvecs, function(v) all(v %in% c(I3_GRID, 9L))))))

## ---- full-data operator (package path) ----------------------------------------------------------------------
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
f$calibrate(method = "normal", combine = "stouffer")
M <- f$mi; f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L); A <- f$operator; f$release(); rm(f)
P <- rz_lazy_P(A, opt$alpha); rm(A)
H_att <- { q <- rz_stationary(P); outer(q, q, "+") }
mi_rm <- (rowSums(M) - diag(M)) / (G - 1); H_mi <- outer(mi_rm, mi_rm, "+"); rm(M)
stopifnot(max(abs(rank(mi_rm, ties.method = "first") - rank(dec$H_mi, ties.method = "first"))) == 0)

## ---- truths (as robust_z.R) + jackknife streams (seed 20261016) -----------------------------------------------
source("analysis/regulondb.R")
if (opt$dataset == "m3d") {
  to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
  rtab <- build_regulon_table("data/RegulonDBExtract", to_bn); operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) { ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck); list(pairs = comembership_pairs(reg, operon_of, genes)) })
  set.seed(opt$seed)
  draw_half <- function(gin) lapply(seq_len(opt$jack), function(r) sort(sample(gin, floor(length(gin) / 2))))
  halves <- list(); for (e in c("SC", "all")) halves[[e]] <- draw_half(bench[[e]]$pairs$genes_in)
  base_file <- "results/baselines_m3d/metrics.csv"; base_map <- c(SC = "SC", all = "all"); rz_map <- c(SC = "SC", all = "all")
} else {
  h5 <- "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"; set.seed(20261015L)
  rd <- function(p) rhdf5::h5read(h5, p)
  cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field)); if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL); if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
  pert <- cat_of("/obs", "gene_transcript"); pval <- as.numeric(rd("/obs/energy_test_p_value"))
  K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
  A_half <- sort(K[seq_len(floor(length(K) / 2))]); B_half <- setdiff(K, A_half); rm(rd, cat_of, pert, pval, K)
  tr <- list(DoRothEA_AB = dorothea_sets(genes), CORUM = corum_sets(genes), KD_response = response_sets(genes, path = h5, use = B_half),
             CollecTRI = collectri_sets(genes))
  cr <- tr$CORUM; cx <- utils::read.csv("analysis/corum_exploratory/corum_complexes.csv", stringsAsFactors = FALSE)
  sz <- lengths(cr$members); nm <- names(cr$members); rp <- nm %in% cx$complex[cx$group == "RP-dominated"]
  if (opt$quick == 0L) stopifnot(sum(rp) == 8L, all(sz == cx$size[match(nm, cx$complex)]))
  sub_of <- function(keep) list(members = cr$members[keep], class = cr$class[keep])
  mk <- function(reg) list(pairs = comembership_pairs(reg, none, genes))
  bench <- list(DoRothEA_AB = mk(tr$DoRothEA_AB), CORUM_all = mk(cr), KD_response = mk(tr$KD_response), CORUM_noRP = mk(sub_of(!rp)), CollecTRI = mk(tr$CollecTRI))
  set.seed(opt$seed)                           # robust_z.R / baselines.R stream order; the A-C draw only advances the stream
  draw_half <- function(gin) lapply(seq_len(opt$jack), function(r) sort(sample(gin, floor(length(gin) / 2))))
  gin_of <- function(reg) sort(unique(unlist(reg$members)))
  halves <- list(); halves$DoRothEA_AB <- draw_half(bench$DoRothEA_AB$pairs$genes_in)
  invisible(draw_half(gin_of(dorothea_sets(genes, c("A", "B", "C"))))); halves$CollecTRI <- draw_half(bench$CollecTRI$pairs$genes_in)
  halves$CORUM_all <- draw_half(bench$CORUM_all$pairs$genes_in); halves$KD_response <- draw_half(bench$KD_response$pairs$genes_in)
  for (r in seq_len(opt$jack_reg)) { k <- logical(G); k[sample(G, G - round(0.2 * G))] <- TRUE }   # advances the stream (d20 gene subsets of robust_z.R)
  halves$CORUM_noRP <- draw_half(bench$CORUM_noRP$pairs$genes_in)
  bench$KD_response <- NULL; halves$KD_response <- NULL
  base_file <- "results/baselines_k562ess_truthgenes_genie3/metrics.csv"; base_map <- c(DoRothEA_AB = "DoRothEA_AB", CollecTRI = "CollecTRI", CORUM_all = "CORUM"); rz_map <- c(DoRothEA_AB = "DoRothEA_AB", CORUM_all = "CORUM_all", CORUM_noRP = "CORUM_noRP")
}
base_rate <- vapply(bench, function(b) mean(b$pairs$label), 0)
for (e in names(bench)) say(sprintf("truth %-12s %8d pairs, base rate %.4f", e, length(bench[[e]]$pairs$label), base_rate[[e]]))
nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))

jack_one <- function(sc, b, keeps) { i <- b$pairs$i; j <- b$pairs$j; lab_ <- b$pairs$label
  unlist(parallel::mclapply(keeps, function(keep) { inh <- logical(G); inh[keep] <- TRUE; ix <- which(inh[i] & inh[j]); ap_fast(sc[ix], lab_[ix]) }, mc.cores = nc)) }
pooled <- function(S) { invisible(gc(verbose = FALSE)); lapply(stats::setNames(names(bench), names(bench)), function(e) { b <- bench[[e]]
  sc <- S[cbind(b$pairs$i, b$pairs$j)]; list(aupr = ap_fast(sc, b$pairs$label), jack = jack_one(sc, b, halves[[e]])) }) }
evaluate <- function(S) { S[!is.finite(S)] <- 0; diag(S) <- 0
  out <- list(raw = pooled(S)); out$hub_att <- pooled(rz_hub_stratify(S, H_att)); invisible(gc(verbose = FALSE))
  out$hub_mi <- pooled(rz_hub_stratify(S, H_mi)); out }

## ---- walk: store P^t at the grid and t = 9; accumulate the chosen PPR and heat-kernel -----------------------------------
tt <- 0:I3_TMAX; w_ppr <- (1 - ppr_a) * ppr_a^tt; w_ppr <- w_ppr / sum(w_ppr); w_hk <- stats::dpois(tt, hk_tau); w_hk <- w_hk / sum(w_hk)
keep_t <- sort(unique(c(I3_GRID, 9L))); Pts <- list(); K_ppr <- matrix(0, G, G); diag(K_ppr) <- w_ppr[1]; K_hk <- matrix(0, G, G); diag(K_hk) <- w_hk[1]
t0 <- Sys.time()
i3_walk(P, I3_TMAX, function(Pt, t) { K_ppr <<- K_ppr + w_ppr[t + 1L] * Pt; K_hk <<- K_hk + w_hk[t + 1L] * Pt
  if (t %in% keep_t) Pts[[as.character(t)]] <<- Pt })
say(sprintf("walk to t = %d done (%.0f s); stored depths: %s", I3_TMAX, as.numeric(difftime(Sys.time(), t0, units = "secs")), paste(keep_t, collapse = ",")))
rm(P)
build_S <- function(tv) { R <- matrix(0, G, G); for (t in unique(tv)) { r <- which(tv == t); R[r, ] <- Pts[[as.character(t)]][r, ] }; i3_symmetrise(R) }
readouts <- c(names(tvecs), "ppr", "hk")
S_of <- function(nm) if (nm == "ppr") i3_symmetrise(K_ppr) else if (nm == "hk") i3_symmetrise(K_hk) else build_S(tvecs[[nm]])

res <- list()
for (nm in c("t9", setdiff(readouts, "t9"))) {
  cf <- file.path(opt$out, "cache", sprintf("score_%s.rds", nm))
  if (file.exists(cf)) { res[[nm]] <- readRDS(cf); next }
  t1 <- Sys.time(); S <- S_of(nm); res[[nm]] <- evaluate(S); rm(S); invisible(gc(verbose = FALSE)); saveRDS(res[[nm]], cf)
  say(sprintf("readout %-12s scored (%.0f s)", nm, as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  if (nm == "t9") {      # ---- equivalence gate 1 and 3 (stop before anything else is reported) ----
    bm <- utils::read.csv(base_file, stringsAsFactors = FALSE); bm <- bm[bm$method == "attention", ]
    d1 <- vapply(names(base_map), function(e) abs(res$t9$raw[[e]]$aupr - bm$aupr[bm$truth == base_map[[e]]]), 0)
    rzp <- utils::read.csv(file.path("results/robust_z", opt$dataset, "pooled_metrics.csv"), stringsAsFactors = FALSE); rzp <- rzp[rzp$variant == "V0" & rzp$kind == "att9", ]
    d3 <- unlist(lapply(names(rz_map), function(e) { r <- rzp[rzp$truth == rz_map[[e]], ]
      c(abs(res$t9$hub_mi[[e]]$aupr / base_rate[[e]] - r$enr_hub_mi), abs(res$t9$hub_att[[e]]$aupr / base_rate[[e]] - r$enr_hub_att), abs(res$t9$raw[[e]]$aupr / base_rate[[e]] - r$enrichment)) }))
    writeLines(c(sprintf("gate 1: t = 9 readout vs %s (attention): max |diff| pooled AUPR %.3e over truths %s; tolerance 1e-9 %s", base_file, max(d1), paste(names(base_map), collapse = ","), if (max(d1) <= 1e-9) "PASS" else "FAIL"),
                 sprintf("gate 3: t = 9 enrichment (raw, hub MI, hub pi) vs results/robust_z/%s V0 att9: max |diff| %.3e over truths %s; tolerance 1e-9 %s", opt$dataset, max(d3), paste(names(rz_map), collapse = ","), if (max(d3) <= 1e-9) "PASS" else "FAIL")),
               file.path(opt$out, "equivalence_scoring.txt"))
    say(sprintf("EQUIVALENCE gate 1: max |diff| %.3e; gate 3: max |diff| %.3e", max(d1), max(d3)))
    if (opt$quick == 0L && (max(d1) > 1e-9 || max(d3) > 1e-9)) stop("equivalence gate FAILED")
  }
}

## ---- tables --------------------------------------------------------------------------------------------------------
long <- do.call(rbind, lapply(names(res), function(nm) do.call(rbind, lapply(names(res[[nm]]), function(md) do.call(rbind, lapply(names(res[[nm]][[md]]), function(e) {
  r <- res[[nm]][[md]][[e]]; data.frame(readout = nm, mode = md, truth = e, aupr = r$aupr, base_rate = base_rate[[e]], enrichment = r$aupr / base_rate[[e]], jack_sd = stats::sd(r$jack)) }))))))
utils::write.csv(long, file.path(opt$out, "scores.csv"), row.names = FALSE)
con <- do.call(rbind, lapply(setdiff(names(res), "t9"), function(nm) do.call(rbind, lapply(names(res[[nm]]), function(md) do.call(rbind, lapply(names(res[[nm]][[md]]), function(e) {
  a <- res[[nm]][[md]][[e]]; b <- res$t9[[md]][[e]]; d <- a$aupr - b$aupr; dj <- a$jack - b$jack; s <- stats::sd(dj)
  data.frame(readout = nm, mode = md, truth = e, aupr = a$aupr, aupr_t9 = b$aupr, diff_aupr = d, diff_lo = d - 1.96 * s, diff_hi = d + 1.96 * s,
             diff_enrichment = d / base_rate[[e]], diff_enr_lo = (d - 1.96 * s) / base_rate[[e]], diff_enr_hi = (d + 1.96 * s) / base_rate[[e]]) }))))))
utils::write.csv(con, file.path(opt$out, "contrasts.csv"), row.names = FALSE)
# literal pass rule on M3D SC (plan section 9); K562 has no SC: the same two lines are reported for DoRothEA A-B and CollecTRI
tr_pass <- if (opt$dataset == "m3d") "SC" else c("DoRothEA_AB", "CollecTRI")
pass <- do.call(rbind, lapply(setdiff(names(res), "t9"), function(nm) do.call(rbind, lapply(tr_pass, function(e) {
  g <- function(md) con$diff_aupr[con$readout == nm & con$truth == e & con$mode == md]
  data.frame(readout = nm, truth = e, hub_mi_diff = g("hub_mi"), hub_pi_diff = g("hub_att"), raw_diff = g("raw"),
             pass_primary_mi = g("hub_mi") >= 0.005 & g("raw") >= -0.005, pass_sensitivity_pi = g("hub_att") >= 0.005 & g("raw") >= -0.005) }))))
utils::write.csv(pass, file.path(opt$out, "pass_fail.csv"), row.names = FALSE)
print(pass, row.names = FALSE)
saveRDS(list(opt = opt, tvecs = tvecs, base_rate = base_rate, tg = tg, ppr_a = ppr_a, hk_tau = hk_tau), file.path(opt$out, "cache", "meta.rds"))
say("scoring done")
