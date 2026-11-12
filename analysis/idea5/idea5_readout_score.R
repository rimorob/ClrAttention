#!/usr/bin/env Rscript
# EXPLORATORY, NOT PREREGISTERED (docs/plans/explore_idea5_plan.md). Learned-logit readout (the project's softmax, lazy
# walk and symmetric walk mass), equivalence gates G1-G3, label-free depth curves, and scoring on the M3D RegulonDB
# truths (SC and all evidence; raw, pi-stratified and MI-row-mean-stratified pooled AUPR; 100-replicate paired
# delete-half gene jackknife of analysis/robust_z/robust_z.R). No p-values.
# Usage (repo root): Rscript analysis/idea5/idea5_readout_score.R --prep results/idea5/prep --models results/idea5/models
#     --out results/idea5/scores [--truth 1] [--workers 8] [--tmax 32] [--jack 100]
args <- commandArgs(TRUE)
opt <- list(prep = "results/idea5/prep", models = "results/idea5/models", out = "results/idea5/scores", truth = 1L,
            workers = 8L, threads = 8L, tmax = 32L, jack = 100L, seed = 20261016L, t = 9L, alpha = 0.5)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); if (!key %in% names(opt)) stop("unknown argument: ", args[k])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[k + 1]) else args[k + 1] }
for (k in c("truth", "workers", "threads", "tmax", "jack", "seed", "t")) opt[[k]] <- as.integer(opt[[k]])
dir.create(file.path(opt$out, "cache"), recursive = TRUE, showWarnings = FALSE)
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = file.path(opt$out, "readout_score.log"), append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/idea5/idea5_lib.R"); source("analysis/robust_z/robust_z_lib.R")
source("analysis/regulondb.R"); source("analysis/regulons.R")   # function definitions only; no truth file is read until the truth section
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))
meta <- jsonlite::fromJSON(file.path(opt$prep, "meta.json")); G <- meta$G; N <- meta$N
genes <- readLines(file.path(opt$prep, "genes.txt")); stopifnot(length(genes) == G)
X <- matrix(readBin(file.path(opt$prep, "X.f64"), "numeric", n = G * N, size = 8L), G, N, dimnames = list(genes, NULL))
split <- utils::read.csv(file.path(opt$prep, "split.csv"))$half
rdZ <- function(nm) idea5_read_mat(file.path(opt$prep, sprintf("Z_%s.f64", nm)), G, size = 8L)
Z <- list(full = rdZ("full"), h1 = rdZ("h1"), h2 = rdZ("h2")); tau0 <- meta$tau0
UT <- which(upper.tri(matrix(FALSE, G, G)))
nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))
sets <- list(full = rep(TRUE, N), h1 = split == 1, h2 = split == 2)

## ---- gates G3 and G1 ----------------------------------------------------------------------------------------
eq_f <- file.path(opt$out, "equivalence.csv")
if (!file.exists(eq_f)) {
  Tm <- idea5_row_temps(Z$full); A1 <- clr:::.softmax_rows(Z$full, keff = 10, cap = 50L); A2 <- idea5_softmax_from_temps(Z$full, Tm)
  g3 <- max(abs(A1 - A2)); rm(A1, A2)
  rows <- list(data.frame(gate = "G3", data = "full", check = "max |softmax weights from ported temperatures - clr:::.softmax_rows|", value = g3))
  for (m in c("full", "h1", "h2")) {
    f <- idea5_fit(X[, sets[[m]]], opt$threads); Cp <- f$clr_scores
    f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L); Ap <- f$operator
    Pp <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(Ap, sparse = TRUE)
    Pt <- as.matrix(Pp); for (s in seq_len(opt$t - 1L)) Pt <- as.matrix(Pt %*% Pp); Sp <- (Pt + t(Pt)) / 2; rm(Pt)
    tag <- sprintf("gate_clr_d32_%s", m)
    D <- idea5_read_mat(file.path(opt$models, paste0(tag, ".D.f32")), G, 4L); gam <- idea5_read_vec(file.path(opt$models, paste0(tag, ".gamma.f64")), G)
    C <- idea5_score_matrix(D, Z[[m]], exp(gam), tau0[[m]])
    A <- rz_operator(C); P <- rz_lazy_P(A); S9 <- rz_walk(P, opt$t)[[as.character(opt$t)]]
    rows[[length(rows) + 1L]] <- data.frame(gate = "G1", data = m, check = c("exported D and gamma are exactly 0 (max |D|, max |gamma|)", "max |Z file - package clr_scores|", "max |C(learned readout, step 0) - package clr_scores|",
      "max |operator - package operator|", sprintf("max |S_%d - package walk mass|", opt$t)),
      value = c(max(abs(D), abs(gam)), max(abs(Z[[m]] - Cp)), max(abs(C - Cp)), max(abs(A - Ap)), max(abs(S9 - Sp))))
    rm(D, C, A, P, S9, Sp, Pp, Ap, Cp); f$release(); rm(f); invisible(gc(verbose = FALSE))
  }
  eq <- do.call(rbind, rows); utils::write.csv(eq, eq_f, row.names = FALSE); print(eq, row.names = FALSE)
  if (any(eq$value[eq$gate == "G1" & !grepl("exactly 0", eq$check)] >= 1e-8) || g3 >= 1e-10 || any(eq$value[grepl("exactly 0", eq$check)] != 0))
    stop("EQUIVALENCE GATE FAILED (G1 < 1e-8, G3 < 1e-10)")
  say("GATES G1 and G3 PASSED")
} else say("equivalence.csv present (gates passed earlier)")

## ---- models ----------------------------------------------------------------------------------------------------
fl <- sub("\\.json$", "", list.files(opt$models, "\\.json$"))
mm <- do.call(rbind, lapply(fl, function(tag) {
  m <- regmatches(tag, regexec("^(clr|scratch)_d([0-9]+)_s([0-9]+)_(h1|h2|full)$", tag))[[1]]
  if (length(m)) return(data.frame(tag = tag, kind = if (m[5] == "full") "refit" else "half", init = m[2], d = as.integer(m[3]), seed = as.integer(m[4]), mode = m[5]))
  m <- regmatches(tag, regexec("^ctrl_scratch_d([0-9]+)_s([0-9]+)$", tag))[[1]]
  if (length(m)) return(data.frame(tag = tag, kind = "control", init = "scratch", d = as.integer(m[2]), seed = as.integer(m[3]), mode = "full"))
  NULL }))
mm <- mm[order(mm$kind, mm$init, mm$d, mm$seed, mm$mode), ]
say(sprintf("models found: %d half, %d refit, %d control", sum(mm$kind == "half"), sum(mm$kind == "refit"), sum(mm$kind == "control")))
get_C <- function(tag, init, mode) {
  D <- idea5_read_mat(file.path(opt$models, paste0(tag, ".D.f32")), G, 4L)
  if (init == "clr") idea5_score_matrix(D, Z[[mode]], exp(idea5_read_vec(file.path(opt$models, paste0(tag, ".gamma.f64")), G)), tau0[[mode]])
  else idea5_score_matrix(D, NULL, NULL, tau0$full)
}

## ---- label-free depth curves (lab_splits.R criterion), half models and half-CLR references -----------------------
cur_f <- file.path(opt$out, "depth_curves.csv")
if (!file.exists(cur_f)) {
  jobs <- c(lapply(which(mm$kind == "half"), function(r) list(label = mm$tag[r], C = function() get_C(mm$tag[r], mm$init[r], mm$mode[r]), val = if (mm$mode[r] == "h1") "h2" else "h1", dir = mm$mode[r])),
            lapply(c("h1", "h2"), function(h) list(label = paste0("CLRhalf_", h), C = function() Z[[h]], val = if (h == "h1") "h2" else "h1", dir = h)))
  say(sprintf("depth curves: %d operators, t = 1..%d, %d workers", length(jobs), opt$tmax, min(nc, length(jobs))))
  curve_one <- function(j) {
    C <- j$C(); A <- rz_operator(C); P <- rz_lazy_P(A); rm(C, A)
    sv <- Z[[j$val]][UT]; target <- as.integer(sv >= stats::quantile(sv, 0.99)); rm(sv)
    crit <- numeric(opt$tmax); Pt <- NULL
    for (t in seq_len(opt$tmax)) { Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P); S <- (Pt + t(Pt)) / 2; crit[t] <- ap_fast(S[UT], target); rm(S) }
    data.frame(label = j$label, dir = j$dir, t = seq_len(opt$tmax), criterion = crit)
  }
  res <- parallel::mclapply(jobs, curve_one, mc.cores = min(nc, length(jobs)), mc.preschedule = FALSE)
  bad <- vapply(res, function(x) !is.data.frame(x), NA); if (any(bad)) stop("depth stage failed: ", paste(unlist(res[bad]), collapse = " | "))
  utils::write.csv(do.call(rbind, res), cur_f, row.names = FALSE); rm(res); invisible(gc(verbose = FALSE))
  say("depth curves written")
}
cur <- utils::read.csv(cur_f, stringsAsFactors = FALSE)
tb <- do.call(rbind, lapply(split(cur, cur$label), function(x) data.frame(label = x$label[1], dir = x$dir[1], t_best = x$t[which.max(x$criterion)], crit_best = max(x$criterion))))
utils::write.csv(tb, file.path(opt$out, "depth_tbest.csv"), row.names = FALSE)
ts_cfg <- do.call(rbind, lapply(split(mm[mm$kind == "half", ], list(mm$init[mm$kind == "half"], mm$d[mm$kind == "half"], mm$seed[mm$kind == "half"]), drop = TRUE), function(x) {
  tb1 <- tb$t_best[match(x$tag[x$mode == "h1"], tb$label)]; tb2 <- tb$t_best[match(x$tag[x$mode == "h2"], tb$label)]
  data.frame(init = x$init[1], d = x$d[1], seed = x$seed[1], t_best_h1 = tb1, t_best_h2 = tb2, t_star = as.integer(round(exp(mean(log(c(tb1, tb2))))))) }))
ts_cfg$at_tmax_edge <- ts_cfg$t_best_h1 == opt$tmax | ts_cfg$t_best_h2 == opt$tmax
utils::write.csv(ts_cfg, file.path(opt$out, "tstar_by_config.csv"), row.names = FALSE)
tclr <- as.integer(round(exp(mean(log(tb$t_best[match(c("CLRhalf_h1", "CLRhalf_h2"), tb$label)]))))); say(sprintf("t* of the half-CLR operators on this lab split: %d (directions %s)", tclr, paste(tb$t_best[match(c("CLRhalf_h1", "CLRhalf_h2"), tb$label)], collapse = ", ")))
if (opt$truth == 0L) { say("truth = 0: stopping before any truth is read (smoke mode)"); quit(save = "no") }

## ---- truths (robust_z.R M3D branch) -------------------------------------------------------------------------------
if (opt$truth == 2L) {     # code smoke test only: random pairs and random labels, no truth file is read
  set.seed(5L)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) { n <- 30000L; i <- sample.int(G - 1L, n, replace = TRUE); j <- i + sample.int(G - 1L, n, replace = TRUE); keep <- j <= G
    i <- i[keep]; j <- j[keep]; u <- !duplicated(paste(i, j)); i <- i[u]; j <- j[u]
    list(reg = list(members = list()), pairs = list(i = i, j = j, label = stats::rbinom(length(i), 1L, 0.1), genes_in = seq_len(G))) })
} else {
  to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
  rtab <- build_regulon_table("data/RegulonDBExtract", to_bn); operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
    ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck); list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes)) })
}
TR <- names(bench); base_rate <- vapply(bench, function(b) mean(b$pairs$label), 0)
for (e in TR) say(sprintf("truth %-4s %5d regulons, %8d pairs, base rate %.4f", e, length(bench[[e]]$reg$members), length(bench[[e]]$pairs$label), base_rate[[e]]))
set.seed(opt$seed)
halves <- list(); for (e in TR) { gin <- bench[[e]]$pairs$genes_in; halves[[e]] <- lapply(seq_len(opt$jack), function(r) sort(sample(gin, floor(length(gin) / 2)))) }
pos_e <- lapply(bench, function(b) as.integer((b$pairs$j - 1) * (b$pairs$j - 2) / 2 + b$pairs$i)); stopifnot(all(unlist(lapply(TR, function(e) UT[pos_e[[e]]] == (bench[[e]]$pairs$j - 1) * G + bench[[e]]$pairs$i))))
lab_e <- lapply(bench, function(b) b$pairs$label)
ix_e <- lapply(TR, function(e) { b <- bench[[e]]; lapply(halves[[e]], function(keep) { inh <- logical(G); inh[keep] <- TRUE; which(inh[b$pairs$i] & inh[b$pairs$j]) }) }); names(ix_e) <- TR

## ---- stratifiers: pi of the full-data CLR attention walk and the row mean of raw MI (both fixed for all methods) ----
P0 <- rz_lazy_P(rz_operator(Z$full)); q_att <- rz_stationary(P0); rm(P0)
mi_rm <- utils::read.csv(file.path(opt$prep, "mi_rowmean_full.csv"))$mi_rowmean
r_ut <- (UT - 1L) %% G + 1L; c_ut <- (UT - 1L) %/% G + 1L
strata <- function(H) { h <- H[r_ut] + H[c_ut]; cut(h, unique(stats::quantile(h, seq(0, 1, length.out = 11L))), include.lowest = TRUE) }
g_att <- strata(q_att); g_mi <- strata(mi_rm); rm(r_ut, c_ut)
pct <- function(x) rank(x) / length(x)
evaluate_ut <- function(sv) {
  out <- list()
  for (scope in c("pooled", "hub_att", "hub_mi")) {
    v <- switch(scope, pooled = sv, hub_att = stats::ave(sv, g_att, FUN = pct), hub_mi = stats::ave(sv, g_mi, FUN = pct))
    out[[scope]] <- lapply(stats::setNames(TR, TR), function(e) { sc <- v[pos_e[[e]]]; lb <- lab_e[[e]]
      list(aupr = ap_fast(sc, lb), jack = vapply(ix_e[[e]], function(ix) ap_fast(sc[ix], lb[ix]), 0)) })
  }
  out
}

## ---- scoring tasks --------------------------------------------------------------------------------------------------
tstar_of <- function(init, d, seed) ts_cfg$t_star[ts_cfg$init == init & ts_cfg$d == d & ts_cfg$seed == seed]
tasks <- list(list(id = "CLR_full", C = function() Z$full, ro = c(pair = NA, att1 = 1L, att9 = 9L)))
for (h in c("h1", "h2")) tasks[[length(tasks) + 1L]] <- list(id = paste0("CLRhalf_", h), C = local({ hh <- h; function() Z[[hh]] }),
  ro = c(pair = NA, att1 = 1L, att9 = 9L, attT = as.integer(tb$t_best[tb$label == paste0("CLRhalf_", h)])))
for (r in seq_len(nrow(mm))) tasks[[length(tasks) + 1L]] <- local({ rr <- r; k <- mm$kind[rr]; ts <- if (k == "control") NA else tstar_of(mm$init[rr], mm$d[rr], mm$seed[rr])
  list(id = mm$tag[rr], C = function() get_C(mm$tag[rr], mm$init[rr], mm$mode[rr]),
       ro = switch(k, refit = c(pair = NA, att1 = 1L, att9 = 9L, attT = ts), half = c(att9 = 9L, attT = ts), control = c(pair = NA, att1 = 1L, att9 = 9L))) })
score_task <- function(tk) {
  f <- file.path(opt$out, "cache", paste0(tk$id, ".rds")); if (file.exists(f)) return("cached")
  t0 <- Sys.time(); C <- tk$C(); res <- list(); ro <- tk$ro
  if ("pair" %in% names(ro)) res$pair <- evaluate_ut(C[UT])
  ts <- unique(as.integer(ro[!is.na(ro)])); if (length(ts)) {
    A <- rz_operator(C); P <- rz_lazy_P(A); rm(A); W <- rz_walk(P, ts); rm(P)
    for (nm in names(ro)[!is.na(ro)]) res[[nm]] <- evaluate_ut(W[[as.character(ro[[nm]])]][UT]); rm(W) }
  res$t <- ro; res$isolated_rows <- { A <- rz_operator(C); mean(diag(A) == 1) }
  tmp <- paste0(f, ".part"); saveRDS(res, tmp); file.rename(tmp, f); rm(C); invisible(gc(verbose = FALSE))
  sprintf("%s %.0f s", tk$id, as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
## G2 first: the full-data CLR through this scoring code against the robust_z V0 cache
say("scoring CLR_full first (gate G2)"); invisible(score_task(tasks[[1]]))
chk <- list(); rzc <- "results/robust_z/m3d/cache"
if (opt$truth == 1L && file.exists(file.path(rzc, "V0__pair.rds"))) {
  mine <- readRDS(file.path(opt$out, "cache", "CLR_full.rds")); ref_p <- readRDS(file.path(rzc, "V0__pair.rds")); ref_a <- readRDS(file.path(rzc, "V0__att9.rds"))
  for (ro in c("pair", "att9")) { ref <- if (ro == "pair") ref_p else ref_a; for (sc in c("pooled", "hub_att", "hub_mi")) for (e in TR)
    chk[[length(chk) + 1L]] <- data.frame(readout = ro, scope = sc, truth = e, d_aupr = abs(mine[[ro]][[sc]][[e]]$aupr - ref[[sc]][[e]]$aupr), d_jack = max(abs(mine[[ro]][[sc]][[e]]$jack - ref[[sc]][[e]]$jack))) }
  chk <- do.call(rbind, chk); utils::write.csv(chk, file.path(opt$out, "gate_G2_vs_robust_z.csv"), row.names = FALSE); print(chk, row.names = FALSE)
  if (any(chk$d_aupr > 1e-10) || any(chk$d_jack > 1e-10)) stop("GATE G2 FAILED: the scoring code does not reproduce results/robust_z/m3d V0")
  say("GATE G2 PASSED (scoring code reproduces robust_z V0 pooled AUPR and the 100 jackknife values to 1e-10)")
} else say("G2 not checked (truth = 2 smoke, or results/robust_z/m3d/cache absent)")
say(sprintf("scoring %d tasks with %d workers", length(tasks) - 1L, nc))
st <- parallel::mclapply(tasks[-1], score_task, mc.cores = nc, mc.preschedule = FALSE)
bad <- vapply(st, function(x) !is.character(x), NA); if (any(bad)) { say("scoring failed for: ", paste(vapply(tasks[-1][bad], function(t) t$id, ""), collapse = ", ")); print(st[bad][1]) }
say("timings: ", paste(unlist(st[!bad]), collapse = "; "))
saveRDS(list(base_rate = base_rate, truths = TR, mm = mm, ts_cfg = ts_cfg, tasks = vapply(tasks, function(t) t$id, ""), opt = opt), file.path(opt$out, "cache", "meta.rds"))
say("done")
