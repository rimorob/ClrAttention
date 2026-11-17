#!/usr/bin/env Rscript
# One cell (N, noise, marginal) of the planted-module simulation, seeds in sequence; each compendium is cached as
# <out>/units/N<N>_<noise>_<marginal>_s<seed>.rds and never recomputed. Spec: docs/plans/sim_planted_choices.md.
# Usage (repo root): Rscript analysis/sim_planted/sim_run.R --N 150 --noise gaussian --marginal identity [--seeds 1:10]
#   [--threads 4] [--workers 4] [--out results/sim_planted] [--gate 1] [--mibench 1]
args <- commandArgs(TRUE)
opt <- list(N = "150", noise = "gaussian", marginal = "identity", seeds = "1:10", threads = "4", workers = "4",
            out = "results/sim_planted", gate = "0", mibench = "0")
for (k in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[k])]] <- args[k + 1]
N <- as.integer(opt$N); seeds <- eval(parse(text = opt$seeds)); threads <- as.integer(opt$threads); nc <- as.integer(opt$workers)
source("analysis/sim_planted/sim_lib.R")
dir.create(file.path(opt$out, "units"), recursive = TRUE, showWarnings = FALSE)
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n") }
vmhwm <- function() { l <- grep("^VmHWM", readLines("/proc/self/status"), value = TRUE); as.numeric(gsub("[^0-9]", "", l)) / 1048576 }
Nidx <- match(N, c(150L, 400L, 800L)); stopifnot(!is.na(Nidx))

## ---- gates (toy sizes) ------------------------------------------------------------------------------------------
if (opt$gate == "1") {
  G0 <- 40L; U0 <- sim_ut(G0); m <- matrix(0, G0, G0); m[U0$UT] <- seq_along(U0$UT)
  stopifnot(all(pair_pos(U0$ri, U0$ci) == seq_along(U0$UT)), all(m[cbind(U0$ri, U0$ci)] == seq_along(U0$UT)))
  set.seed(1); S <- matrix(runif(G0^2), G0, G0); S <- (S + t(S)) / 2; S[sample(G0^2, 50)] <- 0
  h <- stats::rnorm(G0); H <- outer(h, h, "+")
  stopifnot(max(abs(rz_hub_stratify(S, H)[U0$UT] - stratify_vec(S[U0$UT], hub_groups(H[U0$UT])))) < 1e-14)
  set.seed(2); sc <- stats::runif(500); sc[sample(500, 100)] <- 0.5; lb <- stats::rbinom(500, 1, 0.1)
  stopifnot(abs(ap_fast(sc, lb) - pr_summary(sc, lb)[["aupr"]]) < 1e-12)
  say("gates passed (pair index, stratification = rz_hub_stratify, ap_fast = pr_summary aupr)")
}

## ---- one compendium ----------------------------------------------------------------------------------------------
U <- sim_ut(); G <- SIM_G
run_unit <- function(s) {
  tm <- list(); t0 <- Sys.time(); tick <- function(nm) { tm[[nm]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs")); t0 <<- Sys.time(); say(sprintf("s%02d %s %.0f s", s, nm, tm[[nm]])) }
  st <- sim_structure(s); tr <- sim_truth(st, U)
  X <- sim_data(st, N, opt$noise, opt$marginal, Nidx); tick("generate")
  fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = threads); tick("MI")
  fit$calibrate(method = "normal", combine = "stouffer"); M <- fit$mi; C <- fit$clr_scores; tick("CLR")
  fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A <- fit$operator; tick("operator keff10")
  fit$build_operator(alpha = 0.5, softmax_keff = 40, softmax_cap = 200L); A40 <- fit$operator; tick("operator keff40")
  fit$release(); rm(fit)
  mi_rm <- (rowSums(M) - diag(M)) / (G - 1); rm(M)
  P <- sim_P(A); pi <- rz_stationary(P)
  hubpi <- logical(G); hubpi[order(pi, decreasing = TRUE)[seq_len(round(0.01 * G))]] <- TRUE
  hubpl <- logical(G); hubpl[st$hubs] <- TRUE; hubsets <- list(pi = hubpi, planted = hubpl)
  groups <- list(true = hub_groups(st$hubsum[U$ri] + st$hubsum[U$ci]), mi = hub_groups(mi_rm[U$ri] + mi_rm[U$ci]))
  vecs <- list(absr = abs(stats::cor(t(X)))[U$UT], clr = C[U$UT]); mech <- list(); dmech <- list(); dmet <- list(); ident <- NULL
  ainshare <- vapply(st$members, function(m) (sum(A[m, m]) - sum(diag(A)[m])) / length(m), 0)
  ainshare40 <- vapply(st$members, function(m) (sum(A40[m, m]) - sum(diag(A40)[m])) / length(m), 0)
  tick("setup")
  ## main walk, t = 1..72 (descriptive depth grid), one pass
  Pt <- NULL; light <- function(v, meth, t) { m <- sim_metrics(v, tr, U, st, hubsets, full = FALSE); m$method <- meth; m$t <- t; m }
  for (t in 1:72) {
    Pt <- if (is.null(Pt)) as.matrix(P) else step_walk(Pt, P)
    if (t == 1L) {
      p1 <- tca_parts(Pt); vecs$tcaclr_t1 <- stouffer0(p1$Z)[U$UT]
      offd <- A > 0; diag(offd) <- FALSE; ncand <- rowSums(offd); ok <- ncand > 10L
      Zd <- matrix(0, G, G)
      for (i in which(ok)) { cand <- which(offd[i, ]); z <- (C[i, cand] - mean(C[i, cand])) / stats::sd(C[i, cand]); z[z < 0] <- 0; Zd[i, cand] <- z }
      Sd <- stouffer0(Zd); St1 <- stouffer0(p1$Z)
      ident <- list(maxdiff_Z = max(abs(p1$Z[ok, ] - Zd[ok, ])), maxdiff_S = max(abs(St1[ok, ok] - Sd[ok, ok])),
                    rows_le10 = sum(!ok), rows_gt10 = sum(ok)); rm(Zd, Sd, St1, offd)
      mech$t1 <- walk_mech(Pt, st, p1$sigma, p1$nreach); mech$t1$A_inshare <- ainshare; rm(p1)
    }
    if (t %in% c(3L, 6L, 9L, 18L, 36L, 72L)) {
      pp <- tca_parts(Pt); vt <- stouffer0(pp$Z)[U$UT]; va <- sym0(Pt)[U$UT]
      dmech[[length(dmech) + 1L]] <- data.frame(t = t, sigma_median = stats::median(pp$sigma, na.rm = TRUE), cv = stats::sd(pp$sigma, na.rm = TRUE) / mean(pp$sigma, na.rm = TRUE),
                                                 unreach_share = 1 - sum(pp$nreach) / (G * (G - 1)))
      dmet[[length(dmet) + 1L]] <- rbind(light(va, "aclr", t), light(vt, "tcaclr", t))
      if (t == 9L) {
        vecs$aclr <- va; vecs$tcaclr <- vt; mech$t9 <- walk_mech(Pt, st, pp$sigma, pp$nreach)
        lf <- log(pmax(Pt, 1e-300)); lp <- rep(log(pi), each = G)
        vecs$pi_lift <- clr_R(((lf - lp) + t(lf - lp)) / 2)[U$UT]
        lc <- rep(log(colMeans(Pt)), each = G); vecs$colmean_lift <- clr_R(((lf - lc) + t(lf - lc)) / 2)[U$UT]; rm(lf, lp, lc)
        vecs$infl_1.5 <- clr_R(Pt^1.5)[U$UT]; vecs$infl_2 <- clr_R(Pt^2)[U$UT]
        vecs$hubscore <- (pi[U$ri] + pi[U$ci])
      }
      rm(pp, vt, va)
    }
    if (t %% 8L == 0L) invisible(gc(verbose = FALSE))
  }
  rm(Pt); invisible(gc(verbose = FALSE)); tick("walk t1..72 + readouts")
  mech$t9$A_inshare <- ainshare
  ## Sinkhorn-balanced lazy operator, walk to 9
  sk <- sinkhorn_P(P); Pb <- Matrix::Matrix(sk$P, sparse = TRUE); Pt <- as.matrix(Pb); for (t in 2:9) Pt <- step_walk(Pt, Pb)
  pp <- tca_parts(Pt); vecs$sink_tca <- stouffer0(pp$Z)[U$UT]; vecs$sink_raw <- sym0(Pt)[U$UT]; mech$sink <- walk_mech(Pt, st, pp$sigma, pp$nreach)
  mech$sink$sinkhorn_iters <- sk$iters; mech$sink$sinkhorn_resid <- sk$resid; rm(Pt, pp, sk, Pb); invisible(gc(verbose = FALSE)); tick("sinkhorn arm")
  ## wide-row arm
  P40 <- sim_P(A40); Pt <- as.matrix(P40); for (t in 2:9) Pt <- step_walk(Pt, P40)
  pp <- tca_parts(Pt); vecs$tcaclr_k40 <- stouffer0(pp$Z)[U$UT]; vecs$aclr_k40 <- sym0(Pt)[U$UT]; mech$k40 <- walk_mech(Pt, st, pp$sigma, pp$nreach)
  mech$k40$A_inshare <- ainshare40; rm(Pt, pp, P40, A40); invisible(gc(verbose = FALSE)); tick("wide arm")
  ## metrics (parallel over methods)
  meths <- setdiff(names(vecs), "hubscore")
  res <- parallel::mclapply(meths, function(nm) { m <- sim_metrics(vecs[[nm]], tr, U, st, hubsets, groups); m$method <- nm; m }, mc.cores = nc)
  bad <- vapply(res, function(r) inherits(r, "try-error") || !is.data.frame(r), NA); if (any(bad)) stop("metric worker failed: ", paste(meths[bad], collapse = ", "))
  hs <- sim_metrics(vecs$hubscore, tr, U, st, hubsets, full = FALSE); hs$method <- "hubscore"
  metrics <- rbind(do.call(rbind, res), hs); tick("metrics")
  list(meta = list(seed = s, N = N, noise = opt$noise, marginal = opt$marginal, sizes = st$sizes, lam = st$lam, big = st$big, vmhwm_gb = vmhwm()),
       metrics = metrics, depth_metrics = do.call(rbind, dmet), depth_mech = do.call(rbind, dmech),
       mech = mech, ident = ident, timings = tm)
}

if (opt$mibench == "1") {            # time of the MI stage alone at the three sample sizes (pilot only)
  st <- sim_structure(1L); for (nn in c(150L, 400L, 800L)) { X <- sim_data(st, nn, "gaussian", "identity", match(nn, c(150L, 400L, 800L)))
    t0 <- Sys.time(); invisible(clr::bspline_mi(X, bins = "hg", threads = threads)); say(sprintf("MI N=%d threads=%d: %.1f s", nn, threads, as.numeric(difftime(Sys.time(), t0, units = "secs")))) }
}
for (s in seeds) {
  f <- sprintf("%s/units/N%d_%s_%s_s%02d.rds", opt$out, N, opt$noise, opt$marginal, s)
  if (file.exists(f)) { say("cached ", basename(f)); next }
  t1 <- Sys.time()
  r <- tryCatch(run_unit(s), error = function(e) { say("FAILED s", s, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) { r$meta$wall_s <- as.numeric(difftime(Sys.time(), t1, units = "secs")); saveRDS(r, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    say(sprintf("done %s: %.0f s wall, R process peak %.1f GB", basename(f), r$meta$wall_s, r$meta$vmhwm_gb)) }
  invisible(gc(verbose = FALSE))
}
