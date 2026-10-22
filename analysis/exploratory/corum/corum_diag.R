#!/usr/bin/env Rscript
# EXPLORATORY, NOT PREREGISTERED. CORUM-loss diagnosis and cheap fixes for CLR
# attention on the K562 essential truth-gene run (same input matrix and truth
# sets as results/baselines_k562ess_truthgenes_genie3). Width (keff) and
# multi-head variants are deliberately NOT scored here: they are preregistered
# (docs/preregistration_keff.md) and their truths may be touched only once.
#
# Usage: Rscript corum_diag.R [quick_genes]   (quick_genes > 0: bug test on a subset)
args <- commandArgs(trailingOnly = TRUE)
quick <- if (length(args)) as.integer(args[1]) else 0L
setwd("/home/boris/Projects/ClrAttention")
SP <- "/tmp/claude-1000/-home-boris-Projects-ClrAttention/b33914c9-0cfd-485d-8dd4-cf8a2a22b632/scratchpad"
OUT <- file.path(SP, if (quick > 0L) "corum_quick" else "corum")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUT, "run.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R")
NC <- 4L; NJ <- 100L
RUN <- "results/baselines_k562ess_truthgenes_genie3"
T0 <- Sys.time()

## ---- data and truths (exact sets of the truth-gene run) ---------------------------------------
X <- as.matrix(utils::read.delim(file.path(RUN, "expression.tsv"), row.names = 1, check.names = FALSE))
rmt <- utils::read.csv(file.path(RUN, "regulon_members.csv"), stringsAsFactors = FALSE)
if (quick > 0L) {
  rp <- unique(unlist(strsplit(rmt$members[rmt$truth == "CORUM" & grepl("ibosom", rmt$regulon)], ";")))
  set.seed(1); keepg <- union(intersect(rp, rownames(X)), sample(rownames(X), quick))
  X <- X[sort(match(keepg, rownames(X))), ]
}
genes <- rownames(X); G <- length(genes)
say(sprintf("%d genes x %d samples", G, ncol(X)))
cls <- c(DoRothEA_AB = "TF", DoRothEA_ABC = "TF", CollecTRI = "TF", CORUM = "complex", KD_response = "response")
tr <- lapply(names(cls), function(nm) {
  x <- rmt[rmt$truth == nm, ]
  m <- lapply(strsplit(x$members, ";", fixed = TRUE), function(g) { v <- match(g, genes); sort(v[!is.na(v)]) })
  if (quick == 0L && any(lengths(m) != lengths(strsplit(x$members, ";", fixed = TRUE)))) stop("unmatched members in ", nm)
  names(m) <- x$regulon; m <- m[lengths(m) >= 3L]
  list(members = m, class = stats::setNames(rep(cls[[nm]], length(m)), names(m)))
})
names(tr) <- names(cls)
none <- stats::setNames(character(0), character(0))
bench <- lapply(tr, function(r) list(reg = r, pairs = comembership_pairs(r, none, genes)))
for (nm in names(bench)) say(sprintf("truth %-12s %4d sets, %8d pairs, base rate %.4f", nm, length(bench[[nm]]$reg$members),
                                     length(bench[[nm]]$pairs$label), mean(bench[[nm]]$pairs$label)))
set.seed(20261016L)                       # the delete-half gene splits of baselines.R (same RNG use)
halves <- lapply(bench, function(b) lapply(seq_len(NJ), function(r) {
  keep <- sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2)); inh <- logical(G); inh[keep] <- TRUE; inh }))

## ---- scoring -----------------------------------------------------------------------------------
ap_sorted <- function(s, l) { n <- length(s); if (!n || !sum(l)) return(NA_real_)
  last <- c(s[-1L] != s[-n], TRUE); tp <- cumsum(l)[last]; k <- which(last)
  sum(diff(c(0, tp)) * tp / k) / sum(l) }
JT <- c("DoRothEA_AB", "CollecTRI", "CORUM", "KD_response")    # jackknifed truths (ABC: point estimate)
E <- list()
score <- function(name, S) {
  t1 <- Sys.time()
  E[[name]] <<- lapply(stats::setNames(names(bench), names(bench)), function(nm) {
    b <- bench[[nm]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]; sc[!is.finite(sc)] <- 0
    o <- order(sc, decreasing = TRUE, method = "radix"); so <- sc[o]; lo <- b$pairs$label[o]
    out <- list(aupr = ap_sorted(so, lo), jack = rep(NA_real_, NJ))
    if (nm %in% JT) { io <- b$pairs$i[o]; jo <- b$pairs$j[o]
      out$jack <- unlist(parallel::mclapply(halves[[nm]], function(inh) { m <- inh[io] & inh[jo]; ap_sorted(so[m], lo[m]) },
                                            mc.cores = NC)) }
    out })
  say(sprintf("scored %-14s CORUM %.4f (x%.2f)  KD %.4f  DoRoAB %.4f  DoRoABC %.4f  CollecTRI %.4f  [%.0f s]", name,
              E[[name]]$CORUM$aupr, E[[name]]$CORUM$aupr / mean(bench$CORUM$pairs$label), E[[name]]$KD_response$aupr,
              E[[name]]$DoRothEA_AB$aupr, E[[name]]$DoRothEA_ABC$aupr, E[[name]]$CollecTRI$aupr,
              as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  saveRDS(E, file.path(OUT, "E.rds"))
}
summarise <- function() {
  rows <- list()
  for (m in names(E)) for (nm in names(bench)) {
    e <- E[[m]][[nm]]; br <- mean(bench[[nm]]$pairs$label)
    d <- function(ref) { if (is.null(E[[ref]])) return(c(NA, NA, NA)); r <- E[[ref]][[nm]]
      d0 <- e$aupr - r$aupr; s <- stats::sd(e$jack - r$jack); c(d0, d0 - 1.96 * s, d0 + 1.96 * s) }
    da <- d("att_t9"); dc <- d("clr")
    rows[[length(rows) + 1L]] <- data.frame(variant = m, truth = nm, aupr = e$aupr, base_rate = br, ratio = e$aupr / br,
      d_vs_att9 = da[1], att9_lo = da[2], att9_hi = da[3], d_vs_clr = dc[1], clr_lo = dc[2], clr_hi = dc[3])
  }
  res <- do.call(rbind, rows); utils::write.csv(res, file.path(OUT, "metrics_exploratory.csv"), row.names = FALSE); res
}

## ---- references: CLR, |r|, attention at t = 1..9 ----------------------------------------------
core_f <- file.path(OUT, "core.rds")
if (file.exists(core_f)) { core <- readRDS(core_f); C <- core$C; A <- core$A; rm(core); say("core reused")
} else {
  t1 <- Sys.time()
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = NC)$
    calibrate(method = "normal", combine = "stouffer")
  C <- f$clr_scores; dimnames(C) <- NULL
  f$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L)
  A <- f$operator; dimnames(A) <- NULL; f$release(); rm(f); invisible(gc())
  saveRDS(list(C = C, A = A), core_f, compress = FALSE)
  say(sprintf("MI, CLR, softmax operator in %.1f min", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
}
R <- stats::cor(t(X)); R[is.na(R)] <- 0
score("clr", C); score("abs_r", abs(R))
P <- 0.5 * A; diag(P) <- diag(P) + 0.5
St <- list(); Pk <- list(); Pt <- P
for (t in 1:9) {
  if (t > 1L) Pt <- Pt %*% P
  if (t %in% c(1, 2, 3, 4, 6, 9)) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; St[[as.character(t)]] <- S; rm(S) }
  if (t %in% c(1, 2, 9)) Pk[[as.character(t)]] <- Pt
}
rm(Pt); invisible(gc())
for (t in names(St)) score(paste0("att_t", t), St[[t]])
pi <- rep(1 / G, G); eps <- 1e-3
for (it in 1:5000) { q <- (1 - eps) * as.numeric(pi %*% P) + eps / G; if (sum(abs(q - pi)) < 1e-13) break; pi <- q }
say(sprintf("stationary pi: %d iterations; max/median %.0f", it, max(q) / stats::median(q))); pi <- q
res <- summarise()

## ---- diagnostics on CORUM --------------------------------------------------------------------
cm <- bench$CORUM$reg$members; cn <- names(cm); csz <- lengths(cm)
co <- comember_matrix(bench$CORUM$reg, G)           # shares any CORUM complex (diag TRUE)
isRP <- grepl("^RP[LS][0-9]+[A-Z]?[0-9]*$|^RPLP[0-2]$|^RPSA$", genes)
rp_frac <- vapply(cm, function(m) mean(isRP[m]), 0)
grp <- ifelse(rp_frac >= 0.3, "RP-dominated", "other")
say(sprintf("RP genes in data: %d; RP-dominated complexes (>= 30%% RP members): %d: %s", sum(isRP), sum(grp == "RP-dominated"),
            paste(sprintf("%s (%d)", cn[grp == "RP-dominated"], csz[grp == "RP-dominated"]), collapse = "; ")))
utils::write.csv(data.frame(complex = cn, size = csz, rp_frac = rp_frac, group = grp), file.path(OUT, "corum_complexes.csv"), row.names = FALSE)

# per-positive precision and the decomposition of the pooled AP difference by complex
pp <- bench$CORUM$pairs; npos <- sum(pp$label)
key <- (pp$i - 1) * G + pp$j
pos_prec <- function(S) { sc <- S[cbind(pp$i, pp$j)]; sc[!is.finite(sc)] <- 0
  o <- order(sc, decreasing = TRUE, method = "radix"); s <- sc[o]; l <- pp$label[o]; n <- length(s)
  brk <- c(TRUE, s[-1L] != s[-n]); blk <- cumsum(brk); last <- c(s[-1L] != s[-n], TRUE)
  prec <- (cumsum(l)[last] / which(last))[blk]; pr <- numeric(n); pr[o] <- prec
  pct <- rank(sc, ties.method = "average") / n                 # 1 = top of the pooled ranking (ties averaged)
  list(prec = pr, pct = pct) }
# complex membership of each positive pair (fractional weights)
cp <- do.call(rbind, lapply(seq_along(cm), function(k) { m <- cm[[k]]; if (length(m) < 2) return(NULL)
  u <- utils::combn(m, 2); data.frame(k = k, key = (pmin(u[1, ], u[2, ]) - 1) * G + pmax(u[1, ], u[2, ])) }))
cp$p <- match(cp$key, key); cp <- cp[!is.na(cp$p), ]
cp$w <- 1 / as.numeric(table(cp$p)[as.character(cp$p)])
say(sprintf("positive pairs %d; covered by complex pairs %.0f", npos, sum(cp$w)))
decomp <- function(a, b) {     # contribution of each complex to AP(a) - AP(b)
  da <- pos_prec(a)$prec - pos_prec(b)$prec
  v <- tapply(cp$w * da[cp$p], cp$k, sum) / npos; out <- numeric(length(cm)); out[as.integer(names(v))] <- v; out }
dec9 <- decomp(C, St[["9"]]); dec1 <- decomp(C, St[["1"]])
say(sprintf("pooled CORUM AP, CLR - att_t9 = %.4f; RP-dominated complexes carry %.0f%%; top-10 complexes carry %.0f%%",
            sum(dec9), 100 * sum(dec9[grp == "RP-dominated"]) / sum(dec9), 100 * sum(sort(dec9, decreasing = TRUE)[1:10]) / sum(dec9)))
say(sprintf("pooled CORUM AP, CLR - att_t1 = %.4f; RP-dominated complexes carry %.0f%%",
            sum(dec1), 100 * sum(dec1[grp == "RP-dominated"]) / sum(dec1)))
o10 <- order(-dec9)[1:15]
for (k in o10) say(sprintf("   %-60.60s size %3d  RP %.2f  share of loss t9 %5.1f%%  t1 %5.1f%%", cn[k], csz[k], rp_frac[k],
                           100 * dec9[k] / sum(dec9), 100 * dec1[k] / sum(dec1)))
# pooled percentile of positive pairs by group and size
pc <- lapply(list(clr = C, att_t1 = St[["1"]], att_t9 = St[["9"]]), function(S) pos_prec(S)$pct)
pgrp <- tapply(grp[cp$k] == "RP-dominated", cp$p, any)
psz <- tapply(csz[cp$k], cp$p, max)
pos_idx <- as.integer(names(pgrp))
ptab <- data.frame(p = pos_idx, rp = pgrp, size = psz, clr = pc$clr[pos_idx], t1 = pc$att_t1[pos_idx], t9 = pc$att_t9[pos_idx])
ptab$size_bin <- cut(ptab$size, c(2, 4, 9, 24, 49, 200), labels = c("3-4", "5-9", "10-24", "25-49", "50+"))
utils::write.csv(ptab, file.path(OUT, "corum_positive_percentiles.csv"), row.names = FALSE)
for (g in c(TRUE, FALSE)) { x <- ptab[ptab$rp == g, ]
  say(sprintf("positives %-13s n=%6d  median pooled percentile CLR %.4f  t1 %.4f  t9 %.4f; in top 1%%: CLR %.2f t1 %.2f t9 %.2f",
              if (g) "RP-dominated" else "other", nrow(x), stats::median(x$clr), stats::median(x$t1), stats::median(x$t9),
              mean(x$clr > 0.99), mean(x$t1 > 0.99), mean(x$t9 > 0.99))) }
for (b in levels(ptab$size_bin)) { x <- ptab[ptab$size_bin == b & !ptab$rp, ]; if (!nrow(x)) next
  say(sprintf("   non-RP positives, largest complex size %-5s n=%6d  median pct CLR %.4f t1 %.4f t9 %.4f", b, nrow(x),
              stats::median(x$clr), stats::median(x$t1), stats::median(x$t9))) }
# who occupies the top of each pooled ranking
topk <- function(S, k = 5000) { sc <- S[cbind(pp$i, pp$j)]; o <- order(sc, decreasing = TRUE)[1:k]
  c(precision = mean(pp$label[o]), rp_rp = mean(isRP[pp$i[o]] & isRP[pp$j[o]]), rp_any = mean(isRP[pp$i[o]] | isRP[pp$j[o]]),
    pos_rp = mean(pp$label[o] & isRP[pp$i[o]] & isRP[pp$j[o]]), hub = mean((pi[pp$i[o]] + pi[pp$j[o]]) > stats::quantile(pi, 0.99))) }
for (nm in c("clr", "1", "9")) { S <- if (nm == "clr") C else St[[nm]]; v <- topk(S)
  say(sprintf("top 5000 CORUM-universe pairs, %-6s: precision %.3f, RP-RP %.3f (positive RP-RP %.3f), any RP %.3f, a top-1%% pi gene %.3f",
              nm, v[1], v[2], v[4], v[3], v[5])) }

# mass of each member's walk inside its complex, in other co-complex genes, and outside (t = 1, 2, 9)
big <- which(csz >= 5)
mt <- do.call(rbind, lapply(big, function(k) { m <- cm[[k]]
  do.call(rbind, lapply(names(Pk), function(tt) { Pt <- Pk[[tt]]
    self <- Pt[cbind(m, m)]; inm <- rowSums(Pt[m, m, drop = FALSE]) - self
    coc <- rowSums(Pt[m, , drop = FALSE] * co[m, , drop = FALSE]) - self - inm
    data.frame(k = k, complex = cn[k], size = length(m), group = grp[k], t = as.integer(tt), gene = genes[m],
               self = self, in_complex = inm, other_cocomplex = coc, outside = 1 - self - inm - coc) })) }))
utils::write.csv(mt, file.path(OUT, "corum_member_mass.csv"), row.names = FALSE)
mt$sz <- cut(mt$size, c(4, 9, 24, 49, 200), labels = c("5-9", "10-24", "25-49", "50+"))
mt$nonself <- 1 - mt$self
for (tt in c(1, 2, 9)) for (g in c("RP-dominated", "other")) for (b in levels(mt$sz)) {
  x <- mt[mt$t == tt & mt$group == g & mt$sz == b, ]; if (!nrow(x)) next
  say(sprintf("mass t=%d %-12s size %-5s rows %4d | share of non-self mass: in complex %.3f, other co-complex %.3f, outside %.3f | per partner %.4f",
              tt, g, b, nrow(x), stats::median(x$in_complex / x$nonself), stats::median(x$other_cocomplex / x$nonself),
              stats::median(x$outside / x$nonself), stats::median(x$in_complex / (x$size - 1)))) }
# within-complex pair scores vs complex size: attention per pair falls with size, CLR does not
wp <- do.call(rbind, lapply(big, function(k) { m <- cm[[k]]; ut <- upper.tri(diag(length(m)))
  data.frame(k = k, size = length(m), group = grp[k], clr = C[m, m][ut], a1 = (A[m, m] + t(A[m, m]))[ut] / 2,
             s9 = St[["9"]][m, m][ut]) }))
wp$sz <- cut(wp$size, c(4, 9, 24, 49, 200), labels = c("5-9", "10-24", "25-49", "50+"))
for (b in levels(wp$sz)) { x <- wp[wp$sz == b, ]
  say(sprintf("within-complex pairs, size %-5s n=%6d: median CLR %.2f | softmax (A+A')/2 %.4f (zero %.0f%%) | S9 %.5f",
              b, nrow(x), stats::median(x$clr), stats::median(x$a1), 100 * mean(x$a1 == 0), stats::median(x$s9))) }
say(sprintf("Spearman(size, per-pair score) over within-complex pairs: CLR %.2f, softmax %.2f, S9 %.2f",
            stats::cor(wp$size, wp$clr, method = "spearman"), stats::cor(wp$size, wp$a1, method = "spearman"),
            stats::cor(wp$size, wp$s9, method = "spearman")))
# softmax temperature and plateau of near-equal CLR partners, RP vs other genes
tau <- vapply(seq_len(G), function(i) { a <- A[i, ]; s <- which(a > 0); if (length(s) < 3) return(NA_real_)
  j1 <- s[which.max(a[s])]; j0 <- s[which.min(a[s])]; if (C[i, j1] == C[i, j0]) return(NA_real_)
  (C[i, j1] - C[i, j0]) / (log(a[j1]) - log(a[j0])) }, 0)
plateau <- vapply(seq_len(G), function(i) { r <- C[i, ]; r[i] <- -Inf; sum(r >= max(r) - 1) }, 0)
cormem <- sort(unique(unlist(cm)))
say(sprintf("softmax temperature: RP genes median %.3f vs other CORUM genes %.3f vs all %.3f; CLR partners within 1 unit of the row max: RP %.0f, other CORUM %.0f, all %.0f",
            stats::median(tau[isRP], na.rm = TRUE), stats::median(tau[setdiff(cormem, which(isRP))], na.rm = TRUE),
            stats::median(tau, na.rm = TRUE), stats::median(plateau[isRP]), stats::median(plateau[setdiff(cormem, which(isRP))]),
            stats::median(plateau)))
# where do the RP-complex rows send their attention outside the complex? top recipients (t = 9)
ribo <- which(cn == "60S ribosomal subunit, cytoplasmic"); if (!length(ribo)) ribo <- which.max(rp_frac * csz)
m <- cm[[ribo[1]]]; recv <- colSums(Pk[["9"]][m, , drop = FALSE]) / length(m); recv[m] <- 0
ot <- order(-recv)[1:25]
say(sprintf("t=9 mass from %s members to non-members (top 25): %s", cn[ribo[1]],
            paste(sprintf("%s %.3f%s", genes[ot], recv[ot], ifelse(colSums(co[m, ot, drop = FALSE]) > 0, "*", "")), collapse = ", ")))
say("   (* = shares a CORUM complex with at least one member)")
recv1 <- colSums(A[m, , drop = FALSE]) / length(m); recv1[m] <- 0; ot1 <- order(-recv1)[1:15]
inward <- rowSums(A[, m, drop = FALSE]); inward[m] <- 0; oi <- order(-inward)[1:15]
say(sprintf("t=1 softmax mass from members to non-members (top 15): %s", paste(sprintf("%s %.3f", genes[ot1], recv1[ot1]), collapse = ", ")))
say(sprintf("non-members attending most to the complex (sum of A over members, top 15): %s",
            paste(sprintf("%s %.2f%s", genes[oi], inward[oi], ifelse(colSums(co[m, oi, drop = FALSE]) > 0, "*", "")), collapse = ", ")))

## ---- row-local vs pooled ranking ------------------------------------------------------------
# per-member row AP (complex partners vs genes sharing no complex), no comparison across rows
rowap <- function(S) do.call(rbind, lapply(big, function(k) { m <- cm[[k]]
  data.frame(k = k, size = length(m), group = grp[k], gene = m, ap = vapply(m, function(i) {
    cand <- which(!co[i, ] | (seq_len(G) %in% m)); cand <- cand[cand != i]
    sc <- S[i, cand]; lab <- as.integer(cand %in% m); ap_fast(sc, lab) }, 0)) }))
# per-complex AP exactly as the macro view (regulon_ap_table: member pairs vs member-to-non-sharing pairs)
macro <- function(S) { x <- regulon_ap_table(S, list(members = cm[big]), co, none, genes); x$k <- big; x }
RA <- list(); MA <- list()
loc <- function(name, S) { RA[[name]] <<- rowap(S); MA[[name]] <<- macro(S)
  r <- RA[[name]]; mm <- MA[[name]]
  say(sprintf("local %-14s row AP: RP-dominated %.3f, other %.3f | macro AP: 60S %.3f, Ribosome cyto %.3f; RP-dominated mean %.3f, other mean %.3f",
              name, mean(r$ap[r$group == "RP-dominated"]), mean(r$ap[r$group == "other"]),
              mm$ap[match("60S ribosomal subunit, cytoplasmic", mm$regulon)], mm$ap[match("Ribosome, cytoplasmic", mm$regulon)],
              mean(mm$ap[grp[mm$k] == "RP-dominated"], na.rm = TRUE), mean(mm$ap[grp[mm$k] == "other"], na.rm = TRUE))) }
loc("clr", C); loc("att_t1", St[["1"]]); loc("att_t2", St[["2"]]); loc("att_t9", St[["9"]])
saveRDS(list(RA = RA, MA = MA), file.path(OUT, "local.rds"))

## ---- fix (a): normalization after attention (CLR -> attention -> CLR) -------------------------
clr_R <- function(S) { diag(S) <- NA; mu <- rowMeans(S, na.rm = TRUE)
  sdv <- sqrt(rowSums((S - mu)^2, na.rm = TRUE) / (ncol(S) - 2)); Z <- (S - mu) / sdv; Z[is.na(Z) | Z < 0] <- 0
  out <- (Z + t(Z)) / sqrt(2); diag(out) <- 0; out }
chk <- max(abs(clr_R(C) - clr::clr_calibrate(C)))       # sanity: same recipe as the package (up to the diagonal)
say(sprintf("clr_R vs clr_calibrate on a test matrix: max |diff| %.3g (diagonal handling only)", chk))
for (t in c("1", "2", "9")) { S <- clr_R(St[[t]]); score(paste0("a_clr_on_att_t", t), S); if (t == "9") loc("a_clr_on_att_t9", S); rm(S) }
L9 <- log(pmax(Pk[["9"]], 1e-300))
S <- clr_R((L9 + t(L9)) / 2); score("a_clr_on_logatt_t9", S); loc("a_clr_on_logatt_t9", S); rm(S)
## ---- fix (d): readouts without the mass penalty -------------------------------------------------
S <- (L9 + t(L9)) / 2; diag(S) <- 0
score("d_geomean_t9", S); rm(S)                              # mutual attention sqrt(P_ij P_ji), log scale
lp <- log(pi); S <- (L9 + t(L9)) / 2 - outer(lp, lp, "+") / 2   # diffusion density p_t(i,j)/pi_j, symmetrized (log)
score("d_density_t9", S); loc("d_density_t9", S); rm(S, L9); invisible(gc())
res <- summarise()

## ---- fix (b): latent-variable coarse-graining (collapse cliques, attend, expand) ---------------
hc <- stats::hclust(stats::as.dist(1 - R), method = "average")
Zx <- t(scale(t(X)))
coarse <- function(r0) {
  t1 <- Sys.time(); u <- stats::cutree(hc, h = 1 - r0); usz <- tabulate(u)
  Xc <- t(vapply(seq_len(max(u)), function(k) colMeans(Zx[u == k, , drop = FALSE]), numeric(ncol(X))))
  fc <- ClrAttention$new(Xc)$estimate_mi(bins = "hg", transform = "none", threads = NC)$
    calibrate(method = "normal", combine = "stouffer")
  fc$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); Ac <- fc$operator; dimnames(Ac) <- NULL
  fc$release(); rm(fc)
  Pc <- 0.5 * Ac; diag(Pc) <- diag(Pc) + 0.5; Q <- Pc; for (s in 2:9) Q <- Q %*% Pc
  Sc <- (Q + t(Q)) / 2; diag(Sc) <- 0
  Sx <- Sc[u, u]; within <- outer(u, u, "=="); diag(within) <- FALSE
  top <- max(Sx[!within]); Sx[within] <- 2 * top + C[within] / max(C) * top   # same latent variable: ranked first, by CLR
  diag(Sx) <- 0
  say(sprintf("coarse r0=%.1f: %d units (%d multi-gene, largest %d, %d genes in multi-gene units, %d within-unit pairs) in %.1f min",
              r0, max(u), sum(usz > 1), max(usz), sum(usz[u] > 1), sum(within) / 2,
              as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  wl <- within[cbind(pp$i, pp$j)]
  say(sprintf("   within-unit CORUM-universe pairs: %d, precision %.3f (CORUM base rate %.4f); RP genes in units: %d of %d",
              sum(wl), mean(pp$label[wl]), mean(pp$label), sum(isRP & usz[u] > 1), sum(isRP)))
  list(S = Sx, within = within)
}
for (r0 in c(0.6, 0.4, 0.8)) {
  cg <- coarse(r0); tag <- sprintf("%02d", round(10 * r0))
  score(paste0("b_coarse_r", tag), cg$S); if (r0 == 0.6) loc("b_coarse_r06", cg$S)
  if (r0 == 0.6) {
    S <- St[["9"]]; top <- max(S[!cg$within]); S[cg$within] <- 2 * top + C[cg$within] / max(C) * top
    score("b0_boost_on_att_t9", S); loc("b0_boost_on_att_t9", S)
    S <- C; top <- max(S[!cg$within]); S[cg$within] <- 2 * top + C[cg$within]
    score("bC_boost_on_clr", S); rm(S)
  }
  rm(cg); invisible(gc())
  if (as.numeric(difftime(Sys.time(), T0, units = "mins")) > 38) { say("time budget: stopping sensitivity runs"); break }
}
res <- summarise()
saveRDS(list(RA = RA, MA = MA), file.path(OUT, "local.rds"))
say(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
