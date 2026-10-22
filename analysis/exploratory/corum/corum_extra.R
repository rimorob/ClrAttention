#!/usr/bin/env Rscript
# EXPLORATORY, NOT PREREGISTERED. Follow-up to corum_diag.R: (1) CORUM without the
# RP-dominated complexes; (2) the log-domain renormalization at other depths and in
# its bilateral (per-row, then Stouffer) form. Reuses core.rds (CLR, softmax operator).
setwd("/home/boris/Projects/ClrAttention")
SP <- "/tmp/claude-1000/-home-boris-Projects-ClrAttention/b33914c9-0cfd-485d-8dd4-cf8a2a22b632/scratchpad"
OUT <- file.path(SP, "corum"); logf <- file.path(OUT, "extra.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R")
NC <- 4L; NJ <- 100L; T0 <- Sys.time()
RUN <- "results/baselines_k562ess_truthgenes_genie3"
genes <- utils::read.delim(file.path(RUN, "expression.tsv"), colClasses = c("character", rep("NULL", 707)))[[1]]
G <- length(genes)
rmt <- utils::read.csv(file.path(RUN, "regulon_members.csv"), stringsAsFactors = FALSE)
cls <- c(DoRothEA_AB = "TF", DoRothEA_ABC = "TF", CollecTRI = "TF", CORUM = "complex", KD_response = "response")
tr <- lapply(names(cls), function(nm) { x <- rmt[rmt$truth == nm, ]
  m <- lapply(strsplit(x$members, ";", fixed = TRUE), function(g) sort(match(g, genes))); names(m) <- x$regulon
  list(members = m, class = stats::setNames(rep(cls[[nm]], length(m)), names(m))) })
names(tr) <- names(cls)
cc <- utils::read.csv(file.path(OUT, "corum_complexes.csv"), stringsAsFactors = FALSE)
keep <- !(names(tr$CORUM$members) %in% cc$complex[cc$group == "RP-dominated"])
tr$CORUM_noRP <- list(members = tr$CORUM$members[keep], class = tr$CORUM$class[keep])
say(sprintf("CORUM_noRP drops %d RP-dominated complexes", sum(!keep)))
none <- stats::setNames(character(0), character(0))
bench <- lapply(tr, function(r) list(reg = r, pairs = comembership_pairs(r, none, genes)))
for (nm in names(bench)) say(sprintf("truth %-12s %4d sets, %8d pairs, base rate %.4f", nm, length(bench[[nm]]$reg$members),
                                     length(bench[[nm]]$pairs$label), mean(bench[[nm]]$pairs$label)))
set.seed(20261016L)        # same splits as baselines.R for the five original truths; CORUM_noRP drawn after them
halves <- lapply(bench, function(b) lapply(seq_len(NJ), function(r) {
  keep <- sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2)); inh <- logical(G); inh[keep] <- TRUE; inh }))
ap_sorted <- function(s, l) { n <- length(s); if (!n || !sum(l)) return(NA_real_)
  last <- c(s[-1L] != s[-n], TRUE); tp <- cumsum(l)[last]; k <- which(last)
  sum(diff(c(0, tp)) * tp / k) / sum(l) }
JT <- c("DoRothEA_AB", "CollecTRI", "CORUM", "KD_response", "CORUM_noRP")
E <- list()
score <- function(name, S, truths = names(bench)) {
  t1 <- Sys.time()
  E[[name]] <<- lapply(stats::setNames(truths, truths), function(nm) {
    b <- bench[[nm]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]; sc[!is.finite(sc)] <- 0
    o <- order(sc, decreasing = TRUE, method = "radix"); so <- sc[o]; lo <- b$pairs$label[o]
    out <- list(aupr = ap_sorted(so, lo), jack = rep(NA_real_, NJ))
    if (nm %in% JT) { io <- b$pairs$i[o]; jo <- b$pairs$j[o]
      out$jack <- unlist(parallel::mclapply(halves[[nm]], function(inh) { m <- inh[io] & inh[jo]; ap_sorted(so[m], lo[m]) },
                                            mc.cores = NC)) }
    out })
  say(sprintf("scored %-22s %s [%.0f s]", name, paste(sprintf("%s %.4f", truths, vapply(E[[name]], `[[`, 0, "aupr")), collapse = "  "),
              as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  saveRDS(E, file.path(OUT, "E_extra.rds"))
}
core <- readRDS(file.path(OUT, "core.rds")); A <- core$A; rm(core)
rowz <- function(S) { diag(S) <- NA; mu <- rowMeans(S, na.rm = TRUE)
  sdv <- sqrt(rowSums((S - mu)^2, na.rm = TRUE) / (ncol(S) - 2)); Z <- (S - mu) / sdv; Z[is.na(Z) | Z < 0] <- 0; Z }
P <- 0.5 * A; diag(P) <- diag(P) + 0.5; Pt <- P; rm(A)
for (t in 1:9) {
  if (t > 1L) Pt <- Pt %*% P
  if (t %in% c(2, 4, 9)) { z <- Pt < 1e-300; diag(z) <- FALSE
    say(sprintf("t=%d: off-diagonal entries below the 1e-300 floor %.2f%%; rows with any %.1f%%; median row min log P %.1f",
                t, 100 * mean(z), 100 * mean(rowSums(z) > 0), stats::median(apply(log(pmax(Pt, 1e-300)) + diag(Inf, G), 1, min)))) }
}
eps <- 1e-3; q <- solve(diag(G) - (1 - eps) * t(as.matrix(P)), rep(eps / G, G)); q <- q / sum(q); H <- outer(q, q, "+")
hub_stratify <- function(S, H, k = 10L) { S[!is.finite(S)] <- 0; ut <- which(upper.tri(S)); h <- H[ut]
  g <- cut(h, unique(stats::quantile(h, seq(0, 1, length.out = k + 1L))), include.lowest = TRUE)
  r <- stats::ave(S[ut], g, FUN = function(x) rank(x) / length(x))
  R <- matrix(0, nrow(S), ncol(S)); R[ut] <- r; R + t(R) }
Z <- rowz(log(pmax(Pt, 1e-300))); S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; rm(Z, Pt); invisible(gc())
score("a_clr_bilateral_logP_t9_hubstrat", hub_stratify(S, H))
Ef <- readRDS(file.path(OUT, "E_followup.rds"))
rows <- list(); m <- "a_clr_bilateral_logP_t9_hubstrat"
for (nm in names(E[[m]])) { e <- E[[m]][[nm]]; br <- mean(bench[[nm]]$pairs$label)
  d <- function(ref) { r <- Ef[[ref]][[nm]]; d0 <- e$aupr - r$aupr; s <- stats::sd(e$jack - r$jack); c(d0, d0 - 1.96 * s, d0 + 1.96 * s) }
  da <- d("att_t9_hubstrat"); dc <- d("clr_hubstrat")
  rows[[length(rows) + 1L]] <- data.frame(variant = m, truth = nm, aupr = e$aupr, base_rate = br, ratio = e$aupr / br,
    d_vs_att9 = da[1], att9_lo = da[2], att9_hi = da[3], d_vs_clr = dc[1], clr_lo = dc[2], clr_hi = dc[3]) }
utils::write.csv(do.call(rbind, rows), file.path(OUT, "metrics_extra.csv"), row.names = FALSE)
say(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
