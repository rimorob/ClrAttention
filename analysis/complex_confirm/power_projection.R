#!/usr/bin/env Rscript
# Power projection for docs/plans/complex_confirm_design.md. USES ONLY EXISTING OUTPUTS of spent runs
# (results/clr_logpt_test/truth/{m3d,ecomac,k562ess,rpe1}/: jack_replicates.csv, per_complex_ap.csv) and the
# truth STRUCTURE of the candidate yeast truth (docs/complex_confirm_selections/yeast_complexes_sizes.csv).
# Computes no score of any readout on any candidate truth; reads no expression data and no network.
#
# Model. The paired delete-half gene jackknife SD of the pooled AUPR difference d is driven by which complexes
# (and how many of their pairs) a half keeps. Surrogate: d_b ~ sum_k w_k(b) delta_k, where w_k(b) is complex
# k's share of the positive pairs kept in half b (m_k ~ Binomial(n_k, 1/2) genes kept, m_k(m_k - 1)/2 pairs) and
# delta_k is a per-complex effect. Two versions:
#   (iid)  delta_k iid with SD tau; g = E SD_b(sum w_k(b) delta_k) for tau = 1; tau_anchor = SE_obs / g_anchor.
#   (emp)  delta_k = the observed per-complex AP difference (a3r - attention) of the anchor itself;
#          c_anchor = SE_obs / SD_b(...)  (scale factor between per-complex AP and pooled-AUPR units).
# Projection for the yeast truth: SE = tau x g_yeast (iid), and SE = c x SD_b with delta_k resampled from the
# anchors' per-complex differences in the same size class (emp). Ranges over anchors are reported.
# Usage: Rscript analysis/complex_confirm/power_projection.R
set.seed(20261113L)
TR <- "results/clr_logpt_test/truth"; OUT <- "docs/complex_confirm_selections"
anchors <- c(m3d = "complexes", ecomac = "complexes", k562ess = "CORUM_merged", rpe1 = "CORUM_merged")
B <- 400L; ND <- 200L

obs_se <- function(ds, truth, a = "a3r_t9", b = "attention_t9") {
  j <- utils::read.csv(file.path(TR, ds, "jack_replicates.csv")); j <- j[j$truth == truth, ]
  x <- j$aupr[j$method == a][order(j$rep[j$method == a])]; y <- j$aupr[j$method == b][order(j$rep[j$method == b])]
  stats::sd(x - y)
}
pc_of <- function(ds) {
  p <- utils::read.csv(file.path(TR, ds, "per_complex_ap.csv"))
  a <- p[p$method == "a3r_t9", ]; t <- p[p$method == "attention_t9", ]; cl <- p[p$method == "clr", ]
  t <- t[match(a$regulon, t$regulon), ]; cl <- cl[match(a$regulon, cl$regulon), ]
  d <- a$ap - t$ap; d[!is.finite(d)] <- 0; dc <- a$ap - cl$ap; dc[!is.finite(dc)] <- 0
  data.frame(complex = a$regulon, n = a$size, npos = a$n_pos, d = d, dc = dc)
}
# SD over B halves of sum_k w_k(b) delta_k; delta matrix K x ND (each column one draw of effects)
half_sd <- function(n, delta) {
  K <- length(n); delta <- as.matrix(delta)
  W <- vapply(seq_len(B), function(b) { m <- stats::rbinom(K, n, 0.5); p <- m * (m - 1) / 2; p / sum(p) }, numeric(K))
  S <- crossprod(W, delta)                                   # B x ND
  mean(apply(S, 2, stats::sd))
}
size_class <- function(n) cut(n, c(2, 4, 9, 24, 49, Inf), labels = c("3-4", "5-9", "10-24", "25-49", ">=50"))

rows <- list(); pool <- list()
for (ds in names(anchors)) {
  pc <- pc_of(ds); se <- obs_se(ds, anchors[[ds]]); se_c <- obs_se(ds, anchors[[ds]], "a3r_t9", "clr")
  g <- half_sd(pc$n, matrix(stats::rnorm(nrow(pc) * ND), nrow(pc)))
  sd_emp <- half_sd(pc$n, pc$d); sd_emp_c <- half_sd(pc$n, pc$dc)
  w <- pc$npos / sum(pc$npos)
  rows[[ds]] <- data.frame(anchor = ds, truth = anchors[[ds]], n_complexes = nrow(pc), n_pos_sum = sum(pc$npos),
    n_eff = 1 / sum(w^2), largest_share = max(w), se_obs_vs_att = se, se_obs_vs_clr = se_c,
    g_iid = g, tau = se / g, sd_emp = sd_emp, c_emp = se / sd_emp, c_emp_clr = se_c / sd_emp_c,
    sd_perc_d = stats::sd(pc$d), sd_perc_dc = stats::sd(pc$dc))
  pc$cls <- size_class(pc$n); pc$anchor <- ds; pool[[ds]] <- pc
}
A <- do.call(rbind, rows); pool <- do.call(rbind, pool)
utils::write.csv(A, file.path(OUT, "power_anchors.csv"), row.names = FALSE)
print(A, row.names = FALSE, digits = 4)

## ---- projection for each candidate yeast truth variant ----
sz <- utils::read.csv(file.path(OUT, "yeast_complexes_sizes.csv"))
proj <- list()
for (ds in unique(sz$dataset)) for (v in c("main", "noribo")) {
  y <- sz[sz$dataset == ds, ]; if (v == "noribo") y <- y[!y$ribo_dominated, ]
  n <- y$n_genes; w <- y$pos_weight / sum(y$pos_weight)
  g <- half_sd(n, matrix(stats::rnorm(length(n) * ND), length(n)))
  cls <- size_class(n)
  # emp: per-complex differences resampled within size class from the anchor pool (all anchors), ND draws
  draw_delta <- function(col) vapply(seq_len(ND), function(r) vapply(seq_along(n), function(k) {
    cand <- pool[[col]][pool$cls == cls[k]]; if (!length(cand)) 0 else sample(cand, 1) }, 0), numeric(length(n)))
  sd_e <- half_sd(n, draw_delta("d")); sd_ec <- half_sd(n, draw_delta("dc"))
  for (a in seq_len(nrow(A))) proj[[length(proj) + 1]] <- data.frame(dataset = ds, variant = v,
    n_complexes = length(n), n_eff = 1 / sum(w^2), largest_share = max(w), anchor = A$anchor[a],
    se_iid = A$tau[a] * g, se_emp = A$c_emp[a] * sd_e, se_emp_vs_clr = A$c_emp_clr[a] * sd_ec, g_iid = g)
}
P <- do.call(rbind, proj)
pw <- function(delta, se) 1 - stats::pnorm(stats::qnorm(0.975) - delta / se)
P$mde80 <- (stats::qnorm(0.975) + stats::qnorm(0.80)) * pmax(P$se_iid, P$se_emp)
P$power_0.025_iid <- pw(0.025, P$se_iid); P$power_0.025_emp <- pw(0.025, P$se_emp)
utils::write.csv(P, file.path(OUT, "power_projection.csv"), row.names = FALSE)
print(P, row.names = FALSE, digits = 3)
S <- do.call(rbind, lapply(split(P, paste(P$dataset, P$variant)), function(x) data.frame(dataset = x$dataset[1], variant = x$variant[1],
  se_min = min(x$se_iid, x$se_emp), se_median = stats::median(c(x$se_iid, x$se_emp)), se_max = max(x$se_iid, x$se_emp),
  se_vs_clr_median = stats::median(x$se_emp_vs_clr))))
S$mde80_median <- (stats::qnorm(0.975) + stats::qnorm(0.8)) * S$se_median
S$mde80_max <- (stats::qnorm(0.975) + stats::qnorm(0.8)) * S$se_max
for (eff in c(0.015, 0.025, 0.035, 0.05)) { S[[sprintf("power_%.3f_median", eff)]] <- pw(eff, S$se_median); S[[sprintf("power_%.3f_worst", eff)]] <- pw(eff, S$se_max) }
utils::write.csv(S, file.path(OUT, "power_summary.csv"), row.names = FALSE)
print(S, row.names = FALSE, digits = 3)
