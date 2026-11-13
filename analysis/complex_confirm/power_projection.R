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
# anchors' per-complex differences in the same size class (emp). Ranges over anchors are reported as an envelope
# only (the two surrogates disagree in direction across anchors and are not calibrated).
# The precision numbers used in the design come from the k rule at the end (power_summary.csv).
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
  out <- data.frame(complex = a$regulon, n = a$size, npos = a$n_pos, d = d, dc = dc)
  out[out$npos > 0, ]   # rev 2 (review M3): E. coli complexes whose pairs are all same-operon have no positives; dropped
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
utils::write.csv(P, file.path(OUT, "power_projection.csv"), row.names = FALSE)
print(P, row.names = FALSE, digits = 3)
S <- do.call(rbind, lapply(split(P, paste(P$dataset, P$variant)), function(x) data.frame(dataset = x$dataset[1], variant = x$variant[1],
  surrogate_se_min = min(x$se_iid, x$se_emp), surrogate_se_max = max(x$se_iid, x$se_emp))))
utils::write.csv(S, file.path(OUT, "power_surrogate_envelope.csv"), row.names = FALSE)
print(S, row.names = FALSE, digits = 3)

## ---- k rule and precision (rev 2; review M1, M3, R6) ----
# k = SE_J * sqrt(n_eff) per anchor (n_eff overlap-counted, as per_complex_ap.csv); the anchors are about two
# independent pairs (M3D/EcoMAC share a truth, K562/RPE1 share CORUM). Two-point power law through the pair means.
# Experiment-to-truth variance ratios are net4's own, from the DREAM5 Y-M test on the same matrix (scale-free):
# group level SE_B / SE_J and series level SD_series / SE_J (results/dream5/truth_net4/).
cnt <- utils::read.csv(file.path(OUT, "yeast_complexes_counts.csv"))
ym <- utils::read.csv("results/dream5/truth_net4/verdict.csv"); ym <- ym[ym$endpoint == "YM", ]
ss <- utils::read.csv("results/dream5/truth_net4/supplementary_series_jackknife.csv"); ss <- ss[ss$endpoint == "YM", ]
r_group <- ym$SE_B / ym$SE_J; r_series <- ss$SD_series / ym$SE_J
k <- A$se_obs_vs_att * sqrt(A$n_eff); kc <- A$se_obs_vs_clr * sqrt(A$n_eff)
ne_ec <- mean(A$n_eff[1:2]); ne_hu <- mean(A$n_eff[3:4]); se_ec <- mean(A$se_obs_vs_att[1:2]); se_hu <- mean(A$se_obs_vs_att[3:4])
expo <- log(se_ec / se_hu) / log(ne_hu / ne_ec)
REF <- c(min_effect = 0.005, ecoli_mean = mean(c(0.0256587849922917, 0.0285030845859859)))
pw <- function(delta, se) 1 - stats::pnorm(stats::qnorm(0.975) - delta / se)
K <- list()
for (v in c("main", "main_xsub_excl", "main_ribo_joined", "noribo")) {
  ne <- cnt$n_eff_overlap_counted[cnt$dataset == "dream5_net4" & cnt$variant == v]
  se_j <- mean(k) / sqrt(ne); se_j_lo <- min(k) / sqrt(ne); se_j_hi <- max(k) / sqrt(ne)
  se_2pt <- se_hu * (ne_hu / ne)^expo
  for (sc in c("gene_only", "total_group", "total_series")) {
    f <- switch(sc, gene_only = 1, total_group = sqrt(1 + r_group^2), total_series = sqrt(1 + r_series^2))
    K[[length(K) + 1]] <- data.frame(variant = v, n_eff_overlap_counted = ne, interval = sc, ratio_used = switch(sc, gene_only = 0, total_group = r_group, total_series = r_series),
      se = se_j * f, se_lo = se_j_lo * f, se_hi = se_j_hi * f, se_two_point_law = se_2pt * f, half_width_95 = 1.96 * se_j * f,
      power_at_0.010 = pw(0.010, se_j * f), power_at_0.0135 = pw(0.0135, se_j * f), power_at_0.027 = pw(REF[["ecoli_mean"]], se_j * f),
      se_vs_clr = mean(kc) / sqrt(ne) * f)
  }
}
K <- do.call(rbind, K)
attr_rows <- data.frame(quantity = c(paste0("k_", A$anchor), paste0("kCLR_", A$anchor), "two_point_exponent", "ratio_group_net4_YM", "ratio_series_net4_YM", "ref_ecoli_mean"),
  value = c(k, kc, expo, r_group, r_series, REF[["ecoli_mean"]]))
utils::write.csv(K, file.path(OUT, "power_summary.csv"), row.names = FALSE)
utils::write.csv(attr_rows, file.path(OUT, "power_k_rule_inputs.csv"), row.names = FALSE)
print(attr_rows, row.names = FALSE, digits = 4); print(K, row.names = FALSE, digits = 3)
