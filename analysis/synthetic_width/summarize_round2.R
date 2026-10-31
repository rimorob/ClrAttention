#!/usr/bin/env Rscript
# Summaries for sim_round2.R and the mechanical application of the carry-forward
# rule fixed in docs/plans/background_width_exploration.md, section 8.1.
D <- "analysis/synthetic_width/out/round2"
rd <- function(p) do.call(rbind, lapply(list.files(file.path(D, "runs"), p, full.names = TRUE), utils::read.csv))
m <- rd("_metrics.csv$"); g <- rd("_gates.csv$"); o <- rd("_open.csv$")
mets <- c("cpx_x", "small_x", "s5_x", "mid_x", "s20_x", "s80_x", "reg_x", "edge_x")
cell <- function(d) paste(d$gen, sprintf("f%02d", round(100 * d$f)), paste0("b", d$b))
m$cell <- cell(m)
ratios <- do.call(rbind, lapply(split(m, m$cell), function(d) {
  ref <- d[d$op == "K10|att_t9", ]
  do.call(rbind, lapply(split(d, d$op), function(e) { e <- e[match(ref$seed, e$seed), ]
    q <- sapply(mets, function(v) e[[v]] / ref[[v]])
    data.frame(cell = d$cell[1], gen = d$gen[1], f = d$f[1], b = d$b[1], op = e$op[1], n = nrow(e),
      t(sapply(mets, function(v) exp(mean(log(q[, v]))))),
      t(stats::setNames(sapply(mets, function(v) sum(q[, v] > 1)), paste0(mets, "_nb")))) })) }))
utils::write.csv(ratios, file.path(D, "ratios_vs_K10.csv"), row.names = FALSE)
# carry-forward rule (attention-readout operators only)
att_ops <- unique(ratios$op[grepl("\\|att_t9$", ratios$op) & ratios$op != "K10|att_t9"])
rule <- do.call(rbind, lapply(att_ops, function(op) {
  s <- ratios[ratios$op == op & ratios$gen == "std" & ratios$f == 0.10, ]
  r <- ratios[ratios$op == op & ratios$gen == "regheavy", ]
  data.frame(op = op,
    C1_small = all(s$small_x > 1 & s$small_x_nb >= 3) && nrow(s) == 3,
    C2_all = all(s$cpx_x > 1 & s$cpx_x_nb >= 3) && nrow(s) == 3,
    C3_reg = all(c(s$reg_x, r$reg_x) >= 0.95) && nrow(r) == 2,
    C4_edge = all(r$edge_x >= 0.95) && nrow(r) == 2,
    small_b0_b1_b2 = paste(sprintf("%.2f(%d)", s$small_x, s$small_x_nb), collapse = " "),
    all_b0_b1_b2 = paste(sprintf("%.2f(%d)", s$cpx_x, s$cpx_x_nb), collapse = " "),
    reg_b0_b1_b2 = paste(sprintf("%.2f", s$reg_x), collapse = " "),
    regheavy_reg_b0_b1 = paste(sprintf("%.2f", r$reg_x), collapse = " "),
    regheavy_edge_b0_b1 = paste(sprintf("%.2f", r$edge_x), collapse = " ")) }))
rule$carry_forward <- rule$C1_small & rule$C2_all & rule$C3_reg & rule$C4_edge
utils::write.csv(rule, file.path(D, "carry_forward.csv"), row.names = FALSE)
# cosine falsifier: cosine vs a3r on the same operator, f = 10%, b = 0 and 1, small and regulons
cf <- do.call(rbind, lapply(c("GR_RECIP", "GE_RECIP"), function(op) {
  d <- m[m$gen == "std" & m$f == 0.10 & m$b %in% c(0, 1) & m$op %in% paste0(op, c("|cos_t9", "|a3r_t9")), ]
  do.call(rbind, lapply(split(d, d$b), function(e) { cs <- e[grepl("cos", e$op), ]; a3 <- e[grepl("a3r", e$op), ]
    a3 <- a3[match(cs$seed, a3$seed), ]
    data.frame(op = op, b = e$b[1], small_cos_over_a3r = exp(mean(log(cs$small_x / a3$small_x))),
               reg_cos_over_a3r = exp(mean(log(cs$reg_x / a3$reg_x)))) })) }))
cf$cos_survives <- with(cf, ave(small_cos_over_a3r > 1 & reg_cos_over_a3r > 1, op, FUN = all))
utils::write.csv(cf, file.path(D, "cosine_falsifier.csv"), row.names = FALSE)
g$cell <- cell(g)
ga <- stats::aggregate(cbind(f_real, auc_repro, auc_efronE, auc_efron_lfdr, spearman_E_repro, kappa) ~ cell, data = g, FUN = mean)
utils::write.csv(ga, file.path(D, "gate_detection.csv"), row.names = FALSE)
o$cell <- cell(o)
oa <- stats::aggregate(cbind(open_GR, open_GE, E_median, repro_median) ~ cell + role, data = o, FUN = mean)
utils::write.csv(oa, file.path(D, "gate_open_by_role.csv"), row.names = FALSE)
options(width = 250)
print(rule, row.names = FALSE); print(cf, row.names = FALSE); print(ga, row.names = FALSE, digits = 3)
