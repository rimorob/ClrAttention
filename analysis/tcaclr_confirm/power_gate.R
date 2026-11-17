#!/usr/bin/env Rscript
# Power gate arithmetic of docs/plans/tcaclr_confirm_plan.md section 6 gate 2, from the counted n_eff only (no score).
# SE_gene = k / sqrt(n_eff) with k = 0.0417 (net4, docs/plans/tcaclr_prereg_candidates.md C2), SE_total = SE_gene sqrt(1 + rho^2), rho in {0, 0.893, 1.78}.
# Power of the one-sided 0.025 test = Phi(d / SE - 1.96). Gate (Arm A): n_eff >= 27 -> registered sensitivity +0.03; 15 <= n_eff < 27 -> sensitivity stated as +0.04;
# n_eff < 15 -> estimate only.
# Usage: Rscript power_gate.R <n_eff arm A> <n_eff arm B> <outfile.csv>
a <- commandArgs(trailingOnly = TRUE); ne <- c(A = as.numeric(a[1]), B = as.numeric(a[2])); k <- 0.0417
rows <- do.call(rbind, lapply(names(ne), function(arm) do.call(rbind, lapply(c(0, 0.893, 1.78), function(rho) { se <- k * sqrt(1 + rho^2) / sqrt(ne[[arm]])
  data.frame(arm = arm, n_eff = round(ne[[arm]], 2), rho = rho, se = round(se, 4), mde_80 = round(se * (1.96 + 0.8416), 3),
             power_at_0.03 = round(stats::pnorm(0.03 / se - 1.96), 2), power_at_0.05 = round(stats::pnorm(0.05 / se - 1.96), 2)) }))))
rows$gate <- ifelse(rows$arm == "A", ifelse(ne[["A"]] >= 27, "n_eff >= 27: test at +0.03", ifelse(ne[["A"]] >= 15, "15 <= n_eff < 27: test at +0.04", "n_eff < 15: estimate only")), "secondary estimate, no p-value")
utils::write.csv(rows, a[3], row.names = FALSE); print(rows, row.names = FALSE)
