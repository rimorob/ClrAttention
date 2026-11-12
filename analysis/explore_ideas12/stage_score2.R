#!/usr/bin/env Rscript
# EXPLORATORY. Step 5 of docs/plans/explore_ideas12_plan.md: score the idea-2 gated operators (K562, in-sample). Exits cleanly if the
# label-free diagnostic stopped the idea.
args <- commandArgs(TRUE); workers <- if (length(args)) as.integer(args[1]) else 6L
source("analysis/explore_ideas12/ideas12_lib.R"); source("analysis/explore_ideas12/ideas12_tables.R")
out <- file.path(EI_OUT, "k562"); cd <- file.path(out, "cache"); say <- ei_say(file.path(out, "stage_score2.log"))
flag <- readLines(file.path(out, "i2_flag.txt"))[1]; say("diagnostic flag: ", flag)
if (flag == "stop") { say("idea 2 STOPPED by the diagnostic (plan rule); nothing scored"); quit(status = 0) }
core <- readRDS(file.path(cd, "k562_core.rds")); G <- length(core$genes)
q <- core$q_att; H_att <- outer(q, q, "+"); H_mi <- outer(core$mi_rm, core$mi_rm, "+")
A <- as.matrix(core$A); rm(core)
env <- ei_setup_k562(say, workers); stopifnot(env$G == G)
ev <- ei_make_eval(env, H_att, H_mi, say); rm(H_att, H_mi)
rh <- readRDS(file.path(cd, "i2_rho.rds")); hh <- readRDS(file.path(cd, "halves_top20.rds"))
## half-split detector: same top-20 hypergeometric rule between the two halves of each of the 6 splits, soft weights averaged over splits
rho_s <- lapply(1:6, function(si) ei_rho(hh$top20[[si]][[1]], hh$top20[[si]][[2]], G))
rho_half <- rowMeans(do.call(cbind, lapply(rho_s, function(r) r$rho)))
dir.create(file.path(out, "tables"), showWarnings = FALSE); dg <- utils::read.csv(file.path(out, "i2_diagnostic.csv"))
dh <- data.frame(detector = "half_split", n_profiles_Y = NA, n_universe = G, n_outside = 0, frac_rho_lt_0.5 = mean(vapply(rho_s, function(r) mean(r$rho < 0.5), 0)),
  frac_rho_lt_0.5_all_genes = mean(rho_half < 0.5), frac_overlap0 = mean(vapply(rho_s, function(r) mean(r$overlap == 0), 0)), frac_p_le_0.05 = mean(vapply(rho_s, function(r) mean(r$p <= 0.05), 0)),
  median_overlap = stats::median(unlist(lapply(rho_s, function(r) r$overlap))), mean_overlap = mean(unlist(lapply(rho_s, function(r) r$overlap))), chance_overlap = 20 * 20 / (G - 1))
utils::write.csv(rbind(dg, dh), file.path(out, "tables", "idea2_diagnostic_all.csv"), row.names = FALSE); say("half-split detector: mean per-split fraction rho < 0.5 = ", sprintf("%.4f", dh$frac_rho_lt_0.5))
rho_g <- rh$gwps$rho; rho_ghard <- ifelse(rh$gwps$in_universe, as.numeric(rh$gwps$p <= 0.05), 1)
rho_list <- list("G-SOFT" = rho_g, "G-HARD" = rho_ghard, "HALF-SOFT" = rho_half, "R-SOFT" = rh$rpe1$rho)
for (k in 1:5) { set.seed(20261400L + k); rho_list[[sprintf("PERM-%d", k)]] <- sample(rho_g) }
gate_op <- function(rho) { B <- A * outer(rho, rho); diag(B) <- 0; diag(B) <- 1 - rowSums(B); B }
summ <- list()
for (nm in names(rho_list)) { f <- file.path(cd, sprintf("i2_%s.rds", nm)); rho <- rho_list[[nm]]
  summ[[nm]] <- data.frame(readout = nm, frac_rho_lt_0.5 = mean(rho < 0.5), mean_rho = mean(rho), n_rho_eq_1 = sum(rho >= 1), n_rho_eq_0 = sum(rho <= 0))
  if (file.exists(f)) next
  t0 <- Sys.time(); B <- gate_op(rho); stopifnot(max(abs(rowSums(B) - 1)) < 1e-9); S <- ei_attention(B); rm(B)
  r <- ev(S); rm(S); saveRDS(r, f); say(sprintf("%s scored (%.0f s)", nm, as.numeric(difftime(Sys.time(), t0, units = "secs")))); invisible(gc(verbose = FALSE)) }
utils::write.csv(do.call(rbind, summ), file.path(out, "tables", "idea2_rho_summary.csv"), row.names = FALSE)
R <- list(CLR = readRDS(file.path(cd, "ref_CLR.rds")), K10 = readRDS(file.path(cd, "ref_K10.rds")))
for (nm in names(rho_list)) R[[nm]] <- readRDS(file.path(cd, sprintf("i2_%s.rds", nm)))
tb <- ei_tables(R, env$base_rate, file.path(out, "tables"), "idea2", say = say)
pm <- tb$pm; get <- function(r, e) pm$aupr[pm$readout == r & pm$scope == "pooled" & pm$truth == e]
ds <- get("G-SOFT", "CORUM_small") - get("K10", "CORUM_small"); dd <- get("G-SOFT", "DoRothEA_AB") - get("K10", "DoRothEA_AB")
rule <- data.frame(readout = "G-SOFT", dCORUM_small_aupr_vs_K10 = ds, dDoRothEA_AB_aupr_vs_K10 = dd, verdict = if (ds >= 0.005 && dd >= -0.005) "PASS" else "FAIL", diagnostic_flag = flag)
utils::write.csv(rule, file.path(out, "tables", "idea2_passfail.csv"), row.names = FALSE); print(rule); say("done")
