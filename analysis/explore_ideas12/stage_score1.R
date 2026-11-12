#!/usr/bin/env Rscript
# EXPLORATORY. Step 4 of docs/plans/explore_ideas12_plan.md: score the idea-1 consensus co-membership readouts (K562, in-sample).
args <- commandArgs(TRUE); workers <- if (length(args)) as.integer(args[1]) else 6L
source("analysis/explore_ideas12/ideas12_lib.R"); source("analysis/explore_ideas12/ideas12_tables.R")
out <- file.path(EI_OUT, "k562"); cd <- file.path(out, "cache"); say <- ei_say(file.path(out, "stage_score1.log"))
sw <- file.path(out, "sweep"); wt <- utils::read.csv(file.path(sw, "weights.csv")); NG <- nrow(wt)
core <- readRDS(file.path(cd, "k562_core.rds")); G <- length(core$genes)
q <- core$q_att; H_att <- outer(q, q, "+"); H_mi <- outer(core$mi_rm, core$mi_rm, "+"); rm(core)
env <- ei_setup_k562(say, workers); stopifnot(env$G == G)
ev <- ei_make_eval(env, H_att, H_mi, say); rm(H_att, H_mi)
nruns <- c(obs = 60L, null = 60L, full = 5L)
read_cnt <- function(kind, gi) matrix(readBin(file.path(sw, sprintf("counts_%s_g%02d.u8", kind, gi - 1L)), "integer", n = G * G, size = 1L, signed = FALSE), G, G)
forms <- list(R = wt$w_R, U = wt$w_U, G = wt$w_G)
## single-gamma diagnostics: pooled point estimates only (raw scope, no jackknife)
sg_f <- file.path(cd, "single_gamma.rds")
if (!file.exists(sg_f)) { sg <- list(); for (gi in seq_len(NG)) { S <- read_cnt("obs", gi) / nruns[["obs"]]; r <- ev(S, hub = FALSE, jack = FALSE)$pooled
  sg[[gi]] <- data.frame(gamma = wt$gamma[gi], truth = names(r), aupr = vapply(r, function(x) x$aupr, 0)); rm(S); say(sprintf("single gamma %d/%d", gi, NG)) }
  saveRDS(do.call(rbind, sg), sg_f) }
for (kind in c("obs", "full", "null")) {
  fs <- if (kind == "full") c("R", "U") else if (kind == "null") c("U", "R") else c("R", "U", "G")
  need <- fs[!file.exists(file.path(cd, sprintf("%s_%s.rds", kind, fs)))]; if (!length(need)) next
  Sacc <- lapply(need, function(f) matrix(0, G, G)); names(Sacc) <- need
  for (gi in seq_len(NG)) { Cn <- read_cnt(kind, gi) / nruns[[kind]]; for (f in need) Sacc[[f]] <- Sacc[[f]] + forms[[f]][gi] * Cn; rm(Cn) }
  invisible(gc(verbose = FALSE))
  for (f in need) { t0 <- Sys.time(); r <- ev(Sacc[[f]]); saveRDS(r, file.path(cd, sprintf("%s_%s.rds", kind, f)))
    say(sprintf("%s-%s scored (%.0f s)", kind, f, as.numeric(difftime(Sys.time(), t0, units = "secs")))); Sacc[[f]] <- NULL; invisible(gc(verbose = FALSE)) } }
## tables
R <- list(CLR = readRDS(file.path(cd, "ref_CLR.rds")), K10 = readRDS(file.path(cd, "ref_K10.rds")))
nm <- c(obs_R = "CONS-R", obs_U = "CONS-U", obs_G = "CONS-G", full_R = "FULL-R", full_U = "FULL-U", null_U = "NULL-U", null_R = "NULL-R")
for (k in names(nm)) R[[nm[[k]]]] <- readRDS(file.path(cd, paste0(k, ".rds")))
tb <- ei_tables(R, env$base_rate, file.path(out, "tables"), "idea1", say = say)
## literal pass/fail on CONS-R (and the same rule applied to the other weightings, for the disguise flag)
pm <- tb$pm; get <- function(r, e, sc = "pooled") pm$aupr[pm$readout == r & pm$scope == sc & pm$truth == e]
rule <- do.call(rbind, lapply(c("CONS-R", "CONS-U", "CONS-G", "FULL-R", "FULL-U", "NULL-U", "NULL-R"), function(r) {
  enr <- get(r, "CORUM_all") / env$base_rate[["CORUM_all"]]; dd <- get(r, "DoRothEA_AB") - get("K10", "DoRothEA_AB")
  data.frame(readout = r, CORUM_all_enrichment = enr, K10_CORUM_all_enrichment = get("K10", "CORUM_all") / env$base_rate[["CORUM_all"]], dDoRothEA_AB_aupr_vs_K10 = dd,
    verdict = if (enr >= 6.4 && dd >= -0.005) "PASS" else if (enr < get("K10", "CORUM_all") / env$base_rate[["CORUM_all"]] || dd <= -0.01) "FAIL" else "between") }))
utils::write.csv(rule, file.path(out, "tables", "idea1_passfail.csv"), row.names = FALSE); print(rule)
say("done")
