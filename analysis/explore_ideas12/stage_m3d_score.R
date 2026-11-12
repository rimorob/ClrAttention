#!/usr/bin/env Rscript
# EXPLORATORY. Conditional stage: score the M3D lab-split consensus readouts (CONS-R, CONS-U, CONS-G, NULL-U) on RegulonDB SC and all.
source("analysis/explore_ideas12/ideas12_lib.R"); source("analysis/explore_ideas12/ideas12_m3d_lib.R"); source("analysis/explore_ideas12/ideas12_tables.R")
out <- file.path(EI_OUT, "m3d"); cd <- file.path(out, "cache"); say <- ei_say(file.path(out, "stage_m3d_score.log"))
sw <- file.path(out, "sweep"); wt <- utils::read.csv(file.path(sw, "weights.csv")); NG <- nrow(wt)
core <- readRDS(file.path(cd, "m3d_core.rds")); G <- length(core$genes); q <- core$q_att; H_att <- outer(q, q, "+"); H_mi <- outer(core$mi_rm, core$mi_rm, "+"); rm(core)
env <- ei_setup_m3d(say); stopifnot(env$G == G); ev <- ei_make_eval(env, H_att, H_mi, say); rm(H_att, H_mi)
read_cnt <- function(kind, gi) matrix(readBin(file.path(sw, sprintf("counts_%s_g%02d.u8", kind, gi - 1L)), "integer", n = G * G, size = 1L, signed = FALSE), G, G)
forms <- list(R = wt$w_R, U = wt$w_U, G = wt$w_G)
for (kind in c("obs", "null")) { fs <- if (kind == "obs") c("R", "U", "G") else "U"; need <- fs[!file.exists(file.path(cd, sprintf("%s_%s.rds", kind, fs)))]; if (!length(need)) next
  Sacc <- lapply(need, function(f) matrix(0, G, G)); names(Sacc) <- need
  for (gi in seq_len(NG)) { Cn <- read_cnt(kind, gi) / 60; for (f in need) Sacc[[f]] <- Sacc[[f]] + forms[[f]][gi] * Cn; rm(Cn) }
  for (f in need) { r <- ev(Sacc[[f]]); saveRDS(r, file.path(cd, sprintf("%s_%s.rds", kind, f))); say(sprintf("%s-%s scored", kind, f)); Sacc[[f]] <- NULL; invisible(gc(verbose = FALSE)) } }
sg <- list(); for (gi in seq_len(NG)) { r <- ev(read_cnt("obs", gi) / 60, hub = FALSE, jack = FALSE)$pooled; sg[[gi]] <- data.frame(gamma = wt$gamma[gi], truth = names(r), aupr = vapply(r, function(x) x$aupr, 0)) }
saveRDS(do.call(rbind, sg), file.path(cd, "single_gamma.rds"))
R <- list(CLR = readRDS(file.path(cd, "ref_CLR.rds")), K10 = readRDS(file.path(cd, "ref_K10.rds")))
nmv <- c(obs_R = "CONS-R", obs_U = "CONS-U", obs_G = "CONS-G", null_U = "NULL-U"); for (k in names(nmv)) R[[nmv[[k]]]] <- readRDS(file.path(cd, paste0(k, ".rds")))
ei_tables(R, env$base_rate, file.path(out, "tables"), "m3d", say = say); say("done")
