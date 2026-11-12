#!/usr/bin/env Rscript
# EXPLORATORY. Step 1 of docs/plans/explore_ideas12_plan.md: equivalence gate (CLR and K10 attention t=9 on the K562 truths reproduce
# the stored robust_z V0 caches and the baselines metrics.csv to 1e-9) and caches for the later stages. STOPS (non-zero exit) on failure.
args <- commandArgs(TRUE); threads <- if (length(args)) as.integer(args[1]) else 8L
source("analysis/explore_ideas12/ideas12_lib.R")
out <- file.path(EI_OUT, "k562"); dir.create(file.path(out, "cache"), recursive = TRUE, showWarnings = FALSE)
say <- ei_say(file.path(out, "stage0_gate.log"))
env <- ei_setup_k562(say); G <- env$G; X <- env$X
t0 <- Sys.time()
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = threads); f$calibrate(method = "normal", combine = "stouffer")
M <- f$mi; C <- f$clr_scores; f$build_operator(alpha = EI_ALPHA, softmax_keff = 10, softmax_cap = 50L); A <- f$operator
say(sprintf("MI + CLR + operator: %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
P <- Matrix::Diagonal(G, 1 - EI_ALPHA) + EI_ALPHA * Matrix::Matrix(A, sparse = TRUE)
S9 <- ei_attention(A)
q_att <- rz_stationary(P); H_att <- outer(q_att, q_att, "+")
mi_rm <- (rowSums(M) - diag(M)) / (G - 1); H_mi <- outer(mi_rm, mi_rm, "+")
saveRDS(list(M = M, C = C, A = Matrix::Matrix(A, sparse = TRUE), q_att = q_att, mi_rm = mi_rm, genes = env$genes), file.path(out, "cache", "k562_core.rds"))
f$release(); rm(f); invisible(gc(verbose = FALSE))
ev <- ei_make_eval(env, H_att, H_mi, say); rm(H_att, H_mi)
fc <- file.path(out, "cache", "ref_CLR.rds"); fk <- file.path(out, "cache", "ref_K10.rds")
if (file.exists(fc)) r_clr <- readRDS(fc) else { say("scoring CLR"); r_clr <- ev(C); saveRDS(r_clr, fc) }
if (file.exists(fk)) r_k10 <- readRDS(fk) else { say("scoring K10 (attention t=9)"); r_k10 <- ev(S9); saveRDS(r_k10, fk) }
## gate (a): robust_z V0 caches, 7 truths, three scopes, point and all 100 replicates
rz <- "results/robust_z/k562/cache"; tol <- 1e-9; rows <- list()
for (nm in c("CLR", "K10")) { kk <- c(CLR = "V0__pair.rds", K10 = "V0__att9.rds")[[nm]]; ref <- readRDS(file.path(rz, kk)); mine <- if (nm == "CLR") r_clr else r_k10
  for (sc in c("pooled", "hub_att", "hub_mi")) for (e in names(ref[[sc]])) {
    rows[[length(rows) + 1L]] <- data.frame(check = "robust_z V0 cache", readout = nm, scope = sc, truth = e, max_abs_diff_aupr = abs(mine[[sc]][[e]]$aupr - ref[[sc]][[e]]$aupr),
      max_abs_diff_jack = max(abs(mine[[sc]][[e]]$jack - ref[[sc]][[e]]$jack))) } }
## gate (b): baselines metrics.csv (CLR, attention t=9 AUPR)
bm <- utils::read.csv("results/baselines_k562ess_truthgenes_genie3/metrics.csv", stringsAsFactors = FALSE)
for (e in c("DoRothEA_AB", "CORUM", "KD_response", "CollecTRI")) { te <- if (e == "CORUM") "CORUM_all" else e; b <- bm[bm$truth == e, ]
  for (mn in c("CLR", "K10")) { m <- c(CLR = "clr", K10 = "attention")[[mn]]; mine <- if (mn == "CLR") r_clr else r_k10
    rows[[length(rows) + 1L]] <- data.frame(check = "baselines metrics.csv", readout = mn, scope = "pooled", truth = te, max_abs_diff_aupr = abs(mine$pooled[[te]]$aupr - b$aupr[b$method == m]), max_abs_diff_jack = NA_real_) } }
eq <- do.call(rbind, rows); eq$pass <- eq$max_abs_diff_aupr < tol & (is.na(eq$max_abs_diff_jack) | eq$max_abs_diff_jack < tol)
utils::write.csv(eq, file.path(out, "equivalence.csv"), row.names = FALSE)
say(sprintf("equivalence: %d checks, max |dAUPR| %.3g, max |djack| %.3g", nrow(eq), max(eq$max_abs_diff_aupr), max(eq$max_abs_diff_jack, na.rm = TRUE)))
if (!all(eq$pass)) { print(eq[!eq$pass, ]); say("EQUIVALENCE FAILED: STOP"); quit(status = 2) }
say("EQUIVALENCE PASSED (1e-9)"); say("done")
