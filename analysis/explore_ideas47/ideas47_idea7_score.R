#!/usr/bin/env Rscript
# EXPLORATORY, IN-SAMPLE, NOT PREREGISTERED. Idea 7 (soft MCL inflation between walk steps), scoring stage.
# Readouts fixed in docs/plans/explore_ideas47_plan.md section 2; depths from the label-free depth stage
# (ideas47_idea7_depth.R: tstar_by_r.csv, joint_selection.csv, written before any truth is scored).
# Usage: Rscript ideas47_idea7_score.R [--out results/explore_ideas47/idea7] [--threads 8] [--workers 8] [--jack 100] [--quick 0]
args <- commandArgs(TRUE)
opt <- list(out = "results/explore_ideas47/idea7", threads = 8L, workers = 8L, jack = 100L, quick = 0L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "ideas47_lib.R")); source(file.path(script_dir, "ideas47_inflate.R"))
tr <- utils::read.csv(file.path(opt$out, "tstar_by_r.csv")); jt <- utils::read.csv(file.path(opt$out, "joint_selection.csv"))
say(sprintf("depths: %s; joint (r*, t*) = (%.1f, %d)", paste(sprintf("r%.1f: t*=%d", tr$r, tr$t_star_r), collapse = "; "), jt$r_star, jt$t_star))

## ---- gate (d), part 1: r = 1 inflation code reproduces the K10 walk; (a)-(c) on references --------------------------------------
S9r1 <- NULL; inflate_walk(P, 1, 9L, function(t, Q) if (t == 9L) S9r1 <<- sym0(Q))
d_walk <- max(abs(S9r1 - S9)); rm(S9r1)
say(sprintf("gate (d): max|S_9(r = 1 inflation code) - S_9(K10)| = %.3g", d_walk))
utils::write.csv(data.frame(check = "(d) r = 1 inflation walk S_9 vs K10 walk S_9", max_abs_diff = d_walk, pass = d_walk < 1e-12), file.path(opt$out, "gate_d_walk.csv"), row.names = FALSE)
if (!QUICK) stopifnot(d_walk < 1e-12)
res_clr <- score_readout("CLR", C); res_k10 <- score_readout("K10", S9); res_a3r <- score_readout("a3r", a3r(Pt9))
if (!QUICK) e47_gate(res_clr, res_k10, res_a3r) else say("QUICK run: gate (a)-(c) skipped (gene subset)")

## ---- references: reciprocity and degree-scaled (corum_expert2 A5) ---------------------------------------------------------------
score_readout("K10_min", minr(Pt9)); score_readout("K10_deg", degr(Pt9))

## ---- readout table, then scoring along each r's walk ---------------------------------------------------------------------------
key <- function(r, t) sprintf("INF_r%.1f_t%d", r, t)
own <- stats::setNames(tr$t_star_r, sprintf("%.1f", tr$r))
rd <- rbind(data.frame(name = "IN-J", r = jt$r_star, t = jt$t_star, min = FALSE),
            do.call(rbind, lapply(c(1.1, 1.2, 1.5), function(r) data.frame(name = c(sprintf("IN-r%.1f-own", r), sprintf("IN-r%.1f-t9", r)), r = r, t = c(own[[sprintf("%.1f", r)]], 9L), min = FALSE))),
            data.frame(name = "IN-J-min", r = jt$r_star, t = jt$t_star, min = TRUE))
rd$canonical <- ifelse(rd$min, paste0(key(rd$r, rd$t), "_min"), key(rd$r, rd$t))
mp <- rbind(data.frame(name = c("K10", "CLR", "a3r", "K10_min", "K10_deg"), canonical = c("K10", "CLR", "a3r", "K10_min", "K10_deg"), r = NA, t = NA),
            data.frame(name = rd$name, canonical = rd$canonical, r = rd$r, t = rd$t))
utils::write.csv(mp, file.path(opt$out, "readout_map.csv"), row.names = FALSE); print(mp, row.names = FALSE)
for (r in unique(rd$r)) { need <- rd[rd$r == r, ]
  inflate_walk(P, r, max(need$t), function(t, Q) { for (k in which(need$t == t)) { cn <- need$canonical[k]
    score_readout(cn, if (need$min[k]) minr(Q) else sym0(Q)) } })
  invisible(gc(verbose = FALSE)) }
saveRDS(list(base_rate = base_rate, truths = TRUTHS, opt = opt, map = mp), file.path(opt$out, "cache", "meta.rds"))
say("done")
