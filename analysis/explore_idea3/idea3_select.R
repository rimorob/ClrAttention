#!/usr/bin/env Rscript
# EXPLORATORY (docs/plans/explore_idea3_plan.md, sections 4, 6, 7). Selection stage: label-free held-out criterion curves per hub decile
# (rows of the train-half P^t against the test half's top-1% CLR pairs), the global criterion at the grid depths, the PPR and heat-kernel
# controls, the rewired-operator and label-permutation nulls, and the frozen selection (t_d per decile and stratifier, control
# parameters). NO truth is loaded here. Needs results/explore_idea3/<dataset>/deciles.csv from idea3_knee_null.R.
# Usage (repo root): Rscript analysis/explore_idea3/idea3_select.R --dataset m3d|k562 [--nsplit 16] [--workers 6] [--quick 0]
args <- commandArgs(TRUE)
opt <- list(dataset = "m3d", out = NULL, nsplit = 16L, workers = 6L, quick = 0L, depth_seed = 20261014L, perm_seed = 20261113L,
            nperm = 4L, rewire = 1L, rz_curves = NULL)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); if (!key %in% names(opt)) stop("unknown argument: ", args[k])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[k + 1]) else args[k + 1] }
for (k in c("nsplit", "workers", "quick", "depth_seed", "perm_seed", "nperm", "rewire")) opt[[k]] <- as.integer(opt[[k]])
if (is.null(opt$out)) opt$out <- file.path("results/explore_idea3", opt$dataset)
if (is.null(opt$rz_curves)) opt$rz_curves <- file.path("results/robust_z", opt$dataset, "heldout_curves.csv")
dir.create(file.path(opt$out, "cache"), recursive = TRUE, showWarnings = FALSE)
script_dir <- (function() { f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis/explore_idea3" })()
source(file.path(script_dir, "idea3_lib.R")); i3_source_all()
logf <- file.path(opt$out, "select.log"); say <- function(...) i3_say(logf, ...)
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

D <- i3_load(opt$dataset, opt$quick); X <- D$X; lab <- D$lab; genes <- D$genes; G <- length(genes); rm(D)
dec <- utils::read.csv(file.path(opt$out, "deciles.csv"), stringsAsFactors = FALSE); stopifnot(identical(dec$gene, genes))
say(sprintf("%s: %d genes x %d samples (%d labs)", opt$dataset, G, ncol(X), length(unique(lab))))
UT <- which(upper.tri(matrix(FALSE, G, G)))
invisible(i3_rewire_cpp())                              # compile the rewiring routine in the parent, before forking

# label sets: real stratifiers and fixed random permutations of the MI decile labels
set.seed(opt$perm_seed)
labsets <- list(mi = dec$dec_mi, pi = dec$dec_pi)
for (k in seq_len(opt$nperm)) labsets[[paste0("perm", k)]] <- sample(dec$dec_mi)
rewire_sets <- "mi"                                        # label sets evaluated on the rewired operator
idx_of <- function(d) lapply(1:10, function(k) { R <- which(d == k); rr <- rep(R, times = G); cc <- rep(seq_len(G), each = length(R))
  (rr + (cc - 1L) * G)[rr != cc] })
IDX <- lapply(labsets, idx_of)
# control weights over t = 0..TMAX
tt <- 0:I3_TMAX
w_ppr <- t(vapply(I3_PPR, function(a) { w <- (1 - a) * a^tt; w / sum(w) }, numeric(length(tt))))
w_hk  <- t(vapply(I3_HK, function(tau) { w <- stats::dpois(tt, tau); w / sum(w) }, numeric(length(tt))))
W <- rbind(w_ppr, w_hk); ctrl <- data.frame(family = c(rep("ppr", length(I3_PPR)), rep("hk", length(I3_HK))), param = c(I3_PPR, I3_HK))

set.seed(opt$depth_seed); splits <- lapply(seq_len(opt$nsplit), function(s) i3_lab_split(lab))
utils::write.csv(data.frame(sample = colnames(X), do.call(cbind, stats::setNames(splits, sprintf("split_%02d", seq_along(splits))))), file.path(opt$out, "split_assignments.csv"), row.names = FALSE)

run_split <- function(si) {
  h <- splits[[si]]; set.seed(opt$depth_seed + si); t0 <- Sys.time()
  Ms <- lapply(1:2, function(d) { Xm <- X[, h == d]; rng <- apply(Xm, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) Xm[rng <= 0, ] <- Xm[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9)
    bspline_mi(Xm, bins = "hg", spline_order = 3, threads = 1, transform = "none") })
  Cs <- lapply(Ms, rz_scores, v = "V0"); rm(Ms)
  rows <- list(); crow <- list()
  for (dir in 1:2) { te <- 3L - dir
    sv <- Cs[[te]][UT]; tf <- sv >= stats::quantile(sv, 0.99); rm(sv)
    Tm <- matrix(FALSE, G, G); Tm[UT] <- tf; Tm <- Tm | t(Tm); tgt_ut <- as.integer(tf); rm(tf)
    TV <- lapply(IDX, function(il) lapply(il, function(ix) Tm[ix])); rm(Tm)
    A <- rz_operator(Cs[[dir]])
    for (op in c("real", if (opt$rewire > 0L) "rewired")) {
      Aop <- if (op == "real") A else i3_rewire(A, 1000L * si + dir)
      P <- rz_lazy_P(Aop); rm(Aop)
      sets <- if (op == "real") names(labsets) else rewire_sets
      acc <- NULL
      if (op == "real") { acc <- lapply(seq_len(nrow(W)), function(k) { m <- matrix(0, G, G); diag(m) <- W[k, 1]; m }) }
      i3_walk(P, I3_TMAX, function(Pt, t) {
        if (!is.null(acc)) for (k in seq_along(acc)) acc[[k]] <<- acc[[k]] + W[k, t + 1L] * Pt
        if (t %in% I3_GRID) {
          for (ls in sets) for (d in 1:10)
            rows[[length(rows) + 1L]] <<- data.frame(split = si, direction = dir, op = op, labelset = ls, decile = d, t = t,
                                                     criterion = ap_fast(Pt[IDX[[ls]][[d]]], TV[[ls]][[d]]))
          if (op == "real") { S <- (Pt + t(Pt)) / 2
            rows[[length(rows) + 1L]] <<- data.frame(split = si, direction = dir, op = op, labelset = "global", decile = 0L, t = t,
                                                     criterion = ap_fast(S[UT], tgt_ut)); rm(S) }
        } })
      if (op == "real") {
        for (k in seq_along(acc)) { S <- (acc[[k]] + t(acc[[k]])) / 2
          crow[[length(crow) + 1L]] <- data.frame(split = si, direction = dir, family = ctrl$family[k], param = ctrl$param[k], criterion = ap_fast(S[UT], tgt_ut)); rm(S) }
        rm(acc) }
      rm(P); invisible(gc(verbose = FALSE))
    }
    rm(TV, A, tgt_ut); invisible(gc(verbose = FALSE)) }
  list(curves = do.call(rbind, rows), controls = do.call(rbind, crow), secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

todo <- which(!file.exists(file.path(opt$out, "cache", sprintf("sel_split_%02d.rds", seq_len(opt$nsplit)))))
say(sprintf("%d splits to run (%d done), %d workers; label sets: %s; rewired operator: %s", length(todo), opt$nsplit - length(todo),
            min(opt$workers, max(1L, length(todo))), paste(names(labsets), collapse = ","), if (opt$rewire > 0L) "yes" else "no"))
if (length(todo)) {
  invisible(gc(verbose = FALSE))                       # gc before forking (R sets the GC threshold from the parent heap)
  res <- parallel::mclapply(todo, function(si) { r <- run_split(si)
    saveRDS(r, file.path(opt$out, "cache", sprintf("sel_split_%02d.rds.tmp", si)))
    file.rename(file.path(opt$out, "cache", sprintf("sel_split_%02d.rds.tmp", si)), file.path(opt$out, "cache", sprintf("sel_split_%02d.rds", si)))
    r$secs }, mc.cores = min(opt$workers, length(todo)), mc.preschedule = FALSE)
  bad <- vapply(res, function(x) !is.numeric(x), NA)
  if (any(bad)) stop("selection workers failed for splits ", paste(todo[bad], collapse = ","), ": ", paste(unique(vapply(res[bad], function(x) paste(as.character(x), collapse = " "), "")), collapse = " | "))
  say("split seconds: ", paste(sprintf("%.0f", unlist(res)), collapse = " "))
}
S_all <- lapply(seq_len(opt$nsplit), function(si) readRDS(file.path(opt$out, "cache", sprintf("sel_split_%02d.rds", si))))
cur <- do.call(rbind, lapply(S_all, `[[`, "curves")); cct <- do.call(rbind, lapply(S_all, `[[`, "controls")); rm(S_all)
utils::write.csv(cur, file.path(opt$out, "heldout_curves.csv"), row.names = FALSE); utils::write.csv(cct, file.path(opt$out, "heldout_controls.csv"), row.names = FALSE)

## ---- equivalence gate 2: global criterion vs results/robust_z (V0) ------------------------------------
rz <- utils::read.csv(opt$rz_curves, stringsAsFactors = FALSE); rz <- rz[rz$variant == "V0", ]
gl <- cur[cur$labelset == "global" & cur$op == "real", ]
m <- merge(gl, rz, by.x = c("split", "direction", "t"), by.y = c("split", "direction", "t"), suffixes = c(".new", ".rz"))
eqd <- if (nrow(m)) max(abs(m$criterion.new - m$criterion.rz)) else NA_real_
writeLines(sprintf("held-out criterion machinery vs %s (V0): %d split-direction-depth rows compared, max |diff| %.3e, tolerance 1e-9 %s",
                   opt$rz_curves, nrow(m), eqd, if (!is.na(eqd) && eqd <= 1e-9) "PASS" else "FAIL"), file.path(opt$out, "equivalence_criterion.txt"))
say(sprintf("EQUIVALENCE (criterion machinery vs robust_z V0): %d rows, max |diff| %.3e", nrow(m), eqd))
if (opt$quick == 0L && (is.na(eqd) || eqd > 1e-9)) stop("criterion equivalence FAILED")

## ---- frozen selection (plan sections 4 and 6) -------------------------------------------------------------
pos_g <- log2(I3_GRID)
sel <- list(); ncurve <- length(unique(paste(cur$split, cur$direction)))
for (op in unique(cur$op)) for (ls in unique(cur$labelset[cur$op == op])) for (d in sort(unique(cur$decile[cur$op == op & cur$labelset == ls]))) {
  x <- cur[cur$op == op & cur$labelset == ls & cur$decile == d, ]
  mc <- tapply(x$criterion, x$t, mean)[as.character(I3_GRID)]
  se <- tapply(x$criterion, x$t, stats::sd)[as.character(I3_GRID)] / sqrt(length(unique(paste(x$split, x$direction))))
  pk <- i3_pick(as.numeric(mc), pos_g)
  am <- vapply(split(x, paste(x$split, x$direction)), function(z) z$t[which.max(z$criterion)], 0)
  sel[[length(sel) + 1L]] <- data.frame(op = op, labelset = ls, decile = d, t_d = I3_GRID[pk$idx], t_argmax = I3_GRID[pk$argmax], t_lo = I3_GRID[pk$lo], t_hi = I3_GRID[pk$hi],
    width_octaves = pk$width_oct, crit_at_t_d = as.numeric(mc[pk$idx]), crit_max = as.numeric(mc[pk$argmax]), se_at_max = as.numeric(se[pk$argmax]), n_curves = length(am),
    splitdir_argmax_median = stats::median(am), splitdir_argmax_q25 = stats::quantile(am, 0.25, names = FALSE), splitdir_argmax_q75 = stats::quantile(am, 0.75, names = FALSE))
}
sel <- do.call(rbind, sel); utils::write.csv(sel, file.path(opt$out, "selected_depths.csv"), row.names = FALSE)
sc <- list()
for (fam in c("ppr", "hk")) { x <- cct[cct$family == fam, ]; ps <- if (fam == "ppr") I3_PPR else I3_HK
  mc <- tapply(x$criterion, x$param, mean)[as.character(ps)]; pos <- if (fam == "ppr") log2(1 / (1 - ps)) else log2(ps)
  pk <- i3_pick(as.numeric(mc), pos)
  sc[[fam]] <- data.frame(family = fam, param = ps[pk$idx], param_argmax = ps[pk$argmax], param_lo = ps[pk$lo], param_hi = ps[pk$hi], width_octaves = pk$width_oct,
                          crit_at_choice = as.numeric(mc[pk$idx]), crit_max = as.numeric(mc[pk$argmax]), n_curves = length(unique(paste(x$split, x$direction)))) }
gx <- cur[cur$labelset == "global" & cur$op == "real", ]; mc <- tapply(gx$criterion, gx$t, mean)[as.character(I3_GRID)]; pk <- i3_pick(as.numeric(mc), pos_g)
sc[["global_grid"]] <- data.frame(family = "global_grid", param = I3_GRID[pk$idx], param_argmax = I3_GRID[pk$argmax], param_lo = I3_GRID[pk$lo], param_hi = I3_GRID[pk$hi],
                                  width_octaves = pk$width_oct, crit_at_choice = as.numeric(mc[pk$idx]), crit_max = as.numeric(mc[pk$argmax]), n_curves = length(unique(paste(gx$split, gx$direction))))
sc <- do.call(rbind, sc); utils::write.csv(sc, file.path(opt$out, "selected_controls.csv"), row.names = FALSE)
tg <- sc$param[sc$family == "global_grid"]
st <- do.call(rbind, lapply(split(sel, paste(sel$op, sel$labelset)), function(s) { s <- s[order(s$decile), ]; l2 <- log2(s$t_d)
  data.frame(op = s$op[1], labelset = s$labelset[1], sd_log2_t = stats::sd(l2), range_octaves = diff(range(l2)),
             spearman_decile_vs_log2t = suppressWarnings(stats::cor(s$decile, l2, method = "spearman")),
             n_deciles_ge1oct_from_9 = sum(s$t_d <= 4 | s$t_d >= 24), n_deciles_differing_from_global_grid = sum(s$t_d != tg),
             mean_width_octaves = mean(s$width_octaves), global_grid_t = tg) }))
utils::write.csv(st, file.path(opt$out, "selection_stats.csv"), row.names = FALSE)
print(sel[sel$op == "real" & sel$labelset %in% c("mi", "pi"), ], row.names = FALSE); print(sc, row.names = FALSE); print(st, row.names = FALSE)
say("selection done")
