#!/usr/bin/env Rscript
# Figure 7b and 7c (single-cell data, BEELINE).
#
# 7b: per-dataset log2 AUPRC ratio of GENIE3, PIDC and GRNBoost2 relative to CLR
#     and to CLR-NS (7 BEELINE datasets; per dataset the mean over gene sets and
#     primary truths of log2(method / reference), exactly as the summaries in
#     analysis/beeline_baselines.R). Competitors grey, one shape each; counts and
#     means are descriptive (no significance marks).
#     Input: results/beeline_baselines/metrics.csv, summary.csv (the production
#     run with the legacy-early-stopping GRNBoost2; see
#     docs/plans/grnboost2_rerun_change_inventory.md).
# 7c: preregistered learned-attention comparison (docs/scfm_scoring.md).
#     (i)  d_i = CLR-NS minus scPRINT R0, log2 AUPRC ratio, per dataset, with the
#          Hodges-Lehmann estimate and exact interval as a band, zero line and the
#          non-inferiority margin -0.15.
#     (ii) controls: R0 minus the permutation control and R0 minus the ESM2
#          gene-embedding cosine (positive favours R0).
#     (iii) the primary contrast by truth type, cell-type ChIP shown plainly.
#     Input: results/scfm_score/primary_cells.csv, contrasts_summary.csv. Every
#     number printed on the figure is read from these files.
# Usage (repo root):  Rscript analysis/fig7bc.R [--out docs/whitepapers/figures]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(out = "docs/whitepapers/figures")
i <- 1L
while (i <= length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1L]; i <- i + 2L }
source("analysis/plot_theme.R")
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
save_both <- function(p, stem, width, height) {
  f <- file.path(opt$out, paste0(stem, ".png"))
  save_fig(p, f, width = width, height = height)
  ggsave(sub("png$", "pdf", f), p, width = width, height = height, device = grDevices::cairo_pdf,
         bg = ink[["surface"]])
}
fmt <- function(x, d = 2) sprintf("%+.*f", d, x)
ds_order <- c("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L")
fam_grey <- c(GENIE3 = ink[["secondary"]], PIDC = ink[["muted"]], GRNBoost2 = "#b9b8b2")
fam_shape <- c(GENIE3 = 16, PIDC = 15, GRNBoost2 = 17)

# ---------------------------------------------------------------- 7b
m <- read.csv("results/beeline_baselines/metrics.csv", check.names = FALSE)
sm <- read.csv("results/beeline_baselines/summary.csv", check.names = FALSE)
primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
m <- m[m$truth %in% primary, ]
m$key <- paste(m$dataset, m$N, m$truth)
meth <- c(genie3 = "GENIE3", pidc = "PIDC", grnboost2 = "GRNBoost2")
refs <- c(clr = "relative to CLR", clr_ns = "relative to CLR-NS")
d7b <- do.call(rbind, lapply(names(refs), function(rf) {
  r <- m[m$method == rf, ]
  do.call(rbind, lapply(names(meth), function(me) {
    y <- m[m$method == me, ]
    y$lr <- log2(y$auprc_ratio / r$auprc_ratio[match(y$key, r$key)])
    a <- tapply(y$lr, y$dataset, mean)
    data.frame(ref = refs[[rf]], refkey = rf, method = meth[[me]], mkey = me,
               dataset = names(a), lr = as.numeric(a))
  }))
}))
# Check against summary.csv (the source of the printed means and counts).
chk <- do.call(rbind, lapply(split(d7b, list(d7b$refkey, d7b$mkey), drop = TRUE), function(z) {
  s <- sm[sm$reference == z$refkey[1] & sm$truths == "primary" & sm$method == z$mkey[1], ]
  data.frame(refkey = z$refkey[1], mkey = z$mkey[1], mean = mean(z$lr), n_better = sum(z$lr > 0),
             d_mean = mean(z$lr) - s$mean_log2_auprc_vs_ref, d_n = sum(z$lr > 0) - s$datasets_better_auprc,
             wilcox = s$p_wilcoxon_datasets, sign = s$p_sign_datasets)
}))
stopifnot(max(abs(chk$d_mean)) < 1e-9, all(chk$d_n == 0))
print(chk[order(chk$refkey, chk$mkey), c("refkey", "mkey", "mean", "n_better", "wilcox", "sign")], digits = 3)
mean7b <- aggregate(lr ~ ref + method, d7b, mean)
mean7b$n_better <- aggregate(lr ~ ref + method, d7b, function(v) sum(v > 0))$lr
mean7b$dataset <- paste("mean:", mean7b$method)
mean7b$lab <- sprintf("%s  (higher than reference in %d of 7)", fmt(mean7b$lr), mean7b$n_better)
mean_rows <- paste("mean:", names(fam_grey))
d7b$dataset <- factor(d7b$dataset, levels = rev(c(ds_order, mean_rows)))
mean7b$dataset <- factor(mean7b$dataset, levels = levels(d7b$dataset))
d7b$method <- factor(d7b$method, levels = names(fam_grey))
mean7b$method <- factor(mean7b$method, levels = names(fam_grey))
d7b$ref <- factor(d7b$ref, levels = unname(refs))
mean7b$ref <- factor(mean7b$ref, levels = unname(refs))
xmax7b <- 0.95
xr <- c(floor(min(d7b$lr) * 5) / 5, xmax7b)

p7b <- ggplot(d7b, aes(lr, dataset)) +
  geom_hline(yintercept = 3.5, colour = ink[["muted"]], linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = pal[["orange"]], linewidth = 0.7) +
  geom_point(aes(colour = method, shape = method), size = 2.6, alpha = 0.95,
             position = position_dodge(width = 0.55)) +
  geom_point(data = mean7b, aes(lr, dataset, colour = method, shape = method),
             size = 3.6, show.legend = FALSE) +
  geom_text(data = mean7b, aes(x = 0.27, y = dataset, label = lab), inherit.aes = FALSE,
            hjust = 0, size = 2.5, colour = ink[["secondary"]]) +
  facet_wrap(~ref, ncol = 2) +
  scale_colour_manual(values = fam_grey, name = NULL) +
  scale_shape_manual(values = fam_shape, name = NULL) +
  scale_x_continuous(limits = xr, breaks = seq(-0.4, 0.2, 0.2)) +
  scale_y_discrete(limits = levels(d7b$dataset), drop = FALSE, expand = expansion(add = c(0.6, 0.6)),
                   labels = function(l) sub("^mean: ", "mean over 7: ", l)) +
  labs(x = "log2 AUPRC ratio, method relative to reference (0 = same as reference; below 0 = lower AUPRC ratio)",
       y = NULL,
       caption = paste0("Per dataset: mean over gene sets and primary truths of log2(method / reference), 7 BEELINE datasets.\n",
                        "Bottom rows: mean over the seven datasets, with the number of datasets where the method is above the reference (descriptive counts). ",
                        "GRNBoost2: legacy early stopping, production seed.")) +
  theme_clr(10) +
  theme(legend.position = "bottom", panel.spacing.x = unit(1.2, "lines"),
        axis.title.x = element_text(size = rel(0.78)))
save_both(p7b, "fig7b_baselines_vs_clr", 9.5, 5.4)

# ---------------------------------------------------------------- 7c
cs <- read.csv("results/scfm_score/contrasts_summary.csv", check.names = FALSE)
pc <- read.csv("results/scfm_score/primary_cells.csv", check.names = FALSE)
getrow <- function(pattern) {
  r <- cs[grepl(pattern, cs$label, fixed = TRUE), ]
  stopifnot(nrow(r) == 1); r
}
parse_values <- function(v) {
  kv <- strsplit(strsplit(v, "; ", fixed = TRUE)[[1]], "=", fixed = TRUE)
  setNames(as.numeric(sapply(kv, `[`, 2)), sapply(kv, `[`, 1))
}
# per-dataset d_i from primary_cells.csv (CLR-NS = ref, R0 = cand, lr = log2(ref/cand))
di <- tapply(pc$lr, pc$dataset, mean)[ds_order]
prim <- getrow("PRIMARY CLR-NS vs R0")
stopifnot(max(abs(di - parse_values(prim$values)[ds_order])) < 1e-3)
dd <- data.frame(dataset = factor(ds_order, levels = rev(ds_order)), d = as.numeric(di))
margin <- -0.15

a_lab <- data.frame(x = c(margin, 0), lab = c("non-inferiority margin -0.15", "no difference"))
pA <- ggplot(dd) +
  annotate("rect", xmin = prim$L, xmax = prim$U, ymin = -Inf, ymax = Inf, fill = pal[["orange"]], alpha = 0.16) +
  geom_vline(xintercept = 0, colour = ink[["primary"]], linewidth = 0.5) +
  geom_vline(xintercept = margin, colour = ink[["secondary"]], linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = prim$HL, colour = pal[["orange"]], linewidth = 0.9) +
  geom_segment(aes(x = 0, xend = d, y = dataset, yend = dataset), colour = pal[["orange"]], alpha = 0.45, linewidth = 0.5) +
  geom_point(aes(d, dataset), colour = pal[["orange"]], size = 3) +
  geom_label(aes(d, dataset, label = fmt(d, 3)), hjust = -0.2, size = 2.7, colour = ink[["secondary"]],
             fill = ink[["surface"]], border.colour = NA, label.padding = unit(0.12, "lines")) +
  annotate("text", x = margin, y = 7.75, label = "margin -0.15", hjust = 1.08, size = 2.6, colour = ink[["secondary"]]) +
  annotate("text", x = 0, y = 7.75, label = "no difference", hjust = -0.1, size = 2.6, colour = ink[["secondary"]]) +
  annotate("label", x = prim$HL, y = 0.25, fill = ink[["surface"]], border.colour = NA,
           label = sprintf("Hodges-Lehmann %s\nexact interval [%s, %s]", fmt(prim$HL, 3), fmt(prim$L, 3), fmt(prim$U, 3)),
           hjust = 0.5, vjust = 0, size = 2.7, colour = ink[["primary"]], lineheight = 0.95) +
  scale_x_continuous(limits = c(-0.95, 2.75), breaks = seq(0, 2.5, 0.5)) +
  scale_y_discrete(expand = expansion(add = c(1.5, 1.2))) +
  labs(title = "(i) CLR-NS minus scPRINT R0, per dataset",
       x = "log2 AUPRC-ratio difference (positive favours CLR-NS)", y = NULL) +
  theme_clr(10) + theme(panel.grid.major.y = element_blank())

# (ii) controls: R0 minus control (positive favours R0)
ctl <- list(
  list(lab = "Permutation control\n(genes permuted across cells)", pat = "scPRINT R0 minus control ctrl_perm (mean over 3 seeds"),
  list(lab = "ESM2 gene-embedding\ncosine (data-free)", pat = "scPRINT R0 minus ESM2 cosine"))
cd <- do.call(rbind, lapply(ctl, function(z) {
  r <- getrow(z$pat); v <- parse_values(r$values)[ds_order]
  list(sum = data.frame(row = z$lab, HL = r$HL, L = r$L, U = r$U),
       pts = data.frame(row = z$lab, d = as.numeric(v)))
}))
csum <- do.call(rbind, cd[, "sum"]); cpts <- do.call(rbind, cd[, "pts"])
rows_c <- sapply(ctl, `[[`, "lab")
csum$row <- factor(csum$row, levels = rev(rows_c)); cpts$row <- factor(cpts$row, levels = rev(rows_c))
csum$txt <- sprintf("%s [%s, %s]", fmt(csum$HL, 3), fmt(csum$L, 3), fmt(csum$U, 3))
pB <- ggplot() +
  geom_vline(xintercept = 0, colour = ink[["primary"]], linewidth = 0.5) +
  geom_point(data = cpts, aes(d, row), colour = ink[["muted"]], size = 1.4, alpha = 0.6,
             position = position_nudge(y = -0.18)) +
  geom_errorbar(data = csum, aes(xmin = L, xmax = U, y = row), orientation = "y", width = 0.18, colour = ink[["secondary"]], linewidth = 0.8) +
  geom_point(data = csum, aes(HL, row), colour = ink[["secondary"]], size = 3.2, shape = 18) +
  geom_label(data = csum, aes(x = -1.6, y = row, label = txt), hjust = 0, vjust = -1.1, size = 2.7, colour = ink[["primary"]],
             fill = ink[["surface"]], border.colour = NA, label.padding = unit(0.12, "lines")) +
  scale_x_continuous(limits = c(-1.6, 0.6), breaks = seq(-1.5, 0.5, 0.5)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(title = "(ii) Controls: scPRINT R0 minus control",
       x = "log2 AUPRC-ratio difference (positive favours R0)", y = NULL) +
  theme_clr(10) + theme(panel.grid.major.y = element_blank())

# (iii) by truth type
tr <- c(celltype_ChIP = "Cell-type ChIP-seq", nonspecific_ChIP = "Non-specific ChIP-seq", STRING = "STRING")
td <- do.call(rbind, lapply(names(tr), function(k) {
  r <- getrow(sprintf("CLR-NS vs R0, truth %s only", k)); v <- parse_values(r$values)[ds_order]
  x <- pc[pc$truth == k, ]
  list(sum = data.frame(row = tr[[k]], HL = r$HL, L = r$L, U = r$U,
                        n_r0 = sum(x$lr < 0), n_cells = nrow(x)),
       pts = data.frame(row = tr[[k]], d = as.numeric(v)))
}))
tsum <- do.call(rbind, td[, "sum"]); tpts <- do.call(rbind, td[, "pts"])
tsum$row <- factor(tsum$row, levels = rev(unname(tr))); tpts$row <- factor(tpts$row, levels = rev(unname(tr)))
tsum$txt <- sprintf("%s [%s, %s]", fmt(tsum$HL, 3), fmt(tsum$L, 3), fmt(tsum$U, 3))
tsum$sub <- sprintf("scPRINT R0 higher in %d of %d cells", tsum$n_r0, tsum$n_cells)
pC <- ggplot() +
  geom_vline(xintercept = 0, colour = ink[["primary"]], linewidth = 0.5) +
  geom_vline(xintercept = margin, colour = ink[["secondary"]], linewidth = 0.5, linetype = "dashed") +
  geom_point(data = tpts, aes(d, row), colour = pal[["orange"]], size = 1.4, alpha = 0.55,
             position = position_nudge(y = -0.2)) +
  geom_errorbar(data = tsum, aes(xmin = L, xmax = U, y = row), orientation = "y", width = 0.18, colour = pal[["orange"]], linewidth = 0.8) +
  geom_point(data = tsum, aes(HL, row), colour = pal[["orange"]], size = 3.2, shape = 18) +
  geom_label(data = tsum, aes(x = -0.7, y = row, label = paste0(txt, "; ", sub)), hjust = 0, vjust = -0.9, size = 2.6, colour = ink[["primary"]],
             fill = ink[["surface"]], border.colour = NA, label.padding = unit(0.12, "lines")) +
  scale_x_continuous(limits = c(-0.7, 7.6), breaks = 0:5) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(title = "(iii) CLR-NS minus scPRINT R0, by truth type",
       x = "log2 AUPRC-ratio difference (positive favours CLR-NS)", y = NULL) +
  theme_clr(10) + theme(panel.grid.major.y = element_blank())

p7c <- (pA | (pB / pC)) + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(caption = paste0("Hodges-Lehmann estimate and exact interval (92.2% coverage, n = 7 datasets); small dots: per-dataset values. ",
                                   "Dashed line: non-inferiority margin -0.15.\n",
                                   sprintf("(iii) is exploratory (not preregistered); on cell-type ChIP-seq scPRINT R0 scores higher than CLR-NS in %d of %d cells and the interval includes 0.", tsum$n_r0[1], tsum$n_cells[1])),
                  theme = theme(plot.caption = element_text(colour = ink[["muted"]], size = 7.5, hjust = 0),
                                plot.background = element_rect(fill = ink[["surface"]], colour = NA)))
save_both(p7c, "fig7c_scfm_comparison", 11, 6.0)
invisible(NULL)
