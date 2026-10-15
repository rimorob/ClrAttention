#!/usr/bin/env Rscript
# Figure: what each component of completed CLR attention contributes.
# A/B: Shapley waterfalls for co-membership AUPR and regulon coherence, from
# the most performant configuration (top) down to the bare baseline (|r|, no
# component), each step removing one component's Shapley share (95% interval
# as a whisker). C: drop-one losses (full model minus the model without that
# component), which expose interactions that Shapley averaging hides.
# Input: results/shapley/shapley_<evidence>.csv (analysis/shapley.R).
# Usage (repo root):  Rscript analysis/fig_shapley.R [--evidence SC] [--out results/figures]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(evidence = "SC", out = "results/figures")
i <- 1L
while (i <= length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1L]; i <- i + 2L }
source("analysis/plot_theme.R")
wrap <- function(x, w) paste(strwrap(x, width = w), collapse = "\n")
d <- utils::read.csv(file.path("results/shapley", sprintf("shapley_%s.csv", opt$evidence)), check.names = FALSE)
mlab <- c(aupr = "Co-membership AUPR", coherence = "Median regulon coherence")
short <- c("Dependence: mutual information" = "Mutual information (vs |r|)",
           "Binning: common HG bins" = "Common HG bins (vs per-gene)",
           "Context normalization (CLR)" = "Context normalization (CLR)",
           "Softmax with temperature" = "Softmax with temperature",
           "Depth (iterated attention)" = "Depth (iterated attention)",
           "Removal of 2 global PCs" = "Removal of 2 global PCs†")
cols <- stats::setNames(component_cols[names(short)], short)

waterfall <- function(metric) {
  x <- d[d$metric == metric, ]
  x <- x[order(-x$shapley), ]
  start <- x$v_full[1] - c(0, cumsum(x$shapley))[seq_len(nrow(x))]
  x$xmax <- start; x$xmin <- start - x$shapley
  x$lab <- short[x$component]
  steps <- c("Full model", x$lab, "Baseline: |r| only")
  pos <- stats::setNames(rev(seq_along(steps)), steps)      # top to bottom
  x$yn <- pos[x$lab]
  ends <- data.frame(yn = pos[c("Full model", "Baseline: |r| only")], v = c(x$v_full[1], x$v_none[1]))
  ggplot(x) +
    geom_segment(data = ends, aes(x = v, xend = v, y = yn - 0.38, yend = yn + 0.38),
                 colour = ink[["primary"]], linewidth = 0.9) +
    geom_text(data = ends, aes(x = v, y = yn, label = sprintf("%.3f", v)), hjust = -0.25, size = 2.9,
              colour = ink[["primary"]], fontface = "bold") +
    geom_segment(data = data.frame(x0 = c(x$v_full[1], x$xmin[-nrow(x)]), y0 = c(pos[["Full model"]], x$yn[-nrow(x)]),
                                   y1 = x$yn),
                 aes(x = x0, xend = x0, y = y0 - 0.32, yend = y1 + 0.32), colour = ink[["muted"]],
                 linewidth = 0.3, linetype = "22") +
    geom_rect(aes(xmin = pmin(xmin, xmax), xmax = pmax(xmin, xmax), ymin = yn - 0.32, ymax = yn + 0.32,
                  fill = lab), alpha = 0.95) +
    geom_errorbar(aes(y = yn, xmin = xmax - shapley_hi, xmax = xmax - shapley_lo), orientation = "y",
                  width = 0.22, colour = ink[["secondary"]], linewidth = 0.45) +
    geom_text(aes(x = pmin(xmin, xmax, xmax - shapley_hi), y = yn, label = sprintf("%+.3f", shapley)),
              hjust = 1.15, size = 2.9, colour = ink[["secondary"]]) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(breaks = pos, labels = names(pos), expand = expansion(add = 0.5)) +
    scale_x_continuous(expand = expansion(mult = c(0.32, 0.22)), labels = label_number(accuracy = 0.01)) +
    labs(title = mlab[[metric]], x = NULL, y = NULL) +
    theme_clr() + theme(panel.grid.major.y = element_blank())
}
dd <- d
dd$lab <- factor(short[dd$component], levels = rev(short[d$component[d$metric == "aupr"][order(d$shapley[d$metric == "aupr"])]]))
dd$metric_lab <- mlab[dd$metric]
pc <- ggplot(dd, aes(x = drop_one, y = lab, colour = lab)) +
  geom_vline(xintercept = 0, colour = ink[["muted"]], linewidth = 0.4) +
  geom_errorbar(aes(xmin = drop_lo, xmax = drop_hi), orientation = "y", width = 0, linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text_repel(aes(label = sprintf("%+.3f", drop_one)), size = 2.7, colour = ink[["secondary"]],
                  nudge_y = 0.32, direction = "x", segment.colour = NA, seed = 1) +
  facet_wrap(~metric_lab, scales = "free_x") +
  scale_colour_manual(values = cols, guide = "none") +
  labs(title = "Loss when one component is removed from the full model", x = "Full model minus model without the component",
       y = NULL) +
  theme_clr()
p <- (waterfall("aupr") | waterfall("coherence")) / pc +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    title = "What each component of completed CLR attention contributes",
    subtitle = wrap(sprintf("M3D, 466 experiments; RegulonDB %s evidence. Exact Shapley values over all 64 on/off combinations of six components; whiskers 95%% intervals (gene jackknife for AUPR, regulon bootstrap for coherence).",
                       if (opt$evidence == "SC") "strong/confirmed" else "all"), 150),
    caption = wrap("† PC removal was chosen after inspecting Supplementary Note 7 and is not part of the main pipeline. Depth at t = 9; 'off' alternatives: |r|, per-gene Freedman-Diaconis bins, raw dependence, hard top-10 neighbours, no attention, no PC removal.", 175),
    tag_levels = "a", theme = theme_clr())
f <- file.path(opt$out, sprintf("fig_shapley_%s.png", opt$evidence))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
save_fig(p, f, width = 10, height = 7.2)
save_fig(p, sub("png$", "pdf", f), width = 10, height = 7.2)
cat("wrote", f, "\n")
