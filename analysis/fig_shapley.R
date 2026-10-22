#!/usr/bin/env Rscript
# Figure: what each component of completed CLR attention contributes.
# a/b: Shapley waterfalls for co-membership AUPR and median regulon coherence,
# built up from the bare baseline (|r|, no component) to the full model in one
# fixed component order (dependence -> bins -> CLR -> softmax -> depth -> PC
# removal); each bar is that component's Shapley share, whisker its 95%
# interval. c: drop-one losses (full model minus the model without that
# component) next to the Shapley share; a gap between the two marks an
# interaction that Shapley averaging hides. Intervals that cross 0 are grey.
# PC removal is post hoc (chosen after Supplementary Note 7) and drawn hatched.
# The title is derived from the data (two largest main-pipeline Shapley shares
# per metric); the caption lists where the strong/confirmed (SC) and
# all-evidence decompositions disagree (needs both CSVs).
# Input: results/shapley/shapley_<evidence>.csv (analysis/shapley.R).
# Usage (repo root):  Rscript analysis/fig_shapley.R [--evidence SC] [--out results/figures]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(evidence = "SC", out = "results/figures")
i <- 1L
while (i <= length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1L]; i <- i + 2L }
source("analysis/plot_theme.R")
wrap <- function(x, w) paste(strwrap(x, width = w), collapse = "\n")
rd <- function(ev) utils::read.csv(file.path("results/shapley", sprintf("shapley_%s.csv", ev)), check.names = FALSE)
d <- rd(opt$evidence)

comp <- names(component_cols)                     # fixed order, pipeline first, PC removal last
lab <- c("Mutual information instead of |r|",
         "Shared histogram bins for MI (vs per-gene bins)",
         "CLR normalization",
         "Softmax with temperature (vs hard top-10)",
         "Depth: iterated attention (t = 9)",
         "Removal of 2 global PCs\n(post hoc, outside main pipeline)")
word <- c("mutual information", "shared bins", "CLR normalization", "softmax", "depth", "PC removal")
names(lab) <- names(word) <- comp
cols <- stats::setNames(component_cols[comp], lab)
pc <- comp[6]
mlab <- c(aupr = "Co-membership AUPR", coherence = "Median regulon coherence (AUROC)")
dlab <- c(aupr = "Δ co-membership AUPR", coherence = "Δ median regulon coherence (AUROC)")
excl0 <- function(lo, hi) lo > 0 | hi < 0
num <- function(x) ifelse(abs(x) < 5e-4, "0.000", sprintf("%+.3f", x))

# Diagonal hatching inside a rectangle (data units), for the post-hoc component.
hatch <- function(x0, x1, y0, y1, step) {
  w <- step * 1.2; s <- seq(x0 - w, x1, by = step)
  a <- data.frame(x = s, y = y0, xend = s + w, yend = y1)
  clip <- function(px, py, qx, qy) {                # clip segment to x in [x0, x1]
    t0 <- pmax(0, (x0 - px) / (qx - px)); t1 <- pmin(1, (x1 - px) / (qx - px))
    data.frame(x = px + t0 * (qx - px), y = py + t0 * (qy - py),
               xend = px + t1 * (qx - px), yend = py + t1 * (qy - py), ok = t1 > t0)
  }
  h <- clip(a$x, a$y, a$xend, a$yend); h[h$ok, 1:4]
}

waterfall <- function(metric) {
  x <- d[d$metric == metric, ]; x <- x[match(comp, x$component), ]
  x$start <- x$v_none[1] + c(0, cumsum(x$shapley))[seq_len(nrow(x))]
  x$end <- x$start + x$shapley
  x$lab <- factor(lab[x$component], lab)
  steps <- c("Baseline: |r| only", lab, "Full model")
  pos <- stats::setNames(rev(seq_along(steps)), steps)      # top to bottom
  x$yn <- pos[as.character(x$lab)]
  x$ci_col <- ifelse(excl0(x$shapley_lo, x$shapley_hi), ink[["primary"]], ink[["muted"]])
  x$tx <- pmax(x$start, x$end, x$start + x$shapley_hi)
  ends <- data.frame(yn = pos[c("Baseline: |r| only", "Full model")], v = c(x$v_none[1], x$v_full[1]))
  conn <- data.frame(x0 = c(x$v_none[1], x$end), y0 = c(pos[[1]], x$yn), y1 = c(x$yn, pos[["Full model"]]))
  rng <- range(c(ends$v, x$start + x$shapley_lo, x$start + x$shapley_hi))
  pcx <- x[x$component == pc, ]
  hh <- hatch(min(pcx$start, pcx$end), max(pcx$start, pcx$end), pcx$yn - 0.3, pcx$yn + 0.3, diff(rng) / 90)
  main <- x[x$component != pc, ]
  ggplot() +
    geom_segment(data = conn, aes(x = x0, xend = x0, y = y0 - 0.32, yend = y1 + 0.32), colour = ink[["muted"]],
                 linewidth = 0.3, linetype = "22") +
    geom_segment(data = ends, aes(x = v, xend = v, y = yn - 0.38, yend = yn + 0.38),
                 colour = ink[["primary"]], linewidth = 0.9) +
    geom_text(data = ends, aes(x = v, y = yn, label = sprintf("%.3f", v)), hjust = -0.25, size = 2.9,
              colour = ink[["primary"]], fontface = "bold") +
    geom_rect(data = main, aes(xmin = pmin(start, end), xmax = pmax(start, end), ymin = yn - 0.3, ymax = yn + 0.3,
                               fill = lab)) +
    geom_rect(data = pcx, aes(xmin = pmin(start, end), xmax = pmax(start, end), ymin = yn - 0.3, ymax = yn + 0.3),
              fill = grDevices::adjustcolor(pal[["yellow"]], 0.18), colour = pal[["yellow"]], linewidth = 0.4) +
    geom_segment(data = hh, aes(x = x, y = y, xend = xend, yend = yend), colour = pal[["yellow"]], linewidth = 0.35) +
    geom_errorbar(data = x, aes(y = yn, xmin = start + shapley_lo, xmax = start + shapley_hi, colour = ci_col),
                  orientation = "y", width = 0.22, linewidth = 0.45) +
    geom_text(data = x, aes(x = tx, y = yn, label = num(shapley), colour = ci_col), hjust = -0.2, size = 2.8) +
    scale_colour_identity() +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(breaks = pos, labels = names(pos), expand = expansion(add = 0.5)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.16)), labels = label_number(accuracy = 0.01)) +
    labs(title = mlab[[metric]], x = sprintf("%s; bar = Shapley share", mlab[[metric]]), y = NULL) +
    theme_clr() + theme(panel.grid.major.y = element_blank())
}

dropone <- function(metric, ylabels = TRUE) {
  x <- d[d$metric == metric, ]; x <- x[match(comp, x$component), ]
  x$lab <- factor(lab[x$component], rev(lab))
  x$ci_col <- ifelse(excl0(x$drop_lo, x$drop_hi), unname(component_cols[x$component]), ink[["muted"]])
  x$pt_fill <- unname(component_cols[x$component])
  ggplot(x, aes(y = lab)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5,             # post-hoc row, greyed
             fill = grDevices::adjustcolor(pal[["yellow"]], 0.08)) +
    geom_vline(xintercept = 0, colour = ink[["secondary"]], linewidth = 0.4) +
    geom_errorbar(aes(xmin = drop_lo, xmax = drop_hi, colour = ci_col), orientation = "y", width = 0,
                  linewidth = 0.8) +
    geom_point(aes(x = shapley), shape = 23, size = 2.3, fill = ink[["surface"]], colour = ink[["primary"]],
               stroke = 0.5, position = position_nudge(y = -0.22)) +
    geom_point(aes(x = drop_one, fill = pt_fill), shape = 21, size = 2.6, colour = ink[["surface"]], stroke = 0.4) +
    geom_text(aes(x = drop_one, label = num(drop_one)), vjust = -0.9, size = 2.6, colour = ink[["secondary"]]) +
    scale_colour_identity() + scale_fill_identity() +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(title = if (ylabels) "Removal loss (dot, 95% interval) against Shapley share (open diamond)" else " ",
         x = sprintf("%s\n(full model minus model without the component)", dlab[[metric]]),
         y = NULL) +
    theme_clr() + theme(plot.title = element_text(size = rel(0.9)),
                        axis.text.y = if (ylabels) element_text() else element_blank())
}

## ---- data-derived text ---------------------------------------------------------------
top2 <- function(metric) {
  x <- d[d$metric == metric & d$component != pc, ]
  x$component[order(-x$shapley)][1:2]
}
ta <- top2("aupr"); tc <- top2("coherence")
seen <- character(0)
nm <- function(cc) { w <- word[[cc]]; if (cc %in% seen && cc == comp[3]) w <- "CLR"; seen <<- c(seen, cc); w }
title <- sprintf("%s and %s drive the coherence gain; %s and %s drive pooled AUPR",
                 nm(tc[1]), nm(tc[2]), nm(ta[1]), nm(ta[2]))
title <- paste0(toupper(substring(title, 1, 1)), substring(title, 2))
# Interactions: the two largest gaps between drop-one loss and Shapley share (main pipeline).
g <- d[d$component != pc, ]; g$gap <- g$drop_one - g$shapley
g <- g[order(-abs(g$gap)), ][1:2, ]
inter <- paste(sprintf("%s on %s (removal %s vs Shapley %s)", word[g$component],
                       ifelse(g$metric == "aupr", "AUPR", "coherence"), num(g$drop_one), num(g$shapley)),
               collapse = " and ")
inter <- paste0(inter, if (all(g$gap > 0)) ": both matter most when the other components are present" else
  if (all(g$gap < 0)) ": both are partly redundant with the other components" else "")
# SC vs all-evidence disagreements: point estimates of opposite sign, with the 95% interval
# excluding 0 at one level at least.
disagree <- ""
if (all(file.exists(file.path("results/shapley", c("shapley_SC.csv", "shapley_all.csv"))))) {
  a <- rd("SC"); b <- rd("all"); b <- b[match(paste(a$metric, a$component), paste(b$metric, b$component)), ]
  sgn <- function(v, lo, hi) ifelse(lo > 0, 1, ifelse(hi < 0, -1, 0))
  out <- character(0)
  for (k in seq_len(nrow(a))) for (q in c("shapley", "drop_one")) {
    lo <- if (q == "shapley") "shapley_lo" else "drop_lo"; hi <- if (q == "shapley") "shapley_hi" else "drop_hi"
    s1 <- sgn(a[[q]][k], a[[lo]][k], a[[hi]][k]); s2 <- sgn(b[[q]][k], b[[lo]][k], b[[hi]][k])
    if (sign(a[[q]][k]) != sign(b[[q]][k]) && (s1 != 0 || s2 != 0))
      out <- c(out, sprintf("%s %s for %s (SC %s, all %s)", if (q == "shapley") "Shapley" else "drop-one",
                            ifelse(a$metric[k] == "aupr", "AUPR", "coherence"), word[[a$component[k]]],
                            num(a[[q]][k]), num(b[[q]][k])))
  }
  if (length(out)) disagree <- paste0("Where strong/confirmed (SC) and all-evidence results disagree in sign (interval excluding 0 at one level at least): ",
                                      paste(out, collapse = "; "), ". All other estimates agree in sign or are within noise at both levels.")
}
cat("title:", title, "\ninteractions:", inter, "\n", disagree, "\n")

pa <- waterfall("aupr"); pb <- waterfall("coherence") + theme(axis.text.y = element_blank())
pc1 <- dropone("aupr"); pc2 <- dropone("coherence", ylabels = FALSE)
p <- (pa | pb) / (pc1 | pc2) + plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = title,
    subtitle = wrap(sprintf("M3D, 466 experiments; RegulonDB %s evidence. a, b: exact Shapley shares over all 64 on/off combinations of six components, added in pipeline order from the |r| baseline to the full model; whiskers 95%% intervals (gene jackknife for AUPR, regulon bootstrap for coherence), grey where they include 0. c: loss when one component is removed from the full model (dot, interval) against its Shapley share (open diamond); a gap is an interaction. Largest interactions: %s.",
                       if (opt$evidence == "SC") "strong/confirmed" else "all", inter), 175),
    caption = wrap(paste("Hatched: PC removal was chosen after inspecting Supplementary Note 7 and is not part of the main pipeline. 'Off' alternatives: |r|, per-gene Freedman-Diaconis bins, raw dependence, hard top-10 neighbours, no attention, no PC removal.",
                         if (opt$evidence != "SC") disagree else ""), 190),
    tag_levels = list(c("a", "b", "c", "")), theme = theme_clr())
f <- file.path(opt$out, sprintf("fig_shapley_%s.png", opt$evidence))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
h <- if (opt$evidence == "SC") 8 else 8.4
save_fig(p, f, width = 11, height = h)
ggsave(sub("png$", "pdf", f), p, width = 11, height = h, device = grDevices::cairo_pdf, bg = ink[["surface"]])
cat("wrote", f, "\n")
