# Shared ggplot2 styling for every figure in the paper (source this file).
#
# Colours come from a validated categorical palette (colour-blind safe on
# adjacent pairs; for scatter plots use at most the first three slots).
# Each method component keeps ONE colour in every figure, so a reader who
# learns "orange = context normalization" in Fig. 1 can rely on it later.
# Text never takes a series colour; identity is carried by marks plus labels.

suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(patchwork); library(scales) })

pal <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
         magenta = "#e87ba4", green = "#008300", violet = "#6250d6", red = "#e34948")
ink <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#8a8984", grid = "#e7e6e2",
         surface = "#fcfcfb")
# Sequential ramp (one hue, light -> dark) for densities and heatmaps.
seq_blue <- c("#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

# Fixed colour per pipeline component / method family.
component_cols <- c(
  "Dependence: mutual information" = pal[["blue"]],
  "Binning: common HG bins"        = pal[["aqua"]],
  "Context normalization (CLR)"    = pal[["orange"]],
  "Softmax with temperature"       = pal[["violet"]],
  "Depth (iterated attention)"     = pal[["magenta"]],
  "Removal of 2 global PCs"        = pal[["yellow"]])
method_cols <- c(CLR = pal[["orange"]], Attention = pal[["violet"]], Correlation = ink[["muted"]],
                 MI = pal[["blue"]], "Differential expression" = pal[["yellow"]])
platform_cols <- c("Antisense (GPL199)" = pal[["blue"]], "Genome 2.0 (GPL3154)" = pal[["orange"]],
                   "no GEO record" = pal[["aqua"]])

theme_clr <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    theme(
      plot.background = element_rect(fill = ink[["surface"]], colour = NA),
      panel.background = element_rect(fill = ink[["surface"]], colour = NA),
      panel.grid.major = element_line(colour = ink[["grid"]], linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text = element_text(colour = ink[["secondary"]], size = rel(0.85)),
      axis.title = element_text(colour = ink[["secondary"]], size = rel(0.9)),
      axis.ticks = element_line(colour = ink[["grid"]], linewidth = 0.35),
      plot.title = element_text(colour = ink[["primary"]], face = "bold", size = rel(1.05),
                                hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = ink[["secondary"]], size = rel(0.85), hjust = 0,
                                   margin = margin(b = 8)),
      plot.caption = element_text(colour = ink[["muted"]], size = rel(0.75), hjust = 0),
      plot.title.position = "plot", plot.caption.position = "plot",
      strip.text = element_text(colour = ink[["primary"]], face = "bold", hjust = 0,
                                size = rel(0.9), margin = margin(b = 4)),
      legend.title = element_text(colour = ink[["secondary"]], size = rel(0.85)),
      legend.text = element_text(colour = ink[["secondary"]], size = rel(0.8)),
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = rel(1.1), colour = ink[["primary"]]),
      plot.margin = margin(10, 12, 8, 10))
}
save_fig <- function(p, file, width = 7.2, height = 4, dpi = 300) {
  ggsave(file, p, width = width, height = height, dpi = dpi, bg = ink[["surface"]])
  invisible(file)
}
