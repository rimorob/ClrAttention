\name{ClrAttention}
\alias{ClrAttention}
\title{Iterative CLR as attention}
\description{
R6 orchestration of the CLR-as-attention pipeline:
\code{$new(data)$estimate_mi()$calibrate()$build_operator()$diffuse()}.
Stages cache privately; \code{$mi}, \code{$clr_scores}, \code{$operator},
\code{$embedding}, \code{$trajectory} are strict read-only active bindings
that error if the stage has not run. Side-effect methods return
\code{invisible(self)} for chaining.
}
\section{Methods}{
\subsection{\code{$new(data)}}{
Validate and store the genes x samples expression matrix.
}
\subsection{\code{$estimate_mi(bins = "fd", spline_order = 3)}}{
B-spline MI matrix. \code{bins} is \code{"fd"} / \code{"scott"} /
\code{"sturges"} (per-gene adaptive) or a fixed integer.
}
\subsection{\code{$calibrate(method = "normal", combine = "euclidean")}}{
CLR scores. \code{combine} is \code{"euclidean"} (historical) or
\code{"stouffer"} (new); normal method only.
}
\subsection{\code{$build_operator(k = 50, tau = NULL, alpha = 0.5)}}{
Sparse row-stochastic attention operator (top-\code{k} per row, or
threshold \code{tau}).
}
\subsection{\code{$diffuse(steps = 10)}}{
Iterate \code{E <- ((1-alpha) I + alpha A_hat) E}, caching the trajectory.
}
\subsection{\code{$reset_diffusion()}}{Drop the cached trajectory.}
\subsection{\code{$print()}}{Summarize pipeline state.}
}
