\name{find_modules}
\alias{find_modules}
\title{Find modules as connected components of the attention graph}
\usage{
find_modules(operator, min_size = 1L)
}
\arguments{
\item{operator}{square numeric matrix (genes x genes), e.g.
\code{ClrAttention$operator}; dimnames used as gene IDs when present.}
\item{min_size}{drop components smaller than this (default 1 keeps all).}
}
\value{named list of character vectors, largest component first.}
\description{
Symmetrizes a thresholded attention operator (\code{(A + t(A)) / 2 > 0}) and
returns its connected components as gene-ID lists. Toy-grade module finder:
on real data the components of a sparse attention graph are a reasonable
first cut, but richer community detection is a separate research question.
Isolated genes are returned as singletons so that no gene silently
disappears from the evaluation.
}
