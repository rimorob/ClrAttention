\name{ModuleTester}
\alias{ModuleTester}
\title{Align found modules to reference regulons}
\description{
R6 class scoring found gene modules against reference modules (e.g.
RegulonDB regulons). \code{$new(reference)} stores reference modules,
\code{$set_found(modules)} stores the modules under test, \code{$align()}
returns one row per reference module with expected precision/recall/F1,
\code{$summary()} aggregates into headline numbers, and the active bindings
expose the assignment state: \code{$alignment} / \code{$assignment} are the
stable Hungarian headline, \code{$soft_alignment} / \code{$soft_assignment}
the soft uncertainty layer, and \code{$soft_uncertainty} per-regulon dummy
mass and entropy.
}
\section{Methods}{
\describe{
\item{\code{$new(reference)}}{Store and validate reference modules (named
list of gene-ID vectors).}
\item{\code{$set_reference(reference)}}{Replace the reference module set
(clears all cached alignments).}
\item{\code{$set_found(modules)}}{Store the found modules to be scored
(clears all cached alignments).}
\item{\code{$align(method, epsilon)}}{Align reference modules to found
modules. \code{method = "hungarian"} (default): hard one-to-one assignment
maximizing total F1 -- the headline recovery metric. Exact recovery scores
exactly 1, and one-to-one matching prevents a single found "hairball" from
earning credit for several regulons. \code{"soft"}: entropic-regularized
probabilistic assignment (log-domain Sinkhorn) as an uncertainty/ambiguity
layer; reported scores are expectations under each regulon's assignment
distribution. \code{"best"}: independent argmax per regulon (many-to-one
allowed; diagnostic). \code{epsilon} is the temperature for \code{"soft"}:
smaller concentrates mass on the best matches (limit: Hungarian), larger
spreads it. Soft results are cached per epsilon value, so re-running with a
new temperature never silently reuses a stale one. Headline
(\code{"hungarian"}) and uncertainty (\code{"soft"}) state are cached
independently, so requesting one never overwrites the other.}
\item{\code{$summary(method, epsilon)}}{Aggregate one cached alignment into
headline numbers. Defaults to the Hungarian headline. Returns a list: n_ref,
n_found, method, epsilon (NA for hard methods), mean_f1, weighted_f1
(regulon-size-weighted), median_f1, recovery_rate (fraction of regulons with
F1 >= 0.5), mean_precision, mean_recall; for \code{method = "soft"}
additionally mean_dummy_mass, mean_entropy, converged and sinkhorn_iters.}
\item{\code{$print()}}{Summarize tester state. Shows the Hungarian headline
only when it is cached; never triggers a computation and never displays
soft state as the headline.}
}
}
\section{Active bindings}{
\describe{
\item{\code{$reference}}{Validated reference module list.}
\item{\code{$found}}{Validated found-module list.}
\item{\code{$alignment}}{Hungarian headline alignment table (computes and
caches if needed). Stable headline state: never reflects a soft or best
alignment.}
\item{\code{$assignment}}{Hungarian headline assignment matrix, reference x
found, 0/1 (computes and caches if needed).}
\item{\code{$best_alignment}, \code{$best_assignment}}{Independent-best
alignment table / assignment matrix (many-to-one allowed; diagnostic).}
\item{\code{$soft_alignment}}{Soft alignment table at the current soft
epsilon (see \code{$align(method = "soft")}); expectations under each
regulon's assignment distribution.}
\item{\code{$soft_assignment}}{Soft assignment matrix, reference x found,
at the current soft epsilon: row i is regulon i's distribution over real
found modules (rows sum to <= 1; the shortfall is dummy/"unexplained"
mass).}
\item{\code{$soft_uncertainty}}{Per-regulon uncertainty at the current soft
epsilon: data.frame with reference, size_ref, dummy_mass (probability
assigned to dummy "unexplained" modules), entropy (nats; over the full row
distribution \emph{including} dummy mass, so unexplained regulons count as
uncertain), top_found and top_prob.}
}
}
\value{\code{$align()}: data.frame with one row per reference module
(\code{reference}, \code{size_ref}, \code{f1}, \code{precision},
\code{recall}, \code{jaccard}, \code{top_found}, \code{top_prob}); under
\code{"soft"} the scores are expectations under the regulon's assignment
distribution, with unassigned (dummy) mass penalizing the expectation.}
