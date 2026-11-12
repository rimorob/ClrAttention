# EXPLORATORY. Walk with soft MCL inflation (docs/plans/explore_ideas47_plan.md section 2):
# Q_0 = I; Q_s = rownorm((Q_{s-1} P)^r), s = 1..tmax; r = 1 is the K10 walk P^s (up to 1e-16 renormalisation).
# `fun(t, Q)` is called for each t = 1..tmax with the dense row-stochastic Q_t (not symmetrised).
inflate_walk <- function(P, r, tmax, fun) {
  Pd <- as.matrix(P); Q <- NULL
  for (s in seq_len(tmax)) {
    Qn <- if (is.null(Q)) Pd else as.matrix(Q %*% P)
    if (r != 1) Qn <- Qn^r
    Qn <- Qn / rowSums(Qn); Q <- Qn; rm(Qn)
    fun(s, Q) }
  invisible(NULL)
}
