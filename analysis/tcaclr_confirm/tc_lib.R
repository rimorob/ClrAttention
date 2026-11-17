# Scoring and reporting library for the preregistration docs/plans/tcaclr_confirm_plan.md (draft 3a).
# Names: aCLR = CLR with attention (internal: attention), tcaCLR = twice cooked aCLR (internal: a3r).
# Nothing here is tuned; every constant is fixed by the plan (section numbers in the comments). Source from the repo root.
# Consumed-data validation: analysis/tcaclr_confirm/validate_net4.R (DREAM5 net4 yeast, Complex Portal yeast). This library never
# reads HepG2, Kemmeren, CYC2008 or Complex Portal human data by itself; it only reads the files it is pointed at.

suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R")               # ap_fast (the tie convention of the plan)
source("analysis/regulons.R")                # comembership_pairs, comember_matrix, regulon_ap_table
source("analysis/explore_stack/stack_lib.R") # sk_a3r (tcaCLR read-out of P^9), unchanged since d8fdb58
source("analysis/robust_z/robust_z_lib.R")   # rz_hub_stratify

## ---- constants fixed by the plan --------------------------------------------------------------------------------------
TC <- list(
  jack_n = 500L, jack_seed = 20261117L,        # SE_gene (section 5)
  exp_n = 100L, exp_seed = 20261118L,          # SE_exp (section 5)
  cjack_n = 500L, cjack_seed = 20261119L,      # complex-level sensitivity (section 5)
  t = 9L, alpha = 0.5, keff = 10, cap = 50L,   # operator (section 3)
  alpha_test = 0.025, min_effect = 0.005,      # H1 rule (section 5)
  ref_zero = 0, ref_min = 0.005, ref_ecoli = 0.027,   # net4 reference values for the fixed interval sentences
  min_pos_complex = 10L,                       # per-complex AP and breadth test (section 5)
  share_large = 0.80, large_min = 20L,         # headline rule: more than 80% of d carried by complexes of 20 or more genes
  carry_share = 0.80,                          # [k] of the fail sentence: fewest complexes whose contributions reach 80% of d (script definition)
  near_zero = 0.01,                            # H2 labels (section 1)
  neff_test = 27, neff_degraded = 15,          # gate 2 (section 6)
  methods = c("tcaCLR", "aCLR", "CLR", "absr", "GENIE3", "PIDC", "GRNBoost2"),
  comparators = c("aCLR", "CLR", "absr", "GENIE3", "PIDC", "GRNBoost2"),
  refit = c("CLR", "aCLR", "tcaCLR", "absr", "PIDC"),      # methods refitted per experiment replicate (section 5)
  no_refit = c("GENIE3", "GRNBoost2"),                     # gene resampling only; required note below
  readouts = c("aCLR_rank", "tcaCLR_t1", "CLR_P9lin", "infl_r1.5", "infl_r2", "sinkhorn"),   # section 3 (plus rewire_NN, hub_*)
  plan_order = c("main", "main_xsub_excl", "main_ribo_joined", "noribo", "small_3to6", "mid_7to19", "large_ge20",
                 "mid_7to10", "mid_11to19", "pert_neither", "pert_rest"),
  net4_order = c("main", "main_xsub_excl", "main_ribo_joined", "noribo", "nomerge", "size2", "main_paralog_incl",
                 "small_3to6", "mid_7to19", "large_ge20"),
  strata = list(small_3to6 = c(3, 6), mid_7to19 = c(7, 19), large_ge20 = c(20, 1e5), mid_7to10 = c(7, 10), mid_11to19 = c(11, 19)))
TC$display <- c(absr = "|r|")

# fixed text (verbatim from the plan)
TC$note_ge <- "error bars for GENIE3 and GRNBoost2 comparisons include gene resampling only and are narrower than they should be."
TC$claim_a <- paste0("tcaCLR exceeds aCLR in pooled co-membership AUPR on the Complex Portal human complexes in HepG2 Perturb-seq pseudobulk ",
                     "(same laboratory and screen family as the discovery data; cross-cell-type and cross-truth transfer, not cross-technology)")
TC$ecoli_fail <- "the E. coli preregistered test of the same readout did not meet the one-sided 0.025 threshold (M3D p = 0.079, EcoMAC p = 0.051)"
TC$large_clause <- "the gain is carried by complexes of 20 or more genes"
TC$trap <- paste0("tcaCLR minus CLR has a smaller standard error than tcaCLR minus aCLR (on E. coli, k 0.016 to 0.027 against 0.037 to 0.041) ",
                  "and may exclude 0 where tcaCLR minus aCLR does not; this is never a basis for 'tcaCLR beats CLR' (E. coli S3).")
TC$none_map <- stats::setNames(character(0), character(0))

tc_fmt_p <- function(p) if (is.na(p)) "= NA" else if (p < 1e-4) "< 0.0001" else sprintf("= %.4f", p)
tc_name <- function(m) ifelse(m %in% names(TC$display), TC$display[m], m)

tc_say <- function(log, ...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)); cat(m, "\n")
  if (!is.null(log)) cat(m, "\n", file = log, append = TRUE)
}
tc_wcsv <- function(x, path) { tmp <- paste0(path, ".tmp", Sys.getpid()); utils::write.csv(x, tmp, row.names = FALSE); file.rename(tmp, path); invisible(path) }
tc_vmhwm_gb <- function() { l <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  v <- sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", grep("^VmHWM:", l, value = TRUE)); if (length(v)) as.numeric(v) / 1048576 else NA_real_ }
tc_options <- function(args, defaults) {
  opt <- defaults; i <- 1L
  while (i <= length(args)) {
    k <- sub("^--", "", args[i]); if (!k %in% names(opt)) stop("unknown option --", k)
    v <- args[i + 1L]; opt[[k]] <- if (is.character(defaults[[k]])) v else if (is.integer(defaults[[k]])) as.integer(v) else as.numeric(v); i <- i + 2L
  }
  opt
}

## ---- frozen truth files: SHA-256 guard (gate 2: the score stage refuses to run on a mismatch) -----------------------------------
tc_sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
# sha_file: lines "<sha256>  <file name>" (names relative to the directory of sha_file). Returns the check table; stops on any
# missing file, mismatch or malformed line. 'must_cover' names that must appear in the file (members, excluded pairs).
tc_hash_guard <- function(sha_file, must_cover = character(0)) {
  if (!file.exists(sha_file)) stop("hash guard: SHA-256 file not found: ", sha_file)
  l <- readLines(sha_file, warn = FALSE); l <- l[nzchar(trimws(l))]
  p <- strsplit(trimws(l), "[[:space:]]+"); if (!length(p) || any(lengths(p) != 2L)) stop("hash guard: malformed SHA-256 file: ", sha_file)
  tab <- data.frame(file = vapply(p, `[`, "", 2L), expected = tolower(vapply(p, `[`, "", 1L)), stringsAsFactors = FALSE)
  tab$path <- file.path(dirname(sha_file), tab$file); tab$got <- vapply(tab$path, function(f) if (file.exists(f)) tc_sha256(f) else NA_character_, "", USE.NAMES = FALSE)
  tab$ok <- !is.na(tab$got) & tab$got == tab$expected
  miss <- setdiff(basename(must_cover), tab$file); if (length(miss)) stop("hash guard: not covered by the SHA-256 file: ", paste(miss, collapse = ", "))
  if (any(!tab$ok)) stop("hash guard: SHA-256 MISMATCH or missing file (refusing to run): ", paste(tab$file[!tab$ok], collapse = ", "))
  tab
}
tc_read_truth <- function(members, excluded)
  list(mem = utils::read.delim(members, stringsAsFactors = FALSE, quote = "", comment.char = ""),
       excl = utils::read.delim(excluded, stringsAsFactors = FALSE, quote = "", comment.char = ""))

## ---- score files -------------------------------------------------------------------------------------------------------------
# A score directory holds genes.txt (one gene name per line: the gene order of every score file) and one file per method:
#   <name>.scores  int32 header (ndim, dims...) then raw doubles. ndim = 1: the packed upper triangle in column-major order (the
#                  d5_pack / upper.tri convention; length G(G-1)/2); ndim = 2: a full G x G symmetric matrix.
#   <name>.rds     a G x G symmetric matrix (dimnames ignored) or the packed vector.
#   mi_rowmean.scores (optional)  the MI row mean over all measured genes (length G), for the hub-stratified read-out.
# Method names: tcaCLR, aCLR, CLR, absr, GENIE3, PIDC, GRNBoost2 (directed methods already symmetrised by the producer),
# the read-outs aCLR_rank, tcaCLR_t1, CLR_P9lin, infl_r1.5, infl_r2, sinkhorn, rewire_01 ... rewire_10.
tc_read_bin <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  nd <- readBin(con, "integer", 1L, size = 4L); dims <- readBin(con, "integer", nd, size = 4L)
  x <- readBin(con, "double", prod(as.numeric(dims)), size = 8L)
  if (nd > 1L) dim(x) <- dims
  x
}
tc_write_bin <- function(x, path) {
  dims <- if (is.null(dim(x))) length(x) else dim(x); tmp <- paste0(path, ".tmp", Sys.getpid())
  con <- file(tmp, "wb"); writeBin(c(length(dims), dims), con, size = 4L); writeBin(as.numeric(x), con, size = 8L); close(con); file.rename(tmp, path)
}
tc_scores_open <- function(dir) {
  gf <- file.path(dir, "genes.txt"); if (!file.exists(gf)) stop("score directory has no genes.txt: ", dir)
  genes <- readLines(gf, warn = FALSE); genes <- genes[nzchar(genes)]
  if (anyDuplicated(genes)) stop("duplicated gene names in genes.txt")
  f <- list.files(dir, pattern = "\\.(scores|rds)$"); nm <- sub("\\.(scores|rds)$", "", f)
  if (anyDuplicated(nm)) stop("two score files for one method in ", dir)
  list(dir = dir, genes = genes, G = length(genes), files = stats::setNames(file.path(dir, f), nm), methods = nm)
}
# packed column-major upper-triangle vector of a method (validated length)
tc_scores_packed <- function(st, name) {
  p <- st$files[[name]]; if (is.null(p)) stop("no score file for ", name)
  v <- if (grepl("\\.rds$", p)) readRDS(p) else tc_read_bin(p)
  G <- st$G; np <- G * (G - 1) / 2
  if (is.matrix(v)) { if (!identical(dim(v), c(G, G))) stop(name, ": matrix is not G x G"); v <- v[upper.tri(v)] }
  if (length(v) != np) stop(name, ": length ", length(v), " is not G(G-1)/2 = ", np, " for G = ", G)
  v
}
# packed-index of pairs (i < j): column-major upper triangle
tc_pair_index <- function(i, j) { lo <- pmin(i, j); hi <- pmax(i, j); lo + (as.numeric(hi) - 1) * (as.numeric(hi) - 2) / 2 }
tc_unpack <- function(v, G) { S <- matrix(0, G, G); S[upper.tri(S)] <- v; S + t(S) }
tc_pack <- function(S) S[upper.tri(S)]
tc_clean_sym <- function(S) { S[!is.finite(S)] <- 0; diag(S) <- 0; S }

## ---- truth: sets, pair universes (the builders of analysis/complex_confirm/cc_lib.R, generalised to a gene vector) -----------------
tc_key <- function(i, j, G) as.numeric(i) * (G + 1) + as.numeric(j)
tc_excl_keys <- function(excl, scopes, genes) {
  G <- length(genes); e <- excl[excl$scope %in% scopes, ]
  a <- match(e$orf_a, genes); b <- match(e$orf_b, genes); ok <- !is.na(a) & !is.na(b)
  unique(pmin(a, b)[ok] * (G + 1) + pmax(a, b)[ok])
}
tc_sets <- function(mem, variant, genes) {
  x <- mem[mem$variant == variant, ]; if (!nrow(x)) stop("variant not in the member list: ", variant)
  rows <- match(x$orf, genes); if (anyNA(rows)) stop("member genes outside the score universe in variant ", variant, ": ", sum(is.na(rows)))
  lapply(split(rows, x$complex), function(r) sort(unique(r)))
}
tc_pairs <- function(sets, genes, excl_k = NULL) {
  G <- length(genes)
  reg <- list(members = sets, class = stats::setNames(rep("complex", length(sets)), names(sets)))
  p <- comembership_pairs(reg, TC$none_map, genes)
  k <- tc_key(p$i, p$j, G); rm_ <- if (length(excl_k)) k %in% excl_k else rep(FALSE, length(k))
  list(i = p$i[!rm_], j = p$j[!rm_], label = p$label[!rm_], genes_in = p$genes_in, n_removed = sum(rm_))
}
# union-find id of paralog sets, names = ORFs (the 'operon' map of the per-complex AP table)
tc_paralog_map <- function(excl) {
  e <- excl[excl$scope == "all", ]; if (!nrow(e)) return(TC$none_map)
  orf <- unique(c(e$orf_a, e$orf_b)); id <- stats::setNames(seq_along(orf), orf)
  repeat { ch <- FALSE; for (r in seq_len(nrow(e))) { a <- id[[e$orf_a[r]]]; b <- id[[e$orf_b[r]]]; if (a != b) { id[id == max(a, b)] <- min(a, b); ch <- TRUE } }; if (!ch) break }
  stats::setNames(paste0("para", id), names(id))
}
# size stratum of the main truth: positives restricted to complexes of lo..hi genes; pairs that are co-members only of complexes
# outside the stratum leave the universe (cc_stratum)
tc_stratum <- function(main_sets, main_pairs, lo, hi, genes, excl_k) {
  sets_in <- main_sets[lengths(main_sets) >= lo & lengths(main_sets) <= hi]; G <- length(genes)
  ps <- tc_pairs(sets_in, genes, excl_k); pos_k <- tc_key(ps$i, ps$j, G)[ps$label == 1L]
  lab <- as.integer(tc_key(main_pairs$i, main_pairs$j, G) %in% pos_k)
  keep <- !(main_pairs$label == 1L & lab == 0L)
  list(i = main_pairs$i[keep], j = main_pairs$j[keep], label = lab[keep], genes_in = main_pairs$genes_in, n_removed = 0L, sets = sets_in)
}
# split of the main pair universe by whether either gene is perturbed in the panel (plan section 5, secondaries)
tc_pert_split <- function(main_pairs, perturbed_idx) {
  pert <- logical(max(c(main_pairs$i, main_pairs$j, perturbed_idx))); pert[perturbed_idx] <- TRUE
  either <- pert[main_pairs$i] | pert[main_pairs$j]
  mk <- function(keep) list(i = main_pairs$i[keep], j = main_pairs$j[keep], label = main_pairs$label[keep], genes_in = main_pairs$genes_in, n_removed = 0L, sets = main_pairs$sets)
  list(pert_neither = mk(!either), pert_rest = mk(either))
}
tc_truth_struct_hash <- function(p) digest::digest(list(p$i, p$j, p$label), algo = "md5")
# all truths of the chosen table order. Order fixes the gene-jackknife draw order (halves are drawn for all truths in this order).
tc_build_truths <- function(tr, genes, order = c("plan", "net4"), perturbed = NULL) {
  order <- match.arg(order); want <- if (order == "plan") TC$plan_order else TC$net4_order
  mem <- tr$mem; excl <- tr$excl
  k_all <- tc_excl_keys(excl, "all", genes); k_x <- tc_excl_keys(excl, c("all", "main_xsub_excl"), genes)
  S <- function(v) tc_sets(mem, v, genes); T <- list()
  mk <- function(sets, k, role, note = "") { p <- tc_pairs(sets, genes, k); p$sets <- sets; p$role <- role; p$note <- note; p }
  has <- function(v) v %in% mem$variant
  T$main <- mk(S("main"), k_all, "primary")
  if (has("main_xsub_excl")) T$main_xsub_excl <- mk(S("main"), k_x, "secondary", "cross-subunit pairs of both ribosomes removed from the universe")
  if (has("main_ribo_joined")) T$main_ribo_joined <- mk(S("main_ribo_joined"), k_all, "secondary", "40S+60S (and 37S+54S) joined")
  if (has("noribo")) T$noribo <- mk(S("noribo"), k_all, "secondary", "pre-declared ribosome list removed")
  if (order == "net4") {
    if (has("nomerge")) T$nomerge <- mk(S("nomerge"), k_all, "sensitivity")
    if (has("size2")) T$size2 <- mk(S("size2"), k_all, "sensitivity", "minimum size 2")
    T$main_paralog_incl <- mk(S("main"), NULL, "sensitivity", "paralog pairs kept")
  }
  ms <- T$main$sets
  for (st in names(TC$strata)) {
    b <- TC$strata[[st]]
    if (!(st %in% want)) next
    if (!any(lengths(ms) >= b[1] & lengths(ms) <= b[2])) { warning("stratum ", st, " is empty"); next }
    p <- tc_stratum(ms, T$main, b[1], b[2], genes, k_all); p$role <- "stratum"; p$note <- "primary universe; positives restricted to the stratum's complexes"; T[[st]] <- p }
  if (order == "plan" && !is.null(perturbed)) {
    sp <- tc_pert_split(T$main, perturbed)
    for (nm in names(sp)) { sp[[nm]]$role <- "perturbed_split"; sp[[nm]]$note <- "main universe split by whether either gene is perturbed in the panel"; T[[nm]] <- sp[[nm]] } }
  T[intersect(want, names(T))]
}

## ---- average precision with the registered tie convention, and a sort-once version for delete-half replicates ---------------------
# ap_fast is the registered AUPR (ties at the end of their block). A delete-half replicate keeps a subset of the pairs; the order of the
# kept pairs is the restriction of the global order and tie blocks stay contiguous, so one sort per (method, truth) serves all
# replicates. The arithmetic below is the same sequence of operations as ap_fast on the subset (bit-identical; tested).
tc_ap_prep <- function(score, label) {
  score[is.na(score)] <- -Inf
  o <- order(score, decreasing = TRUE, method = "radix"); s <- score[o]
  list(o = o, blk = cumsum(c(TRUE, s[-1L] != s[-length(s)])), lab = label[o])
}
tc_ap_sub <- function(pr, ix) {
  l <- pr$lab[ix]; b <- pr$blk[ix]; n <- length(b)
  if (n < 2L) return(NA_real_)
  last <- c(b[-1L] != b[-n], TRUE); tp <- cumsum(l)[last]; k <- which(last)
  sum(diff(c(0, tp)) * tp / k) / sum(l)
}
# paired delete-half gene jackknife of one method on one truth. halves: list of gene-index vectors. Returns the replicate APs.
tc_jack_method <- function(score, p, halves, G, nc = 1L) {
  pr <- tc_ap_prep(score, p$label); is_ <- p$i[pr$o]; js_ <- p$j[pr$o]
  one <- function(keep) { inh <- logical(G); inh[keep] <- TRUE; tc_ap_sub(pr, which(inh[is_] & inh[js_])) }
  r <- if (nc > 1L) parallel::mclapply(halves, one, mc.cores = nc) else lapply(halves, one)
  if (any(vapply(r, function(x) inherits(x, "try-error") || length(x) != 1L, NA))) stop("a jackknife worker failed")
  unlist(r)
}
tc_draw_halves <- function(truths, n, seed) {   # RNG: set.seed(seed), then all truths in table order (cc_draw_halves)
  set.seed(seed); h <- lapply(truths, function(b) lapply(seq_len(n), function(r) sample(b$genes_in, floor(length(b$genes_in) / 2)))); set.seed(NULL); h
}
tc_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(seed); force(expr)
}
# Experiment halves. n_units units (sorted), 'seq': ONE seed, replicates drawn in sequence (the plan: 100 replicates, seed 20261118);
# 'perrep': seed base + b for replicate b (the net4 convention: 20263000 + b groups, 20264000 + b series). Returns indices into the units.
tc_exp_halves <- function(n_units, n, seed, scheme = c("seq", "perrep")) {
  scheme <- match.arg(scheme)
  if (scheme == "seq") tc_with_seed(seed, lapply(seq_len(n), function(b) sort(sample.int(n_units, floor(n_units / 2)))))
  else lapply(seq_len(n), function(b) tc_with_seed(seed + b, sort(sample.int(n_units, floor(n_units / 2)))))
}

## ---- intervals and tests ------------------------------------------------------------------------------------------------------
tc_skew <- function(x) { m <- mean(x); mean((x - m)^3) / mean((x - m)^2)^1.5 }
tc_interval <- function(d, reps) {
  n_na <- sum(!is.finite(reps)); reps <- reps[is.finite(reps)]   # a replicate with no positive pair has no AP; dropped and counted (never occurs on the net4 truths)
  if (length(reps) < 2L) return(data.frame(d = d, se_gene = NA_real_, gene_lo = NA_real_, gene_hi = NA_real_, pct_lo = NA_real_, pct_hi = NA_real_, skew = NA_real_, n_rep = length(reps), n_na = n_na))
  se <- stats::sd(reps); q <- stats::quantile(reps - mean(reps), c(0.025, 0.975), names = FALSE, type = 7)
  data.frame(d = d, se_gene = se, gene_lo = d - 1.96 * se, gene_hi = d + 1.96 * se, pct_lo = d + q[1], pct_hi = d + q[2], skew = tc_skew(reps), n_rep = length(reps), n_na = n_na)
}
# total SE (primary): sqrt(SE_gene^2 + SE_exp^2). SE_exp NA -> gene only. p one-sided against 0 (normal).
tc_total <- function(d, se_gene, se_exp = NA_real_) {
  se <- if (is.finite(se_exp)) sqrt(se_gene^2 + se_exp^2) else se_gene
  data.frame(se_total = se, total_lo = d - 1.96 * se, total_hi = d + 1.96 * se, z = d / se, p_one_sided = 1 - stats::pnorm(d / se), exp_included = is.finite(se_exp))
}
tc_holm <- function(p, alpha = TC$alpha_test) stats::p.adjust(p, "holm") < alpha

## ---- descriptors (section 5) --------------------------------------------------------------------------------------------------
# credit of every pair to AP: precision at the end of its tie block / n_positives for positives, 0 for negatives (AP = sum of credits)
tc_pair_credit <- function(score, label) {
  score[is.na(score)] <- -Inf; o <- order(score, decreasing = TRUE, method = "radix"); s <- score[o]; l <- label[o]; n <- length(s)
  new <- c(TRUE, s[-1L] != s[-n]); blk <- cumsum(new); last <- c(s[-1L] != s[-n], TRUE)
  prec <- (cumsum(l)[last]) / which(last); cr <- numeric(n); cr[o] <- l * prec[blk] / sum(label); cr
}
# complexes (named index vectors) -> sparse table (pair position, complex, weight): a positive pair shared by m complexes counts 1/m each
tc_pair_complex <- function(sets, p, G) {
  pos <- which(p$label == 1L); pk <- tc_key(p$i[pos], p$j[pos], G); out <- vector("list", length(sets))
  for (c in seq_along(sets)) { m <- sort(sets[[c]]); if (length(m) < 2L) next
    a <- rep(m, times = length(m)); b <- rep(m, each = length(m)); s <- a < b; hit <- match(tc_key(a[s], b[s], G), pk)
    hit <- hit[!is.na(hit)]; if (length(hit)) out[[c]] <- data.frame(pair = pos[hit], cx = c) }
  x <- do.call(rbind, out); x$w <- 1 / as.numeric(table(x$pair)[as.character(x$pair)]); x
}
# per-complex contributions to d = AP(a) - AP(b): d_c = sum over its positive pairs of w * (credit_a - credit_b); sum_c d_c = d
tc_complex_contrib <- function(score_a, score_b, p, sets, G, pc = NULL) {
  if (is.null(pc)) pc <- tc_pair_complex(sets, p, G)
  delta <- tc_pair_credit(score_a, p$label) - tc_pair_credit(score_b, p$label)
  dc <- tapply(pc$w * delta[pc$pair], pc$cx, sum); out <- stats::setNames(numeric(length(sets)), names(sets)); out[as.integer(names(dc))] <- dc
  list(d_c = out, d = sum(delta))
}
# n_eff: inverse Herfindahl index of the shares of the positive pairs by complex (multiply-shared pairs split equally)
tc_neff <- function(pc, n_pos) { s <- tapply(pc$w, pc$cx, sum) / n_pos; 1 / sum(s^2) }
# overlap-counted n_eff: a pair that lies in m complexes counts fully in each (the definition of the anchors' k rule and of the plan's power numbers,
# yeast_complexes_counts.csv column n_eff_overlap_counted; net4 main 15.07). The power gate (section 6) uses this one; the plan does not say which.
tc_neff_overlap <- function(pc) { s <- as.numeric(table(pc$cx)); s <- s / sum(s); 1 / sum(s^2) }
# mandatory descriptors from the per-complex contributions
tc_descriptors <- function(d_c, size_c, d) {
  o <- order(d_c, decreasing = TRUE); cs <- cumsum(d_c[o])
  list(share_top3 = if (d > 0) sum(d_c[o][seq_len(min(3L, length(o)))]) / d else NA_real_,
       top3 = paste(names(d_c)[o][seq_len(min(3L, length(o)))], collapse = ";"),
       share_large = if (d > 0) sum(d_c[size_c >= TC$large_min]) / d else NA_real_,
       k_carry = if (d > 0) unname(which(cs >= TC$carry_share * d)[1]) else NA_integer_)
}
# breadth: per-complex AP (>= 10 positives) of two methods -> counts and the one-sided sign test (a > b); ties dropped
tc_breadth <- function(ap_a, ap_b, alpha = TC$alpha_test) {
  ok <- !is.na(ap_a) & !is.na(ap_b); dd <- ap_a[ok] - ap_b[ok]; ni <- sum(dd > 0); nw <- sum(dd < 0)
  p <- if (ni + nw > 0) stats::binom.test(ni, ni + nw, 0.5, alternative = "greater")$p.value else NA_real_
  data.frame(n_complexes = sum(ok), n_improved = ni, n_worse = nw, n_tied = sum(dd == 0), frac_improved = if (sum(ok)) ni / sum(ok) else NA_real_,
             median_d_ap = if (sum(ok)) stats::median(dd) else NA_real_, mean_d_ap = if (sum(ok)) mean(dd) else NA_real_,
             p_sign_one_sided = p, pass = !is.na(p) && p < alpha)
}
# mean within-complex score (paralog pairs excluded) per complex, and its slope on log size
tc_complex_mean_score <- function(S, sets, genes, opm) {
  op <- opm[genes]
  vapply(sets, function(m) { same <- outer(op[m], op[m], "=="); ut <- upper.tri(same); keep <- ut & (is.na(same) | !same)
    if (!any(keep)) NA_real_ else mean(S[m, m][keep]) }, 0)
}
tc_slope <- function(mean_score, size) { ok <- is.finite(mean_score) & size > 0; if (sum(ok) < 3L) return(c(slope = NA_real_, se = NA_real_))
  f <- summary(stats::lm(mean_score[ok] ~ log(size[ok])))$coefficients; c(slope = f[2, 1], se = f[2, 2]) }

## ---- H2 labels and the script verdict (section 1; no p-values anywhere) ------------------------------------------------------------
tc_h2_labels <- function(d, lo, hi) c(positive = lo > 0, near_zero = abs(d) < TC$near_zero && lo <= 0 && hi >= 0, not_positive = d <= TC$near_zero)
# predicted: large positive, mid near zero, small not positive. d/lo/hi: named by class.
tc_h2_verdict <- function(lab_small, lab_mid, lab_large) {
  held <- c(small = unname(lab_small["not_positive"]), mid = unname(lab_mid["near_zero"]), large = unname(lab_large["positive"]))
  n <- sum(held); verdict <- if (n == 3L) "held" else if (n == 0L) "did not hold" else "partly held"
  list(held = held, n_held = n, verdict = verdict,
       sentence = switch(verdict, held = "H2: the prediction held (large positive, mid near zero, small not positive).",
                         `partly held` = sprintf("H2: the prediction partly held (%d of 3 size classes as predicted: %s).", n, paste(names(held)[held], collapse = ", ")),
                         `did not hold` = "H2: the prediction did not hold."))
}

## ---- H3 (section 1): intersection-union for 'best of the [k] methods run', Holm for individual statements -------------------------------
# tab: data.frame(comparator, d, p) one-sided p of tcaCLR minus comparator. k = comparators run + tcaCLR.
tc_h3 <- function(tab) {
  tab$won_iut <- !is.na(tab$p) & tab$p < TC$alpha_test & tab$d > 0
  tab$holm_p <- stats::p.adjust(tab$p, "holm"); tab$won_holm <- !is.na(tab$holm_p) & tab$holm_p < TC$alpha_test & tab$d > 0
  list(tab = tab, best = all(tab$won_iut), k = nrow(tab) + 1L)
}

## ---- fixed sentences ------------------------------------------------------------------------------------------------------------
# net4 interval sentences (design 'Estimands and reporting rule'), parameterised by the names; with the net4 names they reproduce
# results/complex_confirm/estimate.csv byte for byte (tested). L, U: the total 95% interval of the primary contrast.
tc_interval_sentences <- function(L, U, new = "tcaCLR", ref = "aCLR", context = "complexes", tag = "",
                                  fixed = NULL, kind = "preregistered exploratory estimation; not a confirmatory test") {
  c(S_zero = if (L > TC$ref_zero) sprintf("The total interval excludes no difference: %s exceeded %s on %s (%s).", new, ref, context, kind)
             else if (U < TC$ref_zero) sprintf("The total interval lies below no difference: %s was lower than %s on %s.", new, ref, context)
             else "The total interval includes no difference.",
    S_min = if (L > TC$ref_min) "The whole interval lies above the smallest effect of interest (+0.005)."
            else if (U < TC$ref_min) "Effects at or above the smallest effect of interest (+0.005) are excluded."
            else "The interval includes the smallest effect of interest (+0.005).",
    S_ecoli = if (U < TC$ref_ecoli) "Effects as large as the E. coli estimate (+0.027) are excluded."
              else if (L > TC$ref_ecoli) "The whole interval lies above the E. coli estimate (+0.027)."
              else "The interval includes the E. coli estimate (+0.027).",
    S_fixed = if (is.null(fixed)) "Estimation, not a confirmatory test, a replication or a validation." else fixed)
}
# H1 sentence of Arm A. p, d: the primary contrast; neff_class from the power gate; large_share the share of d in complexes of >= 20 genes.
# Arm B and any non-registered dataset get no claim sentence (tc_claim_sentence is only called for arm A in production mode).
tc_claim_sentence <- function(p, d, neff_class, large_share, claim = TC$claim_a) {
  if (neff_class == "estimate_only")
    return(paste0("The registered test is not run (n_eff below 15: estimate only; the arm is not a test); ", TC$ecoli_fail, ". Estimate: d = ", sprintf("%+.4f", d), "."))
  pass <- is.finite(p) && p < TC$alpha_test && d >= TC$min_effect
  if (pass) paste0(claim, if (is.finite(large_share) && large_share > TC$share_large) paste0("; ", TC$large_clause) else "", "; ", TC$ecoli_fail, ".")
  else paste0("The claim is not made: tcaCLR does not exceed aCLR at one-sided p < 0.025 with d at least +0.005 (d = ", sprintf("%+.4f", d), ", one-sided p ", tc_fmt_p(p), "); ", TC$ecoli_fail, ".")
}
tc_breadth_sentence <- function(br, k_carry, d) {
  if (!is.finite(br$p_sign_one_sided)) return("Breadth test not evaluable (no complex with 10 or more positives and a decided difference).")
  if (br$pass) "tcaCLR improved the median complex."
  else if (is.finite(k_carry)) sprintf("The pooled gain is carried by %d complexes; tcaCLR did not improve the median complex.", k_carry)
  else "The pooled difference is not positive; tcaCLR did not improve the median complex."
}
# secondaries: script-chosen sentence from the total (or, where no refit exists, gene-only) interval; no p-value
tc_secondary_sentence <- function(label, d, lo, hi, basis = "gene + experiment") {
  s <- if (lo > 0) "the interval excludes no difference (tcaCLR exceeded aCLR)" else if (hi < 0) "the interval lies below no difference (tcaCLR was lower than aCLR)" else "the interval includes no difference"
  sprintf("%s: d = %+.4f, %s 95%% interval [%+.4f, %+.4f]; %s.", label, d, basis, lo, hi, s)
}
tc_pert_sentence <- function(d_n, lo_n, hi_n, d_r, lo_r, hi_r) {
  if (lo_n > 0) "Perturbed-gene split: the gain is present in pairs where neither gene is perturbed in the panel."
  else if (lo_r > 0) "Perturbed-gene split: the gain is confined to pairs involving perturbed genes; it is a property of the screen design."
  else "Perturbed-gene split: neither part excludes no difference; no statement on the screen design is made."
}

## ---- refits for the experiment component (CLR, aCLR, tcaCLR, |r|, PIDC) ------------------------------------------------------------
# One fit of CLR, aCLR and tcaCLR on X (measured genes x profiles); the body of analysis/complex_confirm/cc_lib.R::cc_fit_three
# (tested bit-identical). |r| is the absolute Pearson correlation. Returns dense symmetric matrices.
tc_fit_three <- function(X, threads = 8L) {
  G <- nrow(X); alpha <- TC$alpha; tt <- TC$t
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = threads)$calibrate(method = "normal", combine = "stouffer")
  out <- list(CLR = tc_clean_sym(f$clr_scores))
  f$build_operator(alpha = alpha, softmax_keff = TC$keff, softmax_cap = as.integer(TC$cap))
  P <- Matrix::Diagonal(G, 1 - alpha) + alpha * Matrix::Matrix(f$operator, sparse = TRUE)
  Pt <- as.matrix(P); for (s in seq_len(tt - 1L)) Pt <- as.matrix(Pt %*% P)
  f$release()
  out$aCLR <- tc_clean_sym((Pt + t(Pt)) / 2)
  out$tcaCLR <- tc_clean_sym(sk_a3r(Pt))
  rm(Pt, P); invisible(gc(verbose = FALSE)); out
}
tc_absr <- function(X) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; tc_clean_sym(r) }
# PIDC by tools/pidc_chunked.jl (bit-identical to the package). Returns a function (i, j) -> weights for gene indices of X.
tc_pidc_refit <- function(X, threads = 8L, tmpdir = tempdir(), julia = "~/.juliaup/bin/julia", script = "tools/pidc_chunked.jl") {
  G <- nrow(X); tag <- paste0("pidc", Sys.getpid(), "_", as.integer(stats::runif(1) * 1e6))
  fin <- file.path(tmpdir, paste0(tag, ".tsv")); fout <- file.path(tmpdir, paste0(tag, ".f64")); on.exit(unlink(c(fin, fout)))
  utils::write.table(cbind(gene = rownames(X), as.data.frame(X)), fin, sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("gene", paste0("c", seq_len(ncol(X)))))
  rc <- system2("env", c("-i", paste0("HOME=", Sys.getenv("HOME")), "PATH=/usr/bin:/bin", path.expand(julia), "+1.6", "-t", threads, script, fin, fout), stdout = FALSE, stderr = FALSE)
  if (!file.exists(fout) || rc != 0L) stop("PIDC refit failed (rc ", rc, ")")
  v <- readBin(fout, "double", G * (G - 1) / 2, size = 8L)
  function(i, j) { lo <- pmin(i, j); hi <- pmax(i, j); v[(as.numeric(lo) - 1) * G - as.numeric(lo) * (as.numeric(lo) - 1) / 2 + (hi - lo)] }
}
