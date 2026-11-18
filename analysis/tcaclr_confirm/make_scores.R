#!/usr/bin/env Rscript
# Production score files for one arm of docs/plans/tcaclr_confirm_plan.md (frozen draft 3f). Reads an expression matrix (first column the gene
# name, then profiles; the file's SHA-256 is checked against the freeze) and writes, in the format of score_confirm.R, genes.txt plus
# packed column-major upper-triangle score files. No truth is read here.
#
#   Rscript analysis/tcaclr_confirm/make_scores.R --mode fit  --expr FILE --sha FILE --out DIR [--threads 8]
#       CLR, aCLR, tcaCLR, |r| (absr) and the hub variable (mi_rowmean: mean MI to the other genes, as stack_m3d.R)
#   Rscript analysis/tcaclr_confirm/make_scores.R --mode pack --expr FILE --sha FILE --out DIR --work DIR
#       PIDC (pidc.f64 of tools/pidc_chunked.jl), GENIE3 (genie3_W.rds) and GRNBoost2 (grnboost2_edges.tsv) of the work directory, each
#       symmetrised by max as in analysis/baselines.R; methods whose file is not yet there are skipped (a later call adds them)
suppressPackageStartupMessages(source("analysis/tcaclr_confirm/tc_lib.R"))
opt <- tc_options(commandArgs(TRUE), list(mode = "", expr = "", sha = "", out = "", work = "", threads = 8L))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
if (!nzchar(opt$expr) || !nzchar(opt$out)) stop("--expr and --out are required")
hs <- tolower(vapply(strsplit(trimws(readLines(opt$sha, warn = FALSE)), "[[:space:]]+"), `[`, "", 1L)); got <- tc_sha256(opt$expr)
say("expression SHA-256 ", got, if (got %in% hs) " matches the freeze" else " DOES NOT MATCH the freeze")
if (!got %in% hs) stop("the expression file does not match its frozen hash")
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
ex <- utils::read.delim(opt$expr, check.names = FALSE); X <- as.matrix(ex[, -1]); rownames(X) <- ex[[1]]; rm(ex); genes <- rownames(X); G <- length(genes)
stopifnot(!anyDuplicated(genes), all(is.finite(X)))
say(sprintf("%s: %d genes x %d profiles; mode %s", opt$expr, G, ncol(X), opt$mode))
gf <- file.path(opt$out, "genes.txt")
if (file.exists(gf)) { stopifnot(identical(readLines(gf), genes)) } else writeLines(genes, gf)
wr <- function(S, name) { tc_write_bin(tc_pack(S), file.path(opt$out, paste0(name, ".scores"))); say("wrote ", name) }

if (opt$mode == "fit") {
  t0 <- Sys.time(); a <- TC$alpha
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)$calibrate(method = "normal", combine = "stouffer")
  M <- f$mi; mi_rm <- (rowSums(M) - diag(M)) / (G - 1); rm(M)
  clr <- tc_clean_sym(f$clr_scores); wr(clr, "CLR"); rm(clr)
  tc_write_bin(mi_rm, file.path(opt$out, "mi_rowmean.scores")); say("wrote mi_rowmean")
  f$build_operator(alpha = a, softmax_keff = TC$keff, softmax_cap = as.integer(TC$cap))
  P <- Matrix::Diagonal(G, 1 - a) + a * Matrix::Matrix(f$operator, sparse = TRUE); f$release()
  Pt <- as.matrix(P); for (s in seq_len(TC$t - 1L)) Pt <- as.matrix(Pt %*% P); rm(P)
  wr(tc_clean_sym((Pt + t(Pt)) / 2), "aCLR"); wr(tc_clean_sym(sk_a3r(Pt)), "tcaCLR"); rm(Pt); invisible(gc(verbose = FALSE))
  wr(tc_absr(X), "absr")
  say(sprintf("fit done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else if (opt$mode == "pack") {
  w <- opt$work; sym_max <- function(W) pmax(W, t(W))
  f <- file.path(w, "pidc.f64")
  if (file.exists(f)) { v <- readBin(f, "double", G * (G - 1) / 2, size = 8L); stopifnot(length(v) == G * (G - 1) / 2, all(is.finite(v)))
    M <- matrix(0, G, G); M[lower.tri(M)] <- v; rm(v); M <- M + t(M); wr(M, "PIDC"); rm(M); invisible(gc(verbose = FALSE)) }   # row-major upper = column-major lower
  f <- file.path(w, "genie3_W.rds")
  if (file.exists(f)) { W <- readRDS(f); stopifnot(identical(rownames(W), genes), identical(colnames(W), genes)); wr(tc_clean_sym(sym_max(W)), "GENIE3"); rm(W); invisible(gc(verbose = FALSE)) }
  f <- file.path(w, "grnboost2_edges.tsv")
  if (file.exists(f)) { e <- utils::read.delim(f, stringsAsFactors = FALSE); ia <- match(e$TF, genes); ib <- match(e$target, genes); stopifnot(!anyNA(ia), !anyNA(ib))
    M <- matrix(0, G, G); M[cbind(ia, ib)] <- e$importance; rm(e); wr(tc_clean_sym(sym_max(M)), "GRNBoost2"); rm(M); invisible(gc(verbose = FALSE)) }
} else stop("--mode must be fit or pack")
