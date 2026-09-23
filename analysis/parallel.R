# Shared compute setup for the analysis scripts: byte-code JIT, a foreach
# backend on the local cores (leaving `reserve` cores free for the user),
# a RAM cap on the number of workers, OpenMP threads split across workers,
# memory release between tasks, and a stop file for clean interruption.
#
#   source("analysis/parallel.R")
#   pc <- par_start(mem_gb = 2.5, workers = opt$workers)   # registers doParallel
#   res <- foreach(task = tasks, .inorder = FALSE) %dopar% { ...; par_release(); out }
#   par_stop(pc)
#
# Workers are PSOCK processes (fresh R sessions), not forks: forking a process
# that has already started an OpenMP thread pool is unsafe with some runtimes,
# and PSOCK also works on Windows. Every object in the master's global
# environment is exported to the workers once, at start-up, so call
# par_start() after the shared data are loaded.

suppressPackageStartupMessages({
  library(foreach); library(doParallel); library(parallel)
})

# R's byte-code JIT has defaulted to level 3 (compile all closures and
# top-level loops) since R 3.4; set it explicitly in case R_ENABLE_JIT was
# changed in the environment. Installed packages are byte-compiled at
# install time (DESCRIPTION: ByteCompile: true).
invisible(compiler::enableJIT(3L))

.par_total_ram_gb <- function() {
  sys <- Sys.info()[["sysname"]]
  b <- tryCatch(switch(sys,
    Darwin = as.numeric(system("sysctl -n hw.memsize", intern = TRUE)),
    Linux = 1024 * as.numeric(sub("\\D*(\\d+).*", "\\1",
                                  grep("^MemTotal", readLines("/proc/meminfo"), value = TRUE))),
    Windows = as.numeric(gsub("\\D", "", system2("powershell",
      c("-NoProfile", "-Command",
        "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"),
      stdout = TRUE)[1])),
    NA_real_), error = function(e) NA_real_)
  b / 1024^3
}

# Physical cores where the OS reports them (hyper-threads do not add
# floating-point throughput to this workload), otherwise logical cores.
.par_cores <- function() {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(n) || n < 1L) n <- parallel::detectCores(logical = TRUE)
  if (is.na(n) || n < 1L) 1L else as.integer(n)
}

# Hardware threads per physical core (2 with Intel Hyper-Threading / SMT on).
.par_threads_per_core <- function() {
  p <- suppressWarnings(parallel::detectCores(logical = FALSE))
  l <- parallel::detectCores(logical = TRUE)
  if (is.na(p) || is.na(l) || p < 1L) 1L else max(1L, as.integer(l %/% p))
}

par_stop_file <- function() file.path("jobs", "control", "stop")
par_should_stop <- function() file.exists(par_stop_file())

# Start the backend. mem_gb = peak RAM of one task (measured); the worker
# count is min(cores - reserve, (RAM - reserve_ram_gb - master_gb) / mem_gb).
# reserve_ram_gb is the RAM left to the user: environment CLR_RESERVE_RAM_GB,
# else max(6, 15% of RAM) (about 6 GB on a 32 GB Mac, 19 GB on 128 GB).
# `workers` (e.g. from --workers) overrides the automatic choice.
#
# Hyper-threading: by default workers are capped at one per *physical* core.
# With ht = TRUE (or environment CLR_HT=1) the cap is counted in hardware
# threads instead, so when RAM allows more workers than physical cores the
# extra workers share cores as hyper-threads, and when it does not, the MI
# kernel's OpenMP threads use the spare hyper-threads. Whether a second
# thread per core raises throughput for this workload is measured by
# tools/bench_ht.sh rather than assumed.
par_plan <- function(mem_gb, workers = NULL, reserve = 2L, reserve_ram_gb = NULL,
                     master_gb = 1.5, ht = identical(Sys.getenv("CLR_HT"), "1")) {
  cores <- .par_cores()
  tpc <- if (isTRUE(ht)) .par_threads_per_core() else 1L
  usable <- max(1L, cores - as.integer(reserve))
  ram <- .par_total_ram_gb()
  if (is.null(reserve_ram_gb)) {
    e <- suppressWarnings(as.numeric(Sys.getenv("CLR_RESERVE_RAM_GB")))
    reserve_ram_gb <- if (!is.na(e)) e else if (is.na(ram)) 6 else max(6, 0.15 * ram)
  }
  by_ram <- if (is.na(ram)) usable * tpc else
    max(1L, floor((ram - reserve_ram_gb - master_gb) / mem_gb))
  W <- if (!is.null(workers) && !is.na(workers)) as.integer(workers) else
    as.integer(min(usable * tpc, by_ram))
  # OpenMP threads per worker: the usable cores' hardware threads (physical
  # cores only, unless ht) divided among the workers, so workers x threads
  # never exceeds them.
  omp <- max(1L, (usable * tpc) %/% W)
  list(W = W, omp = omp, cores = cores, tpc = tpc, reserve = as.integer(reserve),
       ram = ram, reserve_ram_gb = reserve_ram_gb, mem_gb = mem_gb)
}

par_start <- function(mem_gb, workers = NULL, reserve = 2L, reserve_ram_gb = NULL,
                      master_gb = 1.5, export = ls(globalenv()), say = message,
                      ht = identical(Sys.getenv("CLR_HT"), "1")) {
  pl <- par_plan(mem_gb, workers, reserve, reserve_ram_gb, master_gb, ht)
  W <- pl$W; omp <- pl$omp
  say(sprintf("compute: %d physical cores (%d reserved), %s GB RAM (%.0f GB reserved, %.1f GB/worker) -> %d workers x %d OpenMP threads%s (JIT level %d)",
              pl$cores, pl$reserve, if (is.na(pl$ram)) "?" else format(round(pl$ram)),
              pl$reserve_ram_gb, mem_gb, W, omp,
              if (pl$tpc > 1L) sprintf(" [hyper-threading: %d threads/core]", pl$tpc) else "",
              compiler::enableJIT(-1L)))
  Sys.setenv(OMP_NUM_THREADS = omp, VECLIB_MAXIMUM_THREADS = 1,
             OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1,
             # glibc (Linux): serve every allocation >= 1 MB with mmap so a
             # freed matrix goes straight back to the OS. By default glibc
             # raises this threshold after the first large free, and later
             # 150 MB matrices then come from the heap and stay resident
             # after gc(). Read at process start, so it applies to the
             # workers launched below. No effect on macOS or Windows.
             MALLOC_MMAP_THRESHOLD_ = 1048576, MALLOC_TRIM_THRESHOLD_ = 1048576)
  cl <- parallel::makePSOCKcluster(W, outfile = "")
  parallel::clusterCall(cl, function(omp) {
    Sys.setenv(OMP_NUM_THREADS = omp, VECLIB_MAXIMUM_THREADS = 1,
               OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
    invisible(compiler::enableJIT(3L))
    suppressPackageStartupMessages({ library(clr); library(Matrix) })
    NULL
  }, omp)
  export <- setdiff(export, c(".Random.seed"))
  if (length(export)) parallel::clusterExport(cl, export, envir = globalenv())
  parallel::clusterCall(cl, function(omp) { assign("OMP_THREADS", omp, globalenv()); NULL }, omp)
  doParallel::registerDoParallel(cl)
  assign("OMP_THREADS", omp, globalenv())
  list(cl = cl, workers = W, omp = omp)
}

par_stop <- function(pc) {
  if (!is.null(pc$cl)) try(parallel::stopCluster(pc$cl), silent = TRUE)
  foreach::registerDoSEQ()
  invisible(NULL)
}

# Return freed memory to the OS promptly: R frees large vectors only when the
# garbage collector runs, and a worker that has just dropped several
# 150 MB matrices would otherwise keep them until the next allocation
# pressure triggers a collection.
par_release <- function() invisible(gc(verbose = FALSE, full = TRUE))
