# Run: Rscript benchmarks/benchmark-mad-qc.R
# Repeat with VERYMAD_BENCH_ITERATIONS (default: 1).
# Allocation totals include profiled events of at least 1 MiB; they are not peak
# memory. Values are NA when R lacks memory profiling.

if (!requireNamespace("veryMAD", quietly = TRUE)) {
  stop("Install the current veryMAD package before running this benchmark.")
}

iterations <- suppressWarnings(as.integer(Sys.getenv("VERYMAD_BENCH_ITERATIONS", "1")))
if (length(iterations) != 1L || is.na(iterations) || iterations < 1L) {
  stop("VERYMAD_BENCH_ITERATIONS must be a positive integer.")
}

measure <- function(fun) {
  elapsed <- numeric(iterations)
  allocated_mib <- rep(NA_real_, iterations)
  for (i in seq_len(iterations)) {
    invisible(gc())
    if (isTRUE(capabilities("profmem"))) {
      profile <- tempfile()
      utils::Rprofmem(profile, threshold = 1024^2)
      elapsed[i] <- unname(system.time(invisible(fun()))["elapsed"])
      utils::Rprofmem(NULL)
      sizes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(profile, warn = FALSE))))
      allocated_mib[i] <- sum(sizes, na.rm = TRUE) / 1024^2
      unlink(profile)
    } else {
      elapsed[i] <- unname(system.time(invisible(fun()))["elapsed"])
    }
  }
  c(elapsed_seconds = median(elapsed), allocated_MiB_over_1MiB = median(allocated_mib))
}

make_metadata <- function(n) {
  set.seed(600 + n)
  x <- data.frame(
    nCount_RNA = stats::rpois(n, 5000),
    nFeature_RNA = stats::rpois(n, 1500),
    percent.mt = stats::rbeta(n, 2, 18) * 100,
    nCount_peaks = stats::rpois(n, 6000),
    TSS.enrichment = stats::rgamma(n, 3, scale = 2),
    nucleosome_signal = stats::rgamma(n, 2, scale = 0.7),
    nCount_ADT = stats::rpois(n, 500),
    percent.ribo = stats::rbeta(n, 3, 15) * 100,
    log10GenesPerUMI = stats::rnorm(n, 0.7, 0.08),
    doublet_score = stats::rbeta(n, 2, 12)
  )
  rownames(x) <- paste0("cell_", seq_len(n))
  x
}

combine_reference <- function(flags) {
  if (!nrow(flags)) return(logical())
  vapply(seq_len(nrow(flags)), function(i) {
    x <- unlist(flags[i, ], use.names = FALSE)
    if (any(x %in% TRUE)) TRUE else if (anyNA(x)) NA else FALSE
  }, logical(1))
}

directions <- c(
  nCount_RNA = "lower",
  nFeature_RNA = "lower",
  percent.mt = "upper",
  nCount_peaks = "lower",
  TSS.enrichment = "lower",
  nucleosome_signal = "upper",
  nCount_ADT = "lower",
  percent.ribo = "upper",
  log10GenesPerUMI = "lower",
  doublet_score = "upper"
)
combine_flags <- getFromNamespace(".combine_flags", "veryMAD")
check_flags <- data.frame(
  lower = c(FALSE, TRUE, NA, FALSE),
  upper = c(NA, FALSE, FALSE, FALSE)
)
stopifnot(identical(combine_flags(check_flags), combine_reference(check_flags)))
records <- list()

add_record <- function(n, n_metrics, task, result) {
  records[[length(records) + 1L]] <<- data.frame(
    observations = n,
    metrics = n_metrics,
    task = task,
    elapsed_seconds = unname(result[["elapsed_seconds"]]),
    allocated_MiB_over_1MiB = unname(result[["allocated_MiB_over_1MiB"]])
  )
}

# shortcut: caps the row-wise reference at 100k rows; raise when time and memory permit.
baseline_limit <- 100000L
for (n in c(10000L, 100000L, 500000L, 1000000L)) {
  metadata <- make_metadata(n)
  for (k in c(3L, 5L, 10L)) {
    metrics <- directions[seq_len(k)]
    transformable <- intersect(names(metrics), c(
      "nCount_RNA", "nFeature_RNA", "nCount_peaks", "nCount_ADT"
    ))
    transforms <- stats::setNames(rep("log1p", length(transformable)), transformable)

    for (output in c("annotate", "report")) {
      result <- measure(function() veryMAD::mad_qc(
        metadata,
        metrics = metrics,
        transform = transforms,
        output = output,
        verbose = FALSE
      ))
      add_record(n, k, paste0("mad_qc_", output), result)
    }

    flags <- as.data.frame(matrix(
      sample(c(FALSE, TRUE, NA), n * k, replace = TRUE,
             prob = c(0.96, 0.02, 0.02)),
      nrow = n,
      ncol = k
    ))
    names(flags) <- paste0("metric_", seq_len(k))
    result <- measure(function() combine_flags(flags))
    add_record(n, k, "combined_vectorized", result)

    if (n <= baseline_limit) {
      result <- measure(function() combine_reference(flags))
      add_record(n, k, "combined_rowwise_reference", result)
    }
  }
  rm(metadata)
  invisible(gc())
}

utils::write.csv(do.call(rbind, records), row.names = FALSE)
