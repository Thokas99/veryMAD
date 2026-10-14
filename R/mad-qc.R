#' Annotate MAD-based quality-control flags
#'
#' Rows are observations. For each selected metric, estimate limits around the
#' median using `nmads` scaled MADs, then flag the selected tail. Transform
#' metrics before estimation, fit limits by `batch`, or select reference rows
#' with `subset`. Return flags for every observation.
#'
#' @param data A data frame, numeric matrix, or Seurat object. For Seurat,
#'   metrics must already be present in cell metadata.
#' @param metrics A non-empty named character vector with unique numeric
#'   column names and values `"lower"`, `"upper"`, or `"both"`.
#' @param nmads Positive finite number of scaled MADs from the median used for
#'   each limit. The MAD uses the normal-consistency constant 1.4826.
#' @param transform One method for all metrics or named overrides; unspecified
#'   metrics use `"none"`. Supported methods are `"none"`, `"log1p"`, and
#'   `"log10"`.
#' @param output `"annotate"` adds flags to the input; `"report"` returns flags,
#'   limits, and settings without modifying the input.
#' @param min_n Positive integer minimum number of non-missing reference
#'   observations per metric and batch group (default 5). Smaller groups have
#'   unresolved flags.
#' @param zero_mad Policy when a group's MAD is zero: `"na"` (default) leaves
#'   flags unresolved, `"zero"` sets the limit at the median, and `"error"`
#'   stops.
#' @param overwrite Replace existing QC flag columns when annotating?
#' @param verbose Print a summary?
#' @param batch Optional non-missing vector with one label per observation.
#'   Estimate separate limits within each label.
#' @param subset Optional non-missing logical vector selecting reference
#'   observations for estimation. Apply the resulting limits to every row.
#' @return With `output = "annotate"`, returns the input with per-metric and
#'   combined flags; matrix input returns a data frame and Seurat input remains
#'   a Seurat object. With `output = "report"`, returns a `verymad_qc` list
#'   containing `flags`, `thresholds`, and `settings`. Threshold rows report
#'   limits on transformed and raw scales, status, and used, flagged, and
#'   unresolved counts. A batch adds one row per metric and batch.
#' @export
#' @examples
#' qc_metadata <- data.frame(
#'   sample = paste0("sample_", 1:6),
#'   library_size = c(2.4e6, 2.5e6, 2.6e6, 2.7e6, 0.5e6, 2.5e6),
#'   pct_mito = c(4, 5, 6, 5, 7, 30)
#' )
#' annotated <- mad_qc(
#'   qc_metadata,
#'   metrics = c(library_size = "lower", pct_mito = "upper"),
#'   verbose = FALSE
#' )
#' annotated[, c("sample", "mad_qc_outlier")]
#'
#' report <- mad_qc(
#'   qc_metadata,
#'   metrics = c(library_size = "lower", pct_mito = "upper"),
#'   output = "report",
#'   verbose = FALSE
#' )
#' report$thresholds[, c("metric", "n_used", "n_flagged", "n_unresolved")]
mad_qc <- function(data, metrics, nmads = 3, transform = "none",
                   output = c("annotate", "report"), min_n = 5,
                   zero_mad = c("na", "zero", "error"),
                   overwrite = FALSE, verbose = TRUE, batch = NULL,
                   subset = NULL) {
  output <- match.arg(output)
  zero_mad <- match.arg(zero_mad)
  .validate_scalar_flag(overwrite, "overwrite")
  .validate_scalar_flag(verbose, "verbose")
  .validate_positive(nmads, "nmads")
  .validate_positive(min_n, "min_n")
  if (min_n != as.integer(min_n)) stop("`min_n` must be a positive integer.", call. = FALSE)

  is_seurat <- inherits(data, "Seurat")
  original <- data
  if (is_seurat) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Seurat input requires the `SeuratObject` package.", call. = FALSE)
    }
    data <- data[[]]
  } else if (is.matrix(data)) {
    if (!is.numeric(data)) stop("`data` must be a numeric matrix.", call. = FALSE)
  } else if (!is.data.frame(data)) {
    stop("`data` must be a data frame, numeric matrix, or Seurat object.", call. = FALSE)
  }

  metrics <- .validate_metrics(data, metrics)
  strata <- .validate_qc_strata(batch, subset, nrow(data))
  transforms <- .resolve_transforms(transform, names(metrics))
  targets <- c(paste0(names(metrics), "_mad_outlier"), "mad_qc_outlier")
  if (output == "annotate" && !overwrite) {
    data_names <- if (is.matrix(data)) colnames(data) else names(data)
    conflicts <- intersect(targets, data_names)
    if (length(conflicts)) stop(sprintf("Flag column(s) already exist: %s. Use `overwrite = TRUE`.",
      paste(conflicts, collapse = ", ")), call. = FALSE)
  }

  result <- .qc_engine(data, metrics, transforms, nmads, min_n, zero_mad,
                       annotate = output == "annotate" && !is_seurat,
                       batch = strata$batch, reference = strata$reference)
  if (verbose) .inform_qc(result, nrow(data))
  if (output == "report") return(result)

  if (is_seurat) {
    metadata <- result$flags[setdiff(names(result$flags), "id")]
    names(metadata)[seq_along(metrics)] <- paste0(names(metrics), "_mad_outlier")
    rownames(metadata) <- result$flags$id
    return(SeuratObject::AddMetaData(original, metadata = metadata))
  }
  result$annotated
}

.validate_qc_strata <- function(batch, subset, n) {
  if (!is.null(batch)) {
    if (!is.atomic(batch) || !is.null(dim(batch)) || length(batch) != n || anyNA(batch)) {
      stop("`batch` must be a non-missing vector with one value per observation.", call. = FALSE)
    }
    batch <- as.character(batch)
    if (any(batch == "")) stop("`batch` cannot contain empty labels.", call. = FALSE)
  }
  if (!is.null(subset) && (!is.logical(subset) || !is.null(dim(subset)) || length(subset) != n || anyNA(subset))) {
    stop("`subset` must be a non-missing logical vector with one value per observation.", call. = FALSE)
  }
  list(batch = batch, reference = if (is.null(subset)) rep(TRUE, n) else subset)
}

.validate_scalar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
}

.validate_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) stop(sprintf("`%s` must be positive and finite.", name), call. = FALSE)
}

.validate_metrics <- function(data, metrics) {
  if (!is.character(metrics) || !length(metrics) || is.null(names(metrics)) ||
      anyNA(names(metrics)) || any(names(metrics) == "") || anyDuplicated(names(metrics))) {
    stop("`metrics` must be a non-empty named character vector with unique names.", call. = FALSE)
  }
  if (!all(metrics %in% c("lower", "upper", "both"))) stop("Metric directions must be `lower`, `upper`, or `both`.", call. = FALSE)
  metric_names <- if (is.matrix(data)) colnames(data) else names(data)
  if (is.null(metric_names)) metric_names <- character()
  reserved <- intersect(names(metrics), c("id", "mad_qc_outlier"))
  if (length(reserved)) stop(sprintf("Reserved metric name(s): %s.", paste(reserved, collapse = ", ")), call. = FALSE)
  missing <- setdiff(names(metrics), metric_names)
  if (length(missing)) stop(sprintf("Missing metric column(s): %s.", paste(missing, collapse = ", ")), call. = FALSE)
  bad <- names(metrics)[!vapply(names(metrics), function(nm) {
    value <- if (is.matrix(data)) data[, nm] else data[[nm]]
    is.numeric(value)
  }, logical(1))]
  if (length(bad)) stop(sprintf("QC metric column(s) must be numeric: %s.", paste(bad, collapse = ", ")), call. = FALSE)
  invisible(metrics)
}

.resolve_transforms <- function(transform, metrics) {
  out <- stats::setNames(rep("none", length(metrics)), metrics)
  if (!is.character(transform) || !length(transform) || anyNA(transform)) stop("`transform` must be a supported scalar or named partial override.", call. = FALSE)
  if (length(transform) == 1L && is.null(names(transform))) {
    out[] <- transform
  } else {
    if (is.null(names(transform)) || any(names(transform) == "") || anyDuplicated(names(transform))) stop("Named transformations must have unique metric names.", call. = FALSE)
    unknown <- setdiff(names(transform), metrics)
    if (length(unknown)) stop(sprintf("Transformation names must reference selected metrics: %s.", paste(unknown, collapse = ", ")), call. = FALSE)
    out[names(transform)] <- transform
  }
  out[out == "identity"] <- "none"
  if (!all(out %in% c("none", "log1p", "log10"))) stop("Transformations must be `none`, `log1p`, or `log10`.", call. = FALSE)
  out
}

.transform_metric <- function(x, method, metric) {
  if (any(is.nan(x) | is.infinite(x))) stop(sprintf("Metric `%s` contains non-finite values.", metric), call. = FALSE)
  if (method == "log1p" && any(x < 0, na.rm = TRUE)) stop(sprintf("`log1p` requires non-negative values in `%s`.", metric), call. = FALSE)
  if (method == "log10" && any(x <= 0, na.rm = TRUE)) stop(sprintf("`log10` requires positive values in `%s`.", metric), call. = FALSE)
  switch(method, none = x, log1p = log1p(x), log10 = log10(x))
}

.observation_ids <- function(data) {
  ids <- rownames(data)
  if (is.null(ids) || length(ids) != nrow(data) || anyNA(ids) || any(ids == "") || anyDuplicated(ids)) as.character(seq_len(nrow(data))) else ids
}

.qc_engine <- function(data, metrics, transforms, nmads, min_n, zero_mad,
                       annotate = FALSE, batch = NULL, reference = rep(TRUE, nrow(data))) {
  n <- nrow(data); ids <- .observation_ids(data)
  flags <- data.frame(id = ids, stringsAsFactors = FALSE, check.names = FALSE)
  groups <- if (is.null(batch)) list(seq_len(n)) else
    split(seq_len(n), factor(batch, levels = unique(batch)))
  labels <- names(groups)
  if (!length(groups)) {
    groups <- list(integer())
    labels <- if (is.null(batch)) NULL else NA_character_
  }
  thresholds <- list(); insufficient <- character(); threshold_i <- 0L
  for (i in seq_along(metrics)) {
    metric <- names(metrics)[i]; direction <- unname(metrics[[i]])
    raw <- if (is.matrix(data)) data[, metric] else data[[metric]]
    value <- .transform_metric(raw, transforms[[metric]], metric)
    flag <- rep(NA, n)
    for (g in seq_along(groups)) {
      rows <- groups[[g]]
      reference_values <- value[rows[reference[rows]]]
      n_used <- sum(!is.na(reference_values))
      status <- "ok"; centre <- spread <- lower <- upper <- NA_real_
      if (!any(!is.na(value[rows]))) status <- "all_missing" else if (n_used < min_n) {
        status <- "insufficient_n"
        insufficient <- c(insufficient, if (is.null(batch)) metric else paste0(metric, " [", labels[g], "]"))
      } else {
        centre <- stats::median(reference_values, na.rm = TRUE)
        spread <- stats::mad(reference_values, center = centre, constant = 1.4826, na.rm = TRUE)
        if (spread == 0) {
          status <- "zero_mad"
          if (zero_mad == "error") stop(sprintf("MAD is zero for metric `%s`.", metric), call. = FALSE)
          if (zero_mad == "zero") {
            if (direction %in% c("lower", "both")) lower <- centre
            if (direction %in% c("upper", "both")) upper <- centre
          }
        } else {
          if (direction %in% c("lower", "both")) lower <- centre - nmads * spread
          if (direction %in% c("upper", "both")) upper <- centre + nmads * spread
        }
      }
      if (status == "ok" || (status == "zero_mad" && zero_mad == "zero")) {
        present <- rows[!is.na(value[rows])]
        flag[present] <- switch(direction, lower = value[present] < lower, upper = value[present] > upper, both = value[present] < lower | value[present] > upper)
      }
      threshold_i <- threshold_i + 1L
      thresholds[[threshold_i]] <- data.frame(
        metric = metric, direction = direction, transform = transforms[[metric]],
        median = centre, mad = spread, lower = lower, upper = upper,
        lower_raw = if (is.na(lower)) NA_real_ else switch(transforms[[metric]], none = lower, log1p = expm1(lower), log10 = 10^lower),
        upper_raw = if (is.na(upper)) NA_real_ else switch(transforms[[metric]], none = upper, log1p = expm1(upper), log10 = 10^upper),
        status = status, n_used = n_used, n_flagged = sum(flag[rows], na.rm = TRUE),
        n_unresolved = sum(is.na(flag[rows])), stringsAsFactors = FALSE,
        check.names = FALSE
      )
      if (!is.null(batch)) thresholds[[threshold_i]]$batch <- labels[g]
    }
    flags[[metric]] <- flag
  }
  if (length(insufficient)) warning(sprintf("Insufficient usable observations for metric(s): %s.", paste(insufficient, collapse = ", ")), call. = FALSE)
  flags$mad_qc_outlier <- .combine_flags(flags[names(metrics)])
  settings <- list(
    metrics = metrics, transform = transforms, nmads = nmads, constant = 1.4826,
    min_n = min_n, zero_mad = zero_mad)
  if (any(!reference)) settings$n_reference <- sum(reference)
  out <- list(flags = flags, thresholds = do.call(rbind, thresholds), settings = settings)
  class(out) <- c("verymad_qc", "list")
  if (annotate) {
    annotated <- if (is.matrix(data)) as.data.frame(data, check.names = FALSE,
                                                     stringsAsFactors = FALSE) else data
    for (metric in names(metrics)) annotated[[paste0(metric, "_mad_outlier")]] <- flags[[metric]]
    annotated$mad_qc_outlier <- flags$mad_qc_outlier
    out$annotated <- annotated
  }
  out
}

.combine_flags <- function(flags) {
  x <- as.matrix(flags)
  if (!nrow(x)) return(logical())
  out <- rowSums(x, na.rm = TRUE) > 0
  unresolved <- rowSums(is.na(x)) > 0
  out[!out & unresolved] <- NA
  unname(out)
}

.validate_qc_plot_input <- function(data, metrics, batch) {
  if (inherits(data, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Seurat input requires the `SeuratObject` package.", call. = FALSE)
    }
    data <- data[[]]
  }
  if (!is.data.frame(data)) stop("`data` must be annotated metadata or a Seurat object.", call. = FALSE)
  if (!nrow(data)) stop("`data` must contain at least one observation.", call. = FALSE)
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      any(metrics == "") || anyDuplicated(metrics)) {
    stop("`metrics` must be a non-empty character vector with unique metric names.", call. = FALSE)
  }
  flags <- paste0(metrics, "_mad_outlier")
  missing <- setdiff(c(metrics, flags), names(data))
  if (length(missing)) stop(sprintf("Missing metric or QC flag column(s): %s.", paste(missing, collapse = ", ")), call. = FALSE)
  if (!all(vapply(metrics, function(metric) is.numeric(data[[metric]]), logical(1)))) {
    stop("Selected metric columns must be numeric.", call. = FALSE)
  }
  if (!all(vapply(flags, function(flag) is.logical(data[[flag]]), logical(1)))) {
    stop("Selected QC flag columns must be logical.", call. = FALSE)
  }
  list(data = data, batch = .validate_qc_strata(batch, NULL, nrow(data))$batch)
}

#' Plot MAD QC flags by metric and batch
#'
#' Plot raw metric values by batch. Violins and boxplots show distributions;
#' point color marks QC flags. Each metric has its own panel. Requires
#' `ggplot2`.
#'
#' @param data Annotated data frame or Seurat object returned by [mad_qc()].
#' @param metrics Character vector of annotated numeric metric columns.
#' @param batch Optional non-missing vector with one batch label per
#'   observation. If omitted, plot all observations in one group.
#' @return A `ggplot` object with one facet per metric.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- data.frame(metric = c(1:5, 100), sample = rep(c("A", "B"), each = 3))
#'   x <- mad_qc(x, c(metric = "upper"), batch = x$sample, verbose = FALSE)
#'   plot_mad_qc(x, "metric", batch = x$sample)
#' }
plot_mad_qc <- function(data, metrics, batch = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_mad_qc()` requires the `ggplot2` package.", call. = FALSE)
  }
  plot_input <- .validate_qc_plot_input(data, metrics, batch)
  data <- plot_input$data
  batch <- plot_input$batch
  batch_labels <- if (is.null(batch)) rep("All observations", nrow(data)) else batch
  plots <- lapply(seq_along(metrics), function(i) {
    metric <- metrics[i]
    flag <- data[[paste0(metric, "_mad_outlier")]]
    status <- ifelse(is.na(flag), "Unresolved", ifelse(flag, "Flagged", "Not flagged"))
    data.frame(metric = metric, value = data[[metric]],
               status = factor(status, levels = c("Not flagged", "Flagged", "Unresolved")),
               batch = batch_labels, stringsAsFactors = FALSE)
  })
  plot_data <- do.call(rbind, plots)
  plot_data$batch <- factor(plot_data$batch, levels = unique(plot_data$batch))
  observed_data <- plot_data[is.finite(plot_data$value), , drop = FALSE]
  violin_n <- stats::ave(seq_len(nrow(observed_data)), observed_data$metric,
                         observed_data$batch, FUN = length)
  violin_data <- observed_data[violin_n > 1L, , drop = FALSE]
  ggplot2::ggplot(plot_data, ggplot2::aes(x = batch, y = value)) +
    ggplot2::geom_violin(data = violin_data, fill = "grey94", color = "grey65",
                         linewidth = 0.3, trim = TRUE) +
    ggplot2::geom_boxplot(data = observed_data, outlier.shape = NA, width = 0.14,
                          fill = "white", color = "grey35", linewidth = 0.3) +
    ggplot2::geom_jitter(data = observed_data, ggplot2::aes(color = status),
                         width = 0.15, alpha = 0.15, size = 0.35) +
    ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y") +
    ggplot2::scale_color_manual(values = c("Not flagged" = "#697586",
                                           "Flagged" = "#D55E00",
                                           "Unresolved" = "#0072B2")) +
    ggplot2::labs(x = if (is.null(batch)) NULL else "Batch",
                  y = NULL, color = "QC status") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top",
                   strip.background = ggplot2::element_rect(fill = "grey96", color = NA))
}

#' Plot the relationship between two QC metrics
#'
#' Plot paired raw values with per-batch Spearman rho and finite-pair counts.
#' A matching [mad_qc()] report can add raw-scale cutoff lines. Correlations
#' are descriptive; no p-values are calculated. Requires `ggplot2`.
#'
#' @param data Annotated data frame or Seurat object returned by [mad_qc()].
#' @param x,y Character names of two distinct annotated numeric metrics.
#' @param batch Optional non-missing vector with one batch label per
#'   observation. If omitted, plot all observations in one panel.
#' @param report Optional `verymad_qc` report from the same observations and
#'   metrics. Finite raw-scale limits are drawn in each batch panel. A pooled
#'   report is repeated across panels; a batch-specific report requires matching
#'   `batch` labels.
#' @return A `ggplot` object with one panel per batch.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- data.frame(counts = c(1:5, 20), pct_mito = c(2:6, 30))
#'   x <- mad_qc(x, c(counts = "lower", pct_mito = "upper"), verbose = FALSE)
#'   plot_mad_qc_scatter(x, "counts", "pct_mito")
#' }
plot_mad_qc_scatter <- function(data, x, y, batch = NULL, report = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_mad_qc_scatter()` requires the `ggplot2` package.", call. = FALSE)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      !is.character(y) || length(y) != 1L || is.na(y) || !nzchar(y) || identical(x, y)) {
    stop("`x` and `y` must name two distinct metric columns.", call. = FALSE)
  }
  plot_input <- .validate_qc_plot_input(data, c(x, y), batch)
  data <- plot_input$data
  batch <- plot_input$batch
  flags <- paste0(c(x, y), "_mad_outlier")
  batch_labels <- if (is.null(batch)) rep("All observations", nrow(data)) else batch
  flag_x <- data[[flags[1L]]]
  flag_y <- data[[flags[2L]]]
  status <- ifelse(flag_x %in% TRUE | flag_y %in% TRUE, "Flagged",
                   ifelse(is.na(flag_x) | is.na(flag_y), "Unresolved", "Not flagged"))
  plot_data <- data.frame(
    x = data[[x]], y = data[[y]], batch = batch_labels,
    status = factor(status, levels = c("Not flagged", "Flagged", "Unresolved")),
    stringsAsFactors = FALSE
  )
  plot_data$batch <- factor(plot_data$batch, levels = unique(plot_data$batch))
  thresholds <- NULL
  if (!is.null(report)) {
    if (!inherits(report, "verymad_qc") || !is.list(report) ||
        !is.data.frame(report$flags) || !is.data.frame(report$thresholds) ||
        !"id" %in% names(report$flags)) {
      stop("`report` must be a `verymad_qc` report returned by `mad_qc()`.", call. = FALSE)
    }
    if (!identical(.observation_ids(data), report$flags$id)) {
      stop("`report` must use the same observations and order as `data`.", call. = FALSE)
    }
    if (!all(vapply(c(x, y), function(metric) {
      metric %in% names(report$flags) &&
        identical(data[[paste0(metric, "_mad_outlier")]], report$flags[[metric]])
    }, logical(1)))) {
      stop("`report` flags must match the annotations in `data`.", call. = FALSE)
    }
    if (!all(c("metric", "lower_raw", "upper_raw") %in% names(report$thresholds))) {
      stop("`report` does not contain raw-scale thresholds.", call. = FALSE)
    }
    thresholds <- report$thresholds[report$thresholds$metric %in% c(x, y), , drop = FALSE]
    if (!all(c(x, y) %in% thresholds$metric)) {
      stop("`report` must contain thresholds for both plotted metrics.", call. = FALSE)
    }
    batch_levels <- levels(plot_data$batch)
    if ("batch" %in% names(thresholds)) {
      if (is.null(batch) || !setequal(unique(as.character(thresholds$batch)), batch_levels)) {
        stop("Batch-specific `report` thresholds require matching `batch` labels.", call. = FALSE)
      }
      thresholds$batch <- as.character(thresholds$batch)
    } else {
      thresholds <- do.call(rbind, lapply(batch_levels, function(group) {
        out <- thresholds
        out$batch <- group
        out
      }))
    }
    if (anyDuplicated(thresholds[c("metric", "batch")]) ||
        any(table(factor(thresholds$metric, levels = c(x, y)),
                  factor(thresholds$batch, levels = batch_levels)) != 1L)) {
      stop("`report` must contain one threshold row per plotted metric and batch.", call. = FALSE)
    }
    thresholds$batch <- factor(thresholds$batch, levels = batch_levels)
  }
  point_data <- plot_data[is.finite(plot_data$x) & is.finite(plot_data$y), , drop = FALSE]
  correlations <- do.call(rbind, lapply(levels(plot_data$batch), function(group) {
    group_data <- plot_data[plot_data$batch == group, , drop = FALSE]
    finite <- is.finite(group_data$x) & is.finite(group_data$y)
    pair_x <- group_data$x[finite]
    pair_y <- group_data$y[finite]
    rho <- if (length(pair_x) > 1L && length(unique(pair_x)) > 1L &&
               length(unique(pair_y)) > 1L) {
      stats::cor(pair_x, pair_y, method = "spearman")
    } else {
      NA_real_
    }
    label <- if (is.na(rho)) "rho = NA" else sprintf("rho = %.2f", rho)
    data.frame(batch = group, rho = rho, n = length(pair_x),
               label = sprintf("%s\nn = %d", label, length(pair_x)))
  }))
  correlations$batch <- factor(correlations$batch, levels = levels(plot_data$batch))

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y))
  if (!is.null(thresholds)) {
    cutoffs <- rbind(
      data.frame(metric = thresholds$metric, batch = thresholds$batch,
                 cutoff = thresholds$lower_raw),
      data.frame(metric = thresholds$metric, batch = thresholds$batch,
                 cutoff = thresholds$upper_raw)
    )
    cutoffs <- cutoffs[is.finite(cutoffs$cutoff), , drop = FALSE]
    x_cutoffs <- cutoffs[cutoffs$metric == x, , drop = FALSE]
    y_cutoffs <- cutoffs[cutoffs$metric == y, , drop = FALSE]
    if (nrow(x_cutoffs)) {
      plot <- plot + ggplot2::geom_vline(
        data = x_cutoffs,
        ggplot2::aes(xintercept = cutoff, linetype = "MAD cutoff"),
        inherit.aes = FALSE, color = "grey35", linewidth = 0.4
      )
    }
    if (nrow(y_cutoffs)) {
      plot <- plot + ggplot2::geom_hline(
        data = y_cutoffs,
        ggplot2::aes(yintercept = cutoff, linetype = "MAD cutoff"),
        inherit.aes = FALSE, color = "grey35", linewidth = 0.4
      )
    }
    if (nrow(x_cutoffs) || nrow(y_cutoffs)) {
      plot <- plot + ggplot2::scale_linetype_manual(
        values = c("MAD cutoff" = "dashed"), name = NULL
      )
    }
  }
  plot +
    ggplot2::geom_point(data = point_data, ggplot2::aes(color = status),
                        alpha = 0.15, size = 0.4) +
    ggplot2::geom_label(
      data = correlations,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 3.2,
      fill = "white", color = "grey20", linewidth = 0.1
    ) +
    ggplot2::facet_wrap(ggplot2::vars(batch)) +
    ggplot2::scale_color_manual(values = c("Not flagged" = "#697586",
                                           "Flagged" = "#D55E00",
                                           "Unresolved" = "#0072B2")) +
    ggplot2::labs(x = x, y = y, color = "QC status") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top",
                   strip.background = ggplot2::element_rect(fill = "grey96", color = NA))
}

.inform_qc <- function(result, n) {
  metrics <- names(result$settings$metrics)
  combined <- result$flags$mad_qc_outlier
  unsuccessful <- length(unique(result$thresholds$metric[result$thresholds$status != "ok"]))
  cli::cli_inform("veryMAD \u2022 observation QC")
  n_metrics <- length(metrics)
  cli::cli_inform(sprintf("%d observations checked across %d metric%s", n, n_metrics,
                          if (n_metrics == 1L) "" else "s"))
  cli::cli_inform(sprintf("%d observations flagged by at least one metric; %d unresolved",
                           sum(combined %in% TRUE), sum(is.na(combined))))
  cli::cli_inform(sprintf("%d metric%s without estimated thresholds", unsuccessful,
                          if (unsuccessful == 1L) "" else "s"))
}

#' @method print verymad_qc
#' @export
print.verymad_qc <- function(x, ...) {
  metrics <- names(x$settings$metrics)
  status <- table(x$thresholds$status)
  cat("<verymad_qc>\n")
  cat("Observations:", nrow(x$flags), "\n")
  cat("Metrics:", length(metrics), "\n")
  cat("Observations flagged:", sum(x$flags$mad_qc_outlier %in% TRUE), "\n")
  cat("Observations unresolved:", sum(is.na(x$flags$mad_qc_outlier)), "\n")
  cat("Metrics without estimated thresholds:", length(unique(x$thresholds$metric[x$thresholds$status != "ok"])), "\n")
  if (length(status)) cat("Threshold status:", paste(names(status), status, sep = "=", collapse = ", "), "\n")
  invisible(x)
}

utils::globalVariables("cutoff")
