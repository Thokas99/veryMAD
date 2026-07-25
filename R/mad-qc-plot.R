#' Plot an existing MAD QC report
#'
#' @param qc A tidy report created by [mad_qc()]. Calculation groups are read
#'   from the report attributes so group-specific thresholds remain accurate.
#' @param metrics `NULL` for all metrics or a character vector to plot.
#' @param type A lightweight distribution or observation-index view.
#' @param facet_by `NULL` or metadata columns used only to arrange the plot.
#'   Faceting never changes MAD calculations or thresholds.
#' @param show_thresholds Show finite, nonmissing lower and upper thresholds?
#' @param show_legend Show the pass/outlier legend?
#'
#' @details Distribution views combine lightweight boxplots, observations, and finite
#' calculation-group thresholds. Index views show values and finite threshold
#' trajectories against the stable observation index. Missing or infinite
#' thresholds are omitted. No thresholds are recalculated.
#'
#' @return A ggplot object. The plot is returned without being printed.
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' d <- data.frame(sample = paste0("s", 1:20), batch = rep(c("a", "b"), 10),
#'   library_size = c(rnorm(19, 1e6, 1e5), 1e5))
#' qc <- mad_qc(d, c(library_size = "lower"))
#' plot_mad_qc(qc)
#' plot_mad_qc(qc, facet_by = "batch")
plot_mad_qc <- function(qc, metrics = NULL,
                        type = c("distribution", "index"), facet_by = NULL,
                        show_thresholds = TRUE, show_legend = TRUE) {
  .require_namespace("ggplot2", "`plot_mad_qc()`")
  type <- match.arg(type)
  .validate_qc_report(qc)
  qc <- .add_qc_facet_columns(qc, facet_by)
  .arg_flag(show_thresholds, "show_thresholds")
  .arg_flag(show_legend, "show_legend")
  calculation_groups <- attr(qc, "group_by")
  if (is.null(calculation_groups)) calculation_groups <- character()
  missing_groups <- setdiff(calculation_groups, names(qc))
  if (length(missing_groups)) stop("Stored QC calculation groups are missing from `qc`.", call. = FALSE)
  available <- unique(qc$metric)
  if (is.null(metrics)) metrics <- available
  if (!is.character(metrics) || (nrow(qc) && !length(metrics)) || anyNA(metrics) || any(metrics == "") || anyDuplicated(metrics)) {
    stop("`metrics` must be NULL or a character vector of unique metric names.", call. = FALSE)
  }
  missing_metrics <- setdiff(metrics, available)
  if (length(missing_metrics)) stop(sprintf("Unknown metric(s): %s.", paste(missing_metrics, collapse = ", ")), call. = FALSE)
  plot_data <- qc[qc$metric %in% metrics, , drop = FALSE]
  plot_data$metric <- factor(plot_data$metric, levels = metrics)
  plot_data$threshold_group <- .qc_plot_labels(plot_data, calculation_groups, "All observations")
  plot_data$plot_facet <- .qc_plot_labels(plot_data, facet_by, "All observations")
  plot_data$status <- factor(ifelse(is.na(plot_data$is_outlier), "Not evaluated",
    ifelse(plot_data$is_outlier, "Outlier", "Pass")),
    levels = c("Pass", "Not evaluated", "Outlier"))
  if (!nrow(plot_data)) {
    return(ggplot2::ggplot(plot_data) + ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = "QC value"))
  }
  if (type == "distribution") {
    plot <- .plot_qc_distribution(plot_data, facet_by, show_thresholds)
  } else {
    plot <- .plot_qc_index(plot_data, facet_by, show_thresholds)
  }
  plot +
    ggplot2::scale_colour_manual(values = c(Pass = "#8C8C8C", `Not evaluated` = "#CCCCCC", Outlier = "#D55E00"), drop = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = if (show_legend) "right" else "none") +
    ggplot2::labs(colour = "QC status", y = "QC value")
}

.add_qc_facet_columns <- function(qc, facet_by) {
  if (is.null(facet_by)) return(qc)
  if (!is.character(facet_by) || !length(facet_by) || anyNA(facet_by) ||
      any(facet_by == "") || anyDuplicated(facet_by)) {
    stop("`facet_by` must be NULL or one or more unique metadata column names.", call. = FALSE)
  }
  missing <- setdiff(facet_by, names(qc))
  if (!length(missing)) return(qc)
  metadata <- attr(qc, "observation_metadata")
  if (!is.data.frame(metadata) || !all(missing %in% names(metadata))) {
    stop(sprintf("Missing faceting column(s): %s.", paste(missing, collapse = ", ")), call. = FALSE)
  }
  for (name in missing) qc[[name]] <- metadata[[name]][qc$.obs]
  qc
}

.qc_plot_labels <- function(qc, columns, default) {
  if (!length(columns)) return(factor(rep(default, nrow(qc))))
  values <- lapply(columns, function(name) {
    x <- as.character(qc[[name]])
    x[is.na(x)] <- "<missing>"
    paste0(name, "=", x)
  })
  labels <- do.call(paste, c(values, sep = ", "))
  factor(labels, levels = unique(labels))
}

.plot_qc_distribution <- function(data, facet_by, show_thresholds) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = threshold_group, y = value)) +
    ggplot2::geom_boxplot(width = 0.45, outlier.shape = NA, colour = "#555555", fill = NA) +
    ggplot2::geom_jitter(ggplot2::aes(colour = status), width = 0.12, height = 0, alpha = 0.8) +
    ggplot2::labs(x = if (nlevels(data$threshold_group) > 1L) "MAD reference group" else NULL)
  plot <- .facet_qc_plot(plot, facet_by)
  if (show_thresholds) {
    thresholds <- .qc_plot_thresholds(data, by_observation = FALSE)
    if (nrow(thresholds)) {
      plot <- plot + ggplot2::geom_point(data = thresholds,
        ggplot2::aes(x = threshold_group, y = threshold), inherit.aes = FALSE,
        shape = 95, size = 6, colour = "#333333")
    }
  }
  plot
}

.plot_qc_index <- function(data, facet_by, show_thresholds) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .obs, y = value, colour = status)) +
    ggplot2::geom_point(alpha = 0.85) + ggplot2::labs(x = "Observation index")
  plot <- .facet_qc_plot(plot, facet_by)
  if (show_thresholds) {
    thresholds <- .qc_plot_thresholds(data, by_observation = TRUE)
    if (nrow(thresholds)) {
      plot <- plot + ggplot2::geom_step(data = thresholds,
        ggplot2::aes(x = .obs, y = threshold, group = trajectory),
        inherit.aes = FALSE, colour = "#333333", linewidth = 0.45, alpha = 0.8)
    }
  }
  plot
}

.facet_qc_plot <- function(plot, facet_by) {
  if (length(facet_by)) {
    plot + ggplot2::facet_grid(rows = ggplot2::vars(plot_facet),
      cols = ggplot2::vars(metric), scales = "free_y")
  } else {
    plot + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y")
  }
}

.qc_plot_thresholds <- function(data, by_observation) {
  id <- if (by_observation) seq_len(nrow(data)) else !duplicated(data[c("metric", "threshold_group", "plot_facet", "lower", "upper")])
  base <- data[id, c(".obs", "metric", "threshold_group", "plot_facet", "lower", "upper"), drop = FALSE]
  lower <- data.frame(base[c(".obs", "metric", "threshold_group", "plot_facet")], limit = "lower", threshold = base$lower)
  upper <- data.frame(base[c(".obs", "metric", "threshold_group", "plot_facet")], limit = "upper", threshold = base$upper)
  out <- rbind(lower, upper)
  out <- out[is.finite(out$threshold), , drop = FALSE]
  out$trajectory <- interaction(out$metric, out$threshold_group, out$limit, drop = TRUE)
  out
}
