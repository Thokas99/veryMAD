#' Plot an existing MAD QC report
#'
#' @param qc A tidy report created by [mad_qc()]. Calculation groups are read
#'   from report attributes when group-specific thresholds must be displayed.
#' @param metrics `NULL` for all metrics or a character vector to plot.
#' @param type A distribution or observation-index view.
#' @param facet_by `NULL` or metadata columns used only to arrange the plot.
#'   Faceting never changes MAD calculations or thresholds.
#' @param show_thresholds Show active lower and upper thresholds?
#' @param show_legend Show legends?
#'
#' @details Missing QC values are shown as blue triangles at an artificial
#' panel floor. The report itself remains unchanged. Lower thresholds are
#' dashed and upper thresholds are dot-dashed. No thresholds are recalculated.
#'
#' @return A ggplot object. The plot is returned without being printed.
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' d <- data.frame(library_size = c(rnorm(19, 1e6, 1e5), NA))
#' qc <- mad_qc(d, c(library_size = "lower"))
#' plot_mad_qc(qc)
plot_mad_qc <- function(qc, metrics = NULL,
                        type = c("distribution", "index"), facet_by = NULL,
                        show_thresholds = TRUE, show_legend = TRUE) {
  .require_namespace("ggplot2", "`plot_mad_qc()`")
  type <- match.arg(type)
  .validate_qc_report(qc)
  .arg_flag(show_thresholds, "show_thresholds")
  .arg_flag(show_legend, "show_legend")
  calculation_groups <- attr(qc, "group_by")
  if (is.null(calculation_groups)) calculation_groups <- character()
  qc <- .add_qc_facet_columns(qc, facet_by)
  missing_groups <- setdiff(calculation_groups, names(qc))
  if (length(missing_groups)) stop("Stored QC calculation groups are missing from `qc`.", call. = FALSE)
  available <- unique(qc$metric)
  if (is.null(metrics)) metrics <- available
  if (!is.character(metrics) || (nrow(qc) && !length(metrics)) || anyNA(metrics) ||
      any(metrics == "") || anyDuplicated(metrics)) {
    stop("`metrics` must be NULL or a character vector of unique metric names.", call. = FALSE)
  }
  missing_metrics <- setdiff(metrics, available)
  if (length(missing_metrics)) stop(sprintf("Unknown metric(s): %s.", paste(missing_metrics, collapse = ", ")), call. = FALSE)
  data <- qc[qc$metric %in% metrics, , drop = FALSE]
  data$metric <- factor(data$metric, levels = metrics)
  data$threshold_group <- .qc_plot_labels(data, calculation_groups, "All observations")
  data$plot_facet <- .qc_plot_labels(data, facet_by, "All observations")
  data$x_group <- as.integer(data$threshold_group)
  data$status <- factor(ifelse(is.na(data$is_outlier), "Not evaluated",
    ifelse(data$is_outlier, "Outlier", "Pass")),
    levels = c("Pass", "Not evaluated", "Outlier"))
  data$plot_floor <- .qc_plot_floor(data)
  if (!nrow(data)) {
    return(ggplot2::ggplot(data) + ggplot2::theme_minimal() + ggplot2::labs(x = NULL, y = "QC value"))
  }
  plot <- if (type == "distribution") {
    .plot_qc_distribution(data, facet_by, show_thresholds)
  } else {
    .plot_qc_index(data, facet_by, show_thresholds)
  }
  plot <- plot +
    ggplot2::scale_colour_manual(values = c(Pass = "#009E73", `Not evaluated` = "#0072B2", Outlier = "#D55E00"), drop = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (show_legend) "right" else "none") +
    ggplot2::labs(colour = "QC status", linetype = "Threshold", y = "QC value")
  if (anyNA(data$value)) {
    plot <- plot + ggplot2::labs(caption = "Blue triangle at panel floor = missing/not evaluated")
  }
  plot$data <- data
  plot
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

.qc_plot_floor <- function(data) {
  if (!nrow(data)) return(numeric())
  panels <- interaction(data$metric, data$plot_facet, drop = TRUE)
  floors <- vapply(split(data$value, panels), function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(0)
    span <- diff(range(x))
    min(x) - max(span * 0.08, max(abs(x)) * 0.05, 1e-8)
  }, numeric(1))
  unname(floors[as.character(panels)])
}

.plot_qc_distribution <- function(data, facet_by, show_thresholds) {
  finite <- data[!is.na(data$value), , drop = FALSE]
  missing <- data[is.na(data$value), , drop = FALSE]
  plot <- ggplot2::ggplot() +
    ggplot2::geom_boxplot(data = finite, ggplot2::aes(x = x_group, y = value, group = x_group),
      width = 0.45, outlier.shape = NA, colour = "#555555", fill = NA)
  if (show_thresholds) plot <- .add_distribution_thresholds(plot, data)
  plot <- plot + ggplot2::geom_jitter(data = finite,
    ggplot2::aes(x = x_group, y = value, colour = status), width = 0.12, height = 0, alpha = 0.8)
  if (nrow(missing)) plot <- plot + ggplot2::geom_point(data = missing,
    ggplot2::aes(x = x_group, y = plot_floor, colour = status), shape = 17, size = 2.4)
  plot <- .facet_qc_plot(plot, facet_by)
  plot + ggplot2::scale_x_continuous(breaks = seq_along(levels(data$threshold_group)),
    labels = levels(data$threshold_group)) +
    ggplot2::labs(x = if (nlevels(data$threshold_group) > 1L) "MAD reference group" else NULL)
}

.add_distribution_thresholds <- function(plot, data) {
  thresholds <- .qc_plot_thresholds(data, by_observation = FALSE)
  if (!nrow(thresholds)) return(plot)
  thresholds$x_group <- as.integer(thresholds$threshold_group)
  plot + ggplot2::geom_segment(data = thresholds,
    ggplot2::aes(x = x_group - 0.32, xend = x_group + 0.32,
      y = threshold, yend = threshold, linetype = limit),
    inherit.aes = FALSE, colour = "#333333", linewidth = 0.55) +
    ggplot2::scale_linetype_manual(values = c(`Lower MAD threshold` = "dashed",
      `Upper MAD threshold` = "dotdash"), drop = FALSE)
}

.plot_qc_index <- function(data, facet_by, show_thresholds) {
  finite <- data[!is.na(data$value), , drop = FALSE]
  missing <- data[is.na(data$value), , drop = FALSE]
  plot <- ggplot2::ggplot()
  if (show_thresholds) {
    thresholds <- .qc_plot_thresholds(data, by_observation = TRUE)
    if (nrow(thresholds)) {
      plot <- plot + ggplot2::geom_step(data = thresholds,
        ggplot2::aes(x = .obs, y = threshold, group = trajectory, linetype = limit),
        colour = "#333333", linewidth = 0.5) +
        ggplot2::scale_linetype_manual(values = c(`Lower MAD threshold` = "dashed",
          `Upper MAD threshold` = "dotdash"), drop = FALSE)
    }
  }
  plot <- plot + ggplot2::geom_point(data = finite,
    ggplot2::aes(x = .obs, y = value, colour = status), alpha = 0.85)
  if (nrow(missing)) plot <- plot + ggplot2::geom_point(data = missing,
    ggplot2::aes(x = .obs, y = plot_floor, colour = status), shape = 17, size = 2.4)
  .facet_qc_plot(plot, facet_by) + ggplot2::labs(x = "Observation index")
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
  id <- if (by_observation) seq_len(nrow(data)) else
    !duplicated(data[c("metric", "threshold_group", "plot_facet", "lower", "upper")])
  base <- data[id, c(".obs", "metric", "threshold_group", "plot_facet", "lower", "upper"), drop = FALSE]
  lower <- data.frame(base[c(".obs", "metric", "threshold_group", "plot_facet")],
    limit = "Lower MAD threshold", threshold = base$lower)
  upper <- data.frame(base[c(".obs", "metric", "threshold_group", "plot_facet")],
    limit = "Upper MAD threshold", threshold = base$upper)
  out <- rbind(lower, upper)
  out <- out[is.finite(out$threshold), , drop = FALSE]
  out$limit <- factor(out$limit, levels = c("Lower MAD threshold", "Upper MAD threshold"))
  out$trajectory <- interaction(out$metric, out$threshold_group, out$plot_facet, out$limit, drop = TRUE)
  out
}
