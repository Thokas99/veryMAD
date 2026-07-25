#' Plot an existing MAD QC report
#'
#' @param qc A tidy report created by [mad_qc()].
#' @param metrics `NULL` for all metrics or a character vector to plot.
#' @param type A lightweight distribution or observation-index view.
#' @param group_by `NULL` or one or more grouping column names in `qc`.
#' @param show_thresholds Show finite, nonmissing lower and upper thresholds?
#' @param show_legend Show the pass/outlier legend?
#' Distribution views combine lightweight boxplots, observations, and finite
#' group-specific thresholds. Index views show values and finite threshold
#' trajectories against the stable observation index. Missing or infinite
#' thresholds are omitted. No thresholds are recalculated.
#'
#' @return A ggplot object. The plot is returned without being printed.
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' d <- data.frame(x = c(1, 2, 3, 20))
#' qc <- mad_qc(d, c(x = "upper"))
#' plot_mad_qc(qc)
plot_mad_qc <- function(qc, metrics = NULL,
                        type = c("distribution", "index"), group_by = NULL,
                        show_thresholds = TRUE, show_legend = TRUE) {
  .require_namespace("ggplot2", "`plot_mad_qc()`")
  type <- match.arg(type)
  .validate_qc_report(qc, group_by)
  .arg_flag(show_thresholds, "show_thresholds")
  .arg_flag(show_legend, "show_legend")
  available <- unique(qc$metric)
  if (is.null(metrics)) metrics <- available
  if (!is.character(metrics) || (nrow(qc) && !length(metrics)) || anyNA(metrics) || any(metrics == "") || anyDuplicated(metrics)) {
    stop("`metrics` must be NULL or a character vector of unique metric names.", call. = FALSE)
  }
  missing_metrics <- setdiff(metrics, available)
  if (length(missing_metrics)) stop(sprintf("Unknown metric(s): %s.", paste(missing_metrics, collapse = ", ")), call. = FALSE)
  plot_data <- qc[qc$metric %in% metrics, , drop = FALSE]
  plot_data$metric <- factor(plot_data$metric, levels = metrics)
  plot_data$plot_group <- .qc_plot_groups(plot_data, group_by)
  plot_data$status <- factor(ifelse(is.na(plot_data$is_outlier), "Not evaluated",
    ifelse(plot_data$is_outlier, "Outlier", "Pass")),
    levels = c("Pass", "Not evaluated", "Outlier"))
  if (!nrow(plot_data)) {
    return(ggplot2::ggplot(plot_data) + ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = "QC value"))
  }
  if (type == "distribution") {
    plot <- .plot_qc_distribution(plot_data, show_thresholds)
  } else {
    plot <- .plot_qc_index(plot_data, group_by, show_thresholds)
  }
  plot <- plot +
    ggplot2::scale_colour_manual(values = c(Pass = "#8C8C8C", `Not evaluated` = "#CCCCCC", Outlier = "#D55E00"), drop = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = if (show_legend) "right" else "none") +
    ggplot2::labs(colour = "QC status", y = "QC value")
  plot
}

.qc_plot_groups <- function(qc, group_by) {
  if (!length(group_by)) return(factor(rep("All observations", nrow(qc))))
  values <- lapply(group_by, function(name) {
    x <- as.character(qc[[name]])
    x[is.na(x)] <- "<missing>"
    paste0(name, "=", x)
  })
  labels <- do.call(paste, c(values, sep = ", "))
  factor(labels, levels = unique(labels))
}

.plot_qc_distribution <- function(data, show_thresholds) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = plot_group, y = value)) +
    ggplot2::geom_boxplot(width = 0.45, outlier.shape = NA, colour = "#555555", fill = NA) +
    ggplot2::geom_jitter(ggplot2::aes(colour = status), width = 0.12, height = 0, alpha = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y") +
    ggplot2::labs(x = if (nlevels(data$plot_group) > 1L) "Group" else NULL)
  if (show_thresholds) {
    thresholds <- .qc_plot_thresholds(data, by_observation = FALSE)
    if (nrow(thresholds)) {
      plot <- plot + ggplot2::geom_point(data = thresholds,
        ggplot2::aes(x = plot_group, y = threshold), inherit.aes = FALSE,
        shape = 95, size = 6, colour = "#333333")
    }
  }
  plot
}

.plot_qc_index <- function(data, group_by, show_thresholds) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .obs, y = value, colour = status)) +
    ggplot2::geom_point(alpha = 0.85) + ggplot2::labs(x = "Observation index")
  if (length(group_by)) {
    plot <- plot + ggplot2::facet_grid(rows = ggplot2::vars(plot_group),
      cols = ggplot2::vars(metric), scales = "free_y")
  } else {
    plot <- plot + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y")
  }
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

.qc_plot_thresholds <- function(data, by_observation) {
  id <- if (by_observation) seq_len(nrow(data)) else !duplicated(data[c("metric", "plot_group", "lower", "upper")])
  base <- data[id, c(".obs", "metric", "plot_group", "lower", "upper"), drop = FALSE]
  lower <- data.frame(base[c(".obs", "metric", "plot_group")], limit = "lower", threshold = base$lower)
  upper <- data.frame(base[c(".obs", "metric", "plot_group")], limit = "upper", threshold = base$upper)
  out <- rbind(lower, upper)
  out <- out[is.finite(out$threshold), , drop = FALSE]
  out$trajectory <- interaction(out$metric, out$plot_group, out$limit, drop = TRUE)
  out
}
