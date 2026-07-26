.simulate_bulk_qc_metadata <- function(n = 240L, seed = 123L) {
  set.seed(seed)

  if (n < 40L) {
    stop("`n` must be at least 40 for the example QC distribution.", call. = FALSE)
  }

  out <- data.frame(
    sample_id = sprintf("bulk_%03d", seq_len(n)),
    condition = rep(c("control", "treated"), length.out = n),
    batch = rep(c("batch1", "batch2", "batch3"), length.out = n),
    library_size = round(stats::rlnorm(n, log(4.5e7), 0.32)),
    detected_genes = round(stats::rnorm(n, 18000, 1800)),
    mapping_rate = pmin(0.99, pmax(0.80, stats::rnorm(n, 0.94, 0.025))),
    assigned_rate = pmin(0.95, pmax(0.58, stats::rnorm(n, 0.78, 0.055))),
    rrna_rate = pmin(0.18, pmax(0.002, stats::rlnorm(n, log(0.018), 0.5)))
  )

  poor <- seq_len(8L)
  out$library_size[poor] <- round(seq(3.0e6, 1.2e7, length.out = length(poor)))
  out$detected_genes[poor] <- round(seq(4500, 9500, length.out = length(poor)))
  out$mapping_rate[poor] <- seq(0.52, 0.76, length.out = length(poor))
  out$assigned_rate[poor] <- seq(0.35, 0.56, length.out = length(poor))
  out$rrna_rate[poor] <- seq(0.16, 0.30, length.out = length(poor))

  out$library_size <- pmax(out$library_size, 1L)
  out$detected_genes <- pmax(out$detected_genes, 1L)
  out
}

.simulate_single_cell_qc_metadata <- function(n = 1200L, seed = 123L) {
  set.seed(seed)

  if (n < 80L) {
    stop("`n` must be at least 80 for the example QC distribution.", call. = FALSE)
  }

  out <- data.frame(
    cell_id = sprintf("cell_%04d", seq_len(n)),
    nCount_RNA = round(stats::rlnorm(n, log(9000), 0.55)),
    nFeature_RNA = round(stats::rnorm(n, 3200, 550)),
    percent.mt = pmin(25, pmax(0.2, stats::rlnorm(n, log(4.5), 0.45)))
  )

  low <- seq_len(20L)
  high <- 21L:35L
  mito <- 36L:50L

  out$nCount_RNA[low] <- round(seq(150, 900, length.out = length(low)))
  out$nFeature_RNA[low] <- round(seq(80, 600, length.out = length(low)))
  out$nCount_RNA[high] <- round(seq(55000, 140000, length.out = length(high)))
  out$nFeature_RNA[high] <- round(seq(7000, 12000, length.out = length(high)))
  out$percent.mt[mito] <- seq(18, 35, length.out = length(mito))

  out$nCount_RNA <- pmax(out$nCount_RNA, 1L)
  out$nFeature_RNA <- pmax(out$nFeature_RNA, 1L)
  out
}
