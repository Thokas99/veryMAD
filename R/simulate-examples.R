.simulate_bulk_qc_metadata <- function(n = 250L, seed = 123L) {
  set.seed(seed)
  out <- data.frame(
    sample_id = sprintf("bulk_%03d", seq_len(n)),
    condition = rep(c("control", "treated"), length.out = n),
    library_size = round(stats::rlnorm(n, log(3e7), 0.22)),
    detected_genes = round(stats::rnorm(n, 16000, 1100)),
    mapping_rate = pmin(0.99, pmax(0.65, stats::rnorm(n, 0.92, 0.025))),
    duplication_rate = pmin(0.8, pmax(0.05, stats::rnorm(n, 0.2, 0.045))),
    percent_mitochondrial = pmax(0.5, stats::rnorm(n, 4.5, 1.1))
  )
  out$library_size[1:8] <- round(out$library_size[1:8] * 0.12)
  out$detected_genes[c(1:4, 9:16)] <- round(out$detected_genes[c(1:4, 9:16)] * 0.45)
  out$mapping_rate[c(1:4, 17:24)] <- seq(0.52, 0.69, length.out = 12)
  out$duplication_rate[c(1:4, 25:32)] <- seq(0.5, 0.72, length.out = 12)
  out$percent_mitochondrial[1:4] <- c(15, 17, 19, 21)
  out$library_size[40] <- NA_real_
  out$mapping_rate[41] <- NA_real_
  out$duplication_rate[42] <- NA_real_
  rownames(out) <- out$sample_id
  out
}

.simulate_single_cell_qc_metadata <- function(n = 600L, seed = 123L) {
  set.seed(seed)
  counts <- round(stats::rlnorm(n, log(6000), 0.38))
  out <- data.frame(
    cell_id = sprintf("cell_%04d", seq_len(n)),
    nCount_RNA = counts,
    nFeature_RNA = round(pmax(250, 600 + 18 * sqrt(counts) + stats::rnorm(n, 0, 180))),
    percent.mt = pmax(0.2, stats::rnorm(n, 4.5, 1.8)),
    doublet_score = stats::rbeta(n, 2, 16)
  )
  out$nCount_RNA[1:20] <- round(out$nCount_RNA[1:20] * 0.08)
  out$nFeature_RNA[c(1:8, 21:40)] <- round(out$nFeature_RNA[c(1:8, 21:40)] * 0.2)
  out$percent.mt[c(1:8, 41:60)] <- seq(18, 35, length.out = 28)
  out$doublet_score[c(1:8, 61:80)] <- seq(0.65, 0.98, length.out = 28)
  out$nCount_RNA[100] <- NA_real_
  out$nFeature_RNA[101] <- NA_real_
  out$percent.mt[102] <- NA_real_
  out$doublet_score[103] <- NA_real_
  rownames(out) <- out$cell_id
  out
}
