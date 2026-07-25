.simulate_bulk_qc_metadata <- function(n = 150L, seed = 123L) {
  if (length(n) != 1L || n < 100L) stop("`n` must be at least 100.", call. = FALSE)
  set.seed(seed)
  library_size <- round(pmin(1.2e8, pmax(1e7, stats::rlnorm(n, log(4.5e7), 0.38))))
  depth_millions <- library_size / 1e6
  detected_genes <- round(6500 + 13000 * (1 - exp(-depth_millions / 24)) + stats::rnorm(n, 0, 500))
  batch <- rep(c("Batch1", "Batch2", "Batch3"), length.out = n)
  batch_shift <- c(Batch1 = 0, Batch2 = -0.012, Batch3 = 0.008)[batch]
  duplication_rate <- pmin(0.68, pmax(0.15,
    0.42 - 0.055 * log10(depth_millions) + stats::rnorm(n, 0, 0.065)))
  percent_mitochondrial <- pmin(10, pmax(0.5, 0.5 + 9 * stats::rbeta(n, 2, 5)))
  mapping_rate <- pmin(0.985, pmax(0.76,
    0.92 + batch_shift - 0.0025 * pmax(percent_mitochondrial - 5, 0) + stats::rnorm(n, 0, 0.025)))
  out <- data.frame(
    sample_id = sprintf("bulk_%03d", seq_len(n)),
    condition = rep(c("Control", "Treatment"), length.out = n),
    batch = batch,
    library_size = library_size,
    detected_genes = detected_genes,
    mapping_rate = mapping_rate,
    duplication_rate = duplication_rate,
    percent_mitochondrial = percent_mitochondrial
  )
  out$library_size[1:10] <- round(seq(1e6, 6e6, length.out = 10))
  out$detected_genes[1:10] <- round(seq(3000, 7000, length.out = 10))
  out$mapping_rate[1:10] <- seq(0.4, 0.65, length.out = 10)
  out$duplication_rate[1:10] <- seq(0.78, 0.95, length.out = 10)
  out$percent_mitochondrial[1:10] <- seq(15, 35, length.out = 10)
  out$library_size[11:14] <- seq(2e6, 5e6, length.out = 4)
  out$detected_genes[15:18] <- seq(4000, 8500, length.out = 4)
  out$mapping_rate[19:22] <- seq(0.48, 0.68, length.out = 4)
  out$duplication_rate[23:26] <- seq(0.78, 0.92, length.out = 4)
  out$percent_mitochondrial[27:30] <- seq(14, 30, length.out = 4)
  out$library_size[41] <- NA_real_
  out$detected_genes[42] <- NA_real_
  out$mapping_rate[43] <- NA_real_
  out$duplication_rate[44] <- NA_real_
  out$percent_mitochondrial[45] <- NA_real_
  rownames(out) <- out$sample_id
  out
}

.simulate_single_cell_qc_metadata <- function(n = 1000L, seed = 123L) {
  if (length(n) != 1L || n < 200L) stop("`n` must be at least 200.", call. = FALSE)
  set.seed(seed)
  counts <- round(pmin(20000, pmax(900, stats::rlnorm(n, log(5000), 0.55))))
  features <- round(350 + 4000 * (1 - exp(-counts / 4800)) + stats::rnorm(n, 0, 180))
  features <- pmin(counts, pmax(300, features))
  percent_mt <- pmin(18, pmax(0.5, 0.5 + 15 * stats::rbeta(n, 2, 5)))
  doublet_score <- pmin(0.45, stats::rbeta(n, 2, 18))
  out <- data.frame(
    cell_id = sprintf("cell_%04d", seq_len(n)),
    nCount_RNA = counts,
    nFeature_RNA = features,
    percent.mt = percent_mt,
    doublet_score = doublet_score
  )
  out$nCount_RNA[1:20] <- round(seq(100, 650, length.out = 20))
  out$nFeature_RNA[1:20] <- round(seq(50, 350, length.out = 20))
  out$percent.mt[1:10] <- seq(30, 60, length.out = 10)
  out$doublet_score[1:10] <- seq(0.7, 0.98, length.out = 10)
  out$nFeature_RNA[21:35] <- round(seq(60, 380, length.out = 15))
  out$percent.mt[36:50] <- seq(25, 55, length.out = 15)
  out$doublet_score[51:65] <- seq(0.65, 0.95, length.out = 15)
  out$nCount_RNA[51:65] <- pmin(25000, round(out$nCount_RNA[51:65] * 2.2))
  out$nFeature_RNA[51:65] <- pmin(out$nCount_RNA[51:65], round(out$nFeature_RNA[51:65] * 1.5))
  out$nCount_RNA[101] <- NA_real_
  out$nFeature_RNA[102] <- NA_real_
  out$percent.mt[103] <- NA_real_
  out$doublet_score[104] <- NA_real_
  rownames(out) <- out$cell_id
  out
}
