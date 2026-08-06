qc_data <- function() {
  data.frame(
    lower = c(1, 10, 11, 12, 100),
    upper = c(1, 2, 3, 4, 20),
    both = c(1, 2, 3, 4, 100),
    row.names = paste0("sample", 1:5)
  )
}
