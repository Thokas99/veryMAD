# mad_qc errors on invalid inputs

    Code
      mad_qc(cell_metadata, metrics = c(nCount_RNA = "sideways"))
    Condition
      Error in `method(mad_qc, new_S3_class("data.frame"))`:
      ! `metrics` values must be "lower", "upper", or "both".

---

    Code
      mad_qc(cell_metadata, metrics = c("nCount_RNA"))
    Condition
      Error in `method(mad_qc, new_S3_class("data.frame"))`:
      ! `metrics` must be a named character vector.

---

    Code
      mad_qc(cell_metadata, metrics = c(missing_col = "lower"))
    Condition
      Error in `method(mad_qc, new_S3_class("data.frame"))`:
      ! Column(s) "missing_col" not found in `data`.

