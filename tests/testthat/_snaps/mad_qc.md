# mad_qc errors on invalid inputs

    Code
      mad_qc(cell_metadata, metrics = c(nCount_RNA = "sideways"))
    Condition
      Error in `mad_qc()`:
      ! `metrics` values must be "lower", "upper", or "both".

---

    Code
      mad_qc(cell_metadata, metrics = c("nCount_RNA"))
    Condition
      Error in `mad_qc()`:
      ! `metrics` must be a named character vector.

---

    Code
      mad_qc(cell_metadata, metrics = c(missing_col = "lower"))
    Condition
      Error in `mad_qc()`:
      ! Column(s) "missing_col" not found in `data`.

