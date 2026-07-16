# row_mad errors on non-matrix input

    Code
      row_mad(data.frame(x = 1:5))
    Condition
      Error in `row_mad()`:
      ! `x` must be a numeric matrix.

# row_mad errors on invalid margin

    Code
      row_mad(m, margin = 3)
    Condition
      Error in `row_mad()`:
      ! `margin` must be 1 (rows) or 2 (columns).

