# robust_scale handles zero-MAD rows per `zero_mad`

    Code
      robust_scale(m, zero_mad = "error")
    Condition
      Error in `robust_scale()`:
      ! 1 row(s) have MAD = 0; cannot scale.

# robust_scale errors on non-matrix input

    Code
      robust_scale(data.frame(x = 1:5))
    Condition
      Error in `robust_scale()`:
      ! `x` must be a numeric matrix.

# robust_scale errors on invalid margin

    Code
      robust_scale(m, margin = 3)
    Condition
      Error in `robust_scale()`:
      ! `margin` must be 1 (rows) or 2 (columns).

