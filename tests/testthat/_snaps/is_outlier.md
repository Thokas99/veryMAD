# is_outlier errors on invalid threshold

    Code
      is_outlier(1:5, threshold = -1)
    Condition
      Error in `is_outlier()`:
      ! `threshold` must be one positive number.

---

    Code
      is_outlier(1:5, threshold = c(1, 2))
    Condition
      Error in `is_outlier()`:
      ! `threshold` must be one positive number.

