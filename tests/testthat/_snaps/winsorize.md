# winsorize errors on invalid threshold

    Code
      winsorize(1:5, threshold = -1)
    Condition
      Error in `winsorize()`:
      ! `threshold` must be one positive number.

---

    Code
      winsorize(1:5, threshold = c(1, 2))
    Condition
      Error in `winsorize()`:
      ! `threshold` must be one positive number.

