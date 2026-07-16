# mad_score errors on non-numeric input

    Code
      mad_score("a")
    Condition
      Error in `mad_score()`:
      ! `x` must be a numeric vector.

# mad_score errors on invalid constant

    Code
      mad_score(1:5, constant = -1)
    Condition
      Error in `mad_score()`:
      ! `constant` must be one positive number.

---

    Code
      mad_score(1:5, constant = c(1, 2))
    Condition
      Error in `mad_score()`:
      ! `constant` must be one positive number.

