# Interpreting MAD QC flags

Choose `lower` for unusually small values, `upper` for unusually large
values, and `both` for either tail. Count-like metrics may benefit from
explicit `log1p`; bounded rates are often most interpretable on their
raw scale.

`min_n` is a computational safeguard. A zero MAD is reported explicitly
and can be handled as `"na"`, `"zero"`, or `"error"`. These flags are
statistical heuristics, not biological diagnoses or automatic filtering
decisions.

For stratified analyses, split the metadata explicitly and call
[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md) on
each subset. veryMAD does not infer batches, conditions, clusters, or
causes.
