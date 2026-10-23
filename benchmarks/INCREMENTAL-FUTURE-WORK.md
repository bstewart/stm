# Stochastic inference: remaining work

The finite-corpus block-memoized implementation deliberately targets the
standard `LDAbeta` model with CTM or pooled prevalence. The most useful next
steps are:

1. Move block scheduling and sparse sufficient-statistic reduction into the
   native document engine. Profiling shows the exact posterior/Hessian remains
   the dominant cost, but the R document loop is still the next largest item.
2. Store memoized topic-word blocks in single precision or a native packed
   buffer. The current adaptive 256 MB conservative budget prevents runaway
   memory, but double-precision R matrices leave a roughly 2x storage gain.
3. Add a low-memory Robbins--Monro mode for streaming corpora or cases where
   even two memoized blocks exceed the cache budget. Keep the memoized method
   as the finite-corpus default.
4. Consider checkpointing the best global iterate according to the fixed fresh
   evaluation subset. Fresh subset monitoring and optional periodic
   full-corpus bounds are now implemented, but parameters are not yet rolled
   back when a later epoch deteriorates.
5. Support L1 prevalence by periodically refitting `glmnet` from the cached
   lambda matrix. This is straightforward but should be benchmarked separately
   because it lacks the compact cross-moment update used by the pooled prior.
6. Add SAGE/content inference with stratified blocks and a slow proximal or
   periodic deterministic content M-step. Rare content levels need explicit
   coverage guarantees.
7. Add synchronous native parallelism within each block. Deterministic
   reduction should be retained as a validation mode.
8. Consider approximate document covariance only as a separate opt-in tier.
   The exact `hpbcpp_potri` posterior accounts for most stochastic runtime, so
   this is the largest remaining gain but changes the variational approximation.
