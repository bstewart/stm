# Future document-update optimization work

The current branch removes the largest single-threaded optimizer overhead by
moving BFGS into C++, sharing objective/gradient intermediates, avoiding
unnecessary beta copies, and using the existing Cholesky factor more
efficiently. The remaining opportunities fall into three tiers.

## 1. Exact engineering improvements

These preserve the present document objective, stopping rules, Hessian,
covariance, assignments, and sufficient statistics up to floating-point
ordering.

### Fuse the complete document update

The compiled optimizer currently returns to R before posterior construction.
A single compiled routine could retain the final:

- `exp(eta)` and its normalization;
- topic-word mixture denominators;
- prior residual and prior gradient;
- document counts, total count, and square-root counts;
- allocated BFGS, Hessian, assignment, and covariance workspaces.

It would then construct the Hessian, covariance, assignments, bound, and
document contribution without another R/C++ transition or repeated setup.
Likely gain: 5-15% of document-update time.

### Move the E-step document loop into compiled code

The R loop still subsets beta columns, creates per-document objects, receives
lists from C++, and scatters assignments into global sufficient statistics.
Moving this loop into C++ would allow:

- reusable buffers sized to the largest document;
- cached word indices, counts, token totals, and square-root counts;
- direct selection of the content-specific beta block;
- direct accumulation of covariance, bound, lambda, and topic-word statistics;
- fewer temporary R objects and less garbage collection.

Likely gain: 10-25%, with the larger gains expected for many short documents.
This is the highest-priority remaining single-threaded project.

### Improve curvature and covariance kernels

At large K, exact Hessian construction and covariance inversion dominate. The
remaining exact experiments include:

- construct only the required K-1 Hessian rows and columns;
- verify that symmetric rank-k BLAS is consistently selected for `EB * EB'`;
- avoid full-matrix initialization or symmetrization where LAPACK only needs
  one triangle;
- benchmark Accelerate, OpenBLAS, and other optimized BLAS implementations;
- retain reusable K-by-K storage across documents.

Likely gain: 5-20% at K around 100, with little benefit at small K. An
optimized BLAS may provide a larger hardware-specific gain, but should not be
assumed to help small matrices.

### Exact document-level parallelism

Conditional on the global parameters, documents can be optimized
independently. A practical implementation should use:

1. dynamic scheduling based on prior iteration cost or document shape;
2. worker-local lambda, covariance, bound, and topic-word contributions;
3. sparse or vocabulary-sharded topic-word accumulators rather than a dense
   K-by-V-by-content-level copy per worker;
4. deterministic pairwise reduction when reproducibility is requested.

On the four performance cores of the benchmark machine, a reasonable target
is a 2.5-3.5x E-step wall-clock improvement. This reduces elapsed time rather
than total CPU work. Memory bandwidth and sufficient-statistic reduction are
likely to limit scaling before local inference does.

## 2. Controlled changes to local convergence

These retain the same local objective but may change the EM path slightly.
They should be optional until validated on full fits.

### EM-stage-dependent local tolerance

Early EM iterations could use a looser local tolerance or smaller iteration
cap, tightening as the global relative-bound change falls. A strict final
E-step would be required before returning the model. Safeguards should include
objective improvement and gradient checks, escalation for failed documents,
and periodic strict refreshes.

Potential gain: 10-30% of total fitting time, depending on the number of EM
iterations. Validation must compare the full bound trajectory, final bound,
lambda, Sigma, beta, prevalence effects, and fallback frequency.

### Stable-document skipping or short correction

Documents whose warm-started modes have repeatedly changed very little could
receive one or a few correction steps, with periodic full solves and automatic
fallback when the gradient or objective test fails.

This may help long fits with strong warm starts, but it requires per-document
history and is riskier than a global tolerance schedule.

## 3. Changes to the approximation

These can deliver larger gains at high K but are not drop-in engineering
optimizations. They require statistical validation and probably an explicit
user-facing option or separately named inference method.

- Active-topic inference or sparse assignments.
- Lagged, periodic, diagonal, block, or low-rank document covariance.
- Approximate Hessians or Hessian-vector products.
- Memoized or stochastic global updates.
- Amortized document initialization or inference.
- Batched GPU inference.

Exact covariance is the principal ceiling at K around 100 and above. Material
speedups beyond roughly 1.5-2x in that regime will probably require one of
these approximations or parallel execution.

## Suggested order of work

1. Fuse optimizer and posterior construction into one compiled routine.
2. Move the serial E-step loop and sufficient-statistic accumulation into C++.
3. Benchmark optimized BLAS and exact curvature-kernel variants.
4. Add deterministic document-level parallelism with bounded-memory reduction.
5. Only then test adaptive tolerance schedules.
6. Treat active topics or approximate covariance as a separate high-K track.

## Acceptance criteria

Every exact change should continue to satisfy the existing document-update
tests and add checks for:

- objective, gradient, lambda, covariance, assignments, and document bound;
- aggregate sufficient statistics and the entire EM bound trajectory;
- final beta, Sigma, eta, and prevalence quantities;
- Cholesky fallback and optimizer failure counts;
- cold and warm starts across K, unique-term count, covariance conditioning,
  and content levels;
- memory use as well as elapsed and CPU time.

For parallel work, deterministic and nondeterministic modes should be measured
separately. For controlled or approximate work, report a quality-versus-speed
frontier rather than a single runtime number.
