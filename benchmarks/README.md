# STM optimization benchmarks

## Document update

This benchmark isolates the logistic-normal document update on real document
shapes from `poliblog5k`, with deterministic synthetic topic-word
distributions and correlated prevalence priors. It covers cold and perturbed
warm starts at K = 20, 50, and 100. Each reported timing is the median of three
runs.

Run it from the package root with an R library containing RcppArmadillo:

```sh
Rscript benchmarks/document-update.R
```

Add `--profile` for an R-level profile, `--quick` for a smoke test, or
`--quick-k100` for a larger-K smoke test. The companion
`model-fit-validation.R` script runs a deterministic 12-iteration fit for an
end-to-end comparison between separately installed package versions.

## Benchmark environment

- Apple M3 MacBook Air, 8 cores, 16 GB RAM
- R 4.4.3, arm64
- R reference BLAS (`libRblas`)
- 5,000 documents at K = 20; 3,000 at K = 50; 1,500 at K = 100
- Median unique terms per document: 115-116

## Profile

For the K = 50 cold-start baseline, document optimization accounted for 67.1%
of sampled time and posterior/Hessian construction for 23.0%. The objective
and gradient callbacks accounted for 21.2% and 26.9% of self time,
respectively. This made removal of callback overhead and reuse of accepted
objective intermediates the main target.

## Cumulative ablations

`nocopy` removes the mutable K-by-U beta copy from the gradient and one of the
two beta copies from posterior construction. `potri` additionally constructs
the covariance from the existing Cholesky factor using LAPACK's symmetric
inverse. `fused` additionally moves BFGS orchestration into C++ and reuses the
accepted objective's exponentials, mixture denominators, prior residual, and
prior gradient.

| K | Start | Legacy seconds | `nocopy` speedup | `potri` speedup | `fused` seconds | `fused` speedup |
|---:|:------|---------------:|-----------------:|----------------:|----------------:|-----------------:|
| 20 | cold  | 1.354 | 1.017x | 1.021x | 0.666 | 2.033x |
| 20 | warm  | 0.877 | 1.007x | 1.013x | 0.466 | 1.882x |
| 50 | cold  | 1.713 | 1.000x | 1.034x | 1.087 | 1.576x |
| 50 | warm  | 1.358 | 1.062x | 1.054x | 0.860 | 1.579x |
| 100 | cold | 2.313 | 1.029x | 1.095x | 1.622 | 1.426x |
| 100 | warm | 1.977 | 1.030x | 1.126x | 1.444 | 1.369x |

The covariance optimization matters increasingly at large K. Exact Hessian
formation and covariance inversion become the limiting cost there, so the
exact implementation cannot sustain the approximately 2x K = 20 improvement
at K = 100 without changing the approximation.

Across every ablation case, the cumulative objective differed from the legacy
path by at most 9.4e-10 and the cumulative document bound by at most 3.0e-8.
There were no optimizer failures.

## End-to-end validation

A deterministic K = 20 fit on 1,000 `poliblog5k` documents for 12 EM
iterations was run three times with the untouched master package and three
times with the optimized package:

| Measure | Legacy | Optimized |
|:--------|-------:|----------:|
| Median elapsed time | 2.661 s | 1.344 s |
| Speedup | 1.000x | 1.980x |
| Final bound | -1,360,206 (approximately) | baseline - 6.87e-8 |

Maximum absolute differences at the end of the fit were 8.50e-9 for document
means, 2.48e-10 for the topic covariance, and 9.76e-8 for log topic-word
probabilities.

Potential next steps, expected gains, implementation sketches, and acceptance
criteria are recorded in [FUTURE-WORK.md](FUTURE-WORK.md).

## Spectral initialization

`spectral-initialization.R` benchmarks cumulative initialization ablations on
three deterministic `poliblog5k` cases. `spectral-bound-validation.R` compares
the converged variational bounds reached from the exact reference and final
optimized initializers.

```sh
Rscript benchmarks/spectral-initialization.R 3
Rscript benchmarks/spectral-bound-validation.R
```

The optimized initializer streams a deterministic sparse-sign sketch of the
row-normalized second-order moment. The sketch is used only to choose anchors.
For those anchors, two additional sparse corpus passes construct the exact
original-space RecoverL2 kernels, avoiding the dense V-by-V moment while still
recovering against the unprojected geometry. Recovery problems are solved in
one native batch using shared word-anchor products.

The default sketch dimension is `max(256, K)`, its independent hash seed is 9,
and the recovery tolerance is 1e-6. Models with K < 10, vocabularies no larger
than the sketch, explicit `maxV`, and `recoverEG=FALSE` use an exact
compatibility path. `control=list(SpectralRP=FALSE)` also requests the exact
dense path explicitly.

### Baseline profile

An `Rprof` run on all 5,000 documents at K=20 split 3.54 seconds into 0.98
seconds for the moment, 0.99 seconds for anchor search, and 1.58 seconds for
recovery. `expgrad` accounted for 47.9% of sampled total time (25.6% self), and
`fastAnchor` for 20.8% of total time. This ruled out a recovery-only rewrite:
meeting the 10x target required changing moment construction and anchor
geometry as well as batching recovery.

### Cumulative ablations

Each timing below is the median of three runs on the same benchmark environment
listed above. `native exact` preserves the legacy moment, anchors, and recovery
solution. `projected strict` introduces the streamed 256-dimensional anchor
search and exact original-space recovery at the legacy 1e-7 tolerance.
`projected final` changes only that recovery tolerance to 1e-6.

| Documents | V | K | Legacy seconds | Native exact speedup | Projected strict speedup | Final seconds | Final speedup | Exact-anchor overlap |
|----------:|--:|--:|---------------:|---------------------:|-------------------------:|--------------:|--------------:|---------------------:|
| 300 | 2,621 | 10 | 1.172 | 6.37x | 46.88x | 0.021 | 55.81x | 10/10 |
| 1,000 | 2,632 | 20 | 2.619 | 6.10x | 31.18x | 0.059 | 44.39x | 17/20 |
| 5,000 | 2,632 | 30 | 4.628 | 3.91x | 13.11x | 0.191 | 24.23x | 28/30 |

The native exact beta differed from the R reference by at most 4.72e-16. In
the final path, the D=5,000 phase medians were 0.052 seconds for the streamed
moment, 0.026 seconds for anchors, and 0.114 seconds for recovery. The looser
1e-6 tolerance was retained only after 1e-5 failed the downstream bound gate.

### Converged-bound validation

Both versions were run to the same 1e-5 EM convergence criterion with a cap of
150 iterations. All six fits converged. Positive deltas favor the optimized
initializer.

| Documents | V | K | Exact iterations | Optimized iterations | Exact bound | Optimized bound | Bound delta |
|----------:|--:|--:|-----------------:|---------------------:|------------:|----------------:|------------:|
| 300 | 2,621 | 10 | 69 | 81 | -395,365.8 | -395,290.5 | +75.3 |
| 1,000 | 2,632 | 20 | 64 | 49 | -1,300,857 | -1,300,619 | +237.8 |
| 5,000 | 2,632 | 30 | 24 | 22 | -6,893,038 | -6,892,647 | +391.1 |

The bound comparison uses the exact native initializer, which is numerically
equivalent to the legacy R initializer, as its reference. Thus the comparison
isolates the changed initialization geometry while both fits benefit from the
same optimized document update.

### Sketch width and density

`spectral-sketch-grid.R` varies projection width and sign density while holding
the stochastic block seed fixed. Each projected beta is recovered from exact
original-space kernels and passed as a custom initialization with the same
zero mean, diffuse covariance, and zero document modes used by the spectral
initializer. The default has p=256 and 13 signed entries per word. The grid
compares widths 256, 512, and 1,024; densities .05, .25, and fully dense at
p=256; and the exact unprojected initializer. Eight hash seeds were used per
cell. All 104 projected fits and both exact reference fits converged.

On D=1,000, K=20:

| Sketch | Nonzeros / word | Init seconds | Exact-anchor overlap | Pairwise final Hellinger | Bound SD | Median bound delta vs exact |
|:--|--:|--:|--:|--:|--:|--:|
| p256-s05 (default) | 13 | 0.061 | 0.850 | 0.317 | 520 | -310 |
| p512-s05 | 26 | 0.097 | 0.850 | 0.312 | 436 | -124 |
| p1024-s05 | 52 | 0.176 | 0.875 | 0.293 | 400 | -694 |
| p256-s25 | 64 | 0.072 | 0.838 | 0.339 | 610 | -216 |
| p512-s25 | 128 | 0.125 | 0.888 | 0.302 | 353 | -174 |
| p1024-s25 | 256 | 0.265 | 0.919 | 0.268 | 300 | -89 |
| p256-dense | 256 | 0.276 | 0.844 | 0.323 | 292 | -214 |

The exact initializer took 0.558 seconds. Width plus density gave the strongest
small-case robustness: p1024-s25 reduced pairwise final-topic Hellinger by 16%,
bound SD by 42%, and per-token bound range by 43% relative to default, while
improving the median bound by 222 points. p512-s05 was the best inexpensive
quality compromise, beating its paired default-seed fit for six of eight seeds
with a median improvement of 363 points.

On D=5,000, K=30:

| Sketch | Nonzeros / word | Init seconds | Exact-anchor overlap | Pairwise final Hellinger | Bound SD | Median bound delta vs exact |
|:--|--:|--:|--:|--:|--:|--:|
| p256-s05 (default) | 13 | 0.200 | 0.892 | 0.216 | 2,103 | +1,835 |
| p512-s05 | 26 | 0.269 | 0.929 | 0.206 | 1,819 | +1,981 |
| p512-s25 | 128 | 0.397 | 0.929 | 0.183 | 915 | +713 |
| p1024-s05 | 52 | 0.497 | 0.975 | 0.133 | 864 | +220 |
| p1024-s25 | 256 | 0.736 | 0.938 | 0.174 | 1,532 | +1,062 |
| p256-dense | 256 | 0.553 | 0.900 | 0.202 | 1,612 | +2,117 |

The exact initializer took 1.305 seconds. The wider sparse sketch was extremely
stable: p1024-s05 reduced final-topic variation by 38%, bound SD by 59%, and
bound range by 53%. It also reproduced 97.5% of exact anchors. However, its
median converged bound moved toward the exact initializer's basin and was 1,616
points below the default sketch's median. The p512-s25 compromise retained a
median 713-point improvement over exact while halving bound dispersion. The
dense p256 sketch had the highest median and best observed large-case bound,
but only modest stability gains and no consistent paired improvement over
default.

Thus the current sketch is tight enough that increasing fidelity clearly
improves reproducibility, but higher fidelity does not monotonically improve
the final optimum. Approximate anchor perturbations can be beneficial because
they explore basins that the exact initializer misses. The expanded grid did
find better best-of-eight optima than the default grid in both cases: +522
versus +303 over exact at D=1,000, and +4,393 versus +3,428 at D=5,000. This
supports mixing sketch geometries in an explicit multistart procedure. It does
not yet support replacing the fixed default universally; p512-s05 is the most
plausible low-cost candidate, while p512-s25 is the more robust compromise.

Run the full grid or just the large-case middle-width extension with:

```sh
Rscript benchmarks/spectral-sketch-grid.R 8 sketch-grid.rds
Rscript benchmarks/spectral-sketch-grid.R 8 middle.rds large-middle sketch-grid.rds
```

Remaining spectral opportunities and release-level validation ideas are in
[SPECTRAL-FUTURE-WORK.md](SPECTRAL-FUTURE-WORK.md).

## Stochastic inference

`stochastic-inference.R` compares converged batch variational EM with the new
finite-corpus stochastic generalized-EM option. `stochastic-ablation.R`
isolates block storage, block size, prevalence-update frequency, and the final
theta pass.

```sh
Rscript benchmarks/stochastic-inference.R 1
Rscript benchmarks/stochastic-ablation.R
```

The public entry point is `stm(..., inference="Incremental")`. The first
corpus pass fills a complete block cache at the initialized globals. Later
passes visit the fixed randomized blocks in shuffled order, subtract each
stale block statistic, add its exact refreshed statistic, and run the same
beta, covariance, and pooled-prevalence parameter maps as batch STM. Topic-word
block matrices store only vocabulary columns observed in that block. This is a
memoized stochastic generalized-EM method rather than the old SVI branch's
Adam update on unconstrained beta and covariance parameters.

### Adaptive defaults

- Batch size is at least `max(128, 8*K, ceiling(N/12))` and is increased when
  the conservative dense-cache upper bound would exceed 256 MB.
- The epoch cap decreases from 60 to 20 as the corpus grows; convergence needs
  at least three epochs and two consecutive relative changes below `emtol`.
- At 10,000 documents and above, a fixed 1,000-document subset receives a
  separate fresh E-step once per epoch. The documents remain in training. By
  default, both this fresh bound and the memoized cache bound must stabilize.
  `incremental.full.eval.every` optionally requests a fresh full-corpus bound,
  but is off by default because it costs a complete E-step.
- Pooled prevalence is refreshed once per epoch from cached `X'lambda`,
  `lambda'lambda`, and `X'X`; the statistic-only solver matches the original
  regression to machine precision.
- A final exact theta-only pass is on below 10,000 documents and off above that
  threshold. It does not update beta, Gamma, or Sigma. Deferred theta can be
  materialized with `updateTheta()`.
- Block membership and order use a deterministic seed derived from the model
  seed. Content covariates, SAGE beta, and L1 prevalence are rejected explicitly
  in this first implementation.
- `incremental.max.updates` can impose an absolute block-update limit. Epoch and
  update limits are cumulative on restart, and convergence histories are
  retained rather than reset.

### Large-scale batch / memoized / Robbins--Monro comparison

The large-scale harness creates a 50,000-document, 20-topic, 5,000-word
synthetic corpus with ten true anchor words per topic and about four million
tokens. All three estimators receive the same deliberately noisy custom beta,
so final-bound differences reflect optimization paths rather than different
initializers. A common fresh full-corpus E-step is run after fitting, and the
result also records Hellinger distance to the true beta, anchor recovery, and
theta correlation after optimal one-to-one topic matching.

The harness installs the current worktree and `origin/svi` into separate
temporary R libraries before running:

    STM_DEPENDENCY_LIBRARY=/path/to/R/library \
    STM_BENCHMARK_DOCUMENTS=50000 \
    STM_BENCHMARK_EPOCHS=20 \
    sh benchmarks/large-scale-compare.sh

The default is a matched 20-pass comparison. The old SVI branch would choose
only about two passes from its adaptive default at this scale, so the harness
deliberately gives Robbins--Monro substantially more opportunity to converge.
Set `STM_BENCHMARK_SEEDS="1 2 3"` for a basin-stability comparison. Generated
corpora, fitted summaries, per-process timing logs, and the combined CSV are
written outside the repository under a temporary work directory unless
`STM_LARGE_WORK_ROOT` is supplied.

One matched 20-pass run on the benchmark environment above produced:

| Method | Fit seconds | Fresh bound / token | Bound delta vs batch | Mean beta Hellinger | Mean / minimum anchor recall | Mean theta correlation |
|:--|--:|--:|--:|--:|--:|--:|
| Batch variational EM | 67.559 | -7.043109 | 0 | 0.2108 | 1.000 / 1.000 | 0.7599 |
| Memoized stochastic | 66.791 | -7.041246 | +7,452 | 0.2343 | 1.000 / 1.000 | 0.7585 |
| `origin/svi` Robbins--Monro | 960.025 | -7.098653 | -222,119 | 0.1855 | 0.995 / 0.900 | 0.7415 |

The memoized method was 1.01x faster than batch and 14.4x faster than the old
Robbins--Monro implementation. It reached the best variational bound and
retained all true anchors. Robbins--Monro estimated beta closest to the DGP by
mean Hellinger distance, but had a substantially worse fitted objective,
slightly worse theta recovery, and missed one anchor in one matched topic.
Thus this run does not show Robbins--Monro finding a better variational basin;
it does show that closeness to the generating beta and maximization of the
finite-sample variational objective need not rank solutions identically.

All three methods reached their 20-pass cap before their strict convergence
test fired, so this is a fixed-compute-path comparison rather than a
time-to-convergence result. The old branch trains on 49,000 documents after
reserving its default 1,000-document holdout, whereas the other two methods
train on all 50,000. The final bound reported for every method nevertheless
comes from the same fresh E-step over all 50,000 documents.

### Profile

On 1,000 poliblog documents at K=20, 89.9% of sampled stochastic fitting time
was in the E-step. Exact posterior/Hessian construction (`hpbcpp_potri`) was
61.8% self time, the document loop 21.2%, compiled mode optimization 4.9%, and
all stochastic global M-step work 5.8%. The final theta-only pass was 1.4% of
profiled time. This confirms that the stochastic orchestration is no longer
the bottleneck; further large gains require a faster or approximate document
covariance calculation.

### Converged-bound validation

Both methods used the same optimized spectral initialization and `emtol=1e-5`.
All six fits converged. Each reported bound comes from a separate full E-step
at the final globals, rather than the stochastic cache-monitoring bound. The
timings are single end-to-end runs on the benchmark environment above and
include initialization and the default final theta pass.

| Case | Batch seconds | Stochastic seconds | Speedup | Batch bound | Stochastic bound | Bound delta | Batch / stochastic epochs |
|:--|--:|--:|--:|--:|--:|--:|--:|
| D=1,000, K=20, CTM | 4.059 | 3.511 | 1.16x | -1,300,608 | -1,300,596 | +12.1 | 49 / 37 |
| D=1,000, K=15, prevalence | 2.634 | 2.497 | 1.05x | -1,311,293 | -1,311,238 | +55.2 | 38 / 30 |
| D=5,000, K=30, CTM | 17.238 | 13.394 | 1.29x | -6,892,589 | -6,892,276 | +312.5 | 22 / 15 |

The stochastic cache used 2.66 MB, 2.31 MB, and 7.42 MB respectively.

### Seed stability

`stochastic-seed-stability.R` separates the two potential randomness sources.
The projected spectral initializer uses a pseudorandom sparse hash, but its
default `rp.seed=9` is fixed and independent of the ordinary model seed. Thus,
under the package defaults, changing `seed` changes the stochastic block
membership and order but not the spectral sketch. Users can vary the sketch
explicitly with `control=list(rp.seed=...)`.

Run the five-seed benchmark with:

```sh
Rscript benchmarks/stochastic-seed-stability.R 5 results.rds
```

The benchmark first verified that initialized beta was bit-identical across
model seeds 1, 91, and 999 with the default sketch. It then matched topics
one-to-one by maximum Hellinger affinity before comparing beta and theta. All
25 reported comparisons converged. On the D=1,000, K=20 CTM case, the
decomposition was:

| Seeds varied | Fresh-bound range | Range per token | Mean Hellinger | Mean top-20 overlap | Mean theta correlation |
|:--|--:|--:|--:|--:|--:|
| Update/block only | 158 | 0.000821 | 0.045 | 19.43 / 20 | 0.9965 |
| Spectral sketch only | 1,575 | 0.008208 | 0.300 | 13.95 / 20 | 0.7450 |

Changing the sketch itself caused most of the final-fit variation. The
initializers from the four non-default sketch seeds were still close after
topic matching (mean Hellinger 0.044-0.084 and mean top-20 overlap 16.65-18.70),
but EM sometimes amplified these differences into distinct local optima.

Paired sketch and update seeds gave:

| Case | Fresh-bound range | Range per token | Mean Hellinger | Mean top-20 overlap | Mean theta correlation |
|:--|--:|--:|--:|--:|--:|
| D=1,000, K=20, CTM | 1,232 | 0.006424 | 0.298 | 14.16 / 20 | 0.7597 |
| D=1,000, K=15, prevalence | 405 | 0.002110 | 0.234 | 15.53 / 20 | 0.8309 |
| D=5,000, K=30, CTM | 5,031 | 0.005003 | 0.191 | 15.38 / 20 | 0.8902 |

Against the converged batch bounds above, all five paired prevalence fits were
better. For the D=1,000 CTM, the default sketch was 12 points better than batch
and the other four sketches were 45-1,220 points worse. For D=5,000, the
default was 313 points better; one alternate sketch was 3,164 points better
and three were 711-1,868 points worse. These five-seed results support keeping
the sketch seed fixed by default for reproducibility. Randomized sketch
multistarts can occasionally find a better basin, but they are not a free
quality improvement and should be an explicit, separately seeded option.

### Can an early bound select the best sketch?

`stochastic-basin-screening.R` tests whether a randomized spectral start that
will ultimately reach a better basin already has a better variational bound at
the beginning. It holds the stochastic block seed fixed and varies 12
prespecified `rp.seed` values. For every sketch it records the exact
full-corpus E-step bound at initialization, a fresh bound after one exact E/M
step, and a fresh bound after convergence. The first stochastic epoch is a
particularly clean checkpoint: it fills the complete cache before applying a
batch-equivalent M-step, and shuffled block replacement starts only in epoch 2.

All 36 final fits converged with a 75-epoch ceiling. The independently recorded
initial bounds were identical between the one-step and full runs. Results were:

| Case | Initial / final Spearman | One-step / final Spearman | Initial winner's final rank | One-step winner's final rank | Initial / one-step selection regret per token |
|:--|--:|--:|--:|--:|--:|
| D=1,000, K=20, CTM | -0.469 | 0.329 | 8 / 12 | 1 / 12 | 0.003423 / 0 |
| D=1,000, K=15, prevalence | 0.070 | 0.035 | 2 / 12 | 9 / 12 | 0.000150 / 0.002721 |
| D=5,000, K=30, CTM | 0.273 | 0.385 | 1 / 12 | 8 / 12 | 0 / 0.002947 |

None of the six rank correlations was statistically distinguishable from zero
in this small experiment. Pairwise ordering accuracy was 0.35-0.61 at
initialization and 0.52-0.64 after one step. Trajectory crossings are therefore
common: selecting the single apparent early winner found the eventual winner
in only two of the six case/checkpoint combinations.

A shortlist screen was more promising. Advancing the top three one-step starts
found the exact best final basin in two cases and missed it by only 28.8 bound
points, or 0.00015 per token, in the prevalence case. This supports a staged
multistart experiment that advances several survivors, but not choosing one
start from its initial or one-step bound. More corpora and topic counts would be
needed before making that heuristic a package default.

Run the benchmark and produce its RDS, CSV, and six-panel diagnostic plot with:

```sh
Rscript benchmarks/stochastic-basin-screening.R 12 basin-screening.rds
```

### Ablations

- Replacing the existing dense `ngroups=7` implementation on the D=1,000,
  K=20 case reduced time from 6.03 to 3.51 seconds (1.72x) and improved the
  fresh bound from -1,300,628 to -1,300,596.
- Compact block statistics are numerically exact. On a synthetic V=50,000,
  K=20 sparse block they reduced beta-statistic storage from 7.63 MB to 2.10 MB
  (3.63x), with no measurable timing change at this small block size. Poliblog's small vocabulary is
  almost saturated within a block, so its memory reduction is only 4%.
- On the D=1,000 CTM case, raising the block from 40 to the adaptive 160
  documents reduced time from 4.275 to 3.522 seconds (1.21x). Its fresh bound
  decreased by 108 points relative to the much noisier 40-document path but
  remained 12 points above converged batch.
- Refreshing pooled prevalence once per epoch rather than every block reduced
  time from 2.877 to 2.458 seconds (1.17x), with a one-point improvement in the
  fresh bound.
- The theta-only pass was 3.06x faster than a full E-step on D=1,000 and its
  document modes were identical to machine precision.

Remaining scale, streaming, content, and native-parallel opportunities are in
[INCREMENTAL-FUTURE-WORK.md](INCREMENTAL-FUTURE-WORK.md).

## Note on names (1.4.0)

The option that these notes call "stochastic" inference shipped as
`stm(..., inference = "Incremental")`, with control arguments named
`incremental.*` (the `max.updates`, `max.blocks`, `eval.seed` and
`eval.patience` controls were dropped) and `updateTheta()` in place of
`updateStochasticTheta()`. The projected spectral sketch shipped as an opt-in
(`control = list(SpectralRP = TRUE)`, default `rp.p = 1024`) rather than the
default, and the exact native path is the default initialization. The scripts in
this directory were updated to the released names but have not been re-run
since.

## Stochastic inference (1.4.0 release design)

`stm(..., inference = "Stochastic")` runs online EM (Cappe and Moulines 2009;
the natural-gradient SVI update of Hoffman et al. 2013 for these M-steps) on
shuffled minibatches with step size `(t + tau)^-kappa`, then hands the last
pass's minibatch statistics to the memoized incremental method as its cache,
so no cache-filling pass is needed. Progress is measured by paired per-block
bound gains (same block, consecutive passes), because minibatch bounds vary
far more across blocks than across passes.

Benchmark corpus: 25,000 documents made from five token-bootstrap replicates
of each poliblog5k document (V = 2,632, 5.03M tokens), K = 20, spectral
initialization, R 4.3.3 with reference BLAS in a 2-core container. "Fresh
bound" is the full-corpus bound at the final parameters, per token.

| Method | Time | Passes | Fresh bound / token | vs batch |
|:--|--:|--:|--:|--:|
| Batch EM | 100.8 s | 20 | -6.88370 | 0 |
| Incremental, plain fill | 84.5 s | 15 | -6.88343 | +0.00027 |
| Incremental, progressive fill (default) | 77.9 s | 14 | -6.88305 | +0.00065 |
| Stochastic, kappa=0.6 tau=10 | 65.8 s | 2.45 + 10 | -6.88460 | -0.00090 |
| Stochastic, kappa=0.6 tau=10, forced 12 memoized epochs | 75.7 s | 2.45 + 12 | -6.88461 | -0.00091 |
| Stochastic, batch 512 | 76.8 s | 2.45 + 12 | -6.88342 | +0.00028 |
| Stochastic, one pass | 68.7 s | 1 + 12 | -6.88399 | -0.00029 |
| Stochastic, kappa=0.7 tau=50 (default) | 66.4-66.7 s | 2.45 + 10 | -6.88244 | +0.00126 |
| Stochastic, kappa=0.7 tau=100 | 72.4 s | 2.69 + 11 | -6.88248 | +0.00122 |
| Stochastic, kappa=0.8 tau=50 | 72.8 s | 2.69 + 11 | -6.88247 | +0.00123 |
| Incremental, 100 blocks of 250 | 78.2 s | 14 | -6.88305 | +0.00065 |
| Stochastic, kappa=0.7 tau=50, no finish (approximate) | 18.7 s | 2.45 | -6.88881 | -0.00511 |
| Stochastic, kappa=0.6 tau=10, no finish (approximate) | 23.4 s | 3.7 | -6.88681 | -0.00311 |

Two conclusions drove the defaults. First, aggressive early steps (the usual
SVI choice of a small `tau`) move the model out of the basin that the spectral
initialization places it in, and the memoized refinement then converges to a
slightly worse optimum however long it runs; a gentle schedule
(`tau = 50`, `kappa = 0.7`, so the first step is 0.064) keeps the basin and
ends above both batch and incremental. Second, the time to convergence of every
method is dominated by the final epochs under the shared 1e-5 relative-change
criterion; the stochastic phase removes the early epochs but not that tail, so
the speed-up over batch is about 1.5x at convergence and about 2.2x to the
point where the cached bound first exceeds batch's converged value. Without
the refinement the stochastic solution is 4-5x faster than batch but 0.05-0.07%
worse in bound (the gentle default learns more slowly, so for the approximate
mode a larger `stochastic.max.epochs` or a smaller `stochastic.tau` helps).

On the 5,000-document poliblog corpus the stochastic phase has no advantage
over the incremental method (both about 20 s against 24.5 s for batch); it is a
large-corpus tool.
