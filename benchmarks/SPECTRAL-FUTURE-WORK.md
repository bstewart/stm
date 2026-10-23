# Future spectral-initialization work

The current implementation changes the dominant storage from a dense V-by-V
moment to V-by-p sketch storage plus V-by-K exact recovery kernels. It exceeds
the initial speed target on the current benchmarks while retaining a fully
exact fallback. The following work could improve robustness or scaling further.

## Automatic sketch diagnostics

The fixed 256-dimensional sketch performed well on the current cases, but
near-tied extreme points can change across corpora. A production diagnostic
could compare anchor sets from two hash seeds, measure residual margins and
the anchor Gram condition number, then increase p or fall back to exact anchor
search when instability is material. The acceptance test should use final
bounds and topic matching, not anchor identity alone.

## Fuse the projected pipeline

The R orchestration currently materializes a logical row filter and subsets the
V-by-p sketch before calling anchor search. A single native entry point could:

- stream the sketch and valid-row map;
- select anchors without copying valid rows;
- construct exact anchor rows and word-anchor kernels;
- solve recovery and reintroduce zero-moment words;
- return beta and diagnostics only.

This should save allocations and a few milliseconds on moderate vocabularies;
the relative benefit will be larger when V and p are large.

## Parallel sparse corpus passes

Sketch construction and the two exact-recovery passes are independent across
documents. Deterministic document blocks with bounded worker-local V-by-p or
V-by-K tiles could reduce wall time on large corpora. Reduction order, BLAS
threading, and memory use need explicit controls. The current single-threaded
path is already fast enough that this matters mainly for much larger D or V.

## Adaptive projection dimension

The current p is `max(256, K)`. A stability-driven schedule could start there,
increase p only when anchor margins or cross-seed agreement are poor, and use a
lower dimension for clearly separated synthetic or very large corpora. This is
preferable to a universal larger sketch because projected moment construction
is linear in p.

The width/density grid confirms the stability benefit but also shows why the
schedule cannot optimize anchor fidelity alone. At D=5,000, p=1,024 and s=.05
reproduced 97.5% of exact anchors and cut final-topic variation by 38%, yet its
median final bound was lower than the noisier default sketch. Approximation can
act as useful basin exploration. An adaptive policy should therefore retain a
small set of deliberately diverse sketches and compare downstream bounds,
rather than always escalating toward the exact anchor set.

## Candidate vocabularies with all-word recovery

Explicit `maxV` still uses the legacy frequency-truncated compatibility path.
A better large-vocabulary option would choose anchor candidates on a screened
vocabulary, then use the existing exact original-space recovery passes for all
words. Candidate scores should combine frequency, document frequency, and
projected extremeness so uncommon but informative anchors are not discarded.

## Broader statistical validation

Before a release, add corpora with long documents, substantially larger
vocabularies, weak anchor separation, content covariates, and synthetic known
topics. Report converged bound, held-out likelihood, semantic coherence,
exclusivity, permutation-matched topic divergence, anchor stability, peak
memory, and total time to convergence across sketch seeds.
