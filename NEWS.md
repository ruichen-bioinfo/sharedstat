# sharedstat 0.1.3

Maintenance release. Exports are unchanged (thirteen functions).

## Fixed

* `reliability_from_blocks()`: the null-variance test uses the relative tolerance already used by
  `reliability_from_replicates()`. Two blocks that differ only by a constant now return `h = NA` with a warning;
  0.1.2 returned `h = 1`.
* `reliability_from_replicates()` returns `centred` and `n_dropped` and keeps its class on the `NA` path, as documented.
* `iv_gain(n_boot > 0)` restores the caller's `.Random.seed`, so the bootstrap no longer consumes or resets the
  caller's random stream.
* `admission_checklist()`: item 3 requires the measured exponent to lie in the reported range; item 4 no longer
  passes when no count is reported.
* `utils` added to Imports (`utils::combn` is used).

## Documentation and tests

* The two reliability forms are reported as "0.916 (centred) and 0.887 (uncentred)".
* README describes what the test suite asserts: values computable from the bundled tables and each estimator's
  input-domain contract.
* Tests: the pass banner is printed at the end; a comparison that could never fail is replaced by a set-equality
  check; two tautological checks removed; regression tests added for each fix above.

# sharedstat 0.1.2

## New

* `naive_gain_bias()` gains `arm_error_cov_ratio`, the two arms' measurement-error covariance divided by the
  variance of the observed reference. The default of zero is the independent-arm-error baseline, under which the
  expectation is `gamma * h` as before; a non-zero value gives `gamma * h + arm_error_cov_ratio` and moves the
  reversal threshold to `(1 - arm_error_cov_ratio) / gamma`. The function now also returns
  `gamma_naive_baseline`, `arm_error_cov_ratio` and `reversal_threshold_h`. Added because the expectation the
  package previously reported holds only under independence, and a shared library between the arms permits a
  non-zero covariance without fixing its sign. The argument is deliberately not called `kappa`, which already
  names an unrelated quantity in `simulate_closed_sum()`.
* `naive_gain_bias()` accepts `h = 0`. Both reliability functions can return zero, so refusing it meant the
  package's own diagnostic output could not be passed to its own bias function at the boundary a real dataset is
  most likely to reach. At `h = 0` the baseline expectation is zero and the general one is the covariance term.
* `reversal_threshold_h` is reported unclamped. Above 1 it means every admissible reliability lies below the
  threshold, so reversal is unavoidable at that gain; below 0 it means none does. Clamping to `[0, 1]` would make
  those two situations print the same number.
* `reliability_from_blocks(b1, b2)` returns the reliability `h` of a reference arm from two replicate
  blocks of one contrast, together with the signal-to-noise ratio and the number of complete pairs.
  Nothing is fitted. Added because the reliability previously had to be written out by hand at every call
  site, including in the README's own example. The function reproduces all 108 reliabilities reported in
  the accompanying paper to 2.9e-15. A null contrast with zero variance is
  reported as `NA` with a warning rather than as a perfectly reliable reference, because it means the same
  block was supplied twice.
* `reliability_from_replicates(x)` takes a matrix with one column per replicate and returns the same `h` for
  any number of replicates, using orthogonal Helmert null contrasts each scaled to the variance of the signal
  contrast. With two columns it equals `reliability_from_blocks()` exactly. It also returns `h_pairwise`, the
  reliability from each pair of replicates taken alone, because the spread among those values shows how much of
  a two-block reliability depends on which two blocks were used. Added for the third dataset in the
  accompanying paper, which has three biological replicates.
* `rms_true(s, n)` accepts a matrix of null contrasts as `n`, averaging the per-contrast mean squares rather
  than concatenating the values, so that each contrast contributes equally.
* `inst/examples/reproduce_table1.R`, a runnable example that measures a reliability, reads off the naive
  estimate at that reliability, and rebuilds the exaggeration factors and the share non-additivity from the
  bundled tables.

## Changed

* `reliability_from_replicates()` now returns `h` and every element of `h_pairwise` in `[0, 1]`, matching
  `reliability_from_blocks()`. A negative plug-in value means the unreproducible variation is at least as large as
  the variation across features and is reported as zero, not as a negative reliability. `snr`, `rms_signal` and
  `rms_null` are unchanged, so how far below the floor an estimate fell is still readable, and the inconsistency
  warning still reports the raw pairwise ratios before flooring. Every reliability reported in the accompanying
  paper is positive, so no published value changes: the 113 arms of the reliability table were recomputed,
  bootstrap included, and came back byte-identical.
* The help page for `admission_checklist()` no longer says which item the accompanying paper's authors
  failed. That belongs to the paper; a user reading a help page should be told what the check catches.

## Fixed

* `print.sharedstat_checklist` was defined and never registered: `NAMESPACE` carried no `S3method` line, so
  the checklist printed as a raw list.

# sharedstat 0.1.1 (unreleased)

Input-domain hardening. No estimator formula, no default and no return value changes on the inputs
0.1.0 accepted; every change below concerns inputs on which 0.1.0 returned a number that could not be
read as a result. The paper's reported values are unaffected and the tests that assert them are
unchanged.

## Inputs that used to return a number and now raise an error

* `count_vs_amplitude()` divides both ratios by the anchor's values, so an anchor with `n_sig = 0` or
  `amplitude = 0` produced `NaN` and `Inf` ratios and an exaggeration column of mixed `NA` and `0`.
  The anchor's count and amplitude must now be positive, `n_sig` and `amplitude` must be finite and
  non-negative, `noise_floor` strictly positive, and condition labels unique.
* `simulate_closed_sum()` sized each tested class as `round(frac * n_genes)`. A fraction rounding to
  zero gave `idx:(idx - 1L)`, which in R is a descending pair rather than an empty selection: the class
  disappeared from the simulated design and two unrelated genes were overwritten. Separately, a large
  `n_dominators` could overlap the tested classes, so the dominator set no longer lay outside every
  class it was supposed to dilute. Both are rejected before anything is drawn, along with negative
  fractions, `n_rep < 1` and `n_dominators` outside `[0, n_genes)`.
* `reference_class_medians()` now rejects missing class labels, which were dropped silently, and a
  reference class with no finite values, which returned a table of `NA` rather than reporting that the
  estimand does not exist.
* `attenuation_factor()` now rejects negative variances instead of returning `NaN` from `sqrt()`.
* `compensation_lower_bound()` now rejects an `r_observed` outside `[-1, 1]`.
* `crossblock_compensation()` now requires three complete cases and non-constant arms, rather than
  returning the `NA` that `cor()` gives.
* `iv_gain()` now validates `min_abs_cor` and `n_boot`, and rejects a reference or instrument with zero
  variance. `rms_true()` validates `min_n`. `naive_gain_bias()` requires scalar `h` and `gamma`; a
  vector `h` previously returned a vector from a function documented as describing one estimator.

## Changed behaviour

* `count_vs_amplitude()` compares the floor spread to `max_floor_spread` inclusively. The
  documentation defines that argument as the largest spread at which the floors may be treated as
  shared, and a strict comparison called a spread of exactly 1.1 unshared at a threshold of 1.1.
* `reference_class_medians()` uses one definition of an observation for both `n` and the median.
  Counting finite values while taking the median with `na.rm = TRUE` left infinities in the median, so
  a class of two infinities reported `n = 0` beside an infinite median.
* `iv_gain()` returns `naive_minus_iv` for the gap between the naive and instrumented estimates. The
  function never sees the true gain, so that gap is a difference between estimators and not a bias;
  the old field name invited the opposite reading. `bias` remains as a deprecated alias and will be
  removed.
* `iv_gain()` returns `se = NA` and `ci = c(NA, NA)` with a warning when every bootstrap resample was
  degenerate, instead of summarising an all-`NA` vector.
* `reference_class_medians()` accepts a named `se_of_median` and matches it to classes by name. An
  unnamed vector is still taken in the order of `sort(unique(class))`, which is silently wrong if the
  caller supplies it in their own order.
* `simulate_closed_sum()` returns `class_sizes`, so the realised design can be compared with the
  requested one, and takes `lib` as an argument.

## Tests

A second block of tests asserts behaviour outside the inputs the paper supplied.

One duplicated assertion was removed. The exaggeration factor at `37_10` was asserted twice: once from
the bundled full-precision table as 30.3, and once from hand-typed amplitudes rounded to three decimals
as 30.4 with a tolerance loose enough to hide the disagreement. The two are both arithmetically correct
for their own inputs, because rounding 0.8902870 and 0.3775597 to 0.890 and 0.378 moves the factor from
30.3386 to 30.3838. A published value computed in two places at two precisions is how a manuscript and
its software drift apart, so the reported factors now have a single source, and the sensitivity to input
precision is asserted in its place.

The test file also no longer re-sources `R/` over the installed package, which had made the suite test
the working tree when run from a checkout and the installed package when run by `R CMD check`.

# sharedstat 0.1.0

First release.

## Contents

Eleven exported functions covering four constructions in which a statistic is built from a shared
measured component: a contrast formed against a shared reference arm, a residual correlated with the
component subtracted from it, a count of significant features read as an ordering of effect size, and
an interaction contrast on compositional counts.

* `naive_gain_bias()` — closed-form bias of the naive gain estimator, `gamma_hat = gamma * h`.
* `iv_gain()` — instrumental-variable gain estimator, reporting first-stage strength and warning when
  the instrument is too weak for the estimate to be trusted.
* `rms_true()` — noise-corrected amplitude from a signal and a null contrast, returning exactly zero
  when the signal does not exceed its own noise.
* `crossblock_compensation()`, `attenuation_factor()`, `compensation_lower_bound()` — cross-block
  estimation of compensatory divergence, its attenuation, and the resulting bound.
* `count_vs_amplitude()` — whether significant-feature counts order effect sizes, with the
  shared-noise-floor threshold as a documented argument.
* `share_nonadditivity()`, `reference_class_medians()` — compositional share non-additivity and class
  medians read against a designated reference class.
* `simulate_closed_sum()` — known-truth simulator in which compositionality emerges from closed-sum
  sampling rather than being added as an offset.
* `admission_checklist()` — five-item check on a planned shared-component analysis, scored.

## Tests

The suite asserts the numeric values reported in the accompanying paper rather than internal
invariants, so a failing test means the package and the paper have diverged. It runs against the
installed package. Regression tests cover the boundary conditions of each estimator, including paired
missingness in `rms_true()`, weak instruments in `iv_gain()`, out-of-range implied correlations in
`compensation_lower_bound()`, and the sign of the compositional parameter in `simulate_closed_sum()`.
