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
