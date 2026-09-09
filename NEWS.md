# sharedstat 0.1.0

First release. This is the version the accompanying paper's results were computed with; the archived
Zenodo record for this tag, not the current state of the repository, is what the paper cites.

## What it contains

Eleven exported functions covering the four constructions analysed in the paper: a statistic formed
against a shared reference arm, a residual correlated with its own component, a count of significant
features read as an ordering of effect sizes, and an interaction contrast on compositional counts.

* `rms_true()` — noise-corrected amplitude from a signal and a null contrast.
* `naive_gain_bias()` — the closed-form bias of the naive gain estimator, `gamma_hat = gamma * h`.
* `iv_gain()` — instrumental-variable gain estimator, with a first-stage strength report.
* `crossblock_compensation()`, `attenuation_factor()`, `compensation_lower_bound()` — cross-block
  estimation of compensatory divergence and the bound that makes it a lower bound rather than a
  point estimate.
* `count_vs_amplitude()` — whether a count of significant features orders effect sizes.
* `share_nonadditivity()`, `reference_class_medians()` — the compositional pair.
* `simulate_closed_sum()` — known-truth simulator in which compositionality emerges from closed-sum
  sampling rather than being added as an offset.
* `admission_checklist()` — the five-item check, scored.

## Notes on the test suite

The tests assert the values printed in the paper rather than internal invariants, so a failing test
means the package and the manuscript have diverged. Six of them pin behaviour that was **wrong** in
development and is now fixed:

* `rms_true()` dropped missing values independently from its two arguments, comparing a signal
  computed on one feature set against a noise floor computed on another. On the test's inputs it
  returned 4.30 where the paired answer is exactly 0.
* `iv_gain()` warned only when the instrument was exactly uninformative, so a weak instrument
  returned a confident wrong number (24.7 for a true gain of 1.3) in silence.
* `compensation_lower_bound()` could return an implied correlation outside [-1, 1] without comment.
* `count_vs_amplitude()` hard-coded the shared-noise-floor threshold, which is a convention rather
  than a result and is now a documented argument.
* `simulate_closed_sum()` carried a dead loop whose four columns were immediately overwritten,
  implying a more elaborate construction than the one actually used.

A test that only asserts current behaviour cannot catch a regression to a defect that has already
happened once, which is why these are written against the old behaviour.
