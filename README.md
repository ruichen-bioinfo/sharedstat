# sharedstat

Diagnose and correct statistics built from a shared component.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![DOI (v0.1.3)](https://img.shields.io/badge/DOI%20v0.1.3-10.5281%2Fzenodo.23070215-blue.svg)](https://doi.org/10.5281/zenodo.23070215)

<!-- The badge is labelled with the version it resolves to, so that a reader can tell which version’s DOI it is. -->

This source tree is version 0.1.3, a maintenance release of 0.1.2; `NEWS.md` lists what changed. Both export
thirteen functions, including `reliability_from_blocks()` and `reliability_from_replicates()`, which version 0.1.0
does not contain. **A reliability computed with 0.1.2 or later cannot be reproduced with v0.1.0**, which exports eleven.

The version DOI for 0.1.3 is 10.5281/zenodo.23070215. The concept DOI 10.5281/zenodo.22693747 always resolves to the
latest version. Version 0.1.2 is 10.5281/zenodo.22946335; version 0.1.0 is 10.5281/zenodo.22693748.

## What it is

An R package for a specific measurement problem in comparative designs: when a statistic is formed from
a difference, a ratio or a residual of two quantities that share a measured component, the shared
measurement error enters the statistic twice, with opposite signs. The resulting estimate is biased by a
factor that depends on how reliably the shared component was measured, and the bias does not shrink as
more features are added.

`sharedstat` measures the size of that bias in your own data, and applies a correction where one exists.

## The quantity it computes

Let a reference be measured with error, `Rhat = R + u`, and let a response arm carry a gain on the same
underlying quantity, `gamma * R`, measured with its own error. Form the contrast `e = response - Rhat`
and regress it on `Rhat`. Writing `h = var(R) / var(Rhat)` for the reliability of the reference,

```
gamma_hat_naive = gamma * h
```

Three consequences follow:

- at `gamma = 1` the estimator returns `h`, so any `h < 1` produces apparent divergence where there is
  none;
- the bias is `O(1)` in the number of features and `O(1/n_rep)` under independent replicate averaging,
  so adding features does not remove it;
- the sign of the inference inverts once `h < 1/gamma`, which classical attenuation on its own cannot do.

Given an instrument `S` correlated with `R` but independent of the measurement errors,
`1 + cov(e, S) / cov(Rhat, S)` recovers `gamma` without needing to know `h`.

## Install

From the source tarball, which works now:

```r
install.packages("sharedstat_0.1.3.tar.gz", repos = NULL, type = "source")
```

Or from GitHub, pinned to a tag so that the version you get is the version you name:

```r
install.packages("remotes")
remotes::install_github("ruichen-bioinfo/sharedstat@v0.1.3")
```

Requires R >= 4.0. No dependencies beyond `stats`; `data.table` is suggested but not required.

## Usage

To run the package on your own replicate blocks, start with `inst/doc/getting_started.md` and its runnable
script, `source(system.file("examples", "getting_started.R", package = "sharedstat"))`.

Measure the reliability of a reference from two replicate blocks, then read off what a naive estimator
will return:

```r
library(sharedstat)

truth <- rnorm(4000)
b1 <- truth + rnorm(4000, sd = 0.5)   # one contrast, measured in replicate block 1
b2 <- truth + rnorm(4000, sd = 0.5)   # and in block 2

r <- reliability_from_blocks(b1, b2)  # h, snr, complete pairs; nothing is fitted
naive_gain_bias(h = r$h, gamma = 1.1)
```

Correct a shared-reference contrast with an instrument, with a report on instrument strength:

```r
set.seed(1)
R    <- rnorm(4000)
Rhat <- R + rnorm(4000, 0, 0.5)
e    <- 1.3 * R + rnorm(4000, 0, 0.5) - Rhat
S    <- 0.8 * R + rnorm(4000, 0, 0.6)

iv_gain(e, Rhat, S, n_boot = 500)
```

Ask whether a count of significant features can be read as an ordering of effect size:

```r
count_vs_amplitude(condition   = c("37_30", "37_10", "42_10", "42_30"),
                   n_sig       = c(930, 13, 9, 18),
                   amplitude   = c(0.890, 0.378, 0.165, 0.116),
                   noise_floor = c(0.7756, 0.6497, 0.6830, 0.5270),
                   anchor      = "37_30")
```

The `floors_shared` attribute is `FALSE` when the noise floors are too dissimilar for the counts to
carry an ordering, and the `exaggeration` column gives the factor by which they would mislead.

## Where sharedstat fits in an analysis

sharedstat does not replace DESeq2, edgeR or limma, which estimate per-gene changes. It checks a statistic built
across genes in which one measured response is used twice, as the reference inside a difference and as the
regression axis, and says whether the direction of that statistic can be read.

1. **Before the experiment.** With `var_ratio` (per-observation noise variance over signal variance) and `n_rep`,
   `naive_gain_bias()` gives h = 1/(1 + var_ratio/n_rep). The reference arm clears the reversal threshold 1/γ when
   n_rep > var_ratio/(γ − 1); for a 10% gain, n_rep > 10·var_ratio.
2. **Responses.** Normalise counts and compute each gene's response in each replicate block (for example log-CPM
   after TMM, treatment minus baseline). Run per-gene differential expression as usual. If a gain is computed from
   fold changes, use unshrunken estimates: shrinkage adds a second attenuation on top of the one caused by noise.
3. **Reliability.** `reliability_from_blocks(b1, b2)` returns h. `naive_gain_bias(h = h, gamma = ...)` compares it
   with the threshold for the gain you care about; `direction_wrong = TRUE` means the naive direction cannot be read.
   If both arms are measured in the same libraries, their errors may be correlated: with the null contrasts
   n_R and n_A (block 2 minus block 1, halved) of the two arms, κ = cov(n_R, n_A)/var(R̂) can be passed as
   `arm_error_cov_ratio` (the package does not compute κ itself).
4. **Corrections, each with its assumption.** naive/h assumes independent arm errors and no error shared between
   blocks. `iv_gain(e, R, S)` needs an instrument S that shares the signal but not the reference arm's error; read
   its strength warning. With two blocks, the split-block gain (axis from block 1, difference from block 2) cannot
   cross 1 if the blocks share no error, so its side of 1 gives the direction, and 1 + (gain − 1)/h₁ recovers the
   size, with h₁ = cov(b1, b2)/var(b1) the reliability of block 1 (no package function yet).
5. **Report** h with the gain, and run `admission_checklist()` for the five checks.

### Worked example

Simulated data with a true gain of 1.10, 5,000 genes and two blocks (`worked_example.R`, sharedstat 0.1.3):

| Estimate | Value | Reading |
|---|---|---|
| naive gain | 0.958 | reads as a decrease, which is wrong |
| h from two blocks | 0.865 | below 0.909, the threshold for a 10% gain; `direction_wrong = TRUE` |
| naive/h | 1.107 | close to 1.10 |
| split-block gain | 1.073 (expected 1.077) | above 1, so the direction is right |
| split-block gain, size-corrected | 1.095 | close to 1.10 |
| instrumented gain | 1.104 | close to 1.10 |

## Using sharedstat on your own data

The true gain is not needed. h is measured from replicates alone: agreement between blocks is signal,
disagreement is noise. The gain enters only as the effect size you care about.

**At-risk range.** Under independent errors the naive gain is γh, so any true increase between 1 and 1/h can read
as a decrease. With h = 0.865, 1/h = 1.156: true increases up to about 16% are at risk. For an observed naive gain
g, g/h is the gain implied under independent errors; if g < 1 but g/h > 1, the direction cannot be read. In the
yeast hybrids of the accompanying paper the median h is 0.89, so the at-risk range is 1 to 1.12.

| You have | You can |
|---|---|
| two or more replicate blocks | the full workflow: h, the at-risk range, the corrections, the split-block direction |
| one replicate | use h values published for similar designs (an upper bound on h from counting noise is planned) |
| only published summary tables | a sensitivity analysis: does the direction flip for plausible h? |
| no data yet | plan replicates with h = 1/(1 + var_ratio/n_rep) |

**Input:** each gene's response of the reference arm in each replicate block (at least two blocks); optionally the
same for the response arm, which gives κ; optionally the effect size that matters. **Output:** h; the at-risk range
1 to 1/h; the naive gain and naive/h; the direction flag; with two blocks, the split-block gain.

**Why sequencer quality scores do not replace replicates.** Base-call quality describes read errors. The noise in
h comes mainly from differences between replicate cultures, libraries and sampling, which only replicates measure.

## What the data must look like

- The same measured response must appear twice in your statistic: once inside a difference, once as the axis
  the difference is regressed on (or compared with). If your axis comes from separate samples, you do not need
  this package.
- At least two replicate blocks: complete, independent runs of the experiment (separate cultures, libraries and
  sequencing), not technical repeats of one library. Technical repeats measure only part of the noise and make h
  look too high.
- The same genes, normalisation and response definition in every block. Compute h on exactly the genes and the
  scale you use for the gain.
- Enough genes for a stable estimate: several hundred per group; with fewer, the interval on h is wide.

## Where the estimates can still be wrong

- Shared error between blocks. h assumes the blocks do not share error. A batch effect, a shared reference sample
  or a common normalisation step that moves both blocks the same way inflates h, and naive/h then under-corrects.
- Two blocks only. h from two blocks has a wide interval. Report the interval; if it includes the reversal
  threshold for the gain you care about, treat the direction as unresolved, not as a finding.
- Correlated errors between the two arms. If both alleles are read from the same library, their errors are
  usually positively correlated (κ > 0). This moves the reversal threshold to (1 − κ)/γ, as `naive_gain_bias()`
  reports when κ is given as `arm_error_cov_ratio`. Measure κ; do not assume zero.
- Noise that differs between genes. h is one number for all genes. Low-count genes are much noisier, so h can
  differ between expression classes; when the classes differ in expression, compute it per class by passing each
  class's genes to `reliability_from_blocks()`.
- Shrunken fold changes. Shrinkage estimators (for example apeglm in DESeq2) pull fold changes toward zero and add
  their own attenuation. Use unshrunken values for the gain.
- The instrument. The instrumented gain is unbiased only if the instrument shares no error with either arm. That
  cannot be checked from the data. Use the strength warning, and report how large a correlation with the reference
  error would reverse the result (the reversing correlation).
- Which allele is the reference. The gain depends on that choice. Report both orientations when they disagree.

## What it does not do

- It does not estimate the true gain from a single noisy measurement. With only the variances and covariance of
  the two arms, the gain is not identifiable; something more is needed (replicates, an instrument, or a second,
  independent measurement).
- It does not replace differential expression. It checks one kind of statistic built on top of it.

## Planned (not in 0.1.3)

- A function for the split-block gain and its size correction.
- An upper bound on h from counting noise alone, for designs with one replicate.
- Gene-level error weights from the count model, so that noisy genes count less.
- Estimators that need neither an instrument nor extra replicates, for data where the signal is far from normal,
  and bounds on the gain under the weakest assumptions. These need new theory and their own validation.

## Exported functions

| function | purpose |
|---|---|
| `reliability_from_blocks()` | the reliability `h` of a reference arm, from two replicate blocks of one contrast |
| `reliability_from_replicates()` | the same reliability from three or more replicates, pooling every pair |
| `naive_gain_bias()` | closed-form bias of the naive gain estimator, `gamma * h` |
| `iv_gain()` | instrumental-variable gain estimator with first-stage strength report |
| `rms_true()` | noise-corrected amplitude from a signal and a null contrast |
| `crossblock_compensation()` | cross-block estimate of compensatory divergence |
| `attenuation_factor()` | attenuation of that estimate |
| `compensation_lower_bound()` | implied bound on the underlying correlation |
| `count_vs_amplitude()` | whether significant-feature counts order effect sizes |
| `share_nonadditivity()` | non-additivity of compositional shares across a contrast |
| `reference_class_medians()` | class medians read against a designated reference class |
| `simulate_closed_sum()` | known-truth simulator with compositionality from closed-sum sampling |
| `admission_checklist()` | five-item check on a planned shared-component analysis |

Thirteen exported functions. `inst/examples/reproduce_table1.R` is the shortest real use: it measures a
reliability, reads off what a naive estimator would report at that reliability, and rebuilds the
exaggeration factors and the share non-additivity from the bundled tables.

See `help(package = "sharedstat")` for full documentation.

## Tests

The test suite asserts the count, amplitude, exaggeration, floor-spread, compensation and share values that the
package can compute from its bundled tables, alongside each estimator's input-domain contract. The reliability, gain
and mixture-test values of the accompanying paper are re-derived by lock scripts in its analysis deposit, not here.

```r
Rscript -e 'library(sharedstat); source("tests/test-all.R")'
```

## Relation to the accompanying paper

The package implements the estimators, the attenuation result and the checks described in a manuscript
on shared-component bias in comparative transcriptomics. The analysis scripts specific to that paper's
two datasets are not included here; they are archived with the data deposit. This repository holds only
what applies to other datasets.

## Citation

Cite the version you used, by its own version DOI, which is on that version’s Zenodo record:

> Chen R. sharedstat: diagnose and correct statistics built from a shared component. Version 0.1.3.
> Zenodo. doi:10.5281/zenodo.23070215.

The accompanying paper’s results were computed with version 0.1.2, doi:10.5281/zenodo.22946335. `CITATION.cff`
carries the machine-readable entry for this tree: version 0.1.3, with its version DOI and the concept DOI. The DOI 10.5281/zenodo.22693748 belongs to v0.1.0.

The two reliability functions exist only from 0.1.2 (see the note at the top), so a reader who fetches v0.1.0 will not
find the functions that produced the published reliability numbers.

## License

MIT. See [LICENSE.md](LICENSE.md).

## A note for anyone editing this package

`NAMESPACE` is maintained **by hand**. It was not generated by roxygen2, so `roxygenise()` declines to touch it
and reports success while writing nothing: a newly documented `@export` will get a `man/*.Rd` page and still be
invisible to a user who installs the package. Add the `export()` line yourself. The test suite asserts that
every public function is in the installed namespace, so a missing line fails the tests rather than shipping.
