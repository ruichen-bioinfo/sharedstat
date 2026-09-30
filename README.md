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

Version matters for this package rather than being a formality: 0.1.0 exports eleven functions and 0.1.2 and later
export thirteen, and the two reliability functions exist only from 0.1.2. A reader who fetches the wrong version will not
find the functions that produced the published reliability numbers.

## License

MIT. See [LICENSE.md](LICENSE.md).

## A note for anyone editing this package

`NAMESPACE` is maintained **by hand**. It was not generated by roxygen2, so `roxygenise()` declines to touch it
and reports success while writing nothing: a newly documented `@export` will get a `man/*.Rd` page and still be
invisible to a user who installs the package. Add the `export()` line yourself. The test suite asserts that
every public function is in the installed namespace, so a missing line fails the tests rather than shipping.
