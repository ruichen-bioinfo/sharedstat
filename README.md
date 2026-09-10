# sharedstat

Diagnose and correct statistics built from a shared component.

[![DOI](https://zenodo.org/badge/1363200240.svg)](https://doi.org/10.5281/zenodo.22693748)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)

Version 0.1.0

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

```r
# install.packages("remotes")
remotes::install_github("ruichen-bioinfo/sharedstat")
```

Requires R >= 4.0. No dependencies beyond `stats`; `data.table` is suggested but not required.

## Usage

Measure the reliability of a reference from two replicate blocks, then read off what a naive estimator
will return:

```r
library(sharedstat)

b1  <- rnorm(4000)
b2  <- b1 + rnorm(4000, 0, 0.5)
sig <- (b1 + b2) / 2
nul <- (b2 - b1) / 2
h   <- 1 - mean(nul^2) / mean(sig^2)

naive_gain_bias(h = h, gamma = 1.1)
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

See `help(package = "sharedstat")` for full documentation.

## Tests

The test suite asserts the numeric values reported in the accompanying paper, so a failure indicates
that the package and the paper have diverged. Small reference tables are bundled in `inst/extdata/`,
so the tests need no external download.

```r
Rscript -e 'library(sharedstat); source("tests/test-all.R")'
```

## Relation to the accompanying paper

The package implements the estimators, the attenuation result and the checks described in a manuscript
on shared-component bias in comparative transcriptomics. The analysis scripts specific to that paper's
two datasets are not included here; they are archived with the data deposit. This repository holds only
what applies to other datasets.

## Citation

`CITATION.cff` carries machine-readable citation metadata. Cite the archived release for the version you
used; the DOI will be added here once the release is deposited.

## License

MIT. See [LICENSE.md](LICENSE.md).
