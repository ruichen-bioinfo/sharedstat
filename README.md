# sharedstat

Diagnose and correct statistics built from a shared component.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)

Comparative designs routinely form a statistic from a difference, a ratio or a residual of two quantities
that share a component. Four such constructions appear constantly in genomics:

| construction | example |
|---|---|
| a contrast taken against a shared reference arm, then regressed on that same reference | one allele's response measured relative to the other's |
| a residual formed by subtracting a component, then correlated with it | a *trans* term obtained as total minus *cis* |
| a count of significant features read as an ordering of effect sizes | "930 genes here, 13 there, so the effect is larger here" |
| an interaction contrast on compositional counts | a four-genotype epistasis contrast on library-normalised data |

Each has an algebraic failure mode older than genomics. This package computes the size of each one
relative to the effect being estimated, and provides the correction where a correction exists.

## The result the package rests on

Let a reference quantity *R* be measured with error, *R̂* = *R* + *u*, and let a response arm carry a
multiplicative gain on the same underlying quantity, γ*R*, measured with its own error. Form the contrast
*e* = (response) − *R̂* and regress it on *R̂*. With the reference arm's **reliability**
*h* = var(*R*)/var(*R̂*), the ordinary least-squares slope is (γ − 1)*h* − (1 − *h*), so

```
gamma_hat_naive = gamma * h
```

Three consequences, which are not the same statement three times:

1. **At γ = 1 the estimator returns *h*.** With no divergence whatsoever, any *h* < 1 produces apparent
   divergence — deterministically, in every realisation, not in a tail.
2. **The bias is O(1) in the number of features and O(1/n_rep) in replicates.** Measuring more genes makes
   the wrong value more precise. Replication helps, but only by raising *h*, and never reaches *h* = 1.
3. **The direction inverts once *h* < 1/γ.** This is what separates the composed bias from classical
   regression dilution alone, which gives γ*h* + (1 − *h*), lies between 1 and γ, and therefore cannot
   cross 1 or reverse a conclusion.

Given an instrument *S* correlated with *R* but independent of the measurement errors,
`1 + cov(e, S) / cov(R̂, S)` recovers γ free of *h*.

## Install

```r
# install.packages("remotes")
remotes::install_github("ruichen-bioinfo/sharedstat")
```

No dependencies beyond base R and `stats`. `data.table` is suggested but not required.

## Quick start

```r
library(sharedstat)

## 1. How large is the bias in MY design? Measure h from two replicate blocks and read it off.
##    h = 1 - var(null contrast) / var(signal contrast); nothing is fitted.
b1 <- rnorm(4000); b2 <- b1 + rnorm(4000, 0, 0.5)      # two blocks of the same quantity
sig <- (b1 + b2) / 2
nul <- (b2 - b1) / 2
h   <- 1 - mean(nul^2) / mean(sig^2)

naive_gain_bias(h = h, gamma = 1.1)
#> $h              reliability of the reference arm
#> $gamma_naive     what the naive estimator will return
#> $bias            gamma_naive - gamma
#> $direction_wrong TRUE if the truth and the estimate fall on opposite sides of 1

## 2. Is my noise-corrected amplitude real, or is it noise?
##    Returns exactly 0 when the signal does not exceed its own noise. Report the zero as a zero.
rms_true(sig, nul)

## 3. Correct a shared-reference contrast, and be told when the instrument is too weak to trust.
set.seed(1)
R <- rnorm(4000); Rhat <- R + rnorm(4000, 0, 0.5)
e <- 1.3 * R + rnorm(4000, 0, 0.5) - Rhat
S <- 0.8 * R + rnorm(4000, 0, 0.6)                     # shares R's signal, not its error
iv_gain(e, Rhat, S, n_boot = 500)
#> $gamma ~ 1.3 ; $gamma_naive ~ 1.3 * h ; $cor_RS ; $instrument_strength ; $se ; $ci

## 4. Does my count of significant features order effect sizes? Only if the noise floors are shared.
count_vs_amplitude(condition   = c("37_30", "37_10", "42_10", "42_30"),
                   n_sig       = c(930, 13, 9, 18),
                   amplitude   = c(0.890, 0.378, 0.165, 0.116),
                   noise_floor = c(0.7756, 0.6497, 0.6830, 0.5270),
                   anchor      = "37_30")
## attr(, "floors_shared") is FALSE here: the floors span 1.47-fold, so the counts must not be
## read as an amplitude ordering. The exaggeration column says by how much they would mislead.

## 5. Score a planned analysis against the five admission checks before running it.
admission_checklist()
```

## Reproducing the paper's numbers

The test suite **is** the reproduction check. It asserts the values printed in the paper rather than
internal invariants, so a failure means the package and the paper have diverged:

```r
# from a checkout of this repository
Rscript -e 'library(sharedstat); source("tests/test-all.R")'
```

It reproduces, among others, the share non-additivity of −0.1572, the reference-class median of −0.3355,
the count-exaggeration factors of 30.4 / 19.2 / 6.7, the 1.47-fold noise-floor spread, and `gamma * h` at
h = 0.5 / 0.7 / 0.9. Small reference tables live in `inst/extdata/`, so no external download is needed.

The analysis scripts that generated the paper's tables are **not** in this package. They are specific to
two datasets and are not reusable code; they are archived with the data deposit instead. This repository
holds only what someone could apply to their own data.

## What this package does not claim

- It ranks **no** method as more accurate than another, and runs no benchmark.
- It describes **no** standard method as biased. The bias analysed here is a property of a *design shape* —
  reusing one noisy measurement as both predictor and contrast component — not of any software.
- A low reliability does **not** imply that a published conclusion is wrong. It implies that the naive
  estimator was not identified, which is a statement about an estimator.
- The dependence itself is old. Classical regression dilution and the coupling between a difference score
  and its own baseline are long established, and the remedy for the latter — forming the contrast against
  the *mean* of the two arms rather than against one of them — predates genomics. What is offered here is
  the closed form for the two effects composed, the inversion threshold, and an implementation.

## Citation

Please cite the software and the accompanying paper. `CITATION.cff` in this repository carries machine-readable
metadata; the archived release DOI is the one to cite for a specific version.

## Licence

MIT. See [LICENSE.md](LICENSE.md).
