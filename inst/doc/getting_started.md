# Getting started: sharedstat on your own replicate blocks

This walk-through takes a feature × replicate-block table, measures the reliability of the reference arm, and reads
off what a naive shared-reference gain would report. The runnable version is `inst/examples/getting_started.R`:

```r
source(system.file("examples", "getting_started.R", package = "sharedstat"))
```

It runs on a simulated table written to a temporary file and read back, so the code path is the one your own file
takes. To use your data, point `path` at your file and delete the block marked `EXAMPLE INPUT`.

## 1. Input

A tab-separated file with one row per feature and one column per replicate block, plus a first column naming the
feature:

| feature | block1 | block2 | block3 |
|---|---|---|---|
| g0001 | 0.41 | 0.77 | 0.18 |
| … | … | … | … |

Every block column holds the same contrast, computed within that block only: for example the log2 fold change of
heat against baseline from DESeq2, edgeR or limma run separately per block, or log2 of treated over control
normalised counts. Two blocks are the minimum. Rows with a missing value in any block are dropped.

## 2. Reliability

```r
rel <- reliability_from_replicates(x)      # x: the block columns as a numeric matrix
rel$h                                      # reliability of the block average
rel$h_pairwise                             # the same from each pair of blocks, as a consistency check
```

With exactly two blocks, `reliability_from_blocks(x[, 1], x[, 2])` returns the same quantity. `h` is the fraction of
the across-feature variance of the block-averaged contrast that reproduces between blocks. It refers to the average
of the blocks you supply, so pass the blocks that form your reference arm. Nothing is fitted.

## 3. What a naive gain would report

```r
naive_gain_bias(h = rel$h, gamma = 1.1)
```

Under independent arm errors the naive estimate of a true gain γ is γ·h, so any γ < 1/h is reported as a loss.
The script prints this for γ from 1.05 to 2. This step needs no instrument, and for most designs it is the whole
answer: it says whether the naive estimator can be read at the size of effect you are looking for.

## 4. Correction, only with an instrument

`iv_gain(e, R_hat, S)` recovers γ without h when a per-feature instrument `S` shares the reference's signal but none
of its measurement error, for example the same response measured in samples processed separately. The correction is
only as good as that claim, which the data cannot check; `PRACTICAL_USE_GUIDE.md` sets out when it is admissible.

## What the example prints

With the simulated input (three blocks, noise SD 0.6 around a unit-variance response) the script reports h ≈ 0.90,
says that true gains below about 1.11 would read as losses, and recovers the example response arm's gain of 1.1
with the instrument while the naive estimate sits near 0.99.
