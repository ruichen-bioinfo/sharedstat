## getting_started.R -- run sharedstat on your own replicate blocks. Base R only.
## Walk-through: inst/doc/getting_started.md. Run with
##   Rscript -e 'source(system.file("examples", "getting_started.R", package = "sharedstat"))'
## To use your own data, set `path` to your file and skip the block marked "EXAMPLE INPUT".
library(sharedstat)

## ---- 1. The input: one row per feature, one column per replicate block ---------------------------------------------
## Each column holds the SAME contrast (for example a heat-versus-baseline log2 fold change) computed separately in one
## replicate block, for the same features in the same order. Row names or a first column identify the feature.
## Any per-block fold changes will do: from DESeq2/edgeR/limma run within each block, or log2(treated/control) of
## normalised counts. Two columns are the minimum; three or more are used pairwise as a check.
path <- tempfile(fileext = ".tsv")
## ---- EXAMPLE INPUT (replace by your own file) ----
set.seed(1)
n <- 3000
truth <- rnorm(n)                                        # the response the blocks share
ex <- data.frame(feature = sprintf("g%04d", seq_len(n)),
                 block1 = truth + rnorm(n, sd = 0.6),
                 block2 = truth + rnorm(n, sd = 0.6),
                 block3 = truth + rnorm(n, sd = 0.6))
write.table(ex, path, sep = "\t", quote = FALSE, row.names = FALSE)
## ---- end of example input ----
d <- read.delim(path, check.names = FALSE)
x <- as.matrix(d[, -1]); rownames(x) <- d[[1]]
stopifnot(ncol(x) >= 2, is.numeric(x))

## ---- 2. Reliability of the reference arm ----------------------------------------------------------------------------
## h is the fraction of the across-feature variance of the block-averaged contrast that reproduces between blocks.
## It refers to the AVERAGE of the blocks supplied, so use the same blocks that form your reference arm.
rel <- if (ncol(x) == 2) reliability_from_blocks(x[, 1], x[, 2]) else reliability_from_replicates(x)
cat(sprintf("h = %.3f from %d blocks over %d features\n", rel$h, ncol(x), nrow(x)))
if (!is.null(rel$h_pairwise)) cat("pairwise h:", sprintf("%.3f", rel$h_pairwise), "\n")

## ---- 3. What a naive shared-reference gain would report at that reliability -----------------------------------------
## A true gain gamma is reported below 1 (as a loss) whenever h < 1/gamma, under independent arm errors.
gam <- c(1.05, 1.1, 1.2, 1.3, 1.5, 2)
tab <- do.call(rbind, lapply(gam, function(g) {
  b <- naive_gain_bias(h = rel$h, gamma = g)
  data.frame(true_gain = g, naive_reads = round(b$gamma_naive, 3), reported_as_loss = b$direction_wrong)
}))
print(tab, row.names = FALSE)
cat(sprintf("gains below %.3f would be reported as losses by the naive estimator\n", 1 / rel$h))

## ---- 4. Only if you have an instrument: the corrected gain -----------------------------------------------------------
## An instrument S is a per-feature quantity that shares the reference's signal but none of its measurement error,
## for example the same response measured in other samples processed separately. Supply e (response arm minus
## reference arm), R_hat (the reference arm as measured, the block average) and S, all over the same features.
## The correction is only as good as the claim that S carries none of the reference's error; see PRACTICAL_USE_GUIDE.md.
R_hat <- rowMeans(x)
resp  <- 1.1 * truth + rnorm(n, sd = 0.6 / sqrt(ncol(x)))   # EXAMPLE: a response arm with a true gain of 1.1
S     <- truth + rnorm(n, sd = 0.8)                          # EXAMPLE: an independent measurement of the same response
fit <- iv_gain(resp - R_hat, R_hat, S, n_boot = 200)
cat(sprintf("naive gain %.3f; instrumented gain %.3f [%.3f, %.3f]; cor(R_hat, S) = %.2f\n",
            fit$gamma_naive, fit$gamma, fit$ci[1], fit$ci[2], fit$cor_RS))
