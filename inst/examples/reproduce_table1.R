## Rebuild the per-hybrid table of the accompanying paper from the data shipped with the package.
## Run with:  Rscript inst/examples/reproduce_table1.R
##
## This is the shortest real use of the package: measure the reliability, read off what a naive estimator
## will report, and see the ordering change. No file outside the package is needed.

library(sharedstat)

## 1. Reliability from two replicate blocks of the same contrast. Nothing is fitted.
set.seed(11)
truth <- rnorm(3000)                       # the reference signal across features
b1 <- truth + rnorm(3000, sd = 0.55)       # the contrast as measured in block 1
b2 <- truth + rnorm(3000, sd = 0.55)       # and in block 2
r <- reliability_from_blocks(b1, b2)
cat(sprintf("reliability h = %.3f  (snr %.2f, %d complete pairs)\n", r$h, r$snr, r$n_pairs))

## 2. What a naive estimator reports, for a true gain of 1.1.
nb <- naive_gain_bias(h = r$h, gamma = 1.1)
cat(sprintf("true gain 1.100 -> naive estimate %.3f ; direction wrong: %s\n",
            nb$gamma_naive, nb$direction_wrong))

## 3. The four heat conditions of the paper: do the significant-feature counts order the effect sizes?
A <- read.delim(system.file("extdata", "four_condition_amplitude.tsv", package = "sharedstat"))
cv <- count_vs_amplitude(A$condition, A$n_sig, A$rms_tru, A$rms_null, anchor = "37_30")
print(cv[, c("condition", "n_sig", "count_ratio", "amplitude_ratio", "exaggeration")], digits = 3)
cat(sprintf("noise floors span %.2f-fold; shared: %s\n",
            attr(cv, "floor_spread"), attr(cv, "floors_shared")))

## 4. The compositional pair: a zero-sum contrast can still carry an offset.
S <- read.delim(system.file("extdata", "library_share.tsv", package = "sharedstat"))
a <- S[S$cond == "37_30", ]
g <- function(k) a$share[a$genotype == k]
cat(sprintf("share non-additivity at the anchor = %.4f\n",
            share_nonadditivity(g("Wildtype"), g("HSF1.KD"), g("MSN24.KO"), g("Double.KDKO"))))
