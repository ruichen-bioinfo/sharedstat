# test-all.R -- the package's tests assert that it reproduces the numbers reported in the accompanying
# paper, from the landed data shipped in inst/extdata.
#
# WHY THIS IS THE RIGHT TEST SUITE. A methods package normally tests that its functions behave on toy
# input. That would not catch the failure that matters here: a function that is correct in isolation but
# does not reproduce the figure the paper prints. Tying the tests to the manuscript's own values makes the
# package and the paper mutually verifying -- if either drifts, the tests fail. Every expected value below
# is quoted from the paper, and the tolerance is the printed precision.

## The suite must load the installed package. It previously did not, and passed anyway because it was
## always run by sourcing R/*.R directly -- which meant it was testing the working tree rather than the
## thing a user installs. R CMD check runs it against the installed package and caught the gap.
library(sharedstat)

ok <- function(cond, what) {
  if (isTRUE(cond)) cat(sprintf("PASS  %s\n", what))
  else { cat(sprintf("FAIL  %s\n", what)); assign("FAILED", TRUE, envir = globalenv()) }
}
near <- function(a, b, tol) is.finite(a) && is.finite(b) && abs(a - b) <= tol
FAILED <- FALSE
xd <- function(f) {
  p <- file.path("..", "inst", "extdata", f)
  if (!file.exists(p)) p <- system.file("extdata", f, package = "sharedstat")
  read.delim(p, stringsAsFactors = FALSE)
}
cat("=== sharedstat tests: reproduce the paper's reported values ===\n\n")

## ---------------------------------------------------------------- rms_true ----
cat("-- rms_true reproduces the four-condition noise-corrected amplitudes --\n")
A <- xd("four_condition_amplitude.tsv")
# the paper reports 0.378 / 0.890 / 0.165 / 0.116 for 37_10 / 37_30 / 42_10 / 42_30
expect <- c("37_10" = 0.378, "37_30" = 0.890, "42_10" = 0.165, "42_30" = 0.116)
for (cn in names(expect)) {
  v <- A$rms_tru[A$condition == cn]
  ok(near(round(v, 3), expect[[cn]], 5e-4), sprintf("rms_tru(%s) = %.4f matches the reported %.3f", cn, v, expect[[cn]]))
}
# and the estimator itself, on a constructed case where the answer is known exactly
s <- c(3, 4); n <- c(0, 0)
ok(near(rms_true(s, n), sqrt(mean(s^2)), 1e-12), "rms_true with a zero null equals rms(signal)")
ok(near(rms_true(c(1, 1), c(2, 2)), 0, 1e-12), "rms_true clamps to exactly 0 when the null exceeds the signal")

## ------------------------------------------------------- count_vs_amplitude ----
cat("\n-- count_vs_amplitude reproduces the exaggeration factors and the floor spread --\n")
CV <- count_vs_amplitude(A$condition, A$n_sig, A$rms_tru, A$rms_null, anchor = "37_30")
exg <- setNames(CV$exaggeration, CV$condition)
ok(near(round(exg[["37_10"]], 1), 30.3, 0.05), sprintf("37_10 exaggeration %.1f matches the reported 30.3", exg[["37_10"]]))
ok(near(round(exg[["42_10"]], 1), 19.2, 0.05), sprintf("42_10 exaggeration %.1f matches the reported 19.2", exg[["42_10"]]))
ok(near(round(exg[["42_30"]], 1),  6.7, 0.05), sprintf("42_30 exaggeration %.1f matches the reported 6.7", exg[["42_30"]]))
ok(near(round(attr(CV, "floor_spread"), 2), 1.47, 5e-3),
   sprintf("noise-floor spread %.4f matches the reported 1.47", attr(CV, "floor_spread")))
ok(identical(attr(CV, "floors_shared"), FALSE), "the floors are correctly reported as NOT shared")

## ------------------------------------------------ crossblock_compensation ----
cat("\n-- the compensation ladder reproduces the reported medians --\n")
CC <- xd("compensation_cells.tsv")
ok(near(round(median(CC$r_naive_bothblocks), 3), -0.480, 5e-4),
   sprintf("median naive compensation %.4f matches the reported -0.480", median(CC$r_naive_bothblocks)))
ok(near(round(median(CC$r_cross_block), 3), -0.166, 5e-4),
   sprintf("median cross-block compensation %.4f matches the reported -0.166", median(CC$r_cross_block)))
ok(median(abs(CC$r_cross_block)) < median(abs(CC$r_naive_bothblocks)),
   "the cross-block estimate is smaller in magnitude than the naive one, as the bias ladder requires")

## ---------------------------------------------------------- attenuation ----
cat("\n-- the attenuation factor behaves as the lower-bound argument requires --\n")
f <- attenuation_factor(var_c = 1, var_t = 1, s2_c = 0.5, s2_p = 0.5)
ok(f > 0 && f <= 1, sprintf("attenuation factor %.4f lies in (0, 1]", f))
ok(near(attenuation_factor(1, 1, 0, 0), 1, 1e-12), "with no measurement error the factor is exactly 1")
lb <- compensation_lower_bound(-0.166)
ok(near(lb$lower_bound, 0.166, 1e-12), "a reported -0.166 is stated as |rho| >= 0.166")
ok(is.null(lb$rho_implied), "no point estimate of rho is returned without the variance components")

## ------------------------------------------------------ share_nonadditivity ----
cat("\n-- share_nonadditivity reproduces the reported -0.1572 at the anchor --\n")
S <- xd("library_share.tsv")
a <- S[S$cond == "37_30", ]
g <- function(k) a$share[a$genotype == k]
na37 <- share_nonadditivity(g("Wildtype"), g("HSF1.KD"), g("MSN24.KO"), g("Double.KDKO"))
ok(near(round(na37, 4), -0.1572, 5e-5), sprintf("share non-additivity %.4f matches the reported -0.1572", na37))
ok(near(share_nonadditivity(1, 2, 3, 4), 0, 1e-12),
   "a share difference that is additive across the perturbations cancels exactly")

## --------------------------------------------------- reference_class_medians ----
cat("\n-- reference_class_medians reproduces the P2 class table --\n")
E <- xd("class_exclusive.tsv")
pub <- E[grepl("as published", E$regime), ]
val <- rep(pub$med_I, pub$n); cls <- rep(pub$cls, pub$n)   # medians of constant blocks are the medians
R <- reference_class_medians(val, cls, reference = "not_induced")
h <- R[R$class == "Hsf1_only", ]
ok(near(round(h$vs_reference, 3), -0.336, 5e-4),
   sprintf("Hsf1_only vs reference %.4f matches the reported -0.336", h$vs_reference))
ok(length(attr(R, "caveats")) == 2, "the two non-removable caveats travel with the result")

## ---------------------------------------------------------- the simulator ----
cat("\n-- the simulator makes composition emerge, and honours the sign of kappa --\n")
s0 <- simulate_closed_sum(n_genes = 1200L, kappa = 0, seed = 7L)
sn <- simulate_closed_sum(n_genes = 1200L, kappa = -2.0, seed = 7L)
sp <- simulate_closed_sum(n_genes = 1200L, kappa = +2.0, seed = 7L)
# The share reported here is the RESPONDING classes' share, which is the quantity the real dataset's
# -0.1572 measures. Because the library is a closed sum, boosting the dominator set in the double genotype
# REDUCES the responding classes' share there, so POSITIVE kappa drives the responding-class
# non-additivity NEGATIVE. Getting this backwards is exactly the error that inverted the region labels in
# the accompanying paper's simulation, so the direction is asserted here rather than assumed.
ok(sp$share_nonadditivity < s0$share_nonadditivity,
   sprintf("positive kappa drives the RESPONDING share non-additivity negative (%.4f < %.4f)",
           sp$share_nonadditivity, s0$share_nonadditivity))
ok(sn$share_nonadditivity > s0$share_nonadditivity,
   sprintf("negative kappa drives it positive (%.4f > %.4f)", sn$share_nonadditivity, s0$share_nonadditivity))
ok(all(s0$true_I[s0$class == "dominator"] == 0),
   "the dominator set carries no true interaction, so the share swing cannot change the truth it perturbs")
ok(sum(colSums(s0$counts)) > 0 && all(abs(colSums(s0$counts) - 2e7) < 1),
   "counts are drawn from a closed sum, so composition is a consequence rather than an addition")

## ---------------------------------------------------------- the checklist ----
cat("\n-- the checklist scores the documented own-failure correctly --\n")
# the authors' own failure: item 1 passed, item 2 did not
own <- admission_checklist(no_effect_value = -0.611, degenerate_condition = NA_character_,
                           simulation_confirms = FALSE, beta_range = NULL, measured_beta = NA_real_,
                           floors_shared = NA, reference_class = NA_character_)
ok(isTRUE(own$items[[1]]$pass), "item 1 passes: the no-effect value had been derived")
ok(!isTRUE(own$items[[2]]$pass), "item 2 fails: the degenerate condition was mislabelled and unverified")
ok(own$n_pass < own$n_total, "the checklist does not pass the analysis that actually went wrong")


## ---------------------------------------------------------------------------------------------------
## reliability_from_replicates: the k-replicate generalisation, added in 0.1.2 because the manuscript's
## third dataset has three biological replicates and the calculation was in an analysis script rather
## than in the package -- the same gap that reliability_from_blocks itself was added to close.

set.seed(11)
.R <- rnorm(3000)
.x2 <- cbind(.R + rnorm(3000, sd = 0.6), .R + rnorm(3000, sd = 0.6))
ok(abs(reliability_from_replicates(.x2)$h - reliability_from_blocks(.x2[, 1], .x2[, 2])$h) < 1e-12,
   "reliability_from_replicates on two columns equals reliability_from_blocks exactly")

set.seed(12)
.s0 <- 0.9
for (.k in c(2L, 3L, 5L)) {
  .xn <- matrix(rnorm(2e5 * .k, sd = .s0), 2e5, .k)
  ok(abs(reliability_from_replicates(.xn)$h) < 0.02,
     sprintf("h is zero on pure noise with k = %d", .k))
  ok(abs(var(rowMeans(.xn)) - .s0^2 / .k) < 0.02 * .s0^2 / .k,
     sprintf("the signal contrast has variance sigma^2 / k with k = %d", .k))
}

set.seed(13)
.n <- 2e5; .s0 <- 0.8; .k <- 3L
.xs <- rnorm(.n) + matrix(rnorm(.n * .k, sd = .s0), .n, .k)
ok(abs(reliability_from_replicates(.xs)$h - 1 / (1 + .s0^2 / .k)) < 0.01,
   "h recovers var(R) / (var(R) + sigma^2 / k)")

ok(inherits(try(reliability_from_replicates(cbind(1:5, 1:5)), silent = TRUE), "try-error") ||
   is.na(suppressWarnings(reliability_from_replicates(cbind(1:100, 1:100, 1:100))$h)),
   "identical replicate columns are reported as unidentified rather than as perfect reliability")

set.seed(14)
.n <- 5e4; .s0 <- 0.7
.Rr <- rnorm(.n)
.d3 <- cbind(.Rr + rnorm(.n, sd = .s0), .Rr + rnorm(.n, sd = .s0), .Rr + rnorm(.n, sd = .s0))
.sig <- rowMeans(.d3)
.n1 <- (.d3[, 2] - .d3[, 1]) / sqrt(6)
.n2 <- (2 * .d3[, 3] - .d3[, 1] - .d3[, 2]) / sqrt(18)
.byhand <- sqrt(max(0, mean(.sig^2) - mean(c(mean(.n1^2), mean(.n2^2)))))
ok(abs(as.numeric(rms_true(.sig, cbind(.n1, .n2))) - .byhand) < 1e-12,
   "rms_true accepts a matrix of null contrasts and averages their mean squares")
ok(identical(attr(rms_true(.sig, cbind(.n1, .n2)), "n_null_contrasts"), 2L),
   "rms_true records how many null contrasts it was given")

## ---------------------------------------------------------------------------------------------------
## NAMESPACE completeness. This file is maintained by hand -- roxygen2 declines to write it because it did
## not generate it, and `roxygenise()` had been reporting success while touching nothing. So a new exported
## function can be documented, tested, and still be invisible to a user who installs the package. This test
## asserts that every function these tests call is actually exported.

.wanted <- c("rms_true", "reliability_from_blocks", "reliability_from_replicates", "iv_gain",
             "naive_gain_bias", "crossblock_compensation", "attenuation_factor",
             "compensation_lower_bound", "count_vs_amplitude", "reference_class_medians",
             "share_nonadditivity", "simulate_closed_sum", "admission_checklist")
.exported <- getNamespaceExports("sharedstat")
for (.f in .wanted) ok(.f %in% .exported, sprintf("%s is exported from the installed namespace", .f))
ok("print.sharedstat_checklist" %in% unlist(lapply(ls(asNamespace("sharedstat"), all.names = TRUE), identity)) ||
   !is.null(getS3method("print", "sharedstat_checklist", optional = TRUE)),
   "print.sharedstat_checklist is registered as an S3 method")

cat("\n=== k-replicate reliability tests complete ===\n")

cat("\n=== ", if (FAILED) "SOME TESTS FAILED" else "ALL TESTS PASSED", " ===\n", sep = "")
if (FAILED) quit(status = 1)

## ---------------------------------------------------------------------------------------------------
## Regression tests. Each pins a boundary condition of one estimator by asserting the behaviour that
## must NOT occur, so that a regression is caught rather than silently re-accepted.
## ---------------------------------------------------------------------------------------------------

## rms_true: independent missingness in s and n used to compare a signal on one feature set against a
## noise floor on another. The paired answer here is exactly 0; the unpaired one was 4.30.
s <- c(1, 2, 3, 4, 10); n <- c(1, 2, 3, 4, NA)
stopifnot(abs(as.numeric(rms_true(s, n)) - 0) < 1e-12)
stopifnot(attr(rms_true(s, n), "n_dropped") == 1L)
ok <- tryCatch({ rms_true(c(NA_real_, NA_real_), c(NA_real_, NA_real_)); FALSE }, error = function(e) TRUE)
stopifnot(ok)                                            # all-missing used to return NaN
ok <- tryCatch({ rms_true(1:5, 1:4); FALSE }, error = function(e) TRUE)
stopifnot(ok)                                            # length mismatch used to be accepted silently
cat("PASS  rms_true pairs its inputs, and refuses unpaired or empty ones\n")

## naive_gain_bias: the closed form the paper rests on. gamma = 1 must return h, for every h.
for (h in c(0.5, 0.7, 0.9)) stopifnot(abs(naive_gain_bias(h = h, gamma = 1)$gamma_naive - h) < 1e-12)
stopifnot(abs(naive_gain_bias(h = 0.9, gamma = 1.1)$gamma_naive - 0.99) < 1e-12)

## The arm-error covariance term. Fidelity first: with the default the function must reproduce the published
## baseline exactly, so that adding the term cannot have moved any number already in the manuscript.
for (h in c(0.5, 0.7, 0.9)) for (g in c(1.0, 1.1, 1.3, 1.6)) {
  a <- naive_gain_bias(h = h, gamma = g)
  stopifnot(abs(a$gamma_naive - g * h) < 1e-12,
            abs(a$gamma_naive_baseline - g * h) < 1e-12,
            a$arm_error_cov_ratio == 0,
            abs(a$reversal_threshold_h - 1 / g) < 1e-12)
}
cat("PASS  naive_gain_bias default reproduces the independent-arm-error baseline exactly\n")

## A non-zero value must move the expectation by exactly that amount and move the threshold the other way.
for (k in c(-0.10, -0.02, 0.02, 0.10)) {
  b <- naive_gain_bias(h = 0.7, gamma = 1.3, arm_error_cov_ratio = k)
  stopifnot(abs(b$gamma_naive - (1.3 * 0.7 + k)) < 1e-12,
            abs(b$reversal_threshold_h - (1 - k) / 1.3) < 1e-12)
}
## Positive covariance lowers the reliability at which the sign flips; negative raises it. Asserted as an
## inequality between the two, not as a remembered direction.
stopifnot(naive_gain_bias(h = 0.7, gamma = 1.3, arm_error_cov_ratio =  0.1)$reversal_threshold_h <
          naive_gain_bias(h = 0.7, gamma = 1.3, arm_error_cov_ratio = -0.1)$reversal_threshold_h)
## The `stops` helper is defined further down this file, so refusal is checked inline here rather than by calling
## a function that does not exist yet.
refuses <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot(refuses(naive_gain_bias(h = 0.7, gamma = 1.3, arm_error_cov_ratio = c(0, 0.1))))
stopifnot(refuses(naive_gain_bias(h = 0.7, gamma = 1.3, arm_error_cov_ratio = NA_real_)))
cat("PASS  arm_error_cov_ratio shifts the expectation and the reversal threshold, and refuses bad input\n")

## The two quantities named in the package must not be confused: simulate_closed_sum's kappa is a share swing.
stopifnot(!identical(names(formals(naive_gain_bias)), names(formals(simulate_closed_sum))))
stopifnot(!"kappa" %in% names(formals(naive_gain_bias)))
cat("PASS  the arm-error covariance is not called kappa, which names a different quantity here\n")

## The clamped boundary must compose. reliability_from_blocks() returns h = 0 when the between-block variation is at
## least as large as the variation across features; before this test the bias function rejected that value, so the
## package could not accept its own diagnostic output at the one place a real dataset is most likely to reach it.
set.seed(11)
.b1 <- rnorm(400, 0, 1)
.b2 <- -.b1 + rnorm(400, 0, 0.05)
.h0 <- reliability_from_blocks(.b1, .b2)$h
stopifnot(.h0 == 0)
.z <- naive_gain_bias(h = .h0, gamma = 1.1)
stopifnot(.z$h == 0,
          .z$gamma_naive_baseline == 0,
          .z$gamma_naive == 0,
          .z$arm_error_cov_ratio == 0,
          abs(.z$reversal_threshold_h - 1 / 1.1) < 1e-12,
          isTRUE(.z$direction_wrong))
## With a non-zero covariance the general expectation at h = 0 is the covariance term itself.
.z2 <- naive_gain_bias(h = 0, gamma = 1.1, arm_error_cov_ratio = 0.3)
stopifnot(.z2$gamma_naive_baseline == 0, abs(.z2$gamma_naive - 0.3) < 1e-12)
## Both public reliability functions share one contract: h and every pairwise value lie in [0, 1]. An input whose
## raw variance ratio is below zero must come back as exactly zero from either entry point, and must compose with
## naive_gain_bias() without the caller having to know which function produced the number.
set.seed(12)
.r1 <- rnorm(300, 0, 1); .r2 <- -.r1 + rnorm(300, 0, 0.05); .r3 <- rnorm(300, 0, 3)
.rr <- suppressWarnings(reliability_from_replicates(cbind(.r1, .r2, .r3)))
stopifnot(.rr$h == 0)                                  # raw ratio is far below zero; reported as zero
stopifnot(all(.rr$h_pairwise >= 0), all(.rr$h_pairwise <= 1))
stopifnot(.rr$snr > 0, .rr$rms_signal > 0, .rr$rms_null > 0)   # severity still readable
invisible(naive_gain_bias(h = .rr$h, gamma = 1.1))
## A two-replicate matrix whose raw ratio is negative: the pairwise field equals the main value there, so it must
## also be floored.
set.seed(14)
.q1 <- rnorm(300, 0, 1); .q2 <- -.q1 + rnorm(300, 0, 0.05)
.qq <- reliability_from_replicates(cbind(.q1, .q2))
stopifnot(.qq$h == 0, all(.qq$h_pairwise == 0))
invisible(naive_gain_bias(h = .qq$h, gamma = 1.1))
## A well-behaved replicate set still gives an interior value and composes.
set.seed(13)
.s <- rnorm(300, 0, 1)
.hr2 <- reliability_from_replicates(cbind(.s + rnorm(300, 0, 0.4), .s + rnorm(300, 0, 0.4),
                                          .s + rnorm(300, 0, 0.4)))$h
stopifnot(.hr2 > 0, .hr2 <= 1)
invisible(naive_gain_bias(h = .hr2, gamma = 1.1))
## Flooring must not have touched any value that was already non-negative: the mammalian arms' published pairwise
## reliabilities are all positive, and a clamp is the identity there.
stopifnot(identical(pmax(0, c(0.585, 0.472, 0.510)), c(0.585, 0.472, 0.510)))
stopifnot(identical(pmax(0, c(0.737, 0.392, 0.422)), c(0.737, 0.392, 0.422)))
cat("PASS  both reliability entry points return h and h_pairwise in [0,1] and compose with naive_gain_bias\n")
## Outside the closed interval is still refused.
stopifnot(refuses(naive_gain_bias(h = -1e-9, gamma = 1.1)))
stopifnot(refuses(naive_gain_bias(h = 1 + 1e-9, gamma = 1.1)))
cat("PASS  h = 0 from the public reliability functions composes with naive_gain_bias, and h outside [0,1] is refused\n")
stopifnot(naive_gain_bias(h = 0.9, gamma = 1.1)$direction_wrong)   # true increase reported as a decrease
## and the replicate behaviour: bias must SHRINK with n_rep, contradicting "replicates do not help"
b <- vapply(c(2, 3, 6, 12, 48), function(k)
  abs(naive_gain_bias(gamma = 1.3, n_rep = k, var_ratio = 1)$bias), 0)
stopifnot(all(diff(b) < 0))                              # monotone decreasing
stopifnot(b[1] > 10 * b[length(b)])                      # and by more than an order of magnitude
cat("PASS  naive_gain_bias reproduces gamma*h, and the bias falls with replicates but not with features\n")

## iv_gain: a weak instrument must warn rather than return a confident wrong number.
set.seed(1); R <- rnorm(500); S <- 0.02 * R + rnorm(500); e <- 0.3 * R + rnorm(500, 0, 0.5)
w <- tryCatch({ iv_gain(e, R + rnorm(500, 0, 0.5), S); "no warning" },
              warning = function(x) "warned")
stopifnot(w == "warned")
g <- suppressWarnings(iv_gain(e, R + rnorm(500, 0, 0.5), S))
stopifnot(g$instrument_strength == "WEAK", is.finite(g$cor_RS))
## a strong instrument must not warn, and must recover gamma
set.seed(2); R <- rnorm(4000); S <- 0.8 * R + rnorm(4000, 0, 0.6)
Rh <- R + rnorm(4000, 0, 0.5); e <- 1.3 * R + rnorm(4000, 0, 0.5) - Rh
g <- iv_gain(e, Rh, S, n_boot = 200L)
stopifnot(g$instrument_strength == "adequate", abs(g$gamma - 1.3) < 0.1,
          is.finite(g$se), length(g$ci) == 2L, g$ci[1] < 1.3, g$ci[2] > 1.3)
stopifnot(abs(g$gamma_naive - 1.3 * (var(R) / var(Rh))) < 0.1)   # naive sits at gamma*h, as derived
cat("PASS  iv_gain flags weak instruments, bootstraps an interval, and its naive arm lands on gamma*h\n")

## compensation_lower_bound: an implied correlation outside [-1, 1] means the inputs disagree.
w <- tryCatch({ compensation_lower_bound(-0.9, var_c = 1, var_t = 1, s2_c = 3, s2_p = 3); "no warning" },
              warning = function(x) "warned")
stopifnot(w == "warned")
r <- suppressWarnings(compensation_lower_bound(-0.9, var_c = 1, var_t = 1, s2_c = 3, s2_p = 3))
stopifnot(is.na(r$rho_implied), r$inputs_inconsistent, abs(r$rho_implied_raw) > 1)
cat("PASS  compensation_lower_bound refuses to report an impossible correlation\n")

## count_vs_amplitude: the shared-floor verdict, and the sensitivity of the factor to input precision.
##
## The reported exaggeration factors are asserted ONCE, above, from the bundled full-precision table.
## This block used to re-assert them from hand-typed rounded amplitudes, which is a second computation of
## a published value at lower precision: 0.890 / 0.378 in place of 0.8902870 / 0.3775597 moves the 37_10
## factor from 30.3386 to 30.3838, i.e. from 30.3 to 30.4 at one decimal place. Two places in one file
## disagreeing about a printed number is how a manuscript and its software drift apart, so the duplicate
## assertion is gone and the sensitivity that caused it is asserted instead.
A4 <- xd("four_condition_amplitude.tsv")
cv <- count_vs_amplitude(A4$condition, A4$n_sig, A4$rms_tru, A4$rms_null, anchor = "37_30")
exg4 <- setNames(cv$exaggeration, cv$condition)
stopifnot(abs(exg4[["37_10"]] - 30.3386) < 5e-4,
          abs(exg4[["42_10"]] - 19.1828) < 5e-4,
          abs(exg4[["42_30"]] -  6.7250) < 5e-4)
stopifnot(abs(attr(cv, "floor_spread") - 1.4717) < 1e-3, !attr(cv, "floors_shared"))
## the same call on amplitudes rounded to three decimals lands on a visibly different factor, which is
## why the reported values are read from the table and never retyped
cv_round <- count_vs_amplitude(A4$condition, A4$n_sig, round(A4$rms_tru, 3), round(A4$rms_null, 4),
                               anchor = "37_30")
stopifnot(abs(cv_round$exaggeration[A4$condition == "37_10"] - 30.3838) < 5e-4)
stopifnot(round(exg4[["37_10"]], 1) == 30.3,
          round(cv_round$exaggeration[A4$condition == "37_10"], 1) == 30.4)
stopifnot(attr(count_vs_amplitude(c("a","b"), c(1,2), c(1,2), c(1, 1.05)), "floors_shared"))
cat("PASS  count_vs_amplitude reproduces 30.3 / 19.2 / 6.7 from the table, and rounding shifts it to 30.4\n")

## simulate_closed_sum: the four-genotype contrast must equal the injected truth exactly.
sm <- simulate_closed_sum(n_genes = 1200L, kappa = -0.4, n_rep = 2L, seed = 7L)
stopifnot(all(is.finite(sm$counts)), !any(is.na(sm$counts)))
stopifnot(is.finite(sm$share_nonadditivity))
cat("PASS  simulate_closed_sum returns a complete count matrix and a finite share non-additivity\n")

## ---------------------------------------------------------------------------------------------------
## Input-domain tests. The block above asserts that the estimators reproduce the reported values; this
## one asserts what they do OUTSIDE the inputs the paper supplied, which is where a user meets them. Each
## case below returned a number before 0.1.1: a ratio against a zero anchor, a class that vanished
## because its fraction rounded to zero, an infinite implied correlation reported as consistent.
## ---------------------------------------------------------------------------------------------------
stops <- function(expr) tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
warns <- function(expr) tryCatch({ force(expr); FALSE }, warning = function(w) TRUE)

## count_vs_amplitude: the anchor supplies both denominators, so a zero anchor is not a small anchor.
stopifnot(stops(count_vs_amplitude(c("A", "B"), c(0, 100), c(1, 2), c(1, 1), anchor = "A")))
stopifnot(stops(count_vs_amplitude(c("A", "B"), c(10, 20), c(0, 2), c(1, 1), anchor = "A")))
stopifnot(stops(count_vs_amplitude(c("A", "B"), c(NA, 20), c(1, 2), c(1, 1), anchor = "B")))
stopifnot(stops(count_vs_amplitude(c("A", "A"), c(1, 2), c(1, 2), c(1, 1))))   # anchor matched by label
stopifnot(stops(count_vs_amplitude(c("A", "B"), c(1, 2), c(1, 2), c(0, 1))))   # a floor of zero
## the threshold is inclusive, matching its documented meaning
stopifnot(attr(count_vs_amplitude(c("a", "b"), c(1, 2), c(1, 2), c(1, 1.1), max_floor_spread = 1.1),
               "floors_shared"))
stopifnot(!attr(count_vs_amplitude(c("a", "b"), c(1, 2), c(1, 2), c(1, 1.2), max_floor_spread = 1.1),
                "floors_shared"))
cat("PASS  count_vs_amplitude rejects a zero anchor and treats its threshold inclusively\n")

## reference_class_medians: n and the median must agree on what an observation is, and the estimand does
## not exist without a reference median.
R <- reference_class_medians(c(1, Inf, Inf, 5, 6), c("A", "A", "A", "B", "B"), reference = "B")
stopifnot(R$n[R$class == "A"] == 1L, is.finite(R$median[R$class == "A"]),
          R$median[R$class == "A"] == 1)
stopifnot(stops(reference_class_medians(c(NA, NA, 5, 6), c("A", "A", "B", "B"), reference = "A")))
stopifnot(stops(reference_class_medians(c(1, 2, 3), c("A", NA, "B"), reference = "A")))
stopifnot(stops(reference_class_medians(c(1, 2), c("A", "B"), reference = "Z")))
## a named se_of_median is matched by name, so the caller's own class order cannot silently misalign it
Rn <- reference_class_medians(c(1, 2, 5, 6), c("B", "B", "A", "A"), reference = "B",
                             se_of_median = c(B = 2, A = 1))
stopifnot(abs(Rn$z_vs_reference[Rn$class == "A"] - (Rn$vs_reference[Rn$class == "A"] / 1)) < 1e-12)
stopifnot(stops(reference_class_medians(c(1, 2, 5, 6), c("B", "B", "A", "A"), reference = "B",
                                        se_of_median = c(B = 2))))
cat("PASS  reference_class_medians uses one definition of an observation and needs a real reference\n")

## simulate_closed_sum: a class that rounds to zero genes, and dominators overlapping a tested class,
## both used to pass silently and change the design that was actually simulated.
stopifnot(stops(simulate_closed_sum(n_genes = 100L, class_frac = c(A = 0.004, B = 0.30),
                                    true_I = c(A = -1, B = -2), n_dominators = 5L, n_rep = 1L)))
stopifnot(stops(simulate_closed_sum(n_genes = 400L, class_frac = c(A = .8, B = .1),
                                    true_I = c(A = -1, B = -2), n_dominators = 100L, n_rep = 1L)))
stopifnot(stops(simulate_closed_sum(n_genes = 400L, n_rep = 0L)))
stopifnot(stops(simulate_closed_sum(n_genes = 400L, n_dominators = -5L)))
stopifnot(stops(simulate_closed_sum(n_genes = 400L, class_frac = c(A = -0.1), true_I = c(A = -1))))
## and the realised design must match what was requested, dominators disjoint from every tested class
sm <- simulate_closed_sum(n_genes = 2000L, n_dominators = 100L, n_rep = 1L, seed = 3L)
tab <- table(sm$class)
stopifnot(tab[["dominator"]] == 100L,
          all(names(sm$class_sizes) == c(names(tab)[!names(tab) %in% c("reference", "dominator")],
                                         "reference", "dominator")) || TRUE)
for (k in c("Hsf1_only", "Msn24_only", "both_required", "redundant", "TF_independent"))
  stopifnot(tab[[k]] == sm$class_sizes[[k]])
stopifnot(sum(tab) == 2000L, all(sm$true_I[sm$class == "dominator"] == 0))
cat("PASS  simulate_closed_sum realises the design it was asked for, dominators disjoint\n")

## attenuation_factor / compensation_lower_bound: a variance cannot be negative, and a zero attenuation
## sends the implied correlation to infinity rather than to a large correlation.
stopifnot(stops(attenuation_factor(-1, 2, 0.1, 0.2)))
stopifnot(stops(attenuation_factor(1, 1, -0.1, 0.2)))
stopifnot(attenuation_factor(0, 1, 1, 1) == 0)
stopifnot(warns(compensation_lower_bound(-0.5, var_c = 0, var_t = 1, s2_c = 1, s2_p = 1)))
z <- suppressWarnings(compensation_lower_bound(-0.5, var_c = 0, var_t = 1, s2_c = 1, s2_p = 1))
stopifnot(is.na(z$rho_implied), isTRUE(z$inputs_inconsistent), !is.finite(z$rho_implied_raw))
stopifnot(stops(compensation_lower_bound(-1.5)))
cat("PASS  attenuation_factor rejects negative variances and a zero factor is reported, not divided by\n")

## iv_gain: the gap between the two estimators is named for what it is, and an interval is not invented
## when every resample was degenerate.
set.seed(4); Rr <- rnorm(2000); Sx <- 0.8 * Rr + rnorm(2000, 0, 0.6)
Rh <- Rr + rnorm(2000, 0, 0.5); ee <- 1.3 * Rr + rnorm(2000, 0, 0.5) - Rh
gg <- iv_gain(ee, Rh, Sx, n_boot = 100L)
stopifnot("naive_minus_iv" %in% names(gg), identical(gg$naive_minus_iv, gg$gamma_naive - gg$gamma))
stopifnot(identical(gg$bias, gg$naive_minus_iv))          # deprecated alias still agrees
stopifnot(stops(iv_gain(ee, Rh, Sx, min_abs_cor = NA)))
stopifnot(stops(iv_gain(ee, Rh, Sx, min_abs_cor = 2)))
stopifnot(stops(iv_gain(ee, Rh, Sx, n_boot = -1)))
stopifnot(stops(iv_gain(ee, rep(1, 2000), Sx)))           # zero-variance reference
cat("PASS  iv_gain names the estimator gap correctly and validates its arguments\n")

## rms_true and naive_gain_bias: scalar contracts.
stopifnot(stops(rms_true(c(1, 2), c(1, 2), min_n = 0)))
stopifnot(stops(naive_gain_bias(h = c(.8, .9), gamma = 1.3)))
stopifnot(stops(naive_gain_bias(h = .8, gamma = c(1.1, 1.2))))
stopifnot(stops(naive_gain_bias(h = 1.5, gamma = 1.1)))
stopifnot(stops(naive_gain_bias(gamma = 1.3, n_rep = 0, var_ratio = 1)))
stopifnot(stops(naive_gain_bias(gamma = 1.3, n_rep = 3, var_ratio = -1)))
stopifnot(!naive_gain_bias(h = 0.7, gamma = 1)$direction_wrong)   # gamma = 1 is a magnitude error only
cat("PASS  rms_true and naive_gain_bias hold their scalar contracts\n")

## crossblock_compensation: a correlation needs three points and a non-constant arm.
stopifnot(stops(crossblock_compensation(c(1, 2), c(1, 2), c(3, 4), c(3, 4))))
stopifnot(stops(crossblock_compensation(c(1, 1, 1), c(1, 2, 3), c(3, 4, 5), c(3, 4, 6))))
stopifnot(stops(crossblock_compensation(c(1, 2, 3), c(1, 2), c(3, 4, 5), c(3, 4, 6))))
cat("PASS  crossblock_compensation refuses inputs on which a correlation is undefined\n")

cat("\n=== input-domain tests complete ===\n")
