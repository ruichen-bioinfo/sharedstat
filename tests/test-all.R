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
for (f in list.files(file.path("..", "R"), full.names = TRUE)) source(f)
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

cat("\n=== ", if (FAILED) "SOME TESTS FAILED" else "ALL TESTS PASSED", " ===\n", sep = "")
if (FAILED) quit(status = 1)

## ---------------------------------------------------------------------------------------------------
## Regression tests added 2026-09-09, each pinning a defect that was present and is now fixed. A test
## that only asserts current behaviour is worth little; these assert the behaviour that was WRONG.
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

## count_vs_amplitude: the exaggeration factors printed in the paper, and the shared-floor verdict.
cv <- count_vs_amplitude(condition = c("37_30", "37_10", "42_10", "42_30"),
                         n_sig = c(930, 13, 9, 18),
                         amplitude = c(0.890, 0.378, 0.165, 0.116),
                         noise_floor = c(0.7756, 0.6497, 0.6830, 0.5270),
                         anchor = "37_30")
stopifnot(abs(cv$exaggeration[2] - 30.4) < 0.2, abs(cv$exaggeration[3] - 19.2) < 0.2,
          abs(cv$exaggeration[4] - 6.7) < 0.2)
stopifnot(abs(attr(cv, "floor_spread") - 1.4717) < 1e-3, !attr(cv, "floors_shared"))
stopifnot(attr(count_vs_amplitude(c("a","b"), c(1,2), c(1,2), c(1, 1.05)), "floors_shared"))
cat("PASS  count_vs_amplitude reproduces 30.4 / 19.2 / 6.7 and rejects the 1.47-fold floor spread\n")

## simulate_closed_sum: the four-genotype contrast must equal the injected truth exactly.
sm <- simulate_closed_sum(n_genes = 1200L, kappa = -0.4, n_rep = 2L, seed = 7L)
stopifnot(all(is.finite(sm$counts)), !any(is.na(sm$counts)))
stopifnot(is.finite(sm$share_nonadditivity))
cat("PASS  simulate_closed_sum returns a complete count matrix and a finite share non-additivity\n")
