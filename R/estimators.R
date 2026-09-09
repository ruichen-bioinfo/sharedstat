#' Noise-corrected amplitude from a signal and a null contrast
#'
#' `rms_true(s, n) = sqrt(max(0, mean(s^2) - mean(n^2)))`. The subtraction removes the variance the
#' null contrast contributes, so the result is an amplitude rather than a dispersion. It is clamped at
#' zero: when the null exceeds the signal the honest answer is "not distinguishable from noise", not a
#' negative amplitude, and a returned zero must be reported as such rather than treated as a small
#' positive number.
#'
#' The subtraction is only meaningful when both means are taken over the **same** features. Dropping
#' missing values independently from `s` and from `n` silently compares a signal computed on one feature
#' set against a noise floor computed on another, and the error is not symmetric: it inflates the
#' amplitude whenever the features missing from `n` are the loud ones. The function therefore pairs the
#' inputs and uses complete cases only. An earlier version of this function did not, and on
#' `s = c(1,2,3,4,10)`, `n = c(1,2,3,4,NA)` it returned 4.30 where the paired answer is exactly 0 --
#' "not distinguishable from noise" reported as a large positive amplitude.
#'
#' In the analysis this package accompanies the hazard could not fire, because there the signal and null
#' contrasts are two linear combinations of the same four vectors and so are missing on exactly the same
#' features. That is a property of that pipeline, not of the estimator, which is why the guard belongs
#' here: the next caller's contrasts may be assembled separately.
#'
#' @param s numeric, the signal contrast (for two replicate blocks, `(b1 + b2) / 2`).
#' @param n numeric, the null contrast built to contain no effect (for two blocks, `(b2 - b1) / 2`).
#' @param min_n minimum number of complete pairs required; below this the answer is not estimable and an
#'   error is raised rather than a number returned.
#' @return A single non-negative number with attributes `n_pairs` and `n_dropped`. Exactly zero means the
#'   signal did not exceed its own noise and must be reported as such, not as a small positive number.
rms_true <- function(s, n, min_n = 2L) {
  s <- as.numeric(s); n <- as.numeric(n)
  if (length(s) != length(n))
    stop(sprintf("s and n must be the same length (got %d and %d): the noise subtraction is defined per feature",
                 length(s), length(n)))
  keep <- is.finite(s) & is.finite(n)
  if (sum(keep) < min_n)
    stop(sprintf("only %d complete (s, n) pairs; need at least %d", sum(keep), min_n))
  out <- sqrt(max(0, mean(s[keep]^2) - mean(n[keep]^2)))
  attr(out, "n_pairs") <- sum(keep)
  attr(out, "n_dropped") <- sum(!keep)
  out
}

#' Analytic bias of the naive gain estimator against a shared reference
#'
#' The bias has a closed form, and writing it down is more useful than tabulating simulations of it.
#' With true reference signal \eqn{R}, measured \eqn{\hat R = R + u}, response arm \eqn{\gamma R}
#' measured with its own error, and \eqn{h = \mathrm{var}(R) / \mathrm{var}(\hat R)}:
#'
#' \deqn{\mathrm{slope} = (\gamma - 1)h - (1 - h), \qquad \hat\gamma_{naive} = 1 + \mathrm{slope} = \gamma h.}
#'
#' Three consequences follow, and they are not the same consequence stated three ways.
#' \itemize{
#'   \item At \eqn{\gamma = 1} -- no divergence whatsoever -- the naive estimator returns \eqn{h}. Any
#'     \eqn{h < 1} therefore produces apparent divergence from nothing.
#'   \item The bias is \eqn{O(1)} in the number of **features**. Measuring more genes makes the wrong
#'     value more precise; it does not move it.
#'   \item The bias is \eqn{O(1/n_{rep})} in the number of **replicates**, because replication shrinks
#'     the reference arm's measurement variance and so raises \eqn{h}. Replicates help, but only through
#'     \eqn{h}, and no finite number of them reaches \eqn{h = 1}.
#' }
#' The second and third points are easy to conflate, and conflating them gives the false statement that
#' replication cannot help at all.
#'
#' @param h reliability of the reference arm, `var(R) / var(R_hat)`, in (0, 1].
#' @param gamma the true gain.
#' @param n_rep optional replicate count; when given with `var_ratio`, `h` is computed rather than taken.
#' @param var_ratio optional `var(u_per_observation) / var(R)`; with `n_rep` gives `h = 1/(1 + var_ratio/n_rep)`.
#' @return A list with `h`, `gamma_naive`, `bias`, and `direction_wrong` (whether the naive estimate and
#'   the truth fall on opposite sides of 1, which is the failure that changes conclusions).
naive_gain_bias <- function(h = NULL, gamma, n_rep = NULL, var_ratio = NULL) {
  if (is.null(h)) {
    if (is.null(n_rep) || is.null(var_ratio))
      stop("supply either h, or both n_rep and var_ratio")
    h <- 1 / (1 + var_ratio / n_rep)
  }
  stopifnot(all(h > 0), all(h <= 1))
  gn <- gamma * h
  list(h = h, gamma_naive = gn, bias = gn - gamma,
       direction_wrong = (gamma > 1 & gn < 1) | (gamma < 1 & gn > 1))
}

#' Instrumental-variable gain estimator for a shared reference arm
#'
#' When two contrasts are formed against the same reference, the reference's error enters both with
#' opposite sign and the ordinary slope of one on the other acquires a bias that does not vanish with
#' sample size. Using a third quantity that shares the reference's signal but not its error as an
#' instrument removes that bias:  `gamma = 1 + cov(e_i, S) / cov(R, S)`.
#'
#' The estimator does not merely rescale: in the design reported in the accompanying paper it
#' **reorders** the units being compared. That is the practically important consequence, because a bias
#' that rescaled everything would leave comparative conclusions intact.
#'
#' @param e numeric, the contrast of interest (response minus shared reference).
#' @param R numeric, the shared reference arm.
#' @param S numeric, the instrument: correlated with the reference's signal, independent of its error.
#' @param naive logical; if `TRUE`, also return the naive OLS slope for comparison.
#' @param min_abs_cor the smallest `|cor(R, S)|` at which the instrument is treated as adequate. Below it a
#'   warning is raised, because the estimator is a ratio of covariances whose denominator is then near zero:
#'   it returns a number rather than failing, and that number can be far from the truth.
#' @param n_boot number of bootstrap resamples used to attach a standard error and a percentile interval.
#'   Zero, the default, skips it. A ratio of covariances has no useful closed-form standard error at small
#'   samples, so an interval is resampled rather than approximated.
#' @param seed random seed for the bootstrap, so a reported interval is reproducible.
#' @return A list with `gamma`, `n`, `cor_RS` and `instrument_strength`; when `naive = TRUE` also
#'   `gamma_naive` and `bias` (their difference); and when `n_boot > 0` also `se`, `ci` and
#'   `n_boot_failed` (resamples in which the instrument was degenerate).
iv_gain <- function(e, R, S, naive = TRUE, min_abs_cor = 0.1, n_boot = 0L, seed = 1L) {
  e <- as.numeric(e); R <- as.numeric(R); S <- as.numeric(S)
  if (length(unique(c(length(e), length(R), length(S)))) != 1L)
    stop("e, R and S must be the same length: the instrument is applied per feature")
  keep <- is.finite(e) & is.finite(R) & is.finite(S)
  e <- e[keep]; R <- R[keep]; S <- S[keep]
  if (length(e) < 3L) stop(sprintf("only %d complete cases", length(e)))
  den <- stats::cov(R, S)
  if (!is.finite(den) || den == 0)
    stop("cov(R, S) is zero or undefined: the instrument is not informative here")
  ## An instrument that is merely WEAK rather than absent is the dangerous case, because the estimator
  ## returns a number instead of failing. cov(R, S) sits in the denominator, so as it approaches zero the
  ## estimate diverges: with cor(R, S) = -0.02 this estimator returns 24.7 for a true gain of 1.3. The
  ## strength of the first stage is therefore reported always and warned about below a stated threshold,
  ## rather than left for the caller to think of.
  rho_RS <- stats::cor(R, S)
  strength <- if (abs(rho_RS) >= min_abs_cor) "adequate" else "WEAK"
  if (strength == "WEAK")
    warning(sprintf(paste("weak instrument: |cor(R, S)| = %.4f < %.2f. The IV estimate is the ratio of two",
                          "covariances and its denominator is near zero, so it is unstable here and may be",
                          "far from gamma. Treat it as uninformative rather than as a corrected value."),
                    abs(rho_RS), min_abs_cor), call. = FALSE)
  g <- 1 + stats::cov(e, S) / den
  out <- list(gamma = g, n = length(e), cor_RS = rho_RS, instrument_strength = strength)
  if (naive) {
    gn <- 1 + stats::cov(e, R) / stats::var(R)
    out$gamma_naive <- gn
    out$bias <- gn - g
  }
  ## A ratio-of-covariances has no useful closed-form standard error at small n, so when an interval is
  ## wanted it is resampled rather than approximated.
  if (n_boot > 0L) {
    set.seed(seed)
    n <- length(e)
    bs <- vapply(seq_len(n_boot), function(b) {
      i <- sample.int(n, n, replace = TRUE)
      d <- stats::cov(R[i], S[i])
      if (!is.finite(d) || d == 0) NA_real_ else 1 + stats::cov(e[i], S[i]) / d
    }, 0)
    out$se <- stats::sd(bs, na.rm = TRUE)
    out$ci <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    out$n_boot_failed <- sum(!is.finite(bs))
  }
  out
}

#' Cross-block estimator of compensatory divergence, and its attenuation
#'
#' The correlation between a component and a residual formed by subtracting it is negative by
#' construction. Estimating the two from *different* replicate blocks removes the shared error that
#' creates most of that artefact. What remains is attenuated by measurement error, in a known direction:
#' the estimate is therefore a **lower bound** in magnitude on the underlying correlation.
#'
#' Two consequences follow and are easy to get wrong. A weaker reported compensation does not license
#' the inference that compensation is weaker, because it may reflect fewer replicates. And compensation
#' magnitudes are **not comparable across designs with different replicate structure**.
#'
#' @param cis_b1,cis_b2 numeric, the component estimated in block 1 and block 2.
#' @param par_b1,par_b2 numeric, the total estimated in block 1 and block 2.
#' @return A list with `r_cross` (the cross-block estimate), `r_naive` (both from one block, for
#'   comparison) and `artefact` (their difference).
crossblock_compensation <- function(cis_b1, cis_b2, par_b1, par_b2) {
  tr_b2 <- par_b2 - cis_b2
  tr_b1 <- par_b1 - cis_b1
  keep <- is.finite(cis_b1) & is.finite(tr_b2) & is.finite(cis_b2) & is.finite(tr_b1)
  r_cross <- stats::cor(cis_b1[keep], tr_b2[keep])
  r_naive <- stats::cor(cis_b1[keep], tr_b1[keep])
  list(r_cross = r_cross, r_naive = r_naive, artefact = r_naive - r_cross, n = sum(keep))
}

#' Attenuation factor for a cross-block correlation
#'
#' `cor(hat c, hat t) = rho * sd(c) * sd(t) / (sqrt(var(c) + s2_c) * sqrt(var(t) + s2_p + s2_c))`.
#' The factor is at most 1, so the observed magnitude is a lower bound on `|rho|`.
#'
#' @param var_c,var_t true variances of the two components.
#' @param s2_c,s2_p measurement-error variances of the component and of the total.
#' @return The multiplicative attenuation factor, in (0, 1].
attenuation_factor <- function(var_c, var_t, s2_c, s2_p) {
  sqrt(var_c) * sqrt(var_t) / (sqrt(var_c + s2_c) * sqrt(var_t + s2_p + s2_c))
}

#' State a reported compensation as the lower bound it is
#'
#' @param r_observed the cross-block estimate.
#' @param ... passed to [attenuation_factor()]; if supplied, the implied `rho` is also returned.
#' @return A list with `lower_bound` (|r_observed|) and, when the variance components are given,
#'   `rho_implied`. Never returns a point estimate of `rho` without them, because there is none.
compensation_lower_bound <- function(r_observed, ...) {
  out <- list(lower_bound = abs(r_observed),
              statement = sprintf("|rho| >= %.4f; a point estimate of rho is not identified without the variance components",
                                  abs(r_observed)))
  a <- list(...)
  if (all(c("var_c", "var_t", "s2_c", "s2_p") %in% names(a))) {
    f <- attenuation_factor(a$var_c, a$var_t, a$s2_c, a$s2_p)
    out$attenuation <- f
    ri <- r_observed / f
    ## Dividing by an attenuation factor can carry the implied correlation past 1, which no correlation
    ## can be. That is not a rounding artefact: it says the supplied variance components are inconsistent
    ## with the observed correlation, so one of the inputs is wrong. Returning the impossible number
    ## silently would let it propagate into a manuscript.
    if (is.finite(ri) && abs(ri) > 1) {
      warning(sprintf(paste("implied rho = %.4f lies outside [-1, 1], so the supplied variance components",
                            "are inconsistent with r_observed = %.4f (attenuation %.4f). Reporting it as",
                            "not identified rather than as a correlation."), ri, r_observed, f),
              call. = FALSE)
      out$rho_implied <- NA_real_
      out$rho_implied_raw <- ri
      out$inputs_inconsistent <- TRUE
    } else {
      out$rho_implied <- ri
      out$inputs_inconsistent <- FALSE
    }
  }
  out
}
