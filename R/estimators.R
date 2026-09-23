#' Noise-corrected amplitude from a signal and a null contrast
#'
#' `rms_true(s, n) = sqrt(max(0, mean(s^2) - mean(n^2)))`. The subtraction removes the variance the
#' null contrast contributes, so the result is an amplitude rather than a dispersion. It is clamped at
#' zero: when the null exceeds the signal the answer is "not distinguishable from noise", not a negative
#' amplitude, and a returned zero should be reported as such rather than treated as a small positive
#' number.
#'
#' The subtraction is only defined when both means are taken over the same features. Dropping missing
#' values independently from `s` and from `n` would compare a signal computed on one feature set against
#' a noise floor computed on another, and the error is not symmetric: it inflates the amplitude whenever
#' the features missing from `n` are the loud ones. The inputs are therefore paired and complete cases
#' used.
#'
#' @param s numeric, the signal contrast (for two replicate blocks, `(b1 + b2) / 2`).
#' @param n The null contrast, carrying error and no signal. A numeric vector for a two-replicate
#'   design, or a matrix or data frame with one column per orthogonal null contrast when there are
#'   three or more replicates; in the matrix case the per-contrast mean squares are averaged.
#' @param min_n minimum number of complete pairs required; below this the quantity is not estimable and
#'   an error is raised rather than a number returned. Must be a single integer of at least 2, since a
#'   mean square over one pair carries no information about dispersion.
#' @return A single non-negative number with attributes `n_pairs` and `n_dropped`. Exactly zero means the
#'   signal did not exceed its own noise.
rms_true <- function(s, n, min_n = 2L) {
  if (length(min_n) != 1L || !is.finite(min_n) || min_n < 2)
    stop("min_n must be a single finite number of at least 2")
  s <- as.numeric(s)
  ## `n` may be a matrix or data frame with one column per orthogonal null contrast. With three or more
  ## replicates there are two or more such contrasts, and they are averaged as mean squares -- one mean square
  ## per contrast, then the mean of those -- rather than concatenated into one long vector. Concatenating would
  ## weight a contrast by how many features it happens to have, which is the same for all of them here but need
  ## not be once features are dropped for missingness.
  if (is.matrix(n) || is.data.frame(n)) {
    nm <- as.matrix(n); storage.mode(nm) <- "double"
    if (nrow(nm) != length(s))
      stop(sprintf("s has %d features and the null matrix has %d rows: the noise subtraction is per feature",
                   length(s), nrow(nm)))
    keep <- is.finite(s) & apply(is.finite(nm), 1L, all)
    if (sum(keep) < min_n)
      stop(sprintf("only %d complete rows; need at least %d", sum(keep), min_n))
    ms_nul <- mean(apply(nm[keep, , drop = FALSE], 2L, function(v) mean(v^2)))
    out <- sqrt(max(0, mean(s[keep]^2) - ms_nul))
    attr(out, "n_pairs") <- sum(keep)
    attr(out, "n_dropped") <- sum(!keep)
    attr(out, "n_null_contrasts") <- ncol(nm)
    return(out)
  }
  n <- as.numeric(n)
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
#' With true reference signal \eqn{R}, measured \eqn{\hat R = R + u}, response arm \eqn{\gamma R}
#' measured with its own error, and \eqn{h = \mathrm{var}(R) / \mathrm{var}(\hat R)}:
#'
#' \deqn{\mathrm{slope} = (\gamma - 1)h - (1 - h), \qquad \hat\gamma_{naive} = 1 + \mathrm{slope} = \gamma h,}
#' which holds when the two arms' measurement errors are independent. Retaining their covariance,
#' \deqn{\mathbb{E}[\hat\gamma_{naive}] = \gamma h + \mathrm{cov}(u, v) / \mathrm{var}(\hat R),}
#' so the first expression is a baseline rather than the general case.
#'
#' Three consequences follow, and they are not the same consequence stated three ways.
#' \itemize{
#'   \item At \eqn{\gamma = 1} the naive estimator returns \eqn{h}, so any \eqn{h < 1} produces apparent
#'     divergence where there is none.
#'   \item The bias is \eqn{O(1)} in the number of features. Measuring more features makes the wrong
#'     value more precise; it does not move it.
#'   \item The bias is \eqn{O(1/n_{rep})} under independent replicate averaging, because averaging
#'     shrinks the reference arm's measurement variance and so raises \eqn{h}. Replicates help, but only
#'     through \eqn{h}, and no finite number of them reaches \eqn{h = 1}.
#' }
#' The second and third points are easy to conflate, and conflating them gives the false statement that
#' replication cannot help at all.
#'
#' @param h reliability of the reference arm, `var(R) / var(R_hat)`, a single number in [0, 1]. Zero is
#'   admissible and is what [reliability_from_blocks()] and [reliability_from_replicates()] return when the
#'   between-block variation is at least as large as the variation across features. It says that no
#'   reproducible signal was resolved above that contrast's own noise, which is an uninformative reference
#'   rather than a negative reliability, and not a statement that the biological signal is absent. At `h = 0`
#'   the baseline expectation is 0 and the general one is `arm_error_cov_ratio`.
#' @param gamma the true gain, a single finite number.
#' @param n_rep optional replicate count; with `var_ratio`, `h` is computed rather than taken.
#' @param arm_error_cov_ratio the two arms' measurement-error covariance divided by the variance of the
#'   observed reference, `cov(u, v) / var(R_hat)`. Zero, the default, is the independent-arm-error baseline;
#'   the returned `gamma_naive` is then `gamma * h`. A non-zero value shifts both the expectation and the
#'   reversal threshold, which becomes `(1 - arm_error_cov_ratio) / gamma`. Sharing a library between the two
#'   arms permits a non-zero value but does not fix its sign, so it is a quantity to report rather than assume.
#'   Note that `simulate_closed_sum` has an unrelated argument named `kappa`; the two are different quantities
#'   and are deliberately not given the same name.
#'
#' @section The reversal threshold at the edges:
#' `reversal_threshold_h` is `(1 - arm_error_cov_ratio) / gamma` and is reported unclamped. A value above 1
#' means every admissible reliability lies below the threshold, so the reversal is unavoidable for that gain;
#' a value below 0 means none does, so it is impossible. Clamping the number to [0, 1] would hide which of
#' those two situations holds.
#' @param var_ratio optional `var(u_per_observation) / var(R)`; with `n_rep` gives `h = 1/(1 + var_ratio/n_rep)`.
#' @return A list with `h`, `gamma_naive`, `bias` (`gamma_naive - gamma`, a true bias here because
#'   `gamma` is supplied) and `direction_wrong` (whether the naive estimate and the truth fall on
#'   opposite sides of 1, which is the failure that changes conclusions).
naive_gain_bias <- function(h = NULL, gamma, n_rep = NULL, var_ratio = NULL, arm_error_cov_ratio = 0) {
  if (length(gamma) != 1L || !is.finite(gamma))
    stop("gamma must be a single finite number")
  if (is.null(h)) {
    if (is.null(n_rep) || is.null(var_ratio))
      stop("supply either h, or both n_rep and var_ratio")
    if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1)
      stop("n_rep must be a single finite number of at least 1")
    if (length(var_ratio) != 1L || !is.finite(var_ratio) || var_ratio < 0)
      stop("var_ratio must be a single finite non-negative number")
    h <- 1 / (1 + var_ratio / n_rep)
  }
  ## h is the single parameter the whole result turns on, so it is required to be scalar: a vectorised h
  ## would return a vector from a function documented to describe one estimator, and the caller would have
  ## no signal that they had asked for something else.
  ## The interval is closed at zero. reliability_from_blocks() and reliability_from_replicates() clamp a negative
  ## variance estimate to h = 0, so rejecting it here would mean this package's own diagnostic output could not be
  ## passed to its own bias function at the least informative boundary a real dataset can produce.
  if (length(h) != 1L || !is.finite(h) || h < 0 || h > 1)
    stop("h must be a single finite number in [0, 1]")
  ## The expectation of the naive estimator is gamma * h only when the two arms' measurement errors are
  ## independent. In general it carries a third term, the arm-error covariance divided by the variance of the
  ## observed reference: E[gamma_hat] = gamma * h + arm_error_cov_ratio. Sharing a library between the arms
  ## permits a non-zero value without fixing its sign, so the default of zero is the independent-arm-error
  ## baseline and is labelled as such rather than presented as the general case.
  if (length(arm_error_cov_ratio) != 1L || !is.finite(arm_error_cov_ratio))
    stop("arm_error_cov_ratio must be a single finite number")
  baseline <- gamma * h
  gn <- baseline + arm_error_cov_ratio
  list(h = h, gamma_naive = gn, gamma_naive_baseline = baseline,
       arm_error_cov_ratio = arm_error_cov_ratio, bias = gn - gamma,
       reversal_threshold_h = (1 - arm_error_cov_ratio) / gamma,
       direction_wrong = (gamma > 1 && gn < 1) || (gamma < 1 && gn > 1))
}

#' Instrumental-variable gain estimator for a shared reference arm
#'
#' When two contrasts are formed against the same reference, the reference's error enters both with
#' opposite sign and the ordinary slope of one on the other acquires a bias that does not vanish with
#' sample size. Using a third quantity that shares the reference's signal but not its error as an
#' instrument removes that bias: `gamma = 1 + cov(e, S) / cov(R, S)`.
#'
#' The estimator requires an instrument `S` with `cor(S, R) != 0` and `S` independent of both measurement
#' errors. Neither condition can be verified from `(e, R, S)` alone; the first stage is reported so that
#' the part which can be checked is visible.
#'
#' @param e numeric, the contrast of interest (response minus shared reference).
#' @param R numeric, the shared reference arm as measured.
#' @param S numeric, the instrument.
#' @param naive logical; if `TRUE`, also return the naive OLS slope for comparison.
#' @param min_abs_cor the smallest `|cor(R, S)|` at which the instrument is treated as adequate, in
#'   \[0, 1\]. Below it a warning is raised, because the estimator is a ratio of covariances whose
#'   denominator is then near zero: it returns a number rather than failing, and that number can be far
#'   from the truth. The default of 0.1 is a convention, not a derived threshold.
#' @param n_boot number of bootstrap resamples used to attach a standard error and a percentile interval.
#'   Zero, the default, skips it. A ratio of covariances has no useful closed-form standard error at small
#'   samples, so an interval is resampled rather than approximated.
#' @param seed random seed for the bootstrap, so a reported interval is reproducible.
#' @return A list with `gamma`, `n`, `cor_RS` and `instrument_strength`; when `naive = TRUE` also
#'   `gamma_naive` and `naive_minus_iv`, their difference. That difference is **not** a bias: the true
#'   gain is unknown here, so it measures how far the naive estimate sits from the instrumented one, and
#'   nothing more. The field `bias` is retained as a deprecated alias of `naive_minus_iv` and will be
#'   removed. When `n_boot > 0` the list also has `se`, `ci` and `n_boot_failed`; if every resample was
#'   degenerate, `se` and `ci` are `NA` and a warning is raised.
iv_gain <- function(e, R, S, naive = TRUE, min_abs_cor = 0.1, n_boot = 0L, seed = 1L) {
  if (length(min_abs_cor) != 1L || !is.finite(min_abs_cor) || min_abs_cor < 0 || min_abs_cor > 1)
    stop("min_abs_cor must be a single finite number in [0, 1]")
  if (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 0)
    stop("n_boot must be a single finite non-negative number")
  e <- as.numeric(e); R <- as.numeric(R); S <- as.numeric(S)
  if (length(unique(c(length(e), length(R), length(S)))) != 1L)
    stop("e, R and S must be the same length: the instrument is applied per feature")
  keep <- is.finite(e) & is.finite(R) & is.finite(S)
  e <- e[keep]; R <- R[keep]; S <- S[keep]
  if (length(e) < 3L) stop(sprintf("only %d complete cases; need at least 3", length(e)))
  if (stats::var(R) == 0 || stats::var(S) == 0)
    stop("R or S has zero variance, so neither the first stage nor the naive slope is defined")
  den <- stats::cov(R, S)
  if (!is.finite(den) || den == 0)
    stop("cov(R, S) is zero or undefined: the instrument is not informative here")
  ## An instrument that is merely WEAK rather than absent is the dangerous case, because the estimator
  ## returns a number instead of failing. cov(R, S) sits in the denominator, so as it approaches zero the
  ## estimate diverges. The strength of the first stage is therefore reported always and warned about
  ## below a stated threshold, rather than left for the caller to think of.
  rho_RS <- stats::cor(R, S)
  strength <- if (abs(rho_RS) >= min_abs_cor) "adequate" else "WEAK"
  if (identical(strength, "WEAK"))
    warning(sprintf(paste("weak instrument: |cor(R, S)| = %.4f < %.2f. The IV estimate is the ratio of two",
                          "covariances and its denominator is near zero, so it is unstable here and may be",
                          "far from gamma. Treat it as uninformative rather than as a corrected value."),
                    abs(rho_RS), min_abs_cor), call. = FALSE)
  g <- 1 + stats::cov(e, S) / den
  out <- list(gamma = g, n = length(e), cor_RS = rho_RS, instrument_strength = strength)
  if (naive) {
    gn <- 1 + stats::cov(e, R) / stats::var(R)
    out$gamma_naive <- gn
    ## Named for what it is. This function never sees the true gain, so the gap between the two estimates
    ## is a difference between estimators and not a bias, and calling it `bias` invited the reading that
    ## the package had reported how wrong the naive estimate was.
    out$naive_minus_iv <- gn - g
    out$bias <- gn - g            # deprecated alias, to be removed
  }
  if (n_boot > 0L) {
    set.seed(seed)
    nn <- length(e)
    bs <- vapply(seq_len(n_boot), function(b) {
      i <- sample.int(nn, nn, replace = TRUE)
      d <- stats::cov(R[i], S[i])
      if (!is.finite(d) || d == 0) NA_real_ else 1 + stats::cov(e[i], S[i]) / d
    }, 0)
    n_bad <- sum(!is.finite(bs))
    out$n_boot_failed <- n_bad
    if (n_bad == n_boot) {
      ## Every resample was degenerate. Returning sd() and quantile() of an all-NA vector would hand back
      ## an estimate with an interval-shaped object beside it that carries no information.
      warning(sprintf(paste("all %d bootstrap resamples were degenerate, so no standard error or interval",
                            "could be formed. The point estimate is returned without one."), n_boot),
              call. = FALSE)
      out$se <- NA_real_
      out$ci <- c(NA_real_, NA_real_)
    } else {
      out$se <- stats::sd(bs, na.rm = TRUE)
      out$ci <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    }
  }
  out
}

#' Cross-block estimator of compensatory divergence, and its attenuation
#'
#' The correlation between a component and a residual formed by subtracting it is negative by
#' construction. Estimating the two from *different* replicate blocks removes the shared error that
#' creates most of that artefact. What remains is attenuated by measurement error, in a known direction:
#' the estimate is therefore a lower bound in magnitude on the underlying correlation.
#'
#' Two consequences are easy to get wrong. A weaker reported compensation does not license the inference
#' that compensation is weaker, because it may reflect fewer replicates. And compensation magnitudes are
#' not comparable across designs with different replicate structure.
#'
#' @param cis_b1,cis_b2 numeric, the component estimated in block 1 and block 2.
#' @param par_b1,par_b2 numeric, the total estimated in block 1 and block 2.
#' @param min_n minimum number of complete cases required; a correlation on fewer than three points is
#'   either undefined or determined by the points themselves.
#' @return A list with `r_cross` (the cross-block estimate), `r_naive` (both from one block, for
#'   comparison), `artefact` (their difference) and `n`.
crossblock_compensation <- function(cis_b1, cis_b2, par_b1, par_b2, min_n = 3L) {
  if (length(unique(c(length(cis_b1), length(cis_b2), length(par_b1), length(par_b2)))) != 1L)
    stop("cis_b1, cis_b2, par_b1 and par_b2 must be the same length")
  tr_b2 <- par_b2 - cis_b2
  tr_b1 <- par_b1 - cis_b1
  keep <- is.finite(cis_b1) & is.finite(tr_b2) & is.finite(cis_b2) & is.finite(tr_b1)
  if (sum(keep) < min_n)
    stop(sprintf("only %d complete cases; need at least %d for a correlation", sum(keep), min_n))
  ## cor() returns NA with a warning on a constant vector. Saying so here identifies which arm is
  ## constant, which the caller cannot recover from an NA.
  vs <- c(cis_b1 = stats::var(cis_b1[keep]), tr_b2 = stats::var(tr_b2[keep]),
          tr_b1 = stats::var(tr_b1[keep]))
  if (any(vs == 0))
    stop(sprintf("zero variance in: %s. A correlation is not defined on a constant arm.",
                 paste(names(vs)[vs == 0], collapse = ", ")))
  r_cross <- stats::cor(cis_b1[keep], tr_b2[keep])
  r_naive <- stats::cor(cis_b1[keep], tr_b1[keep])
  list(r_cross = r_cross, r_naive = r_naive, artefact = r_naive - r_cross, n = sum(keep))
}

#' Attenuation factor for a cross-block correlation
#'
#' `cor(hat c, hat t) = rho * sd(c) * sd(t) / (sqrt(var(c) + s2_c) * sqrt(var(t) + s2_p + s2_c))`.
#' The factor is at most 1, so an observed magnitude is a lower bound on `|rho|`.
#'
#' @param var_c,var_t true variances of the two components; non-negative.
#' @param s2_c,s2_p measurement-error variances of the component and of the total; non-negative.
#' @return The multiplicative attenuation factor, in \[0, 1\]. It is 0 when either true variance is 0,
#'   in which case no correlation is identified and [compensation_lower_bound()] reports that rather
#'   than dividing by it.
attenuation_factor <- function(var_c, var_t, s2_c, s2_p) {
  a <- list(var_c = var_c, var_t = var_t, s2_c = s2_c, s2_p = s2_p)
  ## Variances cannot be negative. Without this the function returns NaN from sqrt() of a negative
  ## number, which propagates into compensation_lower_bound() as an implied correlation of NaN.
  for (nm in names(a)) {
    v <- a[[nm]]
    if (length(v) != 1L || !is.finite(v) || v < 0)
      stop(sprintf("%s must be a single finite non-negative number (a variance)", nm))
  }
  den <- sqrt(var_c + s2_c) * sqrt(var_t + s2_p + s2_c)
  if (den == 0) return(0)
  sqrt(var_c) * sqrt(var_t) / den
}

#' State a reported compensation as the lower bound it is
#'
#' @param r_observed the cross-block estimate, a single number in \[-1, 1\].
#' @param ... passed to [attenuation_factor()]; if the four variance components are supplied, the
#'   implied `rho` is also returned.
#' @return A list with `lower_bound` (|r_observed|) and, when the variance components are given,
#'   `attenuation`, `rho_implied` and `inputs_inconsistent`. Never returns a point estimate of `rho`
#'   without the components, because there is none. `rho_implied` is `NA` whenever it would lie outside
#'   \[-1, 1\] or the attenuation is zero, with `inputs_inconsistent = TRUE` and the unclamped value in
#'   `rho_implied_raw`.
compensation_lower_bound <- function(r_observed, ...) {
  if (length(r_observed) != 1L || !is.finite(r_observed) || abs(r_observed) > 1)
    stop("r_observed must be a single finite number in [-1, 1] (it is a correlation)")
  out <- list(lower_bound = abs(r_observed),
              statement = sprintf("|rho| >= %.4f; a point estimate of rho is not identified without the variance components",
                                  abs(r_observed)))
  a <- list(...)
  if (all(c("var_c", "var_t", "s2_c", "s2_p") %in% names(a))) {
    f <- attenuation_factor(a$var_c, a$var_t, a$s2_c, a$s2_p)
    out$attenuation <- f
    ## Two ways the division fails, and both must be reported rather than returned. A zero attenuation
    ## sends the implied correlation to infinity, and dividing by a small one can carry it past 1. Neither
    ## is a rounding artefact: both say the supplied components are inconsistent with the observed
    ## correlation, so one of the inputs is wrong.
    ri <- if (f == 0) Inf * sign(r_observed) else r_observed / f
    if (!is.finite(ri) || abs(ri) > 1) {
      warning(sprintf(paste("implied rho = %s lies outside [-1, 1] (attenuation %.4f), so the supplied",
                            "variance components are inconsistent with r_observed = %.4f. Reporting it as",
                            "not identified rather than as a correlation."), format(ri), f, r_observed),
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
