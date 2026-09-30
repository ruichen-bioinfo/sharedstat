#' Reliability of a reference arm, from two replicate blocks
#'
#' The parameter that controls the bias in a shared-reference design is the reliability of the reference
#' arm: the share of its between-feature variance that is signal rather than measurement error. It is
#' measurable inside the design, with nothing fitted, whenever the reference has been measured twice.
#'
#' Given the same contrast computed in two replicate blocks, `b1` and `b2`, form
#' `signal = (b1 + b2)/2` and `null = (b2 - b1)/2`. The null contains no signal by construction, so
#'
#' \deqn{h = 1 - \frac{\mathrm{mean}(null^2)}{\mathrm{mean}(signal^2)} = 1 - 1/snr^2,
#'       \qquad snr = \frac{\mathrm{rms}(signal)}{\mathrm{rms}(null)}.}
#'
#' `h` is a variance ratio taken **over the features entering this contrast**. It is not an intraclass
#' correlation between technical replicates and it is not a general reliability of the assay: `h = 0.05`
#' means five per cent of the between-feature variance of this reference contrast is signal, not that a
#' sequencing measurement is five per cent reproducible.
#'
#' @param b1,b2 numeric vectors of the same length: one contrast, computed in each of two replicate
#'   blocks, over the same features in the same order. Non-finite entries are dropped in pairs, because
#'   the two mean squares must be taken over the same features.
#' @param min_n minimum number of complete pairs required.
#' @param center Whether to subtract each contrast's mean over features before the second moment is taken,
#'   making `h` a ratio of variances. `TRUE` by default, because the gain estimators are `1 + cov(e, S) /
#'   var(S)` and the attenuation of a least-squares slope fitted with an intercept is a ratio of centred
#'   moments; a replicate that is uniformly stronger than another shifts the intercept, not the slope, and
#'   must not be charged to the noise. `center = FALSE` gives the uncentred mean-square form. The two agree
#'   when both contrasts have zero mean over features, and can differ materially when they do not: in one
#'   arm of the Keyport Kik data the replicates differ in mean response by 37% of the response amplitude
#'   and the two forms give 0.916 (centred) and 0.887 (uncentred).
#' @return A list with `h`, `snr`, `rms_signal`, `rms_null`, `n_pairs` and `n_dropped`. `h` is clamped at
#'   zero and reported as such: a value of exactly zero means the paired contrast carried no more variance
#'   than its own null, so the reliability is not distinguishable from nothing and the naive estimator is
#'   uninformative rather than slightly biased. `h` is `NA` when the null has zero variance, since the
#'   ratio is then undefined rather than perfect. `h` therefore lies in `[0, 1]`, the same range as
#'   [reliability_from_replicates()], so either function's output can be passed to [naive_gain_bias()] without
#'   the caller checking which one produced it.
#' @seealso [naive_gain_bias()], which takes this `h` and returns what a naive estimator will report.
#' @examples
#' set.seed(1)
#' truth <- rnorm(2000)
#' b1 <- truth + rnorm(2000, sd = 0.6)
#' b2 <- truth + rnorm(2000, sd = 0.6)
#' r <- reliability_from_blocks(b1, b2)
#' naive_gain_bias(h = r$h, gamma = 1.1)
reliability_from_blocks <- function(b1, b2, min_n = 3L, center = TRUE) {
  if (length(min_n) != 1L || !is.finite(min_n) || min_n < 3)
    stop("min_n must be a single finite number of at least 3")
  b1 <- as.numeric(b1); b2 <- as.numeric(b2)
  if (length(b1) != length(b2))
    stop(sprintf("b1 and b2 must be the same length (got %d and %d): the two mean squares are taken over the same features",
                 length(b1), length(b2)))
  keep <- is.finite(b1) & is.finite(b2)
  if (sum(keep) < min_n)
    stop(sprintf("only %d complete pairs; need at least %d", sum(keep), min_n))
  if (length(center) != 1L || !is.logical(center) || is.na(center))
    stop("center must be TRUE or FALSE")
  x <- b1[keep]; y <- b2[keep]
  sig <- (x + y) / 2
  nul <- (y - x) / 2
  ## `center = TRUE` subtracts each contrast's mean over features before the second moment is taken, which makes h
  ## a ratio of VARIANCES. That is the quantity the gain estimators need: they are `1 + cov(e, S) / var(S)`, and
  ## the attenuation of a least-squares slope fitted with an intercept is a ratio of centred moments. A replicate
  ## whose response is uniformly stronger than the other's shifts the intercept, not the slope, and must not be
  ## charged to the noise. `center = FALSE` gives the uncentred mean-square form. The two agree when both
  ## contrasts have zero mean over features and diverge as either mean grows: in one arm of the Keyport Kik data
  ## the replicates differ in mean response by 37% of the response amplitude, and the forms give 0.916 (centred) and 0.887 (uncentred).
  if (center) { sig <- sig - mean(sig); nul <- nul - mean(nul) }
  ms_sig <- mean(sig^2); ms_nul <- mean(nul^2)
  rms_sig <- sqrt(ms_sig); rms_nul <- sqrt(ms_nul)
  ## A null with no variance would give an infinite signal-to-noise ratio and h = 1. That is not a
  ## perfectly reliable reference; it is two blocks that are numerically identical, which in practice
  ## means the same numbers were supplied twice. Report it as undefined rather than as perfect.
  ## 0.1.3: relative tolerance, as in reliability_from_replicates(). With center = TRUE two blocks that differ only
  ## by a constant leave a null of ~1e-33 instead of 0, and an exact test returned h = 1 with no warning.
  if (ms_nul <= ms_sig * .Machine$double.eps * 8 || ms_nul == 0) {
    warning("the null contrast has zero variance: b1 and b2 are identical on every retained feature, ",
            "so the reliability is not identified. Check that two distinct replicate blocks were supplied.",
            call. = FALSE)
    h <- NA_real_; snr <- NA_real_
  } else {
    snr <- rms_sig / rms_nul
    h <- max(0, 1 - ms_nul / ms_sig)
  }
  list(h = h, snr = snr, rms_signal = rms_sig, rms_null = rms_nul,
       centred = center, n_pairs = sum(keep), n_dropped = sum(!keep))
}
