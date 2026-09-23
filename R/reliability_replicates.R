#' Reliability of a reference arm from any number of replicates
#'
#' `reliability_from_blocks()` takes exactly two replicate blocks. With more than two, the same decomposition
#' still applies: one contrast carries the signal and the remaining orthogonal contrasts carry only error. This
#' function forms that decomposition for any number of replicates \eqn{k \ge 2} and returns the same reliability
#' `h`.
#'
#' The signal contrast is the mean over replicates. The \eqn{k-1} null contrasts are the Helmert contrasts,
#' normalised so that every contrast has the same variance \eqn{\sigma^2/k} as the mean does:
#' \deqn{n_j \propto \left(\underbrace{1,\dots,1}_{j}, -j, 0, \dots\right)}
#' Each null has expectation zero because its coefficients sum to zero, and the contrast vectors are mutually
#' orthogonal, so
#' \deqn{h = 1 - \frac{\overline{\mathrm{mean}(n_j^2)}}{\mathrm{mean}(\mathrm{sig}^2)}}
#' estimates the same quantity as the two-block form and reduces to it exactly when \eqn{k = 2}.
#'
#' The null mean square averages the per-contrast mean squares rather than pooling the values, so each null
#' contributes equally regardless of how many features it happens to have.
#'
#' @param x A numeric matrix or data frame with one column per replicate and one row per feature, holding the
#'   same quantity measured in each replicate. At least two columns.
#' @param min_n Minimum number of complete rows required. Default 3.
#' @param center Whether to subtract each contrast's mean over features before the second moment is taken,
#'   making `h` a ratio of variances. `TRUE` by default, because the gain estimators are `1 + cov(e, S) /
#'   var(S)` and the attenuation of a least-squares slope fitted with an intercept is a ratio of centred
#'   moments; a replicate that is uniformly stronger than another shifts the intercept, not the slope, and
#'   must not be charged to the noise. `center = FALSE` gives the uncentred mean-square form. The two agree
#'   when both contrasts have zero mean over features, and can differ materially when they do not: in one
#'   arm of the Keyport Kik data the replicates differ in mean response by 37% of the response amplitude
#'   and the two forms give 0.887 and 0.916.
#'
#' @return A list with the same fields as [reliability_from_blocks()], plus `k`, the number of replicates, and
#'   `h_pairwise`, the reliability from each pair of replicates taken alone. Spread among `h_pairwise` shows how
#'   much of a two-block reliability is a property of which two blocks were used.
#'
#'   `h` and every element of `h_pairwise` lie in `[0, 1]`, the same range as [reliability_from_blocks()]. A
#'   plug-in value below zero means the unreproducible variation is at least as large as the variation across
#'   features; it is reported as `0`, which says that no reproducible signal was resolved above the noise of
#'   that contrast, and not that reliability is negative or that biological signal is absent. How far the
#'   estimate fell below the floor remains readable from `snr`, `rms_signal` and `rms_null`, and the
#'   inconsistency warning reports the raw pairwise ratios before flooring. `h = NA` is reserved for the
#'   undefined case of identical replicate columns, where the null contrasts have no variance at all.
#'
#' @details
#' Why this exists: the bias in a shared-reference design is O(1) in features and O(1/k) in replicates, so
#' raising `k` is the only way to raise `h`. A design with two blocks cannot show that, and comparing `h` at
#' k = 2 against k = 3 on the same data does. `h_pairwise` is returned for exactly that comparison.
#'
#' @examples
#' set.seed(1)
#' R <- rnorm(2000)
#' x <- cbind(R + rnorm(2000, sd = 0.5), R + rnorm(2000, sd = 0.5), R + rnorm(2000, sd = 0.5))
#' r <- reliability_from_replicates(x)
#' r$h
#' r$h_pairwise
#' # with two columns it is the two-block function exactly
#' all.equal(reliability_from_replicates(x[, 1:2])$h, reliability_from_blocks(x[, 1], x[, 2])$h)
#' @seealso [reliability_from_blocks()] for the two-block case
#' @export
reliability_from_replicates <- function(x, min_n = 3L, center = TRUE) {
  if (length(min_n) != 1L || !is.finite(min_n) || min_n < 3)
    stop("min_n must be a single finite number of at least 3")
  if (length(center) != 1L || !is.logical(center) || is.na(center))
    stop("center must be TRUE or FALSE")
  x <- as.matrix(x)
  if (!is.numeric(x)) storage.mode(x) <- "double"
  k <- ncol(x)
  if (is.na(k) || k < 2L)
    stop(sprintf("x must have at least two columns, one per replicate (got %s)",
                 if (is.na(k)) "none" else k))
  keep <- stats::complete.cases(x) & apply(is.finite(x), 1L, all)
  if (sum(keep) < min_n)
    stop(sprintf("only %d complete rows; need at least %d", sum(keep), min_n))
  x <- x[keep, , drop = FALSE]

  sig <- rowMeans(x)
  ## Helmert contrasts, each scaled to variance sigma^2 / k so that the null mean squares are directly
  ## comparable with the signal mean square. For j = 1 .. k-1 the coefficients are j ones, then -j, then zeros;
  ## dividing by sqrt(j * (j + 1) * k / ...) is done via the norm so no algebra is duplicated here.
  nulls <- vector("list", k - 1L)
  for (j in seq_len(k - 1L)) {
    co <- c(rep(1, j), -j, rep(0, k - j - 1L))
    co <- co / sqrt(sum(co^2)) / sqrt(k)   # variance sigma^2 / k, matching var(sig)
    nulls[[j]] <- as.vector(x %*% co)
  }
  ## `center = TRUE` removes each contrast's mean over features, making h a ratio of variances -- the quantity the
  ## gain estimators need, since a least-squares slope fitted with an intercept is attenuated by a ratio of centred
  ## moments. Here it also removes a per-replicate constant offset: the Helmert contrasts turn a replicate that is
  ## uniformly stronger than the others into a null with a non-zero mean, which the uncentred form charges to the
  ## noise. On GSE310700 the largest per-replicate offset reached 12% of the response amplitude and h moved by
  ## 0.002; the effect is small there because h compares mean squares, but it is not zero and it is not signed in
  ## a predictable direction. `center = FALSE` gives the uncentred mean-square form.
  if (center) {
    sig <- sig - mean(sig)
    nulls <- lapply(nulls, function(v) v - mean(v))
  }
  ms_sig <- mean(sig^2)
  ms_nul <- mean(vapply(nulls, function(v) mean(v^2), 0))   # average of per-contrast mean squares
  ## The two-block function tests `ms_nul == 0` for identical inputs, and exact equality is right there because
  ## the null is a subtraction of two supplied numbers. Here the nulls come out of a matrix product, so
  ## identical columns give 3e-15 rather than 0 and an exact test lets h = 1 through as "perfect reliability" --
  ## precisely the outcome the two-block function documents as impossible. The comparison is therefore relative
  ## to the signal scale, which is the only scale on which "no null variance" is meaningful.
  if (ms_nul <= ms_sig * .Machine$double.eps * 8 || ms_nul == 0) {
    warning("every null contrast has zero variance: the replicate columns are identical, so the reliability ",
            "is not identified. Check that distinct replicates were supplied.", call. = FALSE)
    return(list(h = NA_real_, snr = NA_real_, rms_signal = sqrt(ms_sig), rms_null = 0,
                n_pairs = nrow(x), k = k, h_pairwise = rep(NA_real_, choose(k, 2))))
  }
  ## The pairwise values must use the same second moment as the k-replicate value, or the spread among them is not
  ## comparable with `h` and the inconsistency warning below would be computed on a different quantity.
  hp <- if (k == 2L) NULL else {
    cb <- utils::combn(k, 2L)
    vapply(seq_len(ncol(cb)), function(i) {
      a <- x[, cb[1L, i]]; b <- x[, cb[2L, i]]
      s <- (a + b) / 2; n <- (b - a) / 2
      if (center) { s <- s - mean(s); n <- n - mean(n) }
      1 - mean(n^2) / mean(s^2)
    }, 0)
  }
  ## The raw pairwise ratios are kept for the inconsistency warning below, which is computed on their spread. That
  ## spread is the severity signal, and clamping before measuring it would suppress the warning in exactly the cases
  ## it exists for.
  hp_raw <- hp
  ## Unequal replicate noise. The estimator assumes the replicates are exchangeable. When one is markedly noisier
  ## than the others the k-replicate value uses it anyway, and can fall below the best pair -- which is not a
  ## defect but is easy to misread as one, so it is announced. The threshold is a spread of 0.2 among the pairwise
  ## values; on the control arm of GSE310700 the spread is 0.345 and this fires.
  if (!is.null(hp_raw) && length(hp_raw) > 1L && diff(range(hp_raw)) > 0.2)
    warning(sprintf(paste0("one replicate is inconsistent with the others; the k-replicate value uses all of ",
                           "them. Raw pairwise ratios span %.3f to %.3f (spread %.3f); the reported pairwise ",
                           "reliabilities are these values floored at zero."),
                    min(hp_raw), max(hp_raw), diff(range(hp_raw))), call. = FALSE)
  ## Reliability is a share of variance and cannot be negative. A negative plug-in value means the unreproducible
  ## variation is at least as large as the variation across features, which is reported as zero; snr, rms_signal and
  ## rms_null are returned unchanged so how far below the floor the estimate fell is still visible. This makes the
  ## contract identical to reliability_from_blocks(), so both public entry points can be passed to
  ## naive_gain_bias() without the caller checking which one produced the number.
  h_main <- max(0, 1 - ms_nul / ms_sig)
  out <- list(h = h_main, snr = sqrt(ms_sig / ms_nul),
              rms_signal = sqrt(ms_sig), rms_null = sqrt(ms_nul),
              n_pairs = nrow(x), k = k,
              h_pairwise = if (is.null(hp)) h_main else pmax(0, hp))
  class(out) <- c("sharedstat_reliability", "list")
  out
}
