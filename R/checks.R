#' Is a count of significant features an ordering of effect sizes?
#'
#' It is only when the detection thresholds are comparable. This returns the two orderings side by side
#' and the factor by which the count exaggerates each difference, so that the answer is a number rather
#' than an impression. The exaggeration factor is `(amplitude ratio) / (count ratio)` against a named
#' anchor.
#'
#' @param condition character, condition labels.
#' @param n_sig significant-feature counts; non-negative, and strictly positive at the anchor, since the
#'   anchor's count is a denominator.
#' @param amplitude numeric, a noise-corrected amplitude per condition (see [rms_true()]); strictly
#'   positive at the anchor for the same reason.
#' @param noise_floor numeric, the noise floor per condition; strictly positive.
#' @param anchor the condition to express ratios against.
#' @param max_floor_spread the largest ratio between the highest and lowest noise floor at which the
#'   floors may be treated as shared. The comparison is inclusive, so a spread exactly equal to this
#'   value counts as shared. The default of 1.1 is a convention, not a result: it says the floors must
#'   agree to within ten per cent. It is exposed because the right value depends on how much amplitude
#'   distortion the caller is willing to tolerate, and because a hard-coded threshold invites the reader
#'   to assume it was derived.
#' @return A data frame with count and amplitude ratios and the exaggeration factor, plus attributes
#'   `floor_spread`, `floors_shared`, `max_floor_spread` and `anchor`. When `floors_shared` is `FALSE`
#'   the count should not be read as an amplitude ordering, and the exaggeration column says by how much
#'   it would mislead.
count_vs_amplitude <- function(condition, n_sig, amplitude, noise_floor, anchor = condition[1],
                               max_floor_spread = 1.1) {
  stopifnot(length(condition) == length(n_sig), length(n_sig) == length(amplitude),
            length(amplitude) == length(noise_floor))
  if (anyDuplicated(condition))
    stop("condition labels must be unique: the anchor is matched by label")
  if (!anchor %in% condition)
    stop(sprintf("anchor '%s' is not one of the conditions", as.character(anchor)))
  if (length(max_floor_spread) != 1L || !is.finite(max_floor_spread) || max_floor_spread < 1)
    stop("max_floor_spread must be a single finite number of at least 1")
  n_sig <- as.numeric(n_sig); amplitude <- as.numeric(amplitude); noise_floor <- as.numeric(noise_floor)
  if (any(!is.finite(n_sig)) || any(n_sig < 0))
    stop("n_sig must be finite and non-negative")
  if (any(!is.finite(amplitude)) || any(amplitude < 0))
    stop("amplitude must be finite and non-negative")
  if (any(!is.finite(noise_floor)) || any(noise_floor <= 0))
    stop("noise_floor must be finite and strictly positive")
  i <- match(anchor, condition)
  ## Both ratios are taken against the anchor, so a zero anchor makes every ratio Inf or NaN. Left
  ## unchecked the exaggeration column came back as a mixture of NA and 0 that looked like a result.
  if (n_sig[i] <= 0)
    stop(sprintf("the anchor '%s' has n_sig = %g; the anchor's count is the denominator of every count ratio",
                 as.character(anchor), n_sig[i]))
  if (amplitude[i] <= 0)
    stop(sprintf("the anchor '%s' has amplitude = %g; the anchor's amplitude is the denominator of every amplitude ratio",
                 as.character(anchor), amplitude[i]))
  cr <- n_sig / n_sig[i]
  ar <- amplitude / amplitude[i]
  out <- data.frame(condition = condition, n_sig = n_sig, amplitude = amplitude,
                    noise_floor = noise_floor, count_ratio = cr, amplitude_ratio = ar,
                    exaggeration = ifelse(cr > 0, ar / cr, NA_real_),
                    stringsAsFactors = FALSE)
  sp <- max(noise_floor) / min(noise_floor)
  attr(out, "floor_spread") <- sp
  ## Inclusive, to match the documented meaning of max_floor_spread as the largest spread AT WHICH the
  ## floors may be treated as shared. A strict comparison called a spread of exactly 1.1 unshared at a
  ## threshold of 1.1.
  attr(out, "floors_shared") <- sp <= max_floor_spread
  attr(out, "max_floor_spread") <- max_floor_spread
  attr(out, "anchor") <- anchor
  out
}

#' Non-additivity of a compositional share across two perturbations
#'
#' For a zero-sum interaction contrast `(+1, -1, -1, +1)` over wild type, two single perturbations and
#' the double, any share difference that is additive across the perturbations cancels exactly. What
#' survives is the non-additive part, which is the quantity that predicts an offset.
#'
#' @param wt,single1,single2,double numeric shares (or any additive quantity) in the four genotypes.
#' @return The contrast `double - single1 - single2 + wt`. Zero means an additive share difference,
#'   which cannot produce an offset however large the individual differences are.
share_nonadditivity <- function(wt, single1, single2, double) {
  double - single1 - single2 + wt
}

#' Class medians read against a reference class rather than against zero
#'
#' On compositional data an interaction contrast can carry a non-zero genome-wide median even though its
#' coefficients sum to zero. Reading a class median "against zero" is then a hidden assumption. This
#' returns each class median as a difference from a named reference class, which is the reading that is
#' stable across normalisation regimes.
#'
#' Two limitations are not removable and are returned in the `caveats` attribute. The prescription
#' recovers differences, not absolute values: if the reference class carries real signal it propagates
#' into every estimate. And no within-dataset analysis distinguishes "the reference class has real
#' biology" from "the reference class retains residual composition".
#'
#' @param value numeric, one value per feature. Non-finite values are dropped.
#' @param class character or factor, class assignment; must be mutually exclusive and free of missing
#'   values, since a missing class would be dropped silently and change every denominator.
#' @param reference the class to read against. It must have at least one finite value: without a
#'   reference median the estimand does not exist, so this is an error rather than a row of `NA`.
#' @param se_of_median optional numeric, used to form z values. If named, it is matched to the classes by
#'   name; if unnamed it is taken in the order of `sort(unique(class))`. Naming it is recommended, since
#'   an unnamed vector in the caller's own class order is silently wrong.
#' @return A data frame with `n`, `median`, `vs_reference` and (when available) `z_vs_reference`.
reference_class_medians <- function(value, class, reference, se_of_median = NULL) {
  stopifnot(length(value) == length(class))
  if (anyNA(class))
    stop("class contains missing values; a missing class would be dropped silently and change every count")
  cl_chr <- as.character(class)
  if (!reference %in% cl_chr)
    stop(sprintf("reference class '%s' does not appear in class", as.character(reference)))
  ## One pass with split() instead of two scans of the class vector per class: the original form is
  ## O(n * K) and this is O(n + K). It also removes a subtler problem -- `class == k` on a factor with
  ## unused levels silently yields empty groups, whereas split() over the character vector groups exactly
  ## the classes present.
  cl <- sort(unique(cl_chr))
  grp <- split(as.numeric(value), factor(cl_chr, levels = cl))
  ## n and the median must use ONE definition of an observation. Counting finite values while taking the
  ## median with na.rm = TRUE only removes NA and NaN, leaving Inf in the median: a class of two infinities
  ## then reported n = 0 alongside an infinite median.
  med <- vapply(grp, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) }, 0)
  n   <- vapply(grp, function(v) sum(is.finite(v)), 0L)
  if (n[[reference]] == 0L)
    stop(sprintf("reference class '%s' has no finite values, so there is no reference median to read against",
                 as.character(reference)))
  ref <- med[[reference]]
  out <- data.frame(class = cl, n = as.integer(n), median = as.numeric(med),
                    vs_reference = as.numeric(med) - ref, stringsAsFactors = FALSE)
  if (!is.null(se_of_median)) {
    if (!is.null(names(se_of_median))) {
      miss <- setdiff(cl, names(se_of_median))
      if (length(miss))
        stop(sprintf("se_of_median is named but has no entry for: %s", paste(miss, collapse = ", ")))
      se <- as.numeric(se_of_median[cl])
    } else {
      if (length(se_of_median) != length(cl))
        stop(sprintf("se_of_median has %d entries for %d classes; supply one per class, or name it",
                     length(se_of_median), length(cl)))
      se <- as.numeric(se_of_median)
    }
    out$z_vs_reference <- out$vs_reference / se
  }
  attr(out, "reference") <- reference
  attr(out, "reference_median") <- ref
  attr(out, "caveats") <- c(
    "recovers differences, not absolute values: a reference class carrying real signal propagates into every estimate",
    "no within-dataset analysis separates 'reference has real biology' from 'reference retains residual composition'")
  out
}
