#' Is a count of significant features an ordering of effect sizes?
#'
#' It is only when the noise floors are shared. This returns the two orderings side by side and the
#' factor by which the count exaggerates each difference, so that the answer is a number rather than an
#' impression. The exaggeration factor is `(amplitude ratio) / (count ratio)` against a named anchor.
#'
#' @param condition character, condition labels.
#' @param n_sig integer, significant-feature counts.
#' @param amplitude numeric, a noise-corrected amplitude per condition (see [rms_true()]).
#' @param noise_floor numeric, the noise floor per condition.
#' @param anchor the condition to express ratios against.
#' @param max_floor_spread the largest ratio between the highest and lowest noise floor at which the
#'   floors may be treated as shared. The default of 1.1 is a convention, not a result: it says the
#'   floors must agree to within ten per cent. It is exposed because the right value depends on how much
#'   amplitude distortion the caller is willing to tolerate, and because a hard-coded threshold invites
#'   the reader to assume it was derived. For reference, the dataset in the accompanying paper has a
#'   spread of 1.47 and fails at any threshold below that.
#' @return A data frame with count and amplitude ratios and the exaggeration factor, plus attributes
#'   `floor_spread` and `floors_shared`. When `floors_shared` is `FALSE` the count must not be read as
#'   an amplitude ordering, and the exaggeration column says by how much it would mislead.
count_vs_amplitude <- function(condition, n_sig, amplitude, noise_floor, anchor = condition[1],
                               max_floor_spread = 1.1) {
  stopifnot(length(condition) == length(n_sig), length(n_sig) == length(amplitude),
            length(amplitude) == length(noise_floor), anchor %in% condition)
  i <- match(anchor, condition)
  cr <- n_sig / n_sig[i]
  ar <- amplitude / amplitude[i]
  out <- data.frame(condition = condition, n_sig = n_sig, amplitude = amplitude,
                    noise_floor = noise_floor, count_ratio = cr, amplitude_ratio = ar,
                    exaggeration = ifelse(cr > 0, ar / cr, NA_real_),
                    stringsAsFactors = FALSE)
  stopifnot(all(is.finite(noise_floor)), all(noise_floor > 0), max_floor_spread >= 1)
  sp <- max(noise_floor) / min(noise_floor)
  attr(out, "floor_spread") <- sp
  attr(out, "floors_shared") <- sp < max_floor_spread
  attr(out, "max_floor_spread") <- max_floor_spread
  attr(out, "anchor") <- anchor
  out
}

#' Non-additivity of a compositional share across two perturbations
#'
#' For a zero-sum interaction contrast `(+1, -1, -1, +1)` over wild type, two single perturbations and
#' the double, any share difference that is **additive** across the perturbations cancels exactly. What
#' survives is the non-additive part. This computes it, which is the quantity that predicts an offset.
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
#' returns each class median as a **difference from a named reference class**, which is the reading that
#' is stable across normalisation regimes.
#'
#' Two limitations are not removable and are returned in the `caveats` attribute. The prescription
#' recovers **differences, not absolute values**: if the reference class carries real signal it
#' propagates into every estimate. And no within-dataset analysis distinguishes "the reference class has
#' real biology" from "the reference class retains residual composition".
#'
#' @param value numeric, one value per feature.
#' @param class character or factor, class assignment; must be mutually exclusive.
#' @param reference the class to read against.
#' @param se_of_median optional numeric of the same length as the number of classes, in the order of
#'   `sort(unique(class))`, used to form z values.
#' @return A data frame with `n`, `median`, `vs_reference` and (when available) `z_vs_reference`.
reference_class_medians <- function(value, class, reference, se_of_median = NULL) {
  stopifnot(length(value) == length(class), reference %in% class)
  ## One pass with split() instead of two scans of the class vector per class: the original form is
  ## O(n * K) and this is O(n + K). It matters only for large K, but it also removes a subtler problem --
  ## `class == k` on a factor with unused levels silently yields empty groups, whereas split() over the
  ## character vector groups exactly the classes present.
  cl <- sort(unique(as.character(class)))
  grp <- split(value, factor(as.character(class), levels = cl))
  med <- vapply(grp, function(v) stats::median(v, na.rm = TRUE), 0)
  n   <- vapply(grp, function(v) sum(is.finite(v)), 0L)
  ref <- med[[reference]]
  out <- data.frame(class = cl, n = as.integer(n), median = as.numeric(med),
                    vs_reference = as.numeric(med) - ref, stringsAsFactors = FALSE)
  if (!is.null(se_of_median)) {
    stopifnot(length(se_of_median) == length(cl))
    out$z_vs_reference <- out$vs_reference / se_of_median
  }
  attr(out, "reference") <- reference
  attr(out, "reference_median") <- ref
  attr(out, "caveats") <- c(
    "recovers differences, not absolute values: a reference class carrying real signal propagates into every estimate",
    "no within-dataset analysis separates 'reference has real biology' from 'reference retains residual composition'")
  out
}
