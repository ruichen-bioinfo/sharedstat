#' The five admission checks, as a runnable report
#'
#' Before reporting a statistic built from a shared component: (1) derive its value under no effect
#' algebraically and compare the observation to *that*, not to zero; (2) for each degenerate condition
#' state **which** condition it is, and verify by simulation that the label matches the algebra;
#' (3) if the statistic scales with a nuisance magnitude, report it across the full normalisation
#' exponent and mark the measured one; (4) never read a count as an amplitude ordering without showing
#' the noise floors are shared; (5) for compositional data, state the reference class and prefer the
#' difference to the absolute value.
#'
#' Item 2 exists because the authors of the accompanying paper failed it: the no-effect value had been
#' derived correctly but belonged to a different degenerate condition than the one it was labelled with,
#' and a checklist stopping at item 1 would not have caught that.
#'
#' @param no_effect_value the statistic's value under no effect, or `NA` if not derived.
#' @param degenerate_condition a character label naming which degenerate condition `no_effect_value`
#'   belongs to, and `simulation_confirms` whether a simulation verified the label.
#' @param simulation_confirms logical.
#' @param beta_range numeric of length 2, the normalisation exponent range reported, or `NULL`.
#' @param measured_beta the measured exponent, or `NA`.
#' @param floors_shared logical, from [count_vs_amplitude()], or `NA` if no count is reported.
#' @param reference_class the reference class used, or `NA` if the reading is against zero.
#' @return An object of class `sharedstat_checklist`; print it for the report.
admission_checklist <- function(no_effect_value = NA_real_,
                                degenerate_condition = NA_character_,
                                simulation_confirms = NA,
                                beta_range = NULL, measured_beta = NA_real_,
                                floors_shared = NA, reference_class = NA_character_) {
  items <- list(
    list(n = 1, text = "no-effect value derived algebraically and used as the comparison",
         pass = is.finite(no_effect_value),
         catches = "shared reference arm; residual anticorrelated with its own component"),
    list(n = 2, text = "which degenerate condition, stated and simulation-verified",
         pass = !is.na(degenerate_condition) && isTRUE(simulation_confirms),
         catches = "the authors' own inverted inference: item 1 passed and the conclusion was still wrong"),
    list(n = 3, text = "reported across the normalisation exponent with the measured one marked",
         pass = !is.null(beta_range) && length(beta_range) == 2 && is.finite(measured_beta),
         catches = "amplitude-carried enrichment"),
    list(n = 4, text = "count not read as an amplitude ordering unless the noise floors are shared",
         pass = is.na(floors_shared) || isTRUE(floors_shared),
         catches = "counts exaggerating between-condition differences"),
    list(n = 5, text = "compositional data read against a stated reference class",
         pass = !is.na(reference_class),
         catches = "a zero-sum contrast carrying a non-zero median")
  )
  structure(list(items = items,
                 n_pass = sum(vapply(items, function(x) isTRUE(x$pass), TRUE)),
                 n_total = length(items)),
            class = "sharedstat_checklist")
}

#' @export
print.sharedstat_checklist <- function(x, ...) {
  cat(sprintf("sharedstat admission checklist: %d of %d items satisfied\n\n", x$n_pass, x$n_total))
  for (it in x$items) {
    cat(sprintf(" [%s] %d. %s\n        catches: %s\n",
                if (isTRUE(it$pass)) "x" else " ", it$n, it$text, it$catches))
  }
  if (x$n_pass < x$n_total)
    cat("\nAn unsatisfied item is not a verdict on the analysis; it identifies a check that has not been\nshown to have been made, and the artefact it would have caught.\n")
  invisible(x)
}
