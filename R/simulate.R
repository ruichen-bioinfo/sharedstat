#' Simulate compositional counts in which the composition is not imposed
#'
#' Compositionality emerges here rather than being added: counts are drawn from a closed-sum multinomial,
#' so share competition is a consequence of the sampling. A separate high-expression "dominator" set,
#' belonging to none of the tested classes, carries the share swing, so the tested classes experience
#' dilution only. That separation is what makes the simulation a known-truth one: if the swing were
#' applied to the responding classes' own expression it would change their true interaction at the same
#' time, and every class's apparent bias would equal the injected amount by construction.
#'
#' @param n_genes total genes; at least one per tested class plus the dominators.
#' @param class_frac named numeric, the non-negative fraction of genes in each tested class; the
#'   remainder, less the dominators, is the reference class. Each fraction must be large enough to give
#'   at least one gene, since a class requested but not realised would leave `true_I` describing genes
#'   that do not exist.
#' @param true_I named numeric, the true interaction per class (same names, same order, as `class_frac`).
#' @param reference_true_I the true interaction of the reference class; set non-zero to simulate a
#'   contaminated reference.
#' @param kappa the share swing applied to the dominator set in the double genotype; negative values give
#'   negative share non-additivity, the sign observed in the dataset analysed in the accompanying paper.
#' @param n_dominators,dominator_boost size and log-expression boost of the dominator set. The dominators
#'   occupy the last `n_dominators` genes and must not overlap any tested class.
#' @param n_rep replicates per genotype, at least one.
#' @param lib library size drawn per sample.
#' @param seed random seed.
#' @return A list with the count matrix, the design, the class assignment, the true interaction per gene,
#'   the responding-class share per genotype and its non-additivity, suitable for pushing through any
#'   normalisation pipeline.
simulate_closed_sum <- function(n_genes = 4000L,
                                class_frac = c(Hsf1_only = 0.015, Msn24_only = 0.045,
                                               both_required = 0.0075, redundant = 0.02,
                                               TF_independent = 0.112),
                                true_I = c(Hsf1_only = -0.30, Msn24_only = -1.20,
                                           both_required = -0.50, redundant = -1.50,
                                           TF_independent = -0.60),
                                reference_true_I = 0,
                                kappa = 0,
                                n_dominators = 150L, dominator_boost = 3,
                                n_rep = 3L, lib = 2e7, seed = 1L) {
  ## ---- input validation ----------------------------------------------------------------------------
  ## Previously a single stopifnot() covered the names and the sum of the fractions, which let a negative
  ## fraction, a zero replicate count or a dominator set larger than the genome through to the generator.
  if (!identical(names(class_frac), names(true_I)))
    stop("class_frac and true_I must have the same names in the same order")
  if (length(n_genes) != 1L || !is.finite(n_genes) || n_genes < 1)
    stop("n_genes must be a single finite number of at least 1")
  if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1)
    stop("n_rep must be a single finite number of at least 1")
  if (length(n_dominators) != 1L || !is.finite(n_dominators) || n_dominators < 0 || n_dominators >= n_genes)
    stop("n_dominators must be a single finite number in [0, n_genes)")
  if (any(!is.finite(class_frac)) || any(class_frac < 0))
    stop("class_frac must be finite and non-negative")
  if (any(!is.finite(true_I)))
    stop("true_I must be finite")
  if (sum(class_frac) >= 1)
    stop("sum(class_frac) must be below 1, since the reference class is the remainder")
  n_genes <- as.integer(n_genes); n_rep <- as.integer(n_rep); n_dominators <- as.integer(n_dominators)

  ## The realised size of each class is a rounding of frac * n_genes. A fraction small enough to round to
  ## zero used to produce `idx:(idx + m - 1L)` with m = 0, which in R is a DESCENDING pair of indices
  ## rather than an empty selection: the class silently vanished and two unrelated positions were
  ## overwritten. Both the empty class and the overlap below are therefore rejected before anything is
  ## drawn.
  sizes <- vapply(class_frac, function(x) as.integer(round(x * n_genes)), 0L)
  if (any(sizes < 1L))
    stop(sprintf("class_frac gives no genes to: %s (with n_genes = %d). Raise the fraction or n_genes.",
                 paste(names(sizes)[sizes < 1L], collapse = ", "), n_genes))
  if (sum(sizes) + n_dominators > n_genes)
    stop(sprintf(paste("the tested classes take %d genes and the dominators %d, which exceeds n_genes = %d.",
                       "The dominators must lie outside every tested class, since a dominator inside one",
                       "would carry the share swing into a class whose truth is being measured."),
                 sum(sizes), n_dominators, n_genes))

  set.seed(seed)
  genotypes <- c("WT", "S1", "S2", "D")
  base_lfc  <- stats::rnorm(n_genes, 8, 1.5)
  cls <- rep("reference", n_genes)
  idx <- 1L
  for (k in names(class_frac)) {
    m <- sizes[[k]]
    end <- idx + m - 1L
    cls[idx:end] <- k
    idx <- end + 1L
  }
  dom <- if (n_dominators > 0L) (n_genes - n_dominators + 1L):n_genes else integer(0)
  if (length(dom) && any(cls[dom] != "reference"))
    stop("dominators overlap a tested class; this should have been caught by the size check above")
  cls[dom] <- "dominator"
  base_lfc[dom] <- base_lfc[dom] + dominator_boost
  tI <- ifelse(cls == "reference" | cls == "dominator", reference_true_I, true_I[cls])
  tI[cls == "dominator"] <- 0
  ## The four-genotype contrast is D - S1 - S2 + WT, whose coefficients sum to zero. Putting the whole
  ## interaction on D and leaving the other three at baseline makes that contrast equal tI exactly, so the
  ## simulator's truth needs no separate bookkeeping.
  lam <- matrix(NA_real_, n_genes, length(genotypes), dimnames = list(NULL, genotypes))
  lam[, "WT"] <- base_lfc
  lam[, "S1"] <- base_lfc
  lam[, "S2"] <- base_lfc
  lam[, "D"]  <- base_lfc + tI
  stopifnot(all(abs((lam[, "D"] - lam[, "S1"] - lam[, "S2"] + lam[, "WT"]) - tI) < 1e-12))
  lam[dom, "D"] <- lam[dom, "D"] + kappa                # the share swing, carried by dominators only
  cnt <- matrix(0L, n_genes, length(genotypes) * n_rep)
  cn  <- character(0)
  for (g in genotypes) for (r in seq_len(n_rep)) {
    p <- 2^lam[, g]; p <- p / sum(p)
    draw <- stats::rmultinom(1, lib, p)
    ## as.integer() returns NA above 2^31 - 1 with only a warning, which would put silent holes in the
    ## count matrix. lib is a total across all genes so no single cell should approach the limit, but the
    ## assertion costs nothing and turns a silent corruption into a stop.
    if (any(draw > .Machine$integer.max)) stop("a simulated count exceeds integer range; lower lib")
    cnt[, length(cn) + 1L] <- as.integer(draw)
    cn <- c(cn, paste0(g, "_r", r))
  }
  colnames(cnt) <- cn
  resp <- cls %in% names(class_frac)
  sh <- vapply(genotypes, function(g) {
    cols <- grep(paste0("^", g, "_"), cn); sum(cnt[resp, cols]) / sum(cnt[, cols])
  }, 0)
  list(counts = cnt,
       design = data.frame(sample = cn, genotype = sub("_r[0-9]+$", "", cn), stringsAsFactors = FALSE),
       class = cls, true_I = tI,
       class_sizes = c(sizes, reference = n_genes - sum(sizes) - n_dominators, dominator = n_dominators),
       share = sh,
       share_nonadditivity = share_nonadditivity(sh[["WT"]], sh[["S1"]], sh[["S2"]], sh[["D"]]))
}
