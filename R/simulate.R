#' Simulate compositional counts in which the composition is not imposed
#'
#' The point of the simulator is that compositionality **emerges**: counts are drawn from a closed-sum
#' multinomial, so share competition is a consequence of the sampling rather than an offset added by
#' hand. A separate high-expression "dominator" set, belonging to none of the tested classes, carries the
#' share swing, so the tested classes experience dilution only. That separation matters: an earlier
#' version of this design added the swing to the responding classes' own expression, which changed their
#' true interaction at the same time and made every class's apparent bias exactly equal to the injected
#' amount -- a self-fulfilling result.
#'
#' @param n_genes total genes.
#' @param class_frac named numeric, the fraction of genes in each tested class; the remainder is the
#'   reference class.
#' @param true_I named numeric, the true interaction per class (same names as `class_frac`).
#' @param reference_true_I the true interaction of the reference class; set non-zero to simulate a
#'   contaminated reference.
#' @param kappa the share swing applied to the dominator set in the double genotype; negative values
#'   give negative share non-additivity, which is the sign real data may have.
#' @param n_dominators,dominator_boost size and log-expression boost of the dominator set.
#' @param n_rep replicates per genotype.
#' @param seed random seed.
#' @return A list with the count matrix, the design, the realised share non-additivity, and the true
#'   interaction per gene, suitable for pushing through any normalisation pipeline.
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
                                n_rep = 3L, seed = 1L) {
  stopifnot(identical(names(class_frac), names(true_I)), sum(class_frac) < 1)
  set.seed(seed)
  genotypes <- c("WT", "S1", "S2", "D")
  base_lfc  <- stats::rnorm(n_genes, 8, 1.5)
  cls <- rep("reference", n_genes)
  idx <- 1L
  for (k in names(class_frac)) {
    m <- as.integer(round(class_frac[[k]] * n_genes))
    cls[idx:(idx + m - 1L)] <- k; idx <- idx + m
  }
  dom <- (n_genes - n_dominators + 1L):n_genes          # dominators sit outside every tested class
  cls[dom] <- "dominator"
  base_lfc[dom] <- base_lfc[dom] + dominator_boost
  tI <- ifelse(cls == "reference" | cls == "dominator", reference_true_I, true_I[cls])
  tI[cls == "dominator"] <- 0
  ## The four-genotype contrast is D - S1 - S2 + WT, whose coefficients sum to zero. Putting the whole
  ## interaction on D and leaving the other three at baseline makes that contrast equal tI exactly, so the
  ## simulator's "truth" needs no separate bookkeeping. (An earlier version built lam in a loop and then
  ## overwrote all four columns with these four lines, leaving the loop and its coefficient vector as dead
  ## code that suggested a more elaborate construction than the one actually used.)
  lam <- matrix(NA_real_, n_genes, length(genotypes), dimnames = list(NULL, genotypes))
  lam[, "WT"] <- base_lfc
  lam[, "S1"] <- base_lfc
  lam[, "S2"] <- base_lfc
  lam[, "D"]  <- base_lfc + tI
  stopifnot(all(abs((lam[, "D"] - lam[, "S1"] - lam[, "S2"] + lam[, "WT"]) - tI) < 1e-12))
  lam[dom, "D"] <- lam[dom, "D"] + kappa                # the share swing, carried by dominators only
  cnt <- matrix(0L, n_genes, length(genotypes) * n_rep)
  cn  <- character(0)
  lib <- 2e7
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
       share = sh,
       share_nonadditivity = share_nonadditivity(sh[["WT"]], sh[["S1"]], sh[["S2"]], sh[["D"]]))
}
