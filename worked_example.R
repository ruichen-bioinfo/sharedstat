suppressMessages(library(sharedstat)); set.seed(20261010)
G <- 5000; gamma <- 1.10
R <- rnorm(G, 0, 1)                      # true heat response of the reference allele
delta <- rnorm(G, 0, 0.3)                # gene departures (cis divergence not captured by gamma)
sdu <- 0.55                              # per-block measurement noise
Rh1 <- R + rnorm(G,0,sdu); Rh2 <- R + rnorm(G,0,sdu)               # reference arm, block 1 and 2
A1 <- gamma*R + delta + rnorm(G,0,sdu); A2 <- gamma*R + delta + rnorm(G,0,sdu)
Rh <- (Rh1+Rh2)/2; A <- (A1+A2)/2; e <- A - Rh
naive <- 1 + coef(lm(e ~ Rh))[2]
rb <- reliability_from_blocks(Rh1, Rh2); print(str(rb))
h <- if (is.list(rb)) rb$h else rb
cat(sprintf("true gamma %.3f | naive %.3f | h %.3f\n", gamma, naive, h))
print(naive_gain_bias(h = h, gamma = gamma))
cat(sprintf("naive/h %.3f\n", naive/h))
## split-block: axis from block 1, difference from block 2
e2 <- A2 - Rh2; split <- 1 + coef(lm(e2 ~ Rh1))[2]
h1 <- cov(Rh1, Rh2)/var(Rh1)   # reliability of block 1, estimated from the two blocks (no true R needed)
cat(sprintf("split-block gain %.3f (expected 1+(g-1)h1 = %.3f) ; magnitude 1+(split-1)/h1 = %.3f\n", split, 1+(gamma-1)*h1, 1+(split-1)/h1))
## instrument: the same gene's response measured independently elsewhere
S <- R + rnorm(G,0,0.6)
print(iv_gain(e, Rh, S))
cat("sharedstat version", as.character(packageVersion("sharedstat")), "\n")
