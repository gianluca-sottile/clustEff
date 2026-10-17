# Effects curves: polynomial coefficient functions, G groups of m curves
sim_effects <- function(G = 3, m = 8, seed = 1, p = seq(0.01, 0.99, 0.01)) {
  set.seed(seed)
  th <- list(c(.5, .5, 1, 2), c(-3, .5, 1, .5), c(.3, -.5, -1, 2), c(2.5, .5, 1, 2))
  B <- NULL
  for (g in seq_len(G)) for (r in seq_len(m)) {
    t <- th[[g]] + rnorm(4, 0, 0.3)
    B <- cbind(B, t[1] + t[2] * p + t[3] * p^2 + t[4] * p^3)
  }
  colnames(B) <- paste0("g", rep(seq_len(G), each = m), "_", rep(seq_len(m), G))
  list(B = B, L = B - 0.6, U = B + 0.6, p = p, group = rep(seq_len(G), each = m))
}
# waveforms (Sottile and Adelfio, 2019, Section 4.2)
sim_waves <- function(m = c(6, 6, 4), seed = 1, n = 150, sd = 0.15) {
  set.seed(seed)
  tt <- seq(0, 1, length.out = n)
  f <- list(function(t) sin(3 * pi * t), function(t) cos(3 * pi * t), function(t) sin(3 * pi * t) * cos(pi * t))
  W <- do.call(cbind, lapply(seq_along(m), function(g) sapply(seq_len(m[g]), function(i) f[[g]](tt) + rnorm(n, 0, sd))))
  list(W = W, t = tt, group = rep(seq_along(m), m))
}
same_partition <- function(a, b) {
  tab <- table(a, b)
  all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)
}
# proportion of curves in the majority true group of their cluster
agreement <- function(a, b) sum(apply(table(a, b), 1, max)) / length(a)
