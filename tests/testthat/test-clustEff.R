test_that("legacy options reproduce clustEff 0.3.1", {
  ref <- readRDS(test_path("fixtures", "ref_0.3.1.rds"))
  w <- ref$wave
  o1 <- clustEff(w$W, k = c(2, 5), legacy = TRUE)
  expect_equal(unname(o1$clusters), unname(w$o1$clusters))
  expect_equal(unname(o1$criterion), unname(w$o1$distance))
  expect_equal(o1$mean, w$o1$X.mean, ignore_attr = TRUE)
  o2 <- clustEff(w$W, k = 3, method = "average", legacy = TRUE)
  expect_equal(unname(o2$clusters), unname(w$o2$clusters))
  o3 <- clustEff(w$W, k = c(2, 5), criterion = "height", legacy = TRUE)
  expect_equal(o3$k, w$o3$k)
  e <- ref$eff
  o4 <- clustEff(e$B, e$Bl, e$Bu, k = c(2, 5), criterion = "conf.int", legacy = TRUE)
  expect_equal(unname(o4$clusters), unname(e$o4$clusters))
  expect_equal(o4$mean.lower, e$o4$X.mean.lower, ignore_attr = TRUE)
})

test_that("effects curves: the true groups are recovered", {
  s <- sim_effects(G = 3)
  fit <- clustEff(s$B, s$L, s$U, p = s$p, k = c(2, 6))
  expect_s3_class(fit, "clustEff")
  expect_equal(fit$k, 3)
  expect_gte(agreement(fit$clusters, s$group), 0.9)
  expect_equal(fit$settings$alpha, 0.25)
  expect_equal(dim(fit$mean), c(length(s$p), 3))
  expect_true(is.logical(fit$signif))
  expect_true(all(diff(as.integer(table(fit$clusters))) >= 0))   # numbered by size
  for (cr in c("pi.out", "dist.rel", "silhouette")) {
    f <- clustEff(s$B, s$L, s$U, p = s$p, k = c(2, 6), criterion = cr)
    expect_true(f$k %in% 2:6)
    expect_length(f$criterion, if (cr == "silhouette") 5 else 7)
  }
})

test_that("waveforms with smoothing", {
  s <- sim_waves()
  fit <- clustEff(s$W, p = s$t, k = c(2, 5), smooth = TRUE)
  expect_equal(fit$k, 3)
  expect_true(same_partition(fit$clusters, s$group))
  expect_equal(fit$settings$alpha, 0.5)
  fit2 <- clustEff(s$W, p = s$t, k = 3, smooth = TRUE, stand = TRUE)
  expect_null(fit2$criterion)
})

test_that("a user-supplied dissimilarity", {
  s <- sim_effects(G = 2, m = 4)
  D <- distshape(s$B, p = s$p, alpha = 0.25)
  f1 <- clustEff(s$B, p = s$p, k = 2, diss = D)
  f2 <- clustEff(s$B, p = s$p, k = 2, alpha = 0.25)
  expect_equal(f1$clusters, f2$clusters)
  expect_error(clustEff(s$B, k = 2, diss = dist(1:3)), "between the 'curves'")
})

test_that("methods", {
  s <- sim_effects(G = 3)
  fit <- clustEff(s$B, s$L, s$U, p = s$p, k = c(2, 6))
  expect_output(print(fit), "3 clusters")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.clustEff")
  expect_output(print(sm), "Within-cluster measures")
  expect_equal(nrow(sm$measures), 4)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (ty in c("clusters", "dendrogram", "silhouette", "dissimilarity", "criterion")) expect_silent(plot(fit, ty))
  expect_silent(plot(fit, which = 2))
  expect_error(plot(fit, which = 5), "which")
})

test_that("removed arguments give informative errors", {
  s <- sim_effects(G = 2, m = 4)
  expect_error(clustEff(Beta = s$B), "removed in version 1.0.0")
  expect_error(clustEff(s$B, cut.method = "mindist"), "criterion")
  expect_error(clustEff(s$B, plot = FALSE), "removed")
  expect_error(clustEff(s$B, foo = 1), "unused argument")
  expect_error(extract.object(), "removed in clustEff 1.0.0")
  expect_error(clustEff(s$B, k = c(1, 3)), "start from 2")
  expect_error(clustEff(s$B, criterion = "pi.out", k = c(2, 3)), "requires 'lower'")
})

test_that("extract_curves() and clustEff() on qrcm models", {
  skip_if_not_installed("qrcm")
  set.seed(1)
  n <- 300
  x1 <- runif(n); x2 <- runif(n); x3 <- runif(n)
  u <- runif(n)
  y <- qnorm(u) + (1 + u) * x1 + (1 + u) * x2 - (2 * u^2) * x3
  m <- qrcm::iqr(y ~ x1 + x2 + x3, formula.p = ~ I(p) + I(p^2))
  ec <- extract_curves(m)
  expect_s3_class(ec, "effects_curves")
  expect_equal(colnames(ec$curves), c("x1", "x2", "x3"))
  expect_true(all(ec$lower <= ec$curves & ec$curves <= ec$upper))
  # second derivatives of a quadratic basis: constant in p
  expect_lt(max(apply(ec$d2, 2, function(v) diff(range(v)))), 1e-4)
  expect_output(print(ec), "Effects curves")
  fit <- clustEff(m, k = 2)
  expect_equal(fit$k, 2)
  m2 <- qrcm::iqr(y ~ x1, formula.p = ~ I(p) + I(p^2))
  ec2 <- extract_curves(list(a = m, b = m2), which = "x1")
  expect_equal(colnames(ec2$curves), c("a:x1", "b:x1"))
  expect_error(extract_curves(m, which = "zz"), "does not match")
})

test_that("edge cases: few curves, k = q, criteria at the boundary", {
  s <- sim_effects(G = 2, m = 3)
  f <- clustEff(s$B[, 1:4], p = s$p)
  expect_equal(f$settings$k, c(2, 3))
  f2 <- clustEff(s$B, p = s$p, k = ncol(s$B))
  expect_equal(f2$k, 6)
  expect_null(f2$silhouette)
  f3 <- clustEff(s$B, p = s$p, k = c(2, 5), criterion = "silhouette")
  expect_true(f3$k %in% 2:5)
  f4 <- clustEff(s$B[, 1:2], p = s$p)
  expect_equal(f4$k, 2)
  expect_error(clustEff(s$B, s$L, s$U, c(2, 5)), "not 'k'")
  expect_error(fpcac(s$B, 2), "not 'K'")
})
