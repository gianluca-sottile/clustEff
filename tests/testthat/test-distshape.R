test_that("legacy dissimilarity reproduces clustEff 0.3.1", {
  ref <- readRDS(test_path("fixtures", "ref_0.3.1.rds"))$wave
  expect_equal(as.matrix(distshape(ref$W, alpha = 0.5, legacy = TRUE)), ref$D1, ignore_attr = TRUE)
  expect_equal(as.matrix(distshape(ref$W, alpha = 0.25, step = "shape", legacy = TRUE)), ref$D2,
               ignore_attr = TRUE)
  expect_equal(as.matrix(distshape(ref$W, alpha = 0.25, step = "distance", legacy = TRUE)), ref$D3,
               ignore_attr = TRUE)
})

test_that("the dissimilarity of eq. (1) has the expected properties", {
  s <- sim_effects(G = 2, m = 4)
  D <- distshape(s$B, p = s$p, alpha = 0.25)
  Dm <- as.matrix(D)
  expect_s3_class(D, "dist")
  expect_true(all(Dm >= 0 & Dm <= 1))
  expect_equal(Dm, t(Dm))
  expect_true(all(diag(Dm) == 0))
  expect_equal(attr(D, "Labels"), colnames(s$B))
  # product <= each component: d >= max(d_shape, d_distance)
  Ds <- as.matrix(distshape(s$B, p = s$p, alpha = 0.25, step = "shape"))
  Dd <- as.matrix(distshape(s$B, p = s$p, alpha = 0.25, step = "distance"))
  expect_true(all(Dm >= pmax(Ds, Dd) - 1e-12))
  # mean of the two components
  Dmean <- as.matrix(distshape(s$B, p = s$p, alpha = 0.25, combine = "mean"))
  expect_equal(Dmean, (Ds + Dd) / 2)
})

test_that("shape component: linear curves are concordant, exact derivatives are used", {
  p <- seq(0.01, 0.99, 0.01)
  B <- cbind(a = 1 + p, b = 2 - p, c = p^2, d = -p^2)
  Ds <- as.matrix(distshape(B, p = p, step = "shape"))
  expect_equal(Ds["a", "b"], 0)                # both linear
  expect_equal(Ds["c", "d"], 1)                # convex vs concave
  d2 <- cbind(0, 0, 2, -2)
  Ds2 <- as.matrix(distshape(B, p = p, step = "shape", d = matrix(d2, length(p), 4, byrow = TRUE)))
  expect_equal(Ds, Ds2, ignore_attr = TRUE)
  expect_error(distshape(B, p = p, d = matrix(0, 2, 2)), "same dimensions")
})

test_that("input checks", {
  expect_error(distshape(matrix(1:10, ncol = 1)), "two curves")
  expect_error(distshape(cbind(1:10, 2:11), alpha = 2), "alpha")
  expect_error(distshape(cbind(1:10, 2:11), p = 1:5), "one value for each row")
  expect_error(distshape(Beta = cbind(1:10, 2:11)), "removed in version 1.0.0")
})
