# Reference outputs of clustEff 0.3.1 (CRAN), used by the regression tests of the
# legacy options. Run with clustEff 0.3.1 installed.
# Output: tests/testthat/fixtures/ref_0.3.1.rds
library(clustEff)
stopifnot(packageVersion("clustEff") == "0.3.1")
ref <- list()

## waveforms (Sottile and Adelfio, 2019, Section 4.2), smaller grid
set.seed(1234)
tt <- seq(0, 1, length.out = 200)
sdt <- sqrt(pmax(4 * pmax(tt - 0.2, 0) - 8 * pmax(tt - 0.5, 0) + 4 * pmax(tt - 0.8, 0), 0)) * 0.3
W <- cbind(sapply(1:6, function(i) sin(3 * pi * tt) + rnorm(200, 0, sdt)),
           sapply(1:6, function(i) cos(3 * pi * tt) + rnorm(200, 0, sdt)),
           sapply(1:4, function(i) sin(3 * pi * tt) * cos(pi * tt) + rnorm(200, 0, sdt)))
colnames(W) <- paste0("w", 1:16)
D1 <- distshape(W, alpha = 0.5, step = "both", trace = FALSE)
D2 <- distshape(W, alpha = 0.25, step = "shape", trace = FALSE)
D3 <- distshape(W, alpha = 0.25, step = "distance", trace = FALSE)
o1 <- clustEff(W, k = c(2, 5), plot = FALSE, trace = FALSE)
o2 <- clustEff(W, k = 3, plot = FALSE, trace = FALSE, method = "average")
o3 <- clustEff(W, k = c(2, 5), cut.method = "length", plot = FALSE, trace = FALSE)
ref$wave <- list(W = W, D1 = as.matrix(D1), D2 = as.matrix(D2), D3 = as.matrix(D3),
                 o1 = list(clusters = o1$clusters, k = o1$k, distance = o1$distance, X.mean = o1$X.mean),
                 o2 = list(clusters = o2$clusters, X.mean = o2$X.mean, X.mean.lower = o2$X.mean.lower),
                 o3 = list(clusters = o3$clusters, k = o3$k))

## effects curves (Section 4.1), 3 groups x 4 responses
set.seed(1234)
n <- 300
x <- runif(n, 0, 5)
th <- list(c(.5, .5, 1, 2), c(-.3, .5, 1, .5), c(.3, .5, -1, 2))
Y <- NULL
for (g in 1:3) for (r in 1:4) {
  tg <- th[[g]] + rnorm(4, 0, 0.2)
  u <- runif(n)
  Y <- cbind(Y, (1 + qnorm(u)) + (tg[1] + tg[2] * u + tg[3] * u^2 + tg[4] * u^3) * x)
}
eo <- extract.object(Y, cbind(x = x), formula.p = ~ I(qnorm(p)) + p + I(p^2) + I(p^3))
B <- eo$X[[1]]; Bl <- eo$Xl[[1]]; Bu <- eo$Xr[[1]]
o4 <- clustEff(B, Bl, Bu, k = c(2, 5), cut.method = "conf.int", plot = FALSE, trace = FALSE)
ref$eff <- list(Y = Y, x = x, p = eo$p, B = B, Bl = Bl, Bu = Bu,
                o4 = list(clusters = o4$clusters, k = o4$k, distance = o4$distance,
                          X.mean.lower = o4$X.mean.lower, Signif.interval = o4$Signif.interval))
ref$info <- list(version = "0.3.1", R = R.version.string, date = Sys.Date())
saveRDS(ref, "tests/testthat/fixtures/ref_0.3.1.rds", version = 2)
