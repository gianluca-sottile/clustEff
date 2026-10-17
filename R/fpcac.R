#' Functional Principal Components Analysis of Curves on a B-spline Basis
#'
#' Represents each curve on a B-spline basis by least squares and computes the
#' functional principal components (harmonics) and the scores, as in
#' Ramsay and Silverman (2005, Chapter 8). It is used by [fpcac()].
#'
#' @param curves a numeric matrix with one column for each curve and one row for
#'   each time point.
#' @param t the time points (default: equally spaced on \[0, 1\]).
#' @param nbasis number of B-spline basis functions.
#' @param norder order of the B-splines (degree + 1; default 3, quadratic).
#' @param nharmonics number of principal components.
#' @param lambda ridge penalty on the basis coefficients (default 1e-8, for
#'   numerical stability).
#' @return A list with items `scores` (curves x harmonics), `values` (the
#'   eigenvalues), `varprop` (proportion of variance explained), `coefs` (basis
#'   coefficients of the curves), `harmonics` (basis coefficients of the
#'   harmonics), `mean` (basis coefficients of the mean curve), `fitted` (the
#'   smoothed curves on `t`), `basis` (the basis matrix on `t`) and `gram` (the
#'   Gram matrix of the basis).
#' @references Ramsay, J. O. and Silverman, B. W. (2005). *Functional Data
#'   Analysis*, 2nd ed. Springer.
#' @examples
#' tt <- seq(0, 1, length.out = 100)
#' X <- sapply(1:20, function(i) sin(2 * pi * tt + runif(1)) + rnorm(100, 0, 0.1))
#' f <- fpca_bspline(X, tt, nbasis = 10, nharmonics = 2)
#' f$varprop
#' @export
fpca_bspline <- function(curves, t = NULL, nbasis = 5, norder = 3, nharmonics = 3, lambda = 1e-8) {
  X <- as.matrix(curves)
  n <- nrow(X); N <- ncol(X)
  if (is.null(t)) t <- seq(0, 1, length.out = n)
  if (length(t) != n) stop("'t' must have one value for each row of 'curves'", call. = FALSE)
  if (nbasis < norder) stop("'nbasis' must be at least 'norder'", call. = FALSE)
  if (nharmonics < 1 || nharmonics != round(nharmonics) || nharmonics > min(nbasis, N - 1))
    stop("'nharmonics' must be an integer between 1 and min(nbasis, number of curves - 1)", call. = FALSE)
  rng <- range(t)
  breaks <- seq(rng[1], rng[2], length.out = nbasis - norder + 2)
  knots <- c(rep(rng[1], norder - 1), breaks, rep(rng[2], norder - 1))
  Phi <- splines::splineDesign(knots, t, ord = norder, outer.ok = TRUE)
  # least squares coefficients (nbasis x N)
  PtP <- crossprod(Phi)
  C <- solve(PtP + diag(lambda * max(1, mean(diag(PtP))), nbasis), crossprod(Phi, X))
  # Gram matrix: exact Gauss-Legendre quadrature on each knot interval
  gl <- gauss_legendre(norder + 1)
  W <- matrix(0, nbasis, nbasis)
  for (s in seq_len(length(breaks) - 1)) {
    a <- breaks[s]; b <- breaks[s + 1]
    xs <- (b - a) / 2 * gl$x + (a + b) / 2
    Ps <- splines::splineDesign(knots, xs, ord = norder, outer.ok = TRUE)
    W <- W + crossprod(Ps * ((b - a) / 2 * gl$w), Ps)
  }
  # functional PCA (as fda::pca.fd): eigen-decomposition of R C C' R' / N, with W = R'R
  mC <- rowMeans(C)
  Cc <- C - mC
  R <- chol(W)
  V <- R %*% tcrossprod(Cc) %*% t(R) / N
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  harm <- backsolve(R, e$vectors[, seq_len(nharmonics), drop = FALSE])
  scores <- crossprod(Cc, W %*% harm)
  colnames(scores) <- paste0("PC", seq_len(nharmonics))
  rownames(scores) <- colnames(X)
  values <- pmax(e$values, 0)
  list(scores = scores, values = values[seq_len(nharmonics)],
       varprop = values[seq_len(nharmonics)] / sum(values), coefs = C, harmonics = harm, mean = mC,
       fitted = Phi %*% C, basis = Phi, gram = W, t = t)
}

#' Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
#' @noRd
gauss_legendre <- function(m) {
  if (m == 1) return(list(x = 0, w = 2))
  i <- seq_len(m - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m)
  J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

#' Trimmed k-means (García-Escudero and Gordaliza) on the rows of a matrix.
#' @noRd
trimmed_kmeans <- function(X, K, alpha = 0, nstart = 30, iter.max = 25, init = "kmeans++", verbose = FALSE) {
  n <- nrow(X)
  ntrim <- floor(n * (1 - alpha))
  if (K > ntrim) stop("too many clusters for the number of (untrimmed) curves", call. = FALSE)
  x2 <- rowSums(X^2)
  d2 <- function(C) pmax(outer(x2, rowSums(C^2), "+") - 2 * tcrossprod(X, C), 0)
  seed_centers <- function() {
    if (init == "random") return(X[sample.int(n, K), , drop = FALSE])
    C <- X[sample.int(n, 1), , drop = FALSE]
    while (nrow(C) < K) {
      dmin <- apply(d2(C), 1, min)
      pr <- if (sum(dmin) > 0) dmin / sum(dmin) else NULL
      C <- rbind(C, X[sample.int(n, 1, prob = pr), , drop = FALSE])
    }
    C
  }
  best <- list(obj = Inf)
  for (s in seq_len(nstart)) {
    C <- seed_centers()
    prev <- NULL
    for (it in seq_len(iter.max)) {
      D <- d2(C)
      ind <- max.col(-D, ties.method = "first")
      dist <- D[cbind(seq_len(n), ind)]
      keep <- rank(dist, ties.method = "first") <= ntrim
      if (!is.null(prev) && identical(prev, ind * keep)) break
      prev <- ind * keep
      for (j in seq_len(K)) {
        mem <- keep & ind == j
        if (any(mem)) C[j, ] <- colMeans(X[mem, , drop = FALSE])
        else C[j, ] <- X[which.max(ifelse(keep, dist, -Inf)), ]  # empty cluster: farthest kept point
      }
    }
    D <- d2(C)
    ind <- max.col(-D, ties.method = "first")
    dist <- D[cbind(seq_len(n), ind)]
    obj <- mean(sort(dist)[seq_len(ntrim)])
    if (verbose) message(sprintf("start %3d: objective = %.6g%s", s, obj, if (obj < best$obj) "  <--" else ""))
    if (obj < best$obj) best <- list(obj = obj, centers = C)
  }
  D <- d2(best$centers)
  ind <- max.col(-D, ties.method = "first")
  dist <- D[cbind(seq_len(n), ind)]
  radius <- sort(dist)[ntrim]
  asig <- ifelse(rank(dist, ties.method = "first") <= ntrim, ind, 0L)
  list(clusters = asig, centers = best$centers, obj = best$obj, radius = radius)
}

#' Functional PCA-based Clustering of Curves (FPCAC)
#'
#' Clusters curves (e.g., seismic waveforms) by a trimmed k-means applied to
#' their functional principal component scores (Adelfio et al., 2011).
#'
#' @param curves a numeric matrix with one column for each curve and one row for
#'   each time point.
#' @param t the time points (default: equally spaced on \[0, 1\]).
#' @param K the number of clusters.
#' @param nbasis,norder number and order of the B-spline basis functions.
#' @param nharmonics number of functional principal components.
#' @param alpha trimming proportion: the `N - floor(N * (1 - alpha))` curves
#'   farthest from their centers are not assigned to any cluster (label 0).
#' @param nstart number of random starts.
#' @param iter.max maximum number of concentration steps for each start.
#' @param init initialization of the centers: `"kmeans++"` (default) or `"random"`
#'   (as in 'clustEff' 0.3.x).
#' @param conf.level level of the pointwise bands of the curves of each cluster.
#' @param verbose logical: print the objective function of each start?
#' @param ... arguments of 'clustEff' 0.3.x that were removed (informative error).
#'
#' @details
#' Each curve is represented on a B-spline basis and summarized by the scores of
#' its first `nharmonics` functional principal components (see
#' [fpca_bspline()]). The trimmed k-means of García-Escudero and Gordaliza finds
#' `K` centers minimizing the mean squared distance of the closest
#' \eqn{\lfloor N(1-\alpha)\rfloor} curves to their nearest center. Use
#' [set.seed()] for reproducible results.
#'
#' @return An object of class `"fpcac"`, a list with items `clusters` (0 =
#'   trimmed), `K`, `centers`, `obj.function`, `radius`, `scores`, `fpca` (the
#'   output of [fpca_bspline()]), `curves`, `t`, `mean`, `mean.lower`,
#'   `mean.upper` (mean curve and pointwise bands of each cluster),
#'   `mean.dist` (correlation distance of each curve from the mean of its
#'   cluster) and `call`.
#' @references Adelfio, G., Chiodi, M., D'Alessandro, A. and Luzio, D. (2011).
#'   FPCA algorithm for waveform clustering. *Journal of Communication and
#'   Computer*, 8, 494-502.
#'
#'   García-Escudero, L. A. and Gordaliza, A. (2005). A proposal for robust curve
#'   clustering. *Journal of Classification*, 22, 185-201.
#'   \doi{10.1007/s00357-005-0013-x}
#' @seealso [opt.fpcac()], [clustEff()]
#' @examples
#' set.seed(1)
#' tt <- seq(0, 1, length.out = 100)
#' X <- cbind(sapply(1:10, function(i) sin(3 * pi * tt) + rnorm(100, 0, 0.2)),
#'            sapply(1:10, function(i) cos(3 * pi * tt) + rnorm(100, 0, 0.2)),
#'            sapply(1:2, function(i) rnorm(100, 0, 0.2)))          # two outliers
#' fit <- fpcac(X, t = tt, K = 2, nbasis = 10, alpha = 0.1)
#' fit
#' plot(fit)
#' @export
fpcac <- function(curves, t = NULL, K = 2, nbasis = 5, norder = 3, nharmonics = 3, alpha = 0,
                  nstart = 30, iter.max = 25, init = c("kmeans++", "random"), conf.level = 0.9,
                  verbose = FALSE, ...) {
  check_removed(list(...), "fpcac")
  init <- match.arg(init)
  if (alpha < 0 || alpha >= 1) stop("'alpha' must be in [0, 1)", call. = FALSE)
  X <- as.matrix(curves)
  if (!is.null(t) && length(t) == 1 && nrow(X) > 1)
    stop("'t' must be the time points of the curves; in clustEff 1.0.0 the second argument is 't', ",
         "not 'K': please use K = ", t, call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (is.null(t)) t <- seq(0, 1, length.out = nrow(X))
  fp <- fpca_bspline(X, t, nbasis = nbasis, norder = norder, nharmonics = nharmonics)
  tk <- trimmed_kmeans(fp$scores, K, alpha = alpha, nstart = nstart, iter.max = iter.max, init = init,
                       verbose = verbose)
  # clusters numbered by increasing size (as in clustEff 0.3.x); 0 = trimmed
  asig <- tk$clusters
  tab <- table(factor(asig[asig > 0], levels = seq_len(K)))
  ord <- order(tab)
  asig <- ifelse(asig > 0, match(asig, ord), 0L)
  centers <- tk$centers[ord, , drop = FALSE]
  names(asig) <- colnames(X)
  a <- (1 - conf.level) / 2
  M <- ML <- MU <- matrix(NA_real_, nrow(X), K, dimnames = list(NULL, paste0("cluster", seq_len(K))))
  mdist <- vector("list", K)
  for (j in seq_len(K)) {
    Xj <- X[, asig == j, drop = FALSE]
    if (!ncol(Xj)) next
    M[, j] <- rowMeans(Xj)
    ML[, j] <- apply(Xj, 1, stats::quantile, probs = a, names = FALSE)
    MU[, j] <- apply(Xj, 1, stats::quantile, probs = 1 - a, names = FALSE)
    mdist[[j]] <- sqrt(pmax(2 * (1 - suppressWarnings(c(stats::cor(Xj, M[, j])))), 0))
  }
  out <- list(call = match.call(), clusters = asig, K = K, centers = centers, obj.function = tk$obj,
              radius = tk$radius, scores = fp$scores, fpca = fp, curves = X, t = t,
              mean = M, mean.lower = ML, mean.upper = MU, mean.dist = mdist, conf.level = conf.level,
              alpha = alpha)
  class(out) <- "fpcac"
  out
}

#' @rdname fpcac
#' @param x an object of class `"fpcac"`.
#' @export
print.fpcac <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("FPCAC: ", x$K, " clusters of sizes ", paste(table(factor(x$clusters[x$clusters > 0], levels = seq_len(x$K))),
                                                  collapse = ", "),
      "; trimmed curves: ", sum(x$clusters == 0), "\n", sep = "")
  cat("n. of harmonics: ", ncol(x$scores), " (", format(100 * sum(x$fpca$varprop), digits = 3),
      "% of the variance)\n\n", sep = "")
  invisible(x)
}

#' Summarize an FPCAC Fit
#'
#' @param object an object of class `"fpcac"`.
#' @param x an object of class `"summary.fpcac"`.
#' @param digits number of significant digits.
#' @param ... not used.
#' @return An object of class `"summary.fpcac"` with the clustering table, the
#'   number of trimmed curves, the proportion of variance of each harmonic and,
#'   for each cluster, the average correlation distance of its curves from the
#'   mean curve.
#' @seealso [fpcac()]
#' @export
summary.fpcac <- function(object, ...) {
  tab <- table(factor(object$clusters[object$clusters > 0], levels = seq_len(object$K)))
  avd <- vapply(object$mean.dist, function(v) if (length(v)) mean(v, na.rm = TRUE) else NA_real_, numeric(1))
  avd <- c(avd, Average = sum(avd * tab, na.rm = TRUE) / sum(tab))
  names(avd)[seq_len(object$K)] <- paste0("cluster", seq_len(object$K))
  out <- list(call = object$call, K = object$K, N = ncol(object$curves), table = tab,
              trimmed = sum(object$clusters == 0), varprop = object$fpca$varprop, mean.dist = avd)
  class(out) <- "summary.fpcac"
  out
}

#' @rdname summary.fpcac
#' @export
print.summary.fpcac <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of curves:", x$N, "   n. of clusters:", x$K, "   trimmed curves:", x$trimmed, "\n\n")
  cat("Clustering table:\n"); print(x$table)
  cat("\nProportion of variance of the harmonics:\n"); print(signif(x$varprop, digits))
  cat("\nAverage correlation distance from the cluster mean:\n"); print(signif(x$mean.dist, digits))
  cat("\n")
  invisible(x)
}

#' Plot an FPCAC Fit
#'
#' @param x an object of class `"fpcac"`.
#' @param which the clusters to plot. With a single cluster, its curves are drawn.
#' @param polygon logical: draw the bands as shaded areas?
#' @param col colours of the clusters.
#' @param ... further graphical parameters (`xlab`, `ylab`, `main`, `ylim`, `lwd`).
#' @return No return value, called for its side effects.
#' @seealso [fpcac()]
#' @export
plot.fpcac <- function(x, which, polygon = TRUE, col, ...) {
  L <- list(...)
  K <- x$K
  if (missing(which)) which <- seq_len(K)
  if (any(which < 1 | which > K)) stop(sprintf("'which' must be in 1-%d", K), call. = FALSE)
  if (missing(col)) col <- grDevices::hcl.colors(max(K, 2), "Dark 3")[seq_len(K)]
  col <- rep_len(col, K)
  arg <- function(name, default) if (is.null(L[[name]])) default else L[[name]]
  tt <- x$t
  one <- length(which) == 1
  yl <- arg("ylim", range(if (one) x$curves[, x$clusters == which] else NULL,
                          x$mean.lower[, which], x$mean.upper[, which], na.rm = TRUE))
  plot(range(tt), yl, type = "n", xlab = arg("xlab", "t"), ylab = arg("ylab", "x(t)"),
       main = arg("main", if (one) paste("Cluster", which) else ""))
  for (j in which) {
    if (polygon) graphics::polygon(c(tt, rev(tt)), c(x$mean.lower[, j], rev(x$mean.upper[, j])),
                                   col = grDevices::adjustcolor(col[j], alpha.f = 0.25), border = NA)
    else { lines(tt, x$mean.lower[, j], lty = 2, col = col[j]); lines(tt, x$mean.upper[, j], lty = 2, col = col[j]) }
    if (one) matlines(tt, x$curves[, x$clusters == j, drop = FALSE], lty = 3, col = col[j])
    lines(tt, x$mean[, j], lwd = arg("lwd", 2), col = col[j])
  }
  if (!one) graphics::legend("topleft", legend = paste("cluster", which), col = col[which], lwd = 2, bty = "n")
  invisible(NULL)
}

#' Number of Clusters for FPCAC
#'
#' Fits [fpcac()] for \eqn{K = 1, \dots, k.max} and computes the average
#' silhouette width or the total within-cluster sum of squares.
#'
#' @param curves,t,nbasis,norder,nharmonics,alpha,nstart,iter.max,init arguments of [fpcac()].
#' @param k.max maximum number of clusters.
#' @param criterion `"silhouette"` (default: \eqn{K} maximizing the average
#'   silhouette width, \eqn{K \ge 2}) or `"wss"` (\eqn{K} with the largest
#'   decrease of the within-cluster sum of squares from \eqn{K - 1}, as in
#'   'clustEff' 0.3.x).
#' @param diss an optional dissimilarity between the curves (class `"dist"`). By
#'   default the Euclidean distance between the vectors of scores, consistent with
#'   the trimmed k-means (`'clustEff'` 0.3.x used the correlation distance
#'   \eqn{\sqrt{2(1 - r)}} between the vectors of scores, which is degenerate with
#'   one or two harmonics).
#' @param ... arguments of 'clustEff' 0.3.x that were removed (informative error).
#' @return An object of class `"opt.fpcac"`, a list with items `criterion` (the
#'   values for each \eqn{K}), `K.opt`, `clusters` (matrix of the clusters for
#'   each \eqn{K}), `fits` (the [fpcac()] fits) and `criterion.name`.
#'   Trimmed curves are excluded from the computation of the criteria.
#' @seealso [fpcac()]
#' @examples
#' \donttest{
#' set.seed(1)
#' tt <- seq(0, 1, length.out = 100)
#' X <- cbind(sapply(1:8, function(i) sin(3 * pi * tt) + rnorm(100, 0, 0.2)),
#'            sapply(1:8, function(i) cos(3 * pi * tt) + rnorm(100, 0, 0.2)),
#'            sapply(1:8, function(i) sin(3 * pi * tt) * cos(pi * tt) + rnorm(100, 0, 0.2)))
#' o <- opt.fpcac(X, t = tt, k.max = 5, nbasis = 10)
#' o$K.opt
#' plot(o)
#' }
#' @export
opt.fpcac <- function(curves, t = NULL, k.max = 5, criterion = c("silhouette", "wss"), nbasis = 5,
                      norder = 3, nharmonics = 3, alpha = 0, nstart = 30, iter.max = 25,
                      init = c("kmeans++", "random"), diss = NULL, ...) {
  check_removed(list(...), "opt.fpcac")
  criterion <- match.arg(criterion)
  init <- match.arg(init)
  fits <- lapply(seq_len(k.max), function(K)
    fpcac(curves, t = t, K = K, nbasis = nbasis, norder = norder, nharmonics = nharmonics, alpha = alpha,
          nstart = nstart, iter.max = iter.max, init = init))
  clusters <- vapply(fits, function(f) f$clusters, numeric(ncol(as.matrix(curves))))
  colnames(clusters) <- seq_len(k.max)
  if (is.null(diss)) diss <- stats::dist(fits[[1]]$scores)
  Dm <- as.matrix(diss)
  v <- vapply(seq_len(k.max), function(K) {
    oc <- clusters[, K]; keep <- oc > 0
    if (criterion == "silhouette") {
      s <- safe_silhouette(oc[keep], stats::as.dist(Dm[keep, keep]))
      if (is.null(s)) NA_real_ else mean(s[, 3])
    } else {
      sum(vapply(unique(oc[keep]), function(j) {
        idx <- which(oc == j)
        sum(Dm[idx, idx]^2) / 2 / length(idx)
      }, numeric(1)))
    }
  }, numeric(1))
  names(v) <- seq_len(k.max)
  K.opt <- if (criterion == "silhouette") { if (all(is.na(v))) 1L else as.integer(names(which.max(v))) }
           else if (k.max >= 2) as.integer(names(which.min(diff(v)))) else 1L
  out <- list(criterion = v, K.opt = K.opt, clusters = clusters, fits = fits, criterion.name = criterion)
  class(out) <- "opt.fpcac"
  out
}

#' @rdname opt.fpcac
#' @param x an object of class `"opt.fpcac"`.
#' @export
print.opt.fpcac <- function(x, ...) {
  cat("\nNumber of clusters for FPCAC (criterion: ", x$criterion.name, ")\n\n", sep = "")
  print(signif(x$criterion, 4))
  cat("\nSelected number of clusters:", x$K.opt, "\n\n")
  invisible(x)
}

#' @rdname opt.fpcac
#' @export
plot.opt.fpcac <- function(x, ...) {
  L <- list(...)
  ylab <- if (x$criterion.name == "silhouette") "Average silhouette width" else "Total within sum of squares"
  ks <- as.integer(names(x$criterion))
  plot(ks, x$criterion, type = "b", pch = 19, col = "steelblue",
       xlab = if (is.null(L$xlab)) "Number of clusters K" else L$xlab,
       ylab = if (is.null(L$ylab)) ylab else L$ylab,
       main = if (is.null(L$main)) "Optimal number of clusters" else L$main)
  abline(v = x$K.opt, lty = 2, col = "steelblue")
  invisible(NULL)
}
