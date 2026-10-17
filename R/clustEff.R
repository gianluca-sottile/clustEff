#' Clusters of Effects Curves
#'
#' Hierarchical clustering of curves, such as the effects curves
#' \eqn{\beta_j(p)} of quantile regression coefficients models or curves
#' observed in time, based on the shape-distance dissimilarity of Sottile and
#' Adelfio (2019).
#'
#' @param curves a numeric matrix with one column for each curve and one row
#'   for each value of `p`; or an object returned by [extract_curves()]; or an
#'   object of class `"iqr"` (see [qrcm::iqr()]) or a list of such objects, whose
#'   effects curves are clustered.
#' @param lower,upper optional matrices (same dimensions as `curves`) with the
#'   lower and upper pointwise confidence bands of the curves.
#' @param p the grid on which the curves are evaluated (quantile orders, or
#'   times). By default `(1:n)/n`, or the grid of an `"effects_curves"` object.
#' @param k the number of clusters: a single value, or a range `c(kmin, kmax)`
#'   in which the number of clusters is selected with `criterion`. The default
#'   is `c(2, min(5, q - 1))` for \eqn{q} curves.
#' @param alpha the quantile of the pointwise distances used by [distshape()].
#'   The default is 0.25 when bands are available (effects curves) and 0.5
#'   otherwise (waveforms), as suggested by Sottile and Adelfio (2019).
#' @param step,deriv,combine,shape.tol arguments of [distshape()].
#' @param criterion the criterion used to select the number of clusters when `k`
#'   is a range, see Details (default `"height"`).
#' @param method the agglomeration method of [stats::hclust()] (default
#'   `"ward.D2"`).
#' @param diss an optional dissimilarity (an object of class `"dist"`, e.g.
#'   computed by [distshape()]); if supplied, `alpha`, `step`, `deriv`, `combine`
#'   and `shape.tol` are ignored.
#' @param smooth logical: if `TRUE`, each curve (and its bands) is smoothed with
#'   [stats::smooth.spline()] before computing the dissimilarity (useful for
#'   noisy curves, as the shape component depends on derivatives).
#' @param df degrees of freedom of the smoothing splines (default: chosen by
#'   generalized cross-validation).
#' @param conf.level when the bands are not available, the bands of each cluster
#'   are the pointwise quantiles of order `(1 - conf.level)/2` and
#'   `(1 + conf.level)/2` of its curves.
#' @param stand logical: standardize each curve (center and scale over the grid)?
#' @param legacy logical: reproduce 'clustEff' 0.3.x (dissimilarity with
#'   `legacy = TRUE` in [distshape()], `alpha = 0.5`, criterion `"silhouette"`
#'   unless specified). Remaining differences with 0.3.x: `"height"` (formerly
#'   `cut.method = "length"`) searches only within the range `k`; with
#'   `stand = TRUE` also the bands are standardized; without bands the cluster
#'   bands use `conf.level` (0.3.x always used the quantiles 0.1 and 0.9).
#' @param verbose logical: print the progress?
#' @param ... arguments of 'clustEff' 0.3.x that were removed (they give an
#'   informative error, see [clustEff-defunct]).
#'
#' @details
#' The algorithm (Sottile and Adelfio, 2019, Section 3.2) computes the
#' dissimilarity matrix with [distshape()], builds a dendrogram with
#' [stats::hclust()], selects the number of clusters \eqn{l^*} and computes the
#' mean curve (and mean bands) of each cluster.
#'
#' Criteria for the number of clusters \eqn{l}:
#' * `"height"` (default): \eqn{l^*} corresponds to the largest gap between
#'   consecutive heights of the dendrogram.
#' * `"pi.out"` (effects curves with bands): with \eqn{\bar{LB}_j(p)} and
#'   \eqn{\bar{UB}_j(p)} the mean bands of cluster \eqn{j},
#'   \deqn{\pi_{out}(l) = l^{-1}\sum_{j=1}^l q_j^{-1}\sum_{i \in j}\int
#'   [1 - I\{\bar{LB}_j(p) \le \beta_i(p) \le \bar{UB}_j(p)\}]\,dp}{%
#'   pi.out(l) = mean over clusters of the average proportion of p where a curve is outside the mean bands}
#'   is the average proportion of the curves lying outside the mean bands of their
#'   cluster (Sottile and Adelfio, 2019, Section 3.1). It decreases with \eqn{l};
#'   \eqn{l^*} is the elbow, i.e. the value at which the decrease slows down most.
#' * `"dist.rel"` (waveforms): \eqn{dist_{rel}(l) = l\,\max_j q_j^{-1}\sum_{i\in j}
#'   \int|\bar\beta_j(p) - \beta_i(p)|\,dp} (Section 3.1), and \eqn{l^*} is the first
#'   number of clusters at which \eqn{dist_{rel}} stops decreasing.
#' * `"silhouette"`: \eqn{l^*} maximizes the average silhouette width
#'   (default criterion of 'clustEff' 0.3.x).
#' * `"conf.int"`: criterion `"conf.int"` of 'clustEff' 0.3.x (coverage of the
#'   mean bands weighted by the silhouette widths).
#'
#' The paper selects \eqn{l^*} minimizing \eqn{\pi_{out}(l) - \pi_{out}(l+1)} or
#' \eqn{dist_{rel}(l) - dist_{rel}(l+1)}, and writes \eqn{\pi_{out}} with the
#' indicator of the curves inside the bands and a factor \eqn{l}. In simulations
#' with 2 to 4 groups of effects curves these literal rules rarely found the true
#' number of clusters, while the rules above did better and `"height"` was the
#' most accurate; see the vignette.
#'
#' @return An object of class `"clustEff"`, a list with items:
#' \item{clusters}{the cluster of each curve (clusters are numbered by
#'   increasing size).}
#' \item{k}{the number of clusters.}
#' \item{p, curves, lower, upper}{the input (possibly smoothed or standardized).}
#' \item{mean, mean.lower, mean.upper}{matrices of the mean curves and mean bands
#'   of each cluster.}
#' \item{signif}{logical matrix: are the mean bands of each cluster away from
#'   zero at each p? (`NULL` without bands).}
#' \item{diss, hclust, silhouette}{the dissimilarity, the dendrogram and the
#'   silhouette (see [cluster::silhouette()]).}
#' \item{criterion}{the values of the criterion for each number of clusters.}
#' \item{measures}{a data frame with, for each cluster, its size, the average
#'   area between the curves and the mean curve (eq. 2 of the paper), the
#'   correlation distance \eqn{\rho_{dist} = 1 - } average pairwise correlation
#'   between the curves of the cluster (eq. 3), the average silhouette width and
#'   the average within-cluster dissimilarity.}
#' \item{settings}{the options used.}
#' @references Sottile, G. and Adelfio, G. (2019). Clusters of effects curves in
#'   quantile regression models. *Computational Statistics*, 34, 551-569.
#'   \doi{10.1007/s00180-018-0817-8}
#' @seealso [distshape()], [extract_curves()], [summary.clustEff()],
#'   [plot.clustEff()], [fpcac()].
#' @examples
#' ## Waveforms: three groups of noisy curves
#' set.seed(1234)
#' tt <- seq(0, 1, length.out = 100)
#' W <- cbind(sapply(1:5, function(i) sin(3 * pi * tt) + rnorm(100, 0, 0.1)),
#'            sapply(1:5, function(i) cos(3 * pi * tt) + rnorm(100, 0, 0.1)),
#'            sapply(1:5, function(i) sin(3 * pi * tt) * cos(pi * tt) + rnorm(100, 0, 0.1)))
#' fit <- clustEff(W, p = tt, k = c(2, 5), smooth = TRUE)
#' fit
#' summary(fit)
#' plot(fit)
#'
#' ## Effects curves with confidence bands: three groups of four curves
#' p <- seq(0.01, 0.99, 0.01)
#' B <- cbind(sapply(1:4, function(i) 1 + p + rnorm(1, 0, 0.05)),
#'            sapply(1:4, function(i) -1 + 2 * p^2 + rnorm(1, 0, 0.05)),
#'            sapply(1:4, function(i) 2 - 3 * p^3 + rnorm(1, 0, 0.05)))
#' fit2 <- clustEff(B, lower = B - 0.2, upper = B + 0.2, p = p, k = c(2, 5))
#' fit2$clusters
#' plot(fit2)
#' @export
clustEff <- function(curves, lower = NULL, upper = NULL, p = NULL, k = NULL, alpha = NULL,
                     step = c("both", "shape", "distance"), deriv = 2, combine = c("product", "mean"),
                     shape.tol = 1e-6,
                     criterion = c("height", "pi.out", "dist.rel", "silhouette", "conf.int"),
                     method = "ward.D2", diss = NULL, smooth = FALSE, df = NULL, conf.level = 0.9,
                     stand = FALSE, legacy = FALSE, verbose = FALSE, ...) {
  cl <- match.call()
  check_removed(list(...), "clustEff")
  criterion.missing <- missing(criterion)
  step <- match.arg(step)
  combine <- match.arg(combine)
  criterion <- match.arg(criterion)
  METHODS <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
  method <- METHODS[pmatch(method, METHODS)]
  if (is.na(method)) stop("invalid agglomeration method; use one of: ", paste(METHODS, collapse = ", "),
                          call. = FALSE)
  if (conf.level <= 0 || conf.level >= 1) stop("'conf.level' must be in (0, 1)", call. = FALSE)

  if (!is.null(p) && length(p) %in% 1:2 && NROW(curves) > 2 && !inherits(curves, "effects_curves"))
    stop("'p' must be the grid of the curves; in clustEff 1.0.0 the fourth argument is 'p', ",
         "not 'k': please use k = ", deparse(p), call. = FALSE)
  cv <- as_curves(curves, lower, upper, p)
  B <- cv$curves; L <- cv$lower; U <- cv$upper; pg <- cv$p; d2 <- cv$d2
  n <- nrow(B); q <- ncol(B)
  bands <- !is.null(L)

  if (smooth) {
    sm <- function(M) vapply(seq_len(ncol(M)), function(j) {
      f <- if (is.null(df)) stats::smooth.spline(pg, M[, j]) else stats::smooth.spline(pg, M[, j], df = df)
      stats::predict(f, pg)$y
    }, numeric(n))
    B0 <- B
    B <- sm(B); dimnames(B) <- dimnames(B0)
    if (bands) { L <- sm(L); U <- sm(U); dimnames(L) <- dimnames(U) <- dimnames(B0) }
    d2 <- NULL
  }
  if (stand) {
    mu <- colMeans(B); sdv <- apply(B, 2, stats::sd)
    sdv[sdv == 0] <- 1
    B <- sweep(sweep(B, 2, mu), 2, sdv, "/")
    if (bands) { L <- sweep(sweep(L, 2, mu), 2, sdv, "/"); U <- sweep(sweep(U, 2, mu), 2, sdv, "/") }
    if (!is.null(d2)) d2 <- sweep(d2, 2, sdv, "/")
  }

  # number of clusters
  if (is.null(k)) k <- if (q > 2) c(2, min(5, q - 1)) else 2
  if (!is.numeric(k) || !length(k) %in% 1:2 || any(k < 1) || any(k > q))
    stop(sprintf("'k' must be a number of clusters, or a range c(kmin, kmax), between 1 and %d", q), call. = FALSE)
  k <- sort(unique(round(k)))
  if (length(k) == 2 && k[1] < 2) stop("the range of 'k' must start from 2", call. = FALSE)
  if (legacy) {
    if (is.null(alpha)) alpha <- 0.5
    if (criterion.missing) criterion <- "silhouette"
  }
  if (length(k) == 2 && criterion %in% c("pi.out", "conf.int") && !bands)
    stop(sprintf("criterion = \"%s\" requires 'lower' and 'upper'", criterion), call. = FALSE)
  if (is.null(alpha)) alpha <- if (bands) 0.25 else 0.5

  if (is.null(diss)) {
    if (verbose) message("computing the dissimilarity matrix")
    diss <- distshape(B, p = pg, alpha = alpha, step = step, deriv = deriv, combine = combine,
                      d = if (deriv == 2) d2 else NULL, shape.tol = shape.tol, legacy = legacy)
  } else {
    if (!inherits(diss, "dist")) diss <- stats::as.dist(diss)
    if (attr(diss, "Size") != q) stop("'diss' must be a dissimilarity between the 'curves'", call. = FALSE)
  }
  hc <- stats::hclust(diss, method = method)
  w <- if (legacy) rep.int(1 / n, n) else trap_weights(pg)

  crit <- NULL
  if (length(k) == 1) {
    kopt <- k
  } else {
    ks <- if (criterion %in% c("pi.out", "dist.rel")) max(1, k[1] - 1):min(k[2] + 1, q) else k[1]:k[2]
    parts <- stats::cutree(hc, k = ks)
    if (!is.matrix(parts)) parts <- matrix(parts, ncol = 1)
    val <- vapply(seq_along(ks), function(m) {
      oc <- parts[, m]
      switch(criterion,
             pi.out = crit_pi_out(oc, B, L, U, w),
             dist.rel = crit_dist_rel(oc, B, w),
             silhouette = { s <- safe_silhouette(oc, diss); if (is.null(s)) NA_real_ else mean(s[, 3]) },
             height = { h <- rev(hc$height); if (ks[m] >= 2 && ks[m] <= length(h)) h[ks[m] - 1] - h[ks[m]] else NA_real_ },
             conf.int = crit_conf_int(oc, B, L, U, diss))
    }, numeric(1))
    names(val) <- ks
    crit <- val
    inr <- ks >= k[1] & ks <= k[2]
    kopt <- switch(criterion,
                   pi.out = {
                     # elbow: largest slowdown of the decrease of pi.out
                     dec <- -diff(val)
                     cur <- dec[-length(dec)] - dec[-1]
                     kk <- ks[-c(1, length(ks))]
                     ok <- kk >= k[1] & kk <= k[2] & is.finite(cur)
                     if (any(ok)) kk[ok][which.max(cur[ok])] else k[1]
                   },
                   dist.rel = {
                     # first number of clusters at which dist.rel stops decreasing
                     kr <- ks[inr]
                     up <- which(val[match(kr + 1, ks)] >= val[match(kr, ks)])
                     if (length(up)) kr[up[1]] else k[2]
                   },
                   silhouette = , conf.int = , height = {
                     if (all(is.na(val[inr]))) k[1] else ks[inr][which.max(val[inr])]
                   })
    if (verbose) message("selected number of clusters: ", kopt)
  }

  oc <- stats::cutree(hc, k = kopt)
  # clusters numbered by increasing size (as in clustEff 0.3.x)
  tab <- table(oc); ord <- order(tab)
  oc <- match(oc, as.integer(names(tab))[ord])
  names(oc) <- colnames(B)

  a <- (1 - conf.level) / 2
  M <- ML <- MU <- matrix(NA_real_, n, kopt, dimnames = list(NULL, paste0("cluster", seq_len(kopt))))
  for (j in seq_len(kopt)) {
    Bj <- B[, oc == j, drop = FALSE]
    M[, j] <- rowMeans(Bj)
    if (bands) {
      ML[, j] <- rowMeans(L[, oc == j, drop = FALSE])
      MU[, j] <- rowMeans(U[, oc == j, drop = FALSE])
    } else {
      ML[, j] <- apply(Bj, 1, stats::quantile, probs = a, names = FALSE)
      MU[, j] <- apply(Bj, 1, stats::quantile, probs = 1 - a, names = FALSE)
    }
  }
  signif <- if (bands) (ML > 0 | MU < 0) else NULL
  sil <- safe_silhouette(oc, diss)

  out <- list(call = cl, clusters = oc, k = kopt, p = pg, curves = B, lower = L, upper = U,
              mean = M, mean.lower = ML, mean.upper = MU, signif = signif,
              diss = diss, hclust = hc, silhouette = sil, criterion = crit,
              measures = cluster_measures(oc, B, M, w, diss, sil),
              settings = list(alpha = alpha, step = step, deriv = if (legacy) 1 else deriv,
                              combine = if (legacy) "mean" else combine, criterion = criterion,
                              method = method, k = k, smooth = smooth, stand = stand,
                              legacy = legacy, conf.level = conf.level, bands = bands))
  class(out) <- "clustEff"
  out
}

#' Silhouette, or NULL if not defined (one cluster, or as many clusters as curves).
#' @noRd
safe_silhouette <- function(oc, diss) {
  if (length(unique(oc)) < 2 || length(unique(oc)) >= length(oc)) return(NULL)
  s <- cluster::silhouette(oc, diss)
  if (is.matrix(s)) s else NULL
}

#' Average proportion of curves outside the mean bands of their cluster.
#' @noRd
crit_pi_out <- function(oc, B, L, U, w) {
  vals <- vapply(sort(unique(oc)), function(j) {
    idx <- oc == j
    lb <- rowMeans(L[, idx, drop = FALSE]); ub <- rowMeans(U[, idx, drop = FALSE])
    inside <- (B[, idx, drop = FALSE] >= lb) & (B[, idx, drop = FALSE] <= ub)
    mean(1 - colSums(inside * w))
  }, numeric(1))
  mean(vals)
}

#' l * max over clusters of the mean L1 distance between curves and cluster mean.
#' @noRd
crit_dist_rel <- function(oc, B, w) {
  u <- sort(unique(oc))
  vals <- vapply(u, function(j) {
    Bj <- B[, oc == j, drop = FALSE]
    mean(colSums(abs(Bj - rowMeans(Bj)) * w))
  }, numeric(1))
  length(u) * max(vals)
}

#' Criterion "conf.int" of clustEff 0.3.x.
#' @noRd
crit_conf_int <- function(oc, B, L, U, diss) {
  u <- sort(unique(oc))
  if (length(u) < 2) return(NA_real_)
  tempDist <- vapply(u, function(j) {
    idx <- oc == j
    BR <- rowMeans(U[, idx, drop = FALSE]); BL <- rowMeans(L[, idx, drop = FALSE])
    mean((colMeans(B[, idx, drop = FALSE] - BR < 0) + colMeans(B[, idx, drop = FALSE] - BL >= 0)) / 2)
  }, numeric(1))
  s <- safe_silhouette(oc, diss)
  if (is.null(s)) return(NA_real_)
  sw <- tapply(s[, 3], oc, mean)
  sum(tempDist * sw) / sum(tempDist)
}

#' Within-cluster summaries (Sottile and Adelfio, 2019, eq. 2-3).
#' @noRd
cluster_measures <- function(oc, B, M, w, diss, sil) {
  Dm <- as.matrix(diss)
  k <- ncol(M)
  res <- lapply(seq_len(k), function(j) {
    idx <- which(oc == j)
    Bj <- B[, idx, drop = FALSE]
    area <- mean(colSums(abs(Bj - M[, j]) * w))
    rho <- if (length(idx) > 1) {
      R <- suppressWarnings(stats::cor(Bj))
      1 - mean(R[upper.tri(R)], na.rm = TRUE)
    } else 0
    dis <- if (length(idx) > 1) mean(Dm[idx, idx][upper.tri(Dm[idx, idx])]) else 0
    silw <- if (!is.null(sil)) mean(sil[idx, 3]) else NA_real_
    data.frame(size = length(idx), area = area, rho.dist = rho, silhouette = silw, dissimilarity = dis)
  })
  out <- do.call(rbind, res)
  rownames(out) <- paste0("cluster", seq_len(k))
  out
}

#' @rdname clustEff
#' @param x an object of class `"clustEff"`.
#' @export
print.clustEff <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Clusters of curves: ", x$k, " clusters of sizes ", paste(table(x$clusters), collapse = ", "),
      "\n", sep = "")
  s <- x$settings
  cat("Dissimilarity: ", if (s$legacy) "clustEff 0.3.x (legacy)" else
        paste0("step = \"", s$step, "\", alpha = ", s$alpha, ", derivative = ", s$deriv,
               if (s$step == "both") paste0(", combine = \"", s$combine, "\"")), "\n", sep = "")
  cat("Agglomeration method: ", s$method,
      if (length(s$k) == 2) paste0("   criterion: ", s$criterion), "\n\n", sep = "")
  invisible(x)
}

#' Summarize Clusters of Curves
#'
#' @param object an object of class `"clustEff"`.
#' @param x an object of class `"summary.clustEff"`.
#' @param digits number of significant digits.
#' @param ... not used.
#' @return An object of class `"summary.clustEff"`: a list with the number of
#'   clusters, the number of curves and grid points, the options used, the
#'   clustering table, the within-cluster measures (with their weighted average)
#'   and the criterion used to select the number of clusters.
#' @seealso [clustEff()]
#' @export
summary.clustEff <- function(object, ...) {
  m <- object$measures
  avg <- vapply(m[, -1], function(v) sum(v * m$size, na.rm = TRUE) / sum(m$size[!is.na(v)]), numeric(1))
  m2 <- rbind(m, Average = c(sum(m$size), avg))
  out <- list(call = object$call, k = object$k, n = nrow(object$curves), q = ncol(object$curves),
              settings = object$settings, table = table(object$clusters), measures = m2,
              criterion = object$criterion, members = split(names(object$clusters), object$clusters))
  class(out) <- "summary.clustEff"
  out
}

#' @rdname summary.clustEff
#' @export
print.summary.clustEff <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of curves:", x$q, "   n. of grid points:", x$n, "\n")
  cat("n. of clusters:", x$k, if (length(x$settings$k) == 2) paste0("(selected by \"", x$settings$criterion,
                                                                     "\" in ", x$settings$k[1], "-",
                                                                     x$settings$k[2], ")"), "\n\n")
  cat("Clusters:\n")
  for (j in seq_along(x$members)) cat("  ", j, ": ", paste(x$members[[j]], collapse = ", "), "\n", sep = "")
  cat("\nWithin-cluster measures:\n")
  print(format(x$measures, digits = digits), quote = FALSE)
  if (!is.null(x$criterion)) {
    cat("\nCriterion (", x$settings$criterion, ") by number of clusters:\n", sep = "")
    print(signif(x$criterion, digits))
  }
  cat("\n")
  invisible(x)
}

#' Plot Clusters of Curves
#'
#' @param x an object of class `"clustEff"`.
#' @param type what to plot: `"clusters"` (mean curves with their bands, or the
#'   curves of one cluster if `which` has length one), `"dendrogram"`,
#'   `"silhouette"` (boxplots of the silhouette widths by cluster),
#'   `"dissimilarity"` (boxplots of the within-cluster dissimilarities) or
#'   `"criterion"` (the criterion against the number of clusters).
#' @param which the clusters to plot (for `type = "clusters"`).
#' @param polygon logical: draw the bands as shaded areas (otherwise dashed lines)?
#' @param col colours of the clusters.
#' @param ... further graphical parameters (`xlab`, `ylab`, `main`, `ylim`, `lwd`).
#' @return No return value, called for its side effects.
#' @seealso [clustEff()]
#' @export
plot.clustEff <- function(x, type = c("clusters", "dendrogram", "silhouette", "dissimilarity", "criterion"),
                          which, polygon = TRUE, col, ...) {
  type <- match.arg(type)
  L <- list(...)
  k <- x$k
  if (missing(which)) which <- seq_len(k)
  if (any(which < 1 | which > k)) stop(sprintf("'which' must be in 1-%d", k), call. = FALSE)
  if (missing(col)) col <- grDevices::hcl.colors(max(k, 2), "Dark 3")[seq_len(k)]
  col <- rep_len(col, k)
  arg <- function(name, default) if (is.null(L[[name]])) default else L[[name]]
  switch(type,
    dendrogram = {
      plot(x$hclust, main = arg("main", "Dendrogram"), xlab = arg("xlab", ""), ylab = arg("ylab", "Height"),
           sub = "")
      h <- c(rev(x$hclust$height), 0)
      if (k > 1 && k <= length(h)) abline(h = mean(h[c(k - 1, k)]), col = 2, lty = 2)
    },
    clusters = {
      pg <- x$p
      one <- length(which) == 1
      yl <- arg("ylim", if (one) range(x$curves[, x$clusters == which], x$mean.lower[, which], x$mean.upper[, which])
                        else range(x$mean[, which], x$mean.lower[, which], x$mean.upper[, which]))
      plot(range(pg), yl, type = "n", xlab = arg("xlab", "p"), ylab = arg("ylab", "s(p)"),
           main = arg("main", if (one) paste("Cluster", which) else ""))
      for (j in which) {
        if (polygon) graphics::polygon(c(pg, rev(pg)), c(x$mean.lower[, j], rev(x$mean.upper[, j])),
                                       col = grDevices::adjustcolor(col[j], alpha.f = 0.25), border = NA)
        else { lines(pg, x$mean.lower[, j], lty = 2, col = col[j]); lines(pg, x$mean.upper[, j], lty = 2, col = col[j]) }
        if (one) matlines(pg, x$curves[, x$clusters == j, drop = FALSE], lty = 3, col = col[j])
        lines(pg, x$mean[, j], lwd = arg("lwd", 2), col = col[j])
      }
      abline(h = 0, lty = 3)
      if (!one) graphics::legend("topleft", legend = paste("cluster", which), col = col[which], lwd = 2, bty = "n")
    },
    silhouette = {
      if (is.null(x$silhouette)) stop("only one cluster: no silhouette", call. = FALSE)
      graphics::boxplot(split(x$silhouette[, 3], x$clusters), xlab = arg("xlab", "cluster"),
                        ylab = arg("ylab", "silhouette width"), main = arg("main", "Silhouette within cluster"),
                        col = grDevices::adjustcolor(col, alpha.f = 0.4))
    },
    dissimilarity = {
      Dm <- as.matrix(x$diss)
      vals <- lapply(seq_len(k), function(j) {
        idx <- which(x$clusters == j)
        if (length(idx) > 1) Dm[idx, idx][upper.tri(Dm[idx, idx])] else 0
      })
      graphics::boxplot(vals, names = seq_len(k), ylim = arg("ylim", c(0, 1)), xlab = arg("xlab", "cluster"),
                        ylab = arg("ylab", "dissimilarity"), main = arg("main", "Dissimilarity within cluster"),
                        col = grDevices::adjustcolor(col, alpha.f = 0.4))
    },
    criterion = {
      if (is.null(x$criterion)) stop("the number of clusters was fixed: no criterion to plot", call. = FALSE)
      ks <- as.integer(names(x$criterion))
      plot(ks, x$criterion, type = "b", pch = 19, xlab = arg("xlab", "number of clusters"),
           ylab = arg("ylab", x$settings$criterion), main = arg("main", ""))
      abline(v = x$k, lty = 2, col = 2)
    })
  invisible(NULL)
}

#' Removed arguments of clustEff 0.3.x: informative errors.
#' @noRd
check_removed <- function(dots, fun) {
  if (!length(dots)) return(invisible())
  nm <- names(dots)
  if (is.null(nm) || any(nm == "")) stop("unused unnamed argument(s)", call. = FALSE)
  map <- c(Beta = "curves", Beta.lower = "lower", Beta.upper = "upper", diss.mat = "diss",
           cut.method = "criterion", approx.spline = "smooth", nbasis = "smooth and df",
           plot = "plot()", trace = "verbose", ask = "k (a single value) after looking at plot(fit, \"dendrogram\")",
           X = "curves", fd = "curves (a matrix of curves)", niter = "nstart", Ksteps = "iter.max",
           seed = "set.seed() before the call", disp = "verbose", method = "criterion")
  if (fun == "distshape") map["trace"] <- "nothing (distshape() is silent)"
  old <- intersect(nm, names(map))
  if (length(old))
    stop(sprintf("argument(s) of clustEff 0.3.x removed in version 1.0.0: %s. Use instead: %s (see ?%s)",
                 paste(old, collapse = ", "), paste(map[old], collapse = "; "), fun), call. = FALSE)
  stop("unused argument(s): ", paste(nm, collapse = ", "), call. = FALSE)
}
