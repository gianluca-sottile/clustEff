#' Shape-Distance Dissimilarity between Curves
#'
#' Computes the dissimilarity between curves of Sottile and Adelfio (2019,
#' eq. 1), which combines the concordance of their *shapes* (the signs of their
#' second derivatives) and their *closeness* (their distance relative to the
#' distances between all pairs of curves), pointwise on a grid.
#'
#' @param curves a numeric matrix with one column for each curve and one row
#'   for each value of `p`, or an object returned by [extract_curves()].
#' @param p the grid on which the curves are evaluated (quantile orders, or
#'   times). By default `p = (1:n)/n`, or the grid of an `"effects_curves"` object.
#' @param alpha the order of the quantile of the pointwise distances between all
#'   pairs of curves that defines two curves as close (default 0.5). Sottile and
#'   Adelfio (2019) suggest 0.25 for effects curves and 0.5 for waveforms.
#' @param step which components are used: `"both"` (default), `"shape"` or
#'   `"distance"`.
#' @param deriv the order of the derivative used for the shape component: 2
#'   (default, convexity, as in the paper) or 1 (monotonicity).
#' @param combine how the two components are combined when `step = "both"`:
#'   `"product"` (default, eq. 1 of the paper: two curves are similar at p if
#'   they are concordant in shape *and* close) or `"mean"` (average of the two
#'   proportions, as in 'clustEff' 0.3.x).
#' @param d an optional matrix of derivatives of order `deriv` of the curves (same
#'   dimensions as `curves`). If missing they are computed by interpolating
#'   splines, or taken from an `"effects_curves"` object (second derivatives).
#' @param shape.tol derivatives smaller than `shape.tol * max(abs(d))` in absolute
#'   value are treated as zero, so that, e.g., two linear curves are concordant in
#'   shape.
#' @param legacy logical: if `TRUE`, the dissimilarity of 'clustEff' 0.3.x is
#'   computed (first derivative, `combine = "mean"`, grid `(1:n)/n`, equal weights);
#'   the arguments `p`, `deriv`, `combine`, `d` and `shape.tol` are then ignored.
#' @param ... arguments of 'clustEff' 0.3.x that were removed (they give an
#'   informative error).
#'
#' @details
#' For two curves \eqn{\beta_i} and \eqn{\beta_{i'}},
#' \deqn{d_{shape}(p) = I\{sign(\beta_i''(p)) = sign(\beta_{i'}''(p))\},\qquad
#'       d_{distance}(p) = I\{|\beta_i(p) - \beta_{i'}(p)| \le f_\alpha(p)\},}{%
#' d_shape(p) = I(sign(b_i''(p)) == sign(b_i'''(p))), d_distance(p) = I(|b_i(p) - b_i'(p)| <= f(p)),}
#' where \eqn{f_\alpha(p)} is the \eqn{\alpha}-quantile of the distances
#' \eqn{|\beta_j(p) - \beta_{j'}(p)|} between all pairs of curves at \eqn{p}. The
#' dissimilarity is
#' \deqn{d(i, i') = 1 - \int d_{shape}(p)\, d_{distance}(p)\,dp,}{d(i, i') = 1 - int d_shape(p) d_distance(p) dp,}
#' where the integral (normalized by the length of the interval) is computed with
#' the trapezoidal rule on the grid `p`. The computations are vectorized over the
#' pairs of curves and need \eqn{O(q^2)} memory for \eqn{q} curves.
#'
#' The paper defines the shape concordance as \eqn{sign(\beta_i'') \times
#' sign(\beta_{i'}'') = 1}; here derivatives below `shape.tol` are set to zero
#' and two zero derivatives are concordant, which avoids that, e.g., two linear
#' curves are always discordant. The shape component alone is a Hamming distance
#' between the sign patterns; the distance component, based on a threshold, does
#' not satisfy the triangle inequality in general.
#'
#' @return An object of class `"dist"` (see [stats::dist()]) with attributes
#'   `alpha`, `step`, `deriv`, `combine` and `legacy`.
#' @references Sottile, G. and Adelfio, G. (2019). Clusters of effects curves in
#'   quantile regression models. *Computational Statistics*, 34, 551-569.
#'   \doi{10.1007/s00180-018-0817-8}
#' @seealso [clustEff()]
#' @examples
#' tt <- seq(0, 1, length.out = 100)
#' curves <- cbind(a1 = sin(3 * pi * tt), a2 = sin(3 * pi * tt) + 0.1,
#'                 b1 = cos(3 * pi * tt), b2 = cos(3 * pi * tt) - 0.1)
#' distshape(curves, p = tt)
#' distshape(curves, p = tt, step = "shape")
#' @export
distshape <- function(curves, p = NULL, alpha = 0.5, step = c("both", "shape", "distance"),
                      deriv = 2, combine = c("product", "mean"), d = NULL, shape.tol = 1e-6,
                      legacy = FALSE, ...) {
  check_removed(list(...), "distshape")
  step <- match.arg(step)
  combine <- match.arg(combine)
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1)
    stop("'alpha' must be a number in [0, 1]", call. = FALSE)
  if (!deriv %in% 1:2) stop("'deriv' must be 1 or 2", call. = FALSE)
  cv <- as_curves(curves, p = p)
  B <- cv$curves
  n <- nrow(B); q <- ncol(B)
  if (legacy) {
    pg <- seq_len(n) / n
    w <- rep.int(1 / n, n)
    deriv <- 1; combine <- "mean"
  } else {
    pg <- cv$p
    w <- trap_weights(pg)
  }

  if (step != "distance") {
    D <- if (legacy) curve_derivative(B, pg, 1)
         else if (!is.null(d)) as.matrix(d)
         else if (deriv == 2 && !is.null(cv$d2)) cv$d2
         else curve_derivative(B, pg, deriv)
    if (any(dim(D) != dim(B))) stop("'d' must have the same dimensions as 'curves'", call. = FALSE)
    S <- sign(D)
    if (!legacy) S[abs(D) <= shape.tol * max(abs(D))] <- 0
  }

  up <- upper.tri(diag(q))
  shape <- close <- both <- matrix(0, q, q)
  for (l in seq_len(n)) {
    if (step != "distance") {
      s <- S[l, ]
      A <- if (legacy) outer(s, s) == 1 else outer(s, s, "==")
    }
    if (step != "shape") {
      b <- B[l, ]
      Dl <- abs(outer(b, b, "-"))
      cut <- stats::quantile(Dl[up], alpha, names = FALSE)
      C <- Dl <= cut
    }
    switch(step,
           both = if (combine == "product") both <- both + w[l] * (A & C)
                  else { shape <- shape + w[l] * A; close <- close + w[l] * C },
           shape = shape <- shape + w[l] * A,
           distance = close <- close + w[l] * C)
  }
  sim <- switch(step,
                both = if (combine == "product") both else 0.5 * (shape + close),
                shape = shape, distance = close)
  out <- stats::as.dist(1 - sim)
  attr(out, "Labels") <- colnames(B)
  attr(out, "alpha") <- alpha
  attr(out, "step") <- step
  attr(out, "deriv") <- deriv
  attr(out, "combine") <- combine
  attr(out, "legacy") <- legacy
  out
}
