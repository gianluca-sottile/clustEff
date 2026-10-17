# Input handling: matrices of curves and effects curves of quantile regression
# coefficients models.

#' Effects Curves of Quantile Regression Coefficients Models
#'
#' Extracts the estimated quantile regression coefficient functions
#' \eqn{\hat\beta_j(p)} (the *effects curves*), their pointwise confidence
#' bands and their second derivatives from one or more models fitted with
#' [qrcm::iqr()], in the format used by [clustEff()].
#'
#' @param object an object of class `"iqr"`, or a (possibly named) list of
#'   such objects, for example one model for each of several responses.
#' @param p the grid of quantile orders on which the curves are evaluated.
#' @param which the coefficients to extract: a vector of names or of indices
#'   of the rows of the coefficient matrix. By default all of them except the
#'   intercept.
#' @param level confidence level of the pointwise bands (default 0.95).
#'
#' @details
#' The bands are \eqn{\hat\beta_j(p) \pm z_{(1+level)/2}\,se\{\hat\beta_j(p)\}},
#' with standard errors from `predict(object, type = "beta", se = TRUE)`.
#' The second derivatives \eqn{\hat\beta_j''(p)}, used by the shape component of
#' the dissimilarity of [distshape()], are computed by central differences of the
#' smooth fitted functions (step \eqn{10^{-3}}; exact up to rounding errors when
#' \eqn{\beta(p)} is a polynomial of degree up to 3 in \eqn{p}, and very accurate
#' otherwise, unlike derivatives of interpolating splines of the estimates).
#'
#' With a list of models (e.g., the same covariates for different responses, as in
#' the multivariate case of Sottile and Adelfio, 2019), the curves are named
#' `"model:coefficient"`, where `model` is the name of the list element (or its
#' position).
#'
#' @return An object of class `"effects_curves"`, a list with items `curves`,
#'   `lower`, `upper` and `d2` (matrices with one row for each value of `p` and
#'   one column for each curve), `p`, `level` and `info` (a data frame with the
#'   model and the coefficient of each curve).
#' @references Sottile, G. and Adelfio, G. (2019). Clusters of effects curves in
#'   quantile regression models. *Computational Statistics*, 34, 551-569.
#'   \doi{10.1007/s00180-018-0817-8}
#' @seealso [clustEff()], [qrcm::iqr()]
#' @examples
#' \donttest{
#' if (requireNamespace("qrcm", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   x <- runif(n, 0, 5)
#'   # five responses with three different effects of x
#'   slopes <- list(function(u) 0.5 + u, function(u) 0.5 + u, function(u) -1 + 0.5 * u,
#'                  function(u) -1 + 0.5 * u, function(u) 2 * u^2)
#'   fits <- lapply(slopes, function(b) {
#'     u <- runif(n)
#'     y <- qnorm(u) + b(u) * x
#'     qrcm::iqr(y ~ x, formula.p = ~ I(p) + I(p^2))
#'   })
#'   names(fits) <- paste0("y", 1:5)
#'   ec <- extract_curves(fits, which = "x")
#'   ec
#'   fit <- clustEff(ec, k = 3)
#'   fit$clusters
#' }
#' }
#' @export
extract_curves <- function(object, p = seq(0.01, 0.99, 0.01), which, level = 0.95) {
  if (!requireNamespace("qrcm", quietly = TRUE)) stop("package 'qrcm' is required", call. = FALSE)
  if (inherits(object, "iqr")) object <- list(object)
  if (!is.list(object) || !all(vapply(object, inherits, TRUE, what = "iqr")))
    stop("'object' must be an object of class \"iqr\" or a list of such objects", call. = FALSE)
  if (any(p <= 0 | p >= 1)) stop("0 < p < 1 is required", call. = FALSE)
  if (level <= 0 || level >= 1) stop("'level' must be in (0, 1)", call. = FALSE)
  p <- sort(p)
  z <- stats::qnorm((1 + level) / 2)
  h <- min(1e-3, min(p) / 2, (1 - max(p)) / 2)
  mnames <- names(object)
  if (is.null(mnames)) mnames <- if (length(object) > 1) paste0("model", seq_along(object)) else ""
  mnames[mnames == "" & length(object) > 1] <- paste0("model", base::which(mnames == ""))
  out <- list(curves = NULL, lower = NULL, upper = NULL, d2 = NULL)
  info <- NULL
  for (m in seq_along(object)) {
    obj <- object[[m]]
    rn <- rownames(obj$coefficients)
    sel <- if (missing(which)) setdiff(seq_along(rn), match("(Intercept)", rn, 0L)) else {
      if (is.character(which)) match(which, rn) else as.integer(which)
    }
    if (anyNA(sel) || any(sel < 1 | sel > length(rn)))
      stop(sprintf("'which' does not match the coefficients of model %d (%s)", m, paste(rn, collapse = ", ")),
           call. = FALSE)
    pr <- predict(obj, type = "beta", p = p, se = TRUE)
    pm <- predict(obj, type = "beta", p = p - h, se = FALSE)
    pp <- predict(obj, type = "beta", p = p + h, se = FALSE)
    for (j in sel) {
      b <- pr[[j]]$beta; se <- pr[[j]]$se
      out$curves <- cbind(out$curves, b)
      out$lower <- cbind(out$lower, b - z * se)
      out$upper <- cbind(out$upper, b + z * se)
      out$d2 <- cbind(out$d2, (pp[[j]]$beta - 2 * b + pm[[j]]$beta) / h^2)
      info <- rbind(info, data.frame(model = mnames[m], coefficient = rn[j], stringsAsFactors = FALSE))
    }
  }
  if (is.null(out$curves)) stop("no coefficients selected", call. = FALSE)
  nm <- if (length(object) == 1) info$coefficient else paste(info$model, info$coefficient, sep = ":")
  nm <- make.unique(nm)
  for (a in names(out)) { dimnames(out[[a]]) <- list(NULL, nm) }
  out$p <- p
  out$level <- level
  out$info <- info
  class(out) <- "effects_curves"
  out
}

#' @rdname extract_curves
#' @param x an object of class `"effects_curves"`.
#' @param ... not used.
#' @export
print.effects_curves <- function(x, ...) {
  cat("\nEffects curves:", ncol(x$curves), "curves evaluated at", length(x$p), "values of p in [",
      format(min(x$p)), ",", format(max(x$p)), "]\n")
  cat("Pointwise bands at level", x$level, "\n\n")
  print(x$info, row.names = FALSE)
  cat("\n")
  invisible(x)
}

#' Standardize the input of clustEff() and distshape().
#' @return list(curves, lower, upper, p, d2)
#' @noRd
as_curves <- function(curves, lower = NULL, upper = NULL, p = NULL) {
  d2 <- NULL
  if (inherits(curves, "iqr") || (is.list(curves) && !is.data.frame(curves) && length(curves) &&
                                  all(vapply(curves, inherits, TRUE, what = "iqr")))) {
    curves <- extract_curves(curves)
  }
  if (inherits(curves, "effects_curves")) {
    if (is.null(lower)) lower <- curves$lower
    if (is.null(upper)) upper <- curves$upper
    if (!is.null(p) && !isTRUE(all.equal(as.numeric(p), curves$p)))
      stop("the grid 'p' of an \"effects_curves\" object cannot be changed: use extract_curves(..., p = )",
           call. = FALSE)
    p <- curves$p
    d2 <- curves$d2
    curves <- curves$curves
  }
  curves <- as.matrix(curves)
  if (!is.numeric(curves)) stop("'curves' must be a numeric matrix (one column for each curve)", call. = FALSE)
  if (anyNA(curves)) stop("'curves' contains missing values", call. = FALSE)
  n <- nrow(curves); q <- ncol(curves)
  if (q < 2) stop("at least two curves are needed", call. = FALSE)
  if (n < 4) stop("the curves must be evaluated at four or more points", call. = FALSE)
  if (is.null(colnames(curves))) colnames(curves) <- paste0("X", seq_len(q))
  if (is.null(p)) p <- seq_len(n) / n
  if (length(p) != n) stop("'p' must have one value for each row of 'curves'", call. = FALSE)
  if (is.unsorted(p, strictly = TRUE)) stop("'p' must be strictly increasing", call. = FALSE)
  if (xor(is.null(lower), is.null(upper))) stop("please supply both 'lower' and 'upper'", call. = FALSE)
  if (!is.null(lower)) {
    lower <- as.matrix(lower); upper <- as.matrix(upper)
    if (any(dim(lower) != dim(curves)) || any(dim(upper) != dim(curves)))
      stop("'lower' and 'upper' must have the same dimensions as 'curves'", call. = FALSE)
    if (anyNA(lower) || anyNA(upper)) stop("'lower' and 'upper' cannot contain missing values", call. = FALSE)
    dimnames(lower) <- dimnames(upper) <- dimnames(curves)
  }
  if (!is.null(d2)) dimnames(d2) <- dimnames(curves)
  list(curves = curves, lower = lower, upper = upper, p = p, d2 = d2)
}

#' Trapezoidal weights for the integral over the grid p, normalised to sum to 1
#' (i.e., the integral is divided by the length of the interval).
#' @noRd
trap_weights <- function(p) {
  n <- length(p)
  h <- diff(p)
  w <- c(h, 0) / 2 + c(0, h) / 2
  w / sum(w)
}

#' Derivative of order 'deriv' of each column of curves on the grid p.
#' @noRd
curve_derivative <- function(curves, p, deriv = 2) {
  vapply(seq_len(ncol(curves)), function(j) stats::splinefun(p, curves[, j], method = "fmm")(p, deriv = deriv),
         numeric(nrow(curves)))
}
