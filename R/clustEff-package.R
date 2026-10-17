#' clustEff: Clusters of Effects Curves in Quantile Regression Models
#'
#' Hierarchical clustering of curves based on the shape-distance dissimilarity
#' of Sottile and Adelfio (2019): effects curves of quantile regression
#' coefficients models (univariate or multivariate responses, see
#' [extract_curves()] and [clustEff()]) and curves observed in time. The package
#' also implements the functional principal components clustering (FPCAC) of
#' Adelfio et al. (2011) for waveforms ([fpcac()], [opt.fpcac()]).
#'
#' @references
#' Sottile, G. and Adelfio, G. (2019). Clusters of effects curves in quantile
#' regression models. *Computational Statistics*, 34, 551-569.
#' \doi{10.1007/s00180-018-0817-8}
#'
#' Adelfio, G., Chiodi, M., D'Alessandro, A. and Luzio, D. (2011). FPCA algorithm
#' for waveform clustering. *Journal of Communication and Computer*, 8, 494-502.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom graphics abline lines matlines plot
#' @importFrom stats predict
NULL

#' Defunct Functions and Arguments in clustEff
#'
#' The following were removed in 'clustEff' 1.0.0:
#' * `extract.object()`: use [extract_curves()] on models fitted with
#'   [qrcm::iqr()] (e.g., `extract_curves(lapply(responses, function(y) iqr(y ~ x)))`);
#' * the arguments `Beta`, `Beta.lower`, `Beta.upper`, `diss.mat`, `cut.method`,
#'   `approx.spline`, `nbasis`, `plot`, `trace` and `ask` of [clustEff()]: use
#'   `curves`, `lower`, `upper`, `diss`, `criterion`, `smooth`, `df`, `plot()`,
#'   `verbose`, and a fixed `k` after inspecting `plot(fit, "dendrogram")`;
#' * the arguments `X`, `fd`, `niter`, `Ksteps`, `seed` and `disp` of [fpcac()]
#'   and [opt.fpcac()]: use `curves`, a matrix of curves, `nstart`, `iter.max`,
#'   [set.seed()] and `verbose`; `method` of [opt.fpcac()] is now `criterion`.
#'
#' @param ... ignored.
#' @return No return value: an error is raised.
#' @name clustEff-defunct
#' @keywords internal
NULL

#' @rdname clustEff-defunct
#' @export
extract.object <- function(...) {
  .Defunct("extract_curves", package = "clustEff",
           msg = paste("'extract.object()' was removed in clustEff 1.0.0: fit the models with qrcm::iqr()",
                       "and use 'extract_curves(model)' or 'extract_curves(list_of_models)'"))
}
