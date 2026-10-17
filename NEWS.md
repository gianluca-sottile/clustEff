# clustEff 1.0.0

## Breaking changes

* `distshape()` now computes by default the dissimilarity of Sottile and Adelfio
  (2019, eq. 1): concordance of the signs of the **second** derivatives and
  **product** of the shape and distance components, integrated with the
  trapezoidal rule on the actual grid `p`. Version 0.3.x used the first derivative
  and the mean of the two components on the grid `(1:n)/n`; use `legacy = TRUE`
  (or `deriv = 1`, `combine = "mean"`) to reproduce it.
* `clustEff()`: the number of clusters is selected by default with the largest gap
  between the heights of the dendrogram (`criterion = "height"`); the criteria of the
  paper (`"pi.out"`, `"dist.rel"`) and the silhouette (default of 0.3.x, used with
  `legacy = TRUE`) are options. The default `alpha` is 0.25 with confidence bands
  (effects curves) and 0.5 otherwise, as suggested in the paper.
* 'fda', 'ggplot2' and 'ggpubr' are no longer used: the functional PCA of `fpcac()`
  is computed internally on a B-spline basis (`fpca_bspline()`), graphics use base R;
  'qrcm' is suggested. `library(clustEff)` no longer attaches 'qrcm', 'cluster',
  'fda', 'ggpubr' and 'ggplot2'.
* Removed (with an informative error, see `?"clustEff-defunct"`):

  | Removed in 1.0.0 | Use instead |
  |---|---|
  | `extract.object(Y, X, ...)` | fit the models with `qrcm::iqr()`, then `extract_curves(model)` or `extract_curves(list_of_models)` |
  | `clustEff(Beta, Beta.lower, Beta.upper, ...)` (named) | `clustEff(curves, lower, upper, ...)` (positional calls still work) |
  | `diss.mat`, `cut.method`, `approx.spline`, `nbasis`, `trace`, `ask`, `plot` | `diss`, `criterion`, `smooth`/`df`, `verbose`, a fixed `k`, `plot()` |
  | `fpcac(X, fd, niter, Ksteps, seed, disp)` | `fpcac(curves, t, nstart, iter.max)`, `set.seed()`, `verbose` |
  | `opt.fpcac(..., method)` | `opt.fpcac(..., criterion)` |

## New features

* `clustEff()` accepts a matrix of curves, the output of `extract_curves()`, or
  directly models fitted with `qrcm::iqr()`; the grid `p` can be supplied; exact
  second derivatives are used for effects curves.
* `extract_curves()`: effects curves, pointwise bands (any `level`) and second
  derivatives from one or more `qrcm::iqr()` fits.
* Vectorized dissimilarity with O(q^2) memory (it was O(n q^2)).
* `smooth = TRUE` smooths noisy curves with `smooth.spline()` before computing the
  dissimilarity.
* Within-cluster measures of the paper (area, correlation distance) plus silhouette
  and dissimilarity in `summary()`; new `plot()` types (`"clusters"`,
  `"dendrogram"`, `"silhouette"`, `"dissimilarity"`, `"criterion"`).
* `fpcac()`: vectorized trimmed k-means with k-means++ initialization
  (`init = "random"` as before), handling of empty clusters; `opt.fpcac()` excludes
  the trimmed curves from the criteria and has `print()` and `plot()` methods.

## Bug fixes

* `clustEff()`: wrong check of `Beta.upper`; errors with `k[1] > 2`; `conf.level`
  ignored without bands (quantiles 0.1 and 0.9 were always used); graphical
  parameters are no longer modified.
* `fpcac()`: failure when the objective function exceeded 1e6; wrong evaluation
  grid when an `fd` object was supplied.
* CITATION: authors and reference of the FPCA paper.

# clustEff 0.3.0

* Some bugs fixed.

# clustEff 0.2.0

* New implementation of the FPCAC algorithm and S3 methods for the class `fpcac`.
* `clustEff()` is more stable and faster; changes in `plot.clustEff()` and
  `summary.clustEff()`.
