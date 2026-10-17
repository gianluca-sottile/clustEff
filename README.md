
<!-- README.md is generated from README.Rmd. Please edit that file -->


# clustEff

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/clustEff)](https://CRAN.R-project.org/package=clustEff)
[![R-CMD-check](https://github.com/gianluca-sottile/clustEff/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gianluca-sottile/clustEff/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

'clustEff' clusters curves with the shape-distance dissimilarity of Sottile and
Adelfio (2019): the effects curves of quantile regression coefficients models
(fitted with 'qrcm', with one or more responses) and curves observed in time. It
also implements the FPCAC algorithm (trimmed k-means of functional principal
component scores) for waveforms.

## Installation

```r
install.packages("clustEff")
# development version
remotes::install_github("gianluca-sottile/clustEff")
```

## Example

``` r
library(clustEff)
set.seed(1234)
tt <- seq(0, 1, length.out = 100)
W <- cbind(sapply(1:5, function(i) sin(3 * pi * tt) + rnorm(100, 0, 0.1)),
           sapply(1:5, function(i) cos(3 * pi * tt) + rnorm(100, 0, 0.1)),
           sapply(1:5, function(i) sin(3 * pi * tt) * cos(pi * tt) + rnorm(100, 0, 0.1)))
fit <- clustEff(W, p = tt, k = c(2, 5), smooth = TRUE)
fit
#> 
#> Call:
#> clustEff(curves = W, p = tt, k = c(2, 5), smooth = TRUE)
#> 
#> Clusters of curves: 3 clusters of sizes 5, 5, 5
#> Dissimilarity: step = "both", alpha = 0.5, derivative = 2, combine = "product"
#> Agglomeration method: ward.D2   criterion: height
#> 
fit$clusters
#>  X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 X11 X12 X13 X14 X15 
#>   1   1   1   1   1   2   2   2   2   2   3   3   3   3   3 
```

## References

* Sottile, G. and Adelfio, G. (2019). Clusters of effects curves in quantile regression
  models. *Computational Statistics*, 34, 551-569. <https://doi.org/10.1007/s00180-018-0817-8>
* Adelfio, G., Chiodi, M., D'Alessandro, A. and Luzio, D. (2011). FPCA algorithm for
  waveform clustering. *Journal of Communication and Computer*, 8, 494-502.
