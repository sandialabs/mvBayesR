[![R-CMD-check](https://github.com/sandialabs/mvBayesR/actions/workflows/r.yml/badge.svg)](https://github.com/sandialabs/mvBayesR/actions/workflows/r.yml)
[![CRAN status](https://www.r-pkg.org/badges/version/mvBayes)](https://CRAN.R-project.org/package=mvBayes)


# mvBayes

![](man/figures/logo.png)

An R implementation of the multivariate Bayesian regression (mvBayes) framework. 
Decomposes a multivariate/functional response using a user-specified orthogonal basis decomposition, 
and then models each basis component independently using an arbitrary user-specified (univariate) 
Bayesian regression model. Includes prediction and plotting methods. (Francom et al., 2025 <DOI:10.1137/24M164409>)


## Examples
* [Friedman Example](inst/friedman_demo.R) - An extension of the "Friedman function" to functional response. The Bayesian regression model here is BASS (Bayesian Adaptive Smoothing Splines, see https://github.com/lanl/BASS)


### Installation
------------------------------------------------------------------------------
v1.2.2 is on [CRAN](https://cran.r-project.org/package=mvBayes) and can
be installed as

``` r
install.packages("mvBayes")
```

For a more up to date, but may not be stable version from git
repository.

1.  Download zip or tar.gz of package or clone repository
2.  Install into R (\> 4.3.0)

``` r
library(devtools)
install_github("sandialabs/mvBayesR")
```

------------------------------------------------------------------------------

## References


************

Author: Gavin Q. Collins and J. Derek Tucker
Sandia National Laboratories

