# mvBayes

An R implementation of the multivariate Bayesian regression (mvBayes) framework. Decomposes a multivariate/functional response using a user-specified orthogonal basis decomposition, and then models each basis component independently using an arbitrary user-specified (univariate) Bayesian regression model. Includes prediction and plotting methods.


## Examples
* [Friedman Example](inst/friedman_demo.R) - An extension of the "Friedman function" to functional response. The Bayesian regression model here is BASS (Bayesian Adaptive Smoothing Splines, see https://github.com/lanl/BASS)


### Installation
------------------------------------------------------------------------------
1. Download zip or tar.gz of package or clone repository
2. Install into R (> 4.1.0)

> `install.packages("mvBayes.tar.gz", repos = NULL)`

------------------------------------------------------------------------------

## References


************

Author: Gavin Q. Collins and J. Derek Tucker
Sandia National Laboratories

