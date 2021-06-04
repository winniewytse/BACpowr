
# hcbr

<!-- badges: start -->
<!-- badges: end -->

This package implements the hyrbid classical-Bayesian (HCB) approach to determine the
required number of clusters and cluster size in multilevel randomized trials.
The current package supports cluster randomized trials, while it will be developed
to support different types of multilevel trials. 

## Installation

You can install `hcbr` using the following command:

``` r
devtools::install_github("winniewytse/hcbr")
```

## Example

Suppose that from a pilot or published study, you found that the treatment effect size
is 0.5 with a standard error of 0.2 and the intraclass correlation is 0.1 with a 
standard error of 0.05. In the subsequent study, you expect to recruit 30 clusters (`J`) and would like to determine how large each cluster (`n`) should be. Then you can run the 
following to determine the minimum required cluster size using the HCB approach:

``` r
library(hcbr)
crtJn(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30)
```

If instead you expect the cluster size to be 50 (`n`) and hope to determine the minimum
required number of cluster (`J`), you can run the following:

``` r
crtJn(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 50)
```



