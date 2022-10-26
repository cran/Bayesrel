

<!-- badges: start -->
[![R build status](https://github.com/juliuspf/Bayesrel/workflows/R-CMD-check/badge.svg)](https://github.com/juliuspf/Bayesrel/actions)
<!-- badges: end -->


## Installation

You can install the released version of Bayesrel from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Bayesrel")
```
or install the latest version of Bayesrel from [github] (https://github.com) with the help of the remotes-package:

```r
remotes::install_github("juliuspf/Bayesrel")
```

## Example

### Unidimensional data
This is a basic example which shows you how to compute alpha, lambda2, the glb, and omega for an example real data set:

``` r
library(Bayesrel)
## basic example code
## load example data set from the package
## run the main reliability function
res <- strel(asrm)
## get a full result output
summary(strel)
## return the probability that coefficient alpha is larger than .70
pStrel(res, estimate = "alpha", low.bound = .70)

## get the posterior median of, e.g., alpha instead of the mean:
median(res$Bayes$samp$Bayes_alpha)
```

### Multidimensional data
This is a basic example which shows you how to compute omega_t and omega_h for an example real data set. 
The data follow a second-order factor model with no crossloadings (required):

``` r
library(Bayesrel)
## basic example code
## run the Bayesian omegas, specify 5 group factors
res <- bomegas(upps, n.factors = 5, missing = "listwise")
## get a full result output
summary(res)
## return the probability that coefficient omega_t is larger than .70
pOmegas(res, cutoff.t = .70)
## plot posterior predictive check for the higher-order (second-order) factor model
secoFit(res, upps)
```

In the example above we implicitly assumed that the items of the data set were ordered
so that, with 5 group factors, the first four items load on the first factor, 
items 5-8 load on the second factor and so on. When the data is not organized this way and/or the items 
cannot be distributed among the factors evenly, one can specify a model syntax relating the items 
to the group factors in lavaan style. The item names need to equal the variable names in the data:

``` r
model <- "
  f1 =~ U17_r + U22_r + U29_r + U34_r
  f2 =~ U4 + U14 + U19 + U27
  f3 =~ U6 + U16 + U28 + U48
  f4 =~ U23_r + U31_r + U36_r + U46_r
  f5 =~ U10_r + U20_r + U35_r + U52_r
  "
```

The reliability is then estimated as follows: 

```r
res <- bomegas(upps, n.factors = 5, model = model, missing = "listwise")
```
