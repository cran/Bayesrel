# Bayesrel
R-Package "Bayesrel" provides both Bayesian and frequentist single-test reliability estimates

Functionality for the most common single test reliability estimates is provided: 
    Coefficient alpha, Guttman's lambda-2/-4/-6, greatest lower bound and coefficient omega. 
    The Bayesian estimates are provided with credible intervals. 
    The method for the Bayesian estimates is sampling from the posterior inverse Wishart for all coefficients but omega
    (Murphy, 2007, https://www.seas.harvard.edu/courses/cs281/papers/murphy-2007.pdf).
    Gibbs Sampling from the joint conditional distributions of a single factor model in the case of omega
    (Lee, 2007, https://dx.doi.org/10.1002/9780470024737)
    The Bayesian posterior mcmc samples are provided, which allow for further analysis. 
    The frequentist estimates are provided with bootstrapped confidence intervals. 
    The user can choose between non-parametric or parametric bootstrap. 
    For alpha the interval can also be analytic (Bonett & Wright, 2014, https://dx.doi.org/10.1002/job.1960)
    The frequentist omega can be calculated using a PFA or a CFA. 
    A graphical predictive posterior check can be done for the fit of the single factor model for the Bayesian approach. 
    The CFA fit statistics for the single-factor model are provided for the frequentist approach.
    The often used if-item-dropped statistics can be calculated, too. 
    The Bayesian version of that also offers a plot of the posterior densities when an item would be dropped.
    The package also allows for the calculation of the probability of an estimator being bigger than a given threshold, 
    both for the prior and posterior.
    There is also the functionality to plot the posterior with the prior for a chosen estimator with cutoffs. 
    The PFA method for omega is from Aaron Schlegel (https://www.r-bloggers.com/iterated-principal-factor-method-of-factor-analysis-with-r/). 
    The glb method is adjusted code from the Rcsdp package by Hector Corrada Bravo (https://CRAN.R-project.org/package=Rcsdp) 
    The Lambda4 method is from Tom Benton (https://dx.doi.org/10.1007/978-3-319-07503-7_19). 
    Missing data can be handled both listwise and pairwise for both the Bayesian and frequentist cases.
    
   
   <!-- badges: start -->
[![Travis build status](https://travis-ci.org/juliuspf/Bayesrel.svg?branch=master)](https://travis-ci.org/juliuspf/Bayesrel)
<!-- badges: end -->
