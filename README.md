# Bayesrel
R-Package "Bayesrel" provides both Bayesian and frequentist internal consistency estimates

So far, it provides the most common single test reliability, being: 
    Coefficient Alpha, Guttman's lambda-2/-4/-6, greatest lower bound and Mcdonald's Omega. 
    The Bayesian estimates are provided with credible intervals. 
    The method for the Bayesian estimates is sampling from the posterior inverse wishart for the covariance matrix based measures.
    Gibbs Sampling from the joint conditional distributions of a single factor model in the case of omega.
    The Bayesian posterior mcmc objects are provided, which allow for further analysis. 
    The frequentist estimates are provided with bootstrapped confidence intervals. The user can choose between non-parametric or parametric bootstrap. 
    For Alpha the interval can also be analytic. 
    The frequentist Omega can be calculated using a PFA or a CFA. 
    A graphical predictive posterior check can be done for the single factor model. The often used if-item-dropped statistics can be calculated, too. 
    The package also allows for the calculation of the probability of an estimator being bigger than a given threshold.
    There is also the functionality to plot the posterior against the prior for a chosen estimator with cutoffs. 
    In addition to that a plot for the multiple posteriors for the if-item-dropped cases can be done.
    The PFA method for omega is from Aaron Schlegel (url:https://www.r-bloggers.com/iterated-principal-factor-method-of-factor-analysis-with-r/). 
    The glb.algebraic method is an altered function from Moltner (https://www.rdocumentation.org/packages/psych/versions/1.8.10/topics/glb.algebraic)
    The Lambda4 method is from Tom Benton (doi:10.1007/978-3-319-07503-7_19)
   
