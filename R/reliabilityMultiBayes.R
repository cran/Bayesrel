
#' Estimate Bayesian reliability estimates for multidimensional scales following a second-order factor model
#' @description When supplying a multidimensional data set
#' the function estimates the reliability of the set by means of omega_total
#' and the general factor saturation of the set by means of omega_hierarchical.
#' The data must follow a simple second-order factor model structure with no crossloadings
#' and no error covariances. Otherwise the estimation will fail.
#'
#' The prior distributions of omega_t and omega_h are computed from
#' the prior distributions of the second-order factor model.
#' Specifically, a multivariate normal distribution for the group factor loadings and
#' the factor scores; a normal distribution for the general factor loadings;
#' an inverse gamma distribution for the manifest and latent residuals;
#' an inverse Wishart distribution for the covariance matrix of the latent variables.
#' A Gibbs sampler iteratively draws samples from the conditional posterior distributions
#' of the second-order factor model parameters. The posterior distributions of omega_t
#' and omega_h are computed from the posterior samples of the factor model parameters.
#'
#' The output contains the posterior distributions of omega_t and omega_h,
#' their mean and credible intervals.
#' @param data A matrix or data.frame containing multivariate observations,
#' rows = observations, columns = variables/items
#' @param n.factors A number specifying the number of group factors that the items load on
#' @param model A string that by default NULL (=balanced) distributes the items evenly
#' among the number of group factors. This only works if the items are a multiple of
#' the number of group factors and the items are already grouped in the data set,
#' meaning, e.g., items 1-5 load on one factor, 6-10 on another, and so on.
#' A model file can be specified in lavaan syntax style (f1=~.+.+.) to relate the items
#' to the group factors. The items' names need to equal the column names in the data set,
#' aka the variable names
#' @param n.iter A number for the iterations of the Gibbs Sampler
#' @param n.burnin A number for the burnin in the Gibbs Sampler
#' @param n.chains A number for the chains to run for the MCMC sampling
#' @param thin A number for the thinning of the MCMC samples
#' @param interval A number specifying the credible interval,
#' the interval is the highest posterior desntiy interval (HPD)
#' @param missing A string denoting the missing data handling, can be "impute" or "listwise".
#' With impute the missing data will be estimated during the MCMC sampling
#' as further unknown parameters
#' @param a0 A number for the shape of the prior inverse gamma distribution for the manifest residual variances,
#' by default 2
#' @param b0 A number for the scale of the prior inverse gamma distribution for the manifest residual variances,
#' by default 1
#' @param l0 A number for the mean of the prior normal distribution for the manifest loadings,
#' by default 0, can be a single value or a loading matrix
#' @param A0 A number for scaling the variance of the prior normal distribution for the manifest loadings,
#' by default 1
#' @param c0 A number for the shape of the prior inverse gamma distribution for the latent residual variances,
#' by default 2
#' @param d0 A number for the scale of the prior inverse gamma distribution for the latent residual variances,
#' by default 1
#' @param beta0 A number for the mean of the prior normal distribution for the latent loadings,
#' by default 0, can be a single value or a vector
#' @param B0 A number for scaling the variance of the prior normal distribution for the latent loadings,
#' by default 1
#' @param p0 A number for the shape of the prior inverse gamma distribution for the variance of the g-factor,
#' by default set to q^2-q when q are the number of group factors
#' @param R0 A number for the scale of the prior inverse gamma distribution for the variance of the g-factor,
#' by default set to the number of items
#' @param param.out A logical indicating if loadings and residual variances should be attached to the result,
#' by default FALSE because it saves memory
#' @param callback An empty function for implementing a progressbar call
#' from a higher program (e.g., JASP)
#'
#' @return The posterior means and the highest posterior density intervals for
#' omega_t and omega_h
#'
#' @examples
#' # specify the model syntax relating the group factors to the items:
#' model <- "f1 =~ U17_r + U22_r + U29_r + U34_r
#' f2 =~ U4 + U14 + U19 + U27
#' f3 =~ U6 + U16 + U28 + U48
#' f4 =~ U23_r + U31_r + U36_r + U46_r
#' f5 =~ U10_r + U20_r + U35_r + U52_r"
#' # specifying the model can be omitted if the model structure is simple and balanced,
#' # meaning the an equal number of items load on each factor, and the items relate to the factors
#' # in the same order as they appear in the data.
#' # Note that the iterations are set very low for smoother running examples, you should use
#' # at least the defaults:
#' res <- bomegas(upps, n.factors = 5, model = NULL, n.iter = 200, n.burnin = 50,
#' n.chains = 2, missing = "listwise")
#'

#'
#' @references{
#'
#' \insertRef{Lee2007}{Bayesrel}
#'
#'}
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
bomegas <- function(
  data,
  n.factors,
  model = NULL,
  n.iter = 2e3,
  n.burnin = 200,
  n.chains = 3,
  thin = 1,
  interval = .95,
  missing = "impute",
  a0 = 2, b0 = 1,
  l0 = 0, A0 = 1,
  c0 = 2, d0 = 1,
  beta0 = 0, B0 = 2.5,
  p0 = NULL, R0 = NULL,
  param.out = FALSE,
  callback = function(){}
) {

  # make sure only the data referenced in the model file is used, if a model file is specified
  if (!is.null(model)) {
    mod_opts <- modelSyntaxExtract(model, colnames(data))
    data <- data[, mod_opts$var_names]
  }

  # deal with missings
  listwise <- FALSE
  pairwise <- FALSE
  complete_cases <- nrow(data)
  if (anyNA(data)) {
    if (missing == "listwise") {
      pos <- which(is.na(data), arr.ind = TRUE)[, 1]
      data <- data[-pos, ]
      ncomp <- nrow(data)
      complete_cases <- ncomp
      listwise <- TRUE
    } else { # missing is impute
      pairwise <- TRUE
    }
  }

  data <- scale(data, scale = FALSE)

  if (is.null(p0)) p0 <- n.factors^2 - n.factors
  if (is.null(R0)) R0 <- ncol(data)

  pb <- txtProgressBar(max = n.chains * n.iter, style = 3)
  sum_res <- bomegasMultiOut(data, n.factors, n.iter, n.burnin, thin, n.chains,
                             interval, model, pairwise, a0, b0, l0, A0, c0, d0, beta0, B0, p0, R0,
                             param.out, callback, pbtick = pb)
  close(pb)

  sum_res$complete_cases <- complete_cases
  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$pairwise <- pairwise
  sum_res$listwise <- listwise
  sum_res$interval <- interval
  sum_res$prior_params <- list(a0 = a0, b0 = b0, l0 = l0, A0 = A0, c0 = c0, d0 = d0,
                               beta0 = beta0, B0 = B0, p0 = p0, R0 = R0)
  class(sum_res) <- "bomegas"

  return(sum_res)
}
