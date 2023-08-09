
#' Estimate Bayesian reliability estimates for multidimensional scales following common factor models
#' @description When supplying a multidimensional data set
#' the function estimates the reliability of the set by means of omega_total
#' and the general factor saturation of the set by means of omega_hierarchical.
#' The data may follow multiple multi-factorial factor models:
#' a) the second-order factor model (crossloadings possible),
#' b) the bi-factor model (no crossloadings),
#' c) the correlated factor model (crossloadings possible, only omega_t)
#' Error-covariances are not estimable.
#'
#' The prior distributions of omega_t and omega_h are computed from
#' the prior distributions of the respective factor model parameters.
#' Specifically, normal distributions for the factor loadings and factor scores;
#' an inverse gamma distribution for the manifest and latent residuals;
#' an inverse Wishart distribution for the covariance matrix of the latent variables (correlated model).
#' A Gibbs sampler iteratively draws samples from the conditional posterior distributions
#' of the factor model parameters. The posterior distributions of omega_t
#' and omega_h are computed from the posterior samples of the factor model parameters.
#'
#' The output contains the posterior distributions of omega_t and omega_h (only for second-order and bi-factor model),
#' their mean, and credible intervals. If desired, one may also find the posterior implied covariance matrices,
#' and the posterior factor model parameters in the output.
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
#' @param model.type A string indicating the factor model estimated, by default this is the "second-order" model.
#' Another option is the "bi-factor" model; and eventually "correlated" for the correlated factor model
#' @param n.iter A number for the iterations of the Gibbs Sampler
#' @param n.burnin A number for the burnin in the Gibbs Sampler
#' @param n.chains A number for the chains to run for the MCMC sampling
#' @param thin A number for the thinning of the MCMC samples
#' @param interval A number specifying the credible interval,
#' the interval is the highest posterior density interval (HPD)
#' @param missing A string denoting the missing data handling, can be "impute" or "listwise".
#' With impute the missing data will be estimated during the MCMC sampling
#' as further unknown parameters
#' @param a0 A number for the shape of the prior inverse gamma distribution for the manifest residual variances,
#' when left as NA, the default value is set to 2
#' @param b0 A number for the scale of the prior inverse gamma distribution for the manifest residual variances,
#' when left as NA, the default value is set to 1
#' @param l0 A number for the mean of the prior normal distribution for the manifest loadings,
#' when left as NA, the default value is set to 0; can be a single value or a loading matrix
#' @param A0 A number for scaling the variance of the prior normal distribution for the manifest loadings,
#' when left as NA, the default value is set to 1
#' @param c0 A number for the shape of the prior inverse gamma distribution for the latent residual variances,
#' when left as NA, the default value is set to 2; necessary only for the second-order model
#' @param d0 A number for the scale of the prior inverse gamma distribution for the latent residual variances,
#' when left as NA, the default value is set to 1, necessary only for the second-order model
#' @param beta0 A number for the mean of the prior normal distribution for the g-factor loadings,
#' when left as NA, the default value is set to 0, can be a single value or a vector
#' @param B0 A number for scaling the variance of the prior normal distribution for the g-factor loadings,
#' when left as NA, the default value is set to 2.5
#' @param p0 A number for the shape of the prior inverse gamma distribution for the variance of the g-factor,
#' when left as NA, the default value is set to q^2-q when q are the number of group factors for the second-order
#' and bi-factor model
#' @param R0 A number for the scale of the prior inverse gamma distribution for the variance of the g-factor,
#' when left as NA, the default value is set to the number of items for the second-order
#' and bi-factor model
#' @param param.out A logical indicating if loadings and residual variances should be attached to the result,
#' by default FALSE because it saves memory
#' @param callback An empty function for implementing a progressbar call
#' from a higher program (e.g., JASP)
#' @param disableMcmcCheck A logical disabling the MCMC settings check for running tests and the likes
#'
#'
#' @return The posterior means and the highest posterior density intervals for
#' omega_t and omega_h (for the second-order and bi-factor model)
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
#' n.chains = 2, missing = "listwise", model.type = "second-order")
#'

#'
#' @references{
#'
#' \insertRef{Lee2007}{Bayesrel}
#'
#'}
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats cov2cor
#'
#' @export
bomegas <- function(
  data,
  n.factors = NULL,
  model = NULL,
  model.type = "second-order",
  n.iter = 2e3,
  n.burnin = 200,
  n.chains = 3,
  thin = 1,
  interval = .95,
  missing = "impute",
  a0 = NA,
  b0 = NA,
  l0 = NA,
  A0 = NA,
  c0 = NA,
  d0 = NA,
  beta0 = NA,
  B0 = NA,
  p0 = NA,
  R0 = NA,
  param.out = FALSE,
  callback = function(){},
  disableMcmcCheck = FALSE
) {

  # check mcmc settings
  if (!disableMcmcCheck) {
    checkMcmcSettings(n.iter, n.burnin, n.chains, thin)
  }

  # make sure only the data referenced in the model file is used, if a model file is specified
  if (!is.null(model)) {
    mod_opts <- modelSyntaxExtract(model, colnames(data))
    data <- data[, colnames(data) %in% mod_opts$var_names]
    n.factors <- mod_opts$mod_n.factors
  } else {
    if (is.null(n.factors)) {
      stop("The number of factors needs to be specified when no model file is provided.")
    }
  }

  # check model.type string
  if (!(model.type %in% c("bi-factor", "bifactor", "second-order", "correlated", "secondorder", "hierarchical",
                          "secondOrder", "biFactor"))) {
    stop("model.type invalid; needs to be 'bi-factor', 'second-order', or 'correlated'")
  }

  if (model.type %in% c("bifactor", "biFactor")) model.type <- "bi-factor"
  if (model.type %in% c("secondorder", "hierarchical", "secondOrder")) model.type <- "second-order"


  # deal with missings
  listwise <- FALSE
  impute <- FALSE
  complete_cases <- nrow(data)
  if (anyNA(data)) {
    if (missing == "listwise") {
      pos <- which(is.na(data), arr.ind = TRUE)[, 1]
      data <- data[-pos, ]
      ncomp <- nrow(data)
      complete_cases <- ncomp
      listwise <- TRUE
    } else { # missing is impute
      impute <- TRUE
    }
  }

  data <- scale(data, scale = FALSE)

  # handle the prior hyperparameters
  if (model.type != "correlated") {
    defaults <- list(a0 = 2, b0 = 1, c0 = 2, d0 = 1, l0 = 0, A0 = 1, beta0 = 0, B0 = 2.5,
                     p0 = n.factors^2 - n.factors, R0 = ncol(data))
  } else {
    defaults <- list(a0 = 2, b0 = 1, c0 = 2, d0 = 1, l0 = 0, A0 = 1, beta0 = 0, B0 = 2.5,
                     p0 = n.factors + 2, R0 = matrix(n.factors, n.factors, n.factors))
    diag(defaults[["R0"]]) <- ncol(data)

  }

  set <- list(a0, c0, b0, d0, l0, A0, beta0, B0, p0, R0)
  prior.params <- ifelse(is.na(set), defaults, set)
  names(prior.params) <- names(defaults)

  pb <- txtProgressBar(max = n.chains * n.iter, style = 3)

  sum_res <- bomegasMultiOut(data, n.factors, n.iter, n.burnin, thin, n.chains,
                             interval, model, impute, prior.params,
                             param.out, callback, pbtick = pb, model.type = model.type)
  close(pb)


  sum_res$complete_cases <- complete_cases
  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$impute <- impute
  sum_res$listwise <- listwise
  sum_res$interval <- interval
  sum_res$prior_params <- prior.params
  sum_res$model.type <- model.type

  class(sum_res) <- "bomegas"

  return(sum_res)
}
