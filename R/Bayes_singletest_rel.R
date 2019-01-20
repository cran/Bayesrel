#'
#' calculate single test reliability estimates
#' @description calculate Bayesian and frequentist single test reliability measures.
#' Reported are Bayesian credible intervals (HDI) and frequentist confidence intervals (non parametric or parametric bootstrap).
#' The estimates supported are Cronbach alpha, lambda2/4/6, the glb, and Mcdonald omega.
#'
#' @param x A dataset or covariance matrix
#' @param estimates A character vector containing the estimands, we recommend using lambda4 with only a few items due to the computation time
#' @param interval A number specifying the uncertainty interval
#' @param n.iter A number for the iterations of the Gibbs Sampler
#' @param n.burnin A number for the burnin in the Gibbs Sampler
#' @param boot.n A number for the bootstrap samples
#' @param omega.freq.method A character string for the method of frequentist omega, either pfa or cfa
#' @param omega.fit A logical for calculating the fit of the single factor model
#' @param n.obs A number for the sample observations when a covariance matrix is supplied and the factor model is calculated
#' @param alpha.int.analytic A logical for calculating the alpha confidence interval analytically
#' @param bayes A logical for calculating the Bayesian estimates
#' @param freq A logical for calculating the frequentist estimates
#' @param para.boot A logical for calculating the parametric bootstrap, the default is the non-parametric
#' @param prior.samp A logical for calculating the prior distributions (necessary for plot functions)
#' @param item.dropped A logical for calculating the if-item-dropped statistics
#'
#' @examples
#' summary(strel(cavalini, estimates = "lambda2"))
#' summary(strel(cavalini, estimates = "lambda2", item.dropped = TRUE))
#'
#'
#' @references{
#'   \insertRef{murphy2007}{Bayesrel}
#'   \insertRef{lee2007}{Bayesrel}
#' }
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom graphics arrows axis legend lines par plot text title
#' @importFrom methods is
#' @importFrom stats cov cov2cor density ecdf median qnorm quantile rchisq rgamma rnorm runif sd var
#' @importFrom Rdpack reprompt
#'
#' @export
strel <- function(x, estimates = c("alpha", "lambda2", "glb", "omega"),
               interval = .95, n.iter = 2e3, n.burnin = 50, boot.n = 1000,
               omega.freq.method = "cfa", omega.fit = FALSE,
               n.obs = NULL, alpha.int.analytic = FALSE,
               bayes = TRUE, freq = TRUE, para.boot = FALSE, prior.samp = FALSE,
               item.dropped = FALSE) {

  estimates <- match.arg(estimates, several.ok = T)
  default <- c("alpha", "lambda2", "lambda4", "lambda6", "glb", "omega")
  mat <- match(default, estimates)
  estimates <- estimates[mat]
  estimates <- estimates[!is.na(estimates)]
  p <- NULL
  sum.res <- list()
  sum.res$call <- match.call()

  if (sum(is.na(x)) > 0) {
    return("missing values in data detected, please remove and run again")
  }
  data <- NULL
  sigma <- NULL
  if (ncol(x) == nrow(x)){
    if (is.null(n.obs)) {return("number of observations needs to be specified when entering a covariance matrix")}
    if (sum(x[lower.tri(x)] != t(x)[lower.tri(x)]) > 0) {return("input matrix is not symmetric")}
    if (sum(eigen(x)$values < 0) > 0) {return("input matrix is not positive definite")}
    sigma <- x
    data <- MASS::mvrnorm(n.obs, rep(0, ncol(sigma)), sigma, empirical = TRUE)
  } else{
    data <- scale(x, scale = F)
    sigma <- cov(data)
  }

  # if("glb" %in% estimates){
  #   control <- Rcsdp::csdp.control(printlevel = 0)
  #   write.control.file(control)
  # }
  if (bayes){
    sum.res$bay <- gibbsFun(data, n.iter, n.burnin, estimates, interval, item.dropped, omega.fit)
  }
  if (omega.fit) {omega.freq.method <- "cfa"}
  if(freq){
    if (para.boot){
      sum.res$freq <- freqFun_para(data, boot.n, estimates, interval, omega.freq.method, item.dropped,
                                   alpha.int.analytic, omega.fit)
    } else{
      sum.res$freq <- freqFun_nonpara(data, boot.n, estimates, interval, omega.freq.method, item.dropped,
                                    alpha.int.analytic, omega.fit)
    }
    sum.res$omega.freq.method <- omega.freq.method
  }

  # if("glb" %in% estimates)
  #   unlink("param.csdp")

  if (prior.samp) {
    sum.res$priors <- priorSamp(ncol(data), estimates)
  }

  sum.res$estimates <- estimates
  sum.res$n.iter <- n.iter
  sum.res$n.burnin <- n.burnin
  sum.res$interval <- interval
  sum.res$n.item <- ncol(data)

  class(sum.res) = 'strel'
  return(sum.res)
}
