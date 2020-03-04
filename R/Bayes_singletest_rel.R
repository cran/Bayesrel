#'
#' calculate single test reliability estimates
#' @description calculate Bayesian and frequentist single test reliability measures.
#' Reported are Bayesian credible intervals (HDI) and frequentist confidence intervals (non parametric or parametric bootstrap).
#' The estimates supported are Cronbach alpha, lambda2/4/6, the glb, and Mcdonald omega. Beware of lambda4 with many indicators,
#' the computational effort is considerable
#'
#' @param x A dataset or covariance matrix
#' @param estimates A character vector containing the estimands, we recommend using lambda4 with only a few items due to the computation time
#' @param interval A number specifying the uncertainty interval
#' @param n.iter A number for the iterations of the Gibbs Sampler
#' @param n.burnin A number for the burnin in the Gibbs Sampler
#' @param n.boot A number for the bootstrap samples
#' @param omega.freq.method A character string for the method of frequentist omega, either "pfa" or "cfa"
#' @param n.obs A number for the sample observations when a covariance matrix is supplied and the factor model is calculated
#' @param alpha.int.analytic A logical for calculating the alpha confidence interval analytically
#' @param freq A logical for calculating the frequentist estimates
#' @param Bayes A logical for calculating the Bayesian estimates
#' @param para.boot A logical for calculating the parametric bootstrap, the default is the non-parametric
#' @param item.dropped A logical for calculating the if-item-dropped statistics
#' @param missing A string specifying the way to handle missing data
#'
#' @examples
#' summary(strel(asrm, estimates = "lambda2"))
#' summary(strel(asrm, estimates = "lambda2", item.dropped = TRUE))
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
               interval = .95, n.iter = 1e3, n.burnin = 50, n.boot = 1000,
               omega.freq.method = "cfa",
               n.obs = NULL, alpha.int.analytic = FALSE,
               freq = TRUE, Bayes = TRUE,
               para.boot = FALSE,
               item.dropped = FALSE,
               missing = "listwise") {

  default <- c("alpha", "lambda2", "lambda4", "lambda6", "glb", "omega")
  # estimates <- match.arg(arg = estimates, several.ok = T)
  mat <- match(default, estimates)
  estimates <- estimates[mat]
  estimates <- estimates[!is.na(estimates)]
  p <- NULL
  sum_res <- list()
  sum_res$call <- match.call()

  pairwise <- FALSE
  if (any(is.na(x))) {
    if (missing == "listwise") {
      pos <- which(is.na(x), arr.ind = T)[, 1]
      x <- x[-pos, ]
      ncomp <- nrow(x)
      sum_res$complete <- ncomp
    }
    else if (missing == "pairwise") {
      pairwise = T
      sum_res$miss_pairwise <- T
    } else return("missing values in data detected, please remove and run again")
  }
  data <- NULL
  sigma <- NULL
  if (ncol(x) == nrow(x)){
    if (is.null(n.obs) & "omega" %in% estimates) {
      return("number of observations (n.obs) needs to be specified when entering a covariance matrix")}
    if (sum(x[lower.tri(x)] != t(x)[lower.tri(x)]) > 0) {return("input matrix is not symmetric")}
    if (sum(eigen(x)$values < 0) > 0) {return("input matrix is not positive definite")}
    sigma <- x
    data <- MASS::mvrnorm(n.obs, rep(0, ncol(sigma)), sigma, empirical = TRUE)
  } else{
    data <- scale(x, scale = F)
    # sigma <- cov(data)
  }

  if (omega.freq.method != "cfa" & omega.freq.method != "pfa") {
    return("enter a valid omega method, either 'cfa' or 'pfa'")}

  if (Bayes) {
    sum_res$Bayes <- gibbsFun(data, n.iter, n.burnin, estimates, interval, item.dropped, pairwise)
  }


  if(freq){
    if (para.boot){
      sum_res$freq <- freqFun_para(data, n.boot, estimates, interval, omega.freq.method, item.dropped,
                                   alpha.int.analytic, pairwise)
    } else{
      sum_res$freq <- freqFun_nonpara(data, n.boot, estimates, interval, omega.freq.method, item.dropped,
                                    alpha.int.analytic, pairwise)
    }
    sum_res$omega.freq.method <- omega.freq.method
  }


  sum_res$estimates <- estimates
  sum_res$n.iter <- n.iter
  sum_res$n.burnin <- n.burnin
  sum_res$interval <- interval
  sum_res$data <- data
  sum_res$n.boot <- n.boot

  class(sum_res) = 'strel'
  return(sum_res)
}
