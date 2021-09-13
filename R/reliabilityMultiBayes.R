
#' Estimate reliability estimates for multidimensional scales in the Bayesian framework
#' @description When supplying a data set that is multidimensional
#' the function estimates the reliability of the set by means of omega-total
#' and the general factor saturation of the set by means of omega-hierarchical
#' In the process a higher-order factor model is estimated in the Bayesian framework,
#' and posterior distributions of omega-t and omega-h are obtained from the
#' posterior distributions of the factor model parameters
#' The output contains the posterior distributions of omega-t and omega-h,
#' their mean and credible intervals.
#' @param data A matrix or data.frame containing multivariate observations,
#' rows = observations, columns = variables/items
#' @param n.factors A number specifying the number of group factors that the items load on
#' @param model A string that by default ("balanced") distributes the items evenly
#' among the number of group factors. This only works if the items are a multiple of
#' the number of group factors and the items are already grouped in the data set,
#' meaning, e.g., items 1-5 load on one factor, 6-10 on another, and so on.
#' A model file can be specified in lavaan syntax style (f1=~.+.+.) to relate the items
#' to the group factors. The items can either be named as the columns in the data set
#' or x1, ..., xn, where 1,...,n correspond to the column numbers
#' @param n.iter A number for the iterations of the Gibbs Sampler
#' @param n.burnin A number for the burnin in the Gibbs Sampler
#' @param n.chains A number for the chains to run for the MCMC sampling
#' @param thin A number for the thinning of the MCMC samples
#' @param interval A number specifying the credible interval,
#' the interval is the highest posterior desntiy interval (HPD)
#' @param missing A string denoting the missing data handling, can be "pairwise" or "listwise.
#' With pairwise the missing data will be "imputed" during the MCMC sampling
#' as further unknown parameters
#' @param callback An empty function for implementing a progressbar call
#' from a higher program (e.g., JASP)
#'
#' @examples
#' # note that the iterations are set very low for smoother running examples, you should use
#' # at least the defaults
#' res <- bomegas(upps, n.factors = 5, model = "balanced", n.iter = 200, n.burnin = 50,
#' n.chains = 2, missing = "listwise")
#'
#' # or with specified model syntax relating the group factors to the items:
#' model <- "f1 =~ U17_r + U22_r + U29_r + U34_r
#' f2 =~ U4 + U14 + U19 + U27
#' f3 =~ U6 + U16 + U28 + U48
#' f4 =~ U23_r + U31_r + U36_r + U46_r
#' f5 =~ U10_r + U20_r + U35_r + U52_r"
#' res <- bomegas(upps, n.factors = 5, model = model, n.iter = 200, n.burnin = 50,
#' n.chains = 2, missing = "listwise")
#'
#' @references{
#'
#' \insertRef{Lee2007}{Bayesrel}
#'
#'}
#' @export
bomegas <- function(
  data,
  n.factors,
  model = "balanced",
  n.iter = 2e3,
  n.burnin = 200,
  n.chains = 3,
  thin = 1,
  interval = .95,
  missing = "pairwise",
  callback = function(){}
) {

  listwise <- FALSE
  pairwise <- FALSE
  complete_cases <- nrow(data)
  if (any(is.na(data))) {
    if (missing == "listwise") {
      pos <- which(is.na(data), arr.ind = TRUE)[, 1]
      data <- data[-pos, ]
      ncomp <- nrow(data)
      complete_cases <- ncomp
      listwise <- TRUE
    } else { # missing is pairwise
      pairwise <- TRUE
    }
  }

  data <- scale(data, scale = FALSE)

  sum_res <- bomegasMultiOut(data, n.factors, n.iter, n.burnin, thin, n.chains,
                             interval, model, pairwise, callback)

  sum_res$complete_cases <- complete_cases
  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$pairwise <- pairwise
  sum_res$listwise <- listwise
  sum_res$interval <- interval
  class(sum_res) <- "bomegas"

  return(sum_res)
}
