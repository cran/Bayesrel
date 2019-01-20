#' probability of estimate being bigger than threshold
#' @description
#' takes a mcmc posterior sample of any of the single test reliability estimates
#' and calculates any given probability of the estimate being bigger
#' or smaller than an arbitrary value
#'
#' @param x A strel output object (list)
#' @param estimate A character string indicating what estimate to plot from the strel output object
#' @param low.bound A number for the threshold to be tested against
#'
#' @examples
#' pstrel(strel(cavalini, "lambda2"), "lambda2", .80)
#' @export


pstrel <- function(x, estimate, low.bound){
  posi <- grep(estimate, x$estimates, ignore.case = T)
  samp <- coda::as.mcmc(unlist(x$bay$samp[posi]))
  obj <- ecdf(samp)
  est.prob <- 1 - obj(low.bound)
  return(est.prob)
}

# probic <- function(x, low.bound){
#   obj <- ecdf(x)
#   est.prob <- 1 - obj(low.bound)
#   return(est.prob)
# }
