#' prior and posterior probability of estimate being bigger than threshold
#' @description
#' takes a mcmc posterior sample of any of the single test reliability estimates
#' and calculates the prior and posterior probability of the estimate being bigger
#' or smaller than an arbitrary value (priors are stored in the package)
#'
#' @param x A strel output object (list)
#' @param estimate A character string indicating what estimate to plot from the strel output object
#' @param low.bound A number for the threshold to be tested against
#'
#' @examples
#' p_strel(strel(asrm, "lambda2"), "lambda2", .80)
#' @export


p_strel <- function(x, estimate, low.bound){
  posi1 <- grep(estimate, x$estimates, ignore.case = T)
  samp <- unlist(x$Bayes$samp[posi1])
  obj <- ecdf(samp)
  post_prob <- 1 - obj(low.bound)

  # prior prob
  n.item <- dim(x$Bayes$covsamp)[2]
  prior_all <- priors[[as.character(n.item)]]
  posi2 <- grep(estimate, prior_all, ignore.case = T)
  prior <- prior_all[[posi2]]
  end <- length(prior[["x"]])
  poslow <- end - sum(prior[["x"]] > low.bound)
  prior_prob <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]])
  out <- c(prior_prob, post_prob)
  names(out) <- c("prior_prob", "posterior_prob")
  return(out)
}

# probic <- function(x, low.bound){
#   obj <- ecdf(x)
#   est_prob <- 1 - obj(low.bound)
#   return(est_prob)
# }
