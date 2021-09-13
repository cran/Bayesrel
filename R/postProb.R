#' prior and posterior probability of estimate being bigger than threshold
#' @description
#' takes a mcmc posterior sample of any of the single test reliability estimates
#' and calculates the prior and posterior probability of the estimate being bigger
#' or smaller than an arbitrary value (priors are stored in the package)
#'
#' @param x A strel output object (list)
#' @param estimate A character string indicating what estimate to plot from the strel output object
#' @param low.bound A number for the threshold to be tested against
#' @examples
#' p_strel(strel(asrm, "lambda2", n.chains = 2, n.iter = 100, freq = FALSE), "lambda2", .80)
#' @export
p_strel <- function(x, estimate, low.bound) {

  posi1 <- grep(estimate, x$estimates, ignore.case = TRUE)
  samp <- as.vector(x$Bayes$samp[[posi1]])
  obj <- ecdf(samp)
  post_prob <- 1 - obj(low.bound)

  # prior prob
  n.item <- dim(x$Bayes$covsamp)[3]
  if (n.item > 50) {
    prior <- density(unlist(priorSampUni(n.item, estimate)), from = 0, to = 1, n = 512)
  } else {
    prior_all <- priors[[as.character(n.item)]]
    posi2 <- grep(estimate, prior_all, ignore.case = TRUE)
    prior <- prior_all[[posi2]]
  }
  end <- length(prior[["x"]])
  poslow <- end - sum(prior[["x"]] > low.bound)
  prior_prob <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]])
  out <- c(prior_prob, post_prob)
  names(out) <- c("prior_prob", "posterior_prob")
  return(out)
}


#' prior and posterior probability of omega-t and omega-h being bigger than thresholds
#' @description
#' takes mcmc posterior samples of omega-t and omega-h
#' and calculates the prior and posterior probability of the estimate being bigger
#' or smaller than an arbitrary value
#'
#' @param x A strel output object (list)
#' @param cutoff.t A number indicating the threshold for omega-t
#' @param cutoff.h A number indicating the threshold for omega-h
#'
#' @examples
#' p_omegas(bomegas(upps, n.factors = 5, n.chains = 2, n.iter = 100, n.burnin = 50,
#' missing = "listwise"))
#' @export
p_omegas <- function(x, cutoff.t = .80, cutoff.h = .60) {

  sampt <- as.vector(x$omega_t$chains)
  samph <- as.vector(x$omega_h$chains)
  objt <- ecdf(sampt)
  objh <- ecdf(samph)
  post_prob_t <- 1 - objt(cutoff.t)
  post_prob_h <- 1 - objh(cutoff.h)

  # prior prob
  priors <- omegasPrior(x$k, x$n.factors)
  priort <- ecdf(priors$omt_prior)
  priorh <- ecdf(priors$omh_prior)

  prior_prob_t <- 1 - priort(cutoff.t)
  prior_prob_h <- 1 - priorh(cutoff.h)

  out <- matrix(c(prior_prob_t, prior_prob_h, post_prob_t, post_prob_h, cutoff.t, cutoff.h),
                3, 2, byrow = TRUE)
  colnames(out) <- c("omega_t", "omega_h")
  rownames(out) <- c("prior_prob", "posterior_prob", "cutoff")
  return(out)
}
