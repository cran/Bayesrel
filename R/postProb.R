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
#' pStrel(strel(asrm, "lambda2", n.chains = 2, n.iter = 100, freq = FALSE), "lambda2", .80)
#' @export
pStrel <- function(x, estimate, low.bound) {

  posi1 <- grep(estimate, x$estimates, ignore.case = TRUE)
  samp <- as.vector(x$Bayes$samp[[posi1]])
  obj <- ecdf(samp)
  post_prob <- 1 - obj(low.bound)

  # prior prob
  n.item <- dim(x$Bayes$covsamp)[3]

  if (is.null(x$priors$df0)) {
    x$priors$df0 <- n.item
  }

  prior <- priorSampUni(n.item, estimate, k0 = x$priors$k0, df0 = x$priors$df0, a0 = x$priors$a0, b0 = x$priors$b0,
                        m0 = x$priors$m0)

  end <- length(prior[["x"]])
  poslow <- end - sum(prior[["x"]] > low.bound)
  prior_prob <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]])
  out <- c(prior_prob, post_prob)
  names(out) <- c("prior_prob", "posterior_prob")
  return(out)
}


#' prior and posterior probability of omega_t and omega_h being bigger than thresholds
#' @description
#' takes mcmc posterior samples of omega_t and omega_h
#' and calculates the prior and posterior probability of the estimate being bigger
#' or smaller than an arbitrary value
#'
#' @param x A strel output object (list)
#' @param cutoff.t A number indicating the threshold for omega_t
#' @param cutoff.h A number indicating the threshold for omega_h
#'
#' @examples
#' pOmegas(bomegas(upps, n.factors = 5, n.chains = 2, n.iter = 100, n.burnin = 50,
#' missing = "listwise"))
#' @export
pOmegas <- function(x, cutoff.t = .80, cutoff.h = .60) {

  sampt <- as.vector(x$omega_t$chains)
  samph <- as.vector(x$omega_h$chains)
  objt <- ecdf(sampt)
  objh <- ecdf(samph)
  post_prob_t <- 1 - objt(cutoff.t)
  post_prob_h <- 1 - objh(cutoff.h)

  # prior prob
  priors <- omegasPrior(x$k, x$n.factors, nsamp = 2e3, x$prior_params$a0, x$prior_params$b0,
                        x$prior_params$l0, x$prior_params$A0,
                        x$prior_params$c0, x$prior_params$d0,
                        x$prior_params$beta0, x$prior_params$B0,
                        x$prior_params$p0, x$prior_params$R0)

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
