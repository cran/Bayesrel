
#' Estimate reliability estimates for multidimensional scales in the frequentist framework
#' @description When supplying a data set that is multidimensional
#' the function estimates the reliability of the set by means of omega_total
#' and the general factor saturation of the set by means of omega_hierarchical
#' The procedure entails fitting a hierarchical factor model using a CFA.
#' Both the higher-order (second-order) and the bi-factor model can be used in the CFA.
#' The CFA is fit using lavaan 'Yves Rosseel', <https://CRAN.R-project.org/package=lavaan>.
#' Coefficients omega_t and omega_h can be computed from the factor model parameters.
#'
#' @param data A matrix or data.frame containing multivariate observations,
#' rows = observations, columns = variables/items
#' @param n.factors A number specifying the number of group factors that the items load on
#' @param model A string that by default ("balanced") distributes the items evenly
#' among the number of group factors. This only works if the items are a multiple
#' of the number of group factors and the items are already grouped in the data set,
#' meaning, e.g., items 1-5 load on one factor, 6-10 on another, and so on.
#' A model file can be specified in lavaan syntax style (f1=~.+.+.) to relate the items
#' to the group factors. The items can either be named as the columns in the data set
#' or x1, ..., xn, where 1,...,n correspond to the column numbers
#' @param model.type A string denoting if the model that should be fit is the higher-order or
#' bi-factor model. This comes down to the researcher's theory about the measurement
#' and the model fit
#' @param interval A number specifying the confidence interval, which is Wald-type
#' @param missing A string denoting the missing data handling, can be "fiml" (full information ML) or "listwise".
#' @param fit.measures A logical denoting if fit.measures from the CFA should be computed,
#' the output then contains the chisq statistic, chisq df, chisq p-value, cfi, tli,
#' rmsea, rmsea 90\% ci lower, rmsea 90\% ci upper, rmsea<.05 p-value, aic, bic,
#' unbiased srmr, unbiased srmr 90\% ci lower, unbiased srmr 90\% ci upper, unbiased srmr<.05 p-value
#'
#' @return The point estimates and the Wald-type confidence intervals for
#' omega_t and omega_h
#'
#' @examples
#' res <- omegasCFA(upps, n.factors = 5, model = "balanced", model.type = "bi-factor",
#' missing = "listwise")
#'
#' # or with specified model syntax relating the group factors to the items:
#' model <- "f1 =~ x1+x2+x3+x4
#' f2 =~ x5+x6+x7+x8
#' f3 =~ x9+x10+x11+x12
#' f4 =~ x13+x14+x15+x16
#' f5 =~ x17+x18+x19+x20 "
#' res <- omegasCFA(upps, n.factors = 5, model = model, model.type = "higher-order",
#' missing = "listwise")
#'
#' @export
omegasCFA <- function(
  data,
  n.factors,
  model = "balanced",
  model.type = "higher-order",
  interval = .95,
  missing = "fiml",
  fit.measures = FALSE) {

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
    } else { # missing is pairwise
      pairwise <- TRUE
    }
  }

  data <- scale(data, scale = FALSE)

  sum_res <- omegasCFAMultiOut(data, n.factors, interval, pairwise, model, model.type, fit.measures)

  sum_res$complete_cases <- complete_cases
  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$pairwise <- pairwise
  sum_res$listwise <- listwise
  sum_res$interval <- interval
  class(sum_res) <- "omegasCFA"

  return(sum_res)
}
