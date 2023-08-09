
#' Estimate reliability estimates for multidimensional scales in the frequentist framework
#' @description When supplying a data set that is multidimensional
#' the function estimates the reliability of the set by means of omega_total
#' and the general factor saturation of the set by means of omega_hierarchical
#' The procedure entails fitting a hierarchical factor model using a CFA.
#' The second-order (hierarchical, higher-order), the bi-factor, and the correlated factor model
#' can be used in the CFA. The CFA is performed using lavaan 'Yves Rosseel', <https://CRAN.R-project.org/package=lavaan>.
#' Coefficients omega_t and omega_h (only for second-order and bi-factor model)
#' can be computed from the factor model parameters.
#'
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
#' @param model.type A string denoting if the model that should be fit is the second-order or
#' bi-factor model or the correlated factor model. This comes down to the researcher's theory about the measurement
#' and the model fit.
#' @param interval A number specifying the confidence interval, which is Wald-type
#' @param missing A string denoting the missing data handling, can be "fiml" (full information ML) or "listwise".
#' Specifying "pairwise" will defulat to "fiml"
#' @param fit.measures A logical denoting if fit.measures from the CFA should be computed,
#' the output then contains the chisq statistic, chisq df, chisq p-value, cfi, tli,
#' rmsea, rmsea 90\% ci lower, rmsea 90\% ci upper, rmsea<.05 p-value, aic, bic,
#' unbiased srmr, unbiased srmr 90\% ci lower, unbiased srmr 90\% ci upper, unbiased srmr<.05 p-value
#'
#' @return The point estimates and the Wald-type confidence intervals for
#' omega_t and omega_h (for the second-order and bi-factor model)
#'
#' @examples
#' res <- omegasCFA(upps, n.factors = 5, model = NULL, model.type = "bi-factor",
#' missing = "listwise")
#'
#' # or with specified model syntax relating the group factors to the items:
#' model <- "f1 =~ U17_r + U22_r + U29_r + U34_r
#' f2 =~ U4 + U14 + U19 + U27
#' f3 =~ U6 + U16 + U28 + U48
#' f4 =~ U23_r + U31_r + U36_r + U46_r
#' f5 =~ U10_r + U20_r + U35_r + U52_r"
#' res <- omegasCFA(upps, n.factors = 5, model = model, model.type = "second-order",
#' missing = "listwise")
#'
#' @importFrom utils combn
#'
#' @export
omegasCFA <- function(
  data,
  n.factors = NULL,
  model = NULL,
  model.type = "second-order",
  interval = .95,
  missing = "fiml",
  fit.measures = FALSE
) {

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

  listwise <- FALSE
  fiml <- FALSE
  complete_cases <- nrow(data)
  if (anyNA(data)) {
    if (missing == "listwise") {
      pos <- which(is.na(data), arr.ind = TRUE)[, 1]
      data <- data[-pos, ]
      ncomp <- nrow(data)
      complete_cases <- ncomp
      listwise <- TRUE
    } else { # missing is pairwise
      fiml <- TRUE
    }
  }

  data <- scale(data, scale = FALSE)

  sum_res <- omegasCFAMultiOut(data, n.factors, interval, fiml, model, model.type, fit.measures)

  sum_res$complete_cases <- complete_cases
  sum_res$call <- match.call()
  sum_res$k <- ncol(data)
  sum_res$n.factors <- n.factors
  sum_res$fiml <- fiml
  sum_res$listwise <- listwise
  sum_res$interval <- interval
  sum_res$model.type <- model.type
  class(sum_res) <- "omegasCFA"

  return(sum_res)
}
