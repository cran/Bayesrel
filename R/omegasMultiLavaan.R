
omegaMultiF <- function(data, n.factors, interval, pairwise, model, model.type, fit.measures) {

  k <- ncol(data)
  model_opts <- indexMatrix(model, k, n.factors, colnames(data))

  if (model.type == "higher-order") {

    if (is.null(model)) {
      modfile <- lavMultiFileSeco(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFileSecoSyntax(k, n.factors, model, colnames(data))
    }

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = FALSE, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = FALSE)
    }

    sts <- lavaan::parameterestimates(fit, level = interval, standardized = TRUE)
    gloads <- sts$std.all[1:n.factors]
    lmat <- matrix(0, k, n.factors)
    lmat[model_opts$imat] <- sts$std.all[(n.factors + 1):(n.factors + k)]
    theta <- sts$std.all[(n.factors + k + 2):(n.factors + 2*k + 1)]
    psi <- sts$std.all[(n.factors + 2*k + 2):(2*n.factors + 2*k + 1)]

  } else if (model.type == "bi-factor") { # model.type is bifactor
    if (is.null(model)) {

      modfile <- lavMultiFileBif(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFileBifSyntax(k, n.factors, model, colnames(data))
    }

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = TRUE, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = TRUE, orthogonal = TRUE)
    }

    sts <- lavaan::parameterestimates(fit, level = interval, standardized = TRUE)
    gloads <- sts$std.all[1:k]
    lmat <- matrix(0, k, n.factors)
    lmat[model_opts$imat] <- sts$std.all[(k + 1):(2 * k)]
    psi <- sts$std.all[(2 * k + 2):(2 * k + 1 + n.factors)]
    theta <- sts$std.all[(2 * k + 2 + n.factors):(3 * k + 1 + n.factors)]

  }


  if (fit.measures) {
    modfile$fit.measures <- lavaan::fitmeasures(fit)
    if (pairwise) {
      modfile$srmr.summary <- lavaan::lavResiduals(fit)$summary["total"]
    } else {
      modfile$srmr.summary <- lavaan::lavResiduals(fit)$summary["cov"]
    }
  }


  return(list(omhmean = sts$est[sts$label == "omega_h"], omtmean = sts$est[sts$label == "omega_t"],
              omhlow = sts$ci.lower[sts$label == "omega_h"], omhup = sts$ci.upper[sts$label == "omega_h"],
              omtlow = sts$ci.lower[sts$label == "omega_t"], omtup = sts$ci.upper[sts$label == "omega_t"],
              lambda = lmat, gloads = gloads, theta = theta, psi = psi,
              modfile = modfile))
}
