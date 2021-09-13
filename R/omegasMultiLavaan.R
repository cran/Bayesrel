
omegaMultiF <- function(data, n.factors, interval, pairwise, model, model.type, fit.measures) {

  k <- ncol(data)
  if (model.type == "higher-order") {

    if (model == "balanced") {
      modfile <- lavMultiFileSeco(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFileSecoSyntax(k, n.factors, model, colnames(data))

      if (!modfile$colnames) { # name the variables in the dataset:
        names <- unlist(modfile$names)
        inds <- as.numeric(unlist(regmatches(names, gregexpr("[[:digit:]]+", names))))
        colnames(data)[inds] <- names
      }
    }

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = F)
    }


  } else if (model.type == "bi-factor") { # model.type is bifactor
    if (model == "balanced") {

      modfile <- lavMultiFileBif(k, n.factors)
      colnames(data) <- modfile$names

    } else { # if model syntax is specified
      modfile <- lavMultiFileBifSyntax(k, n.factors, model, colnames(data))

      if (!modfile$colnames) { # name the variables in the dataset:
        names <- unlist(modfile$names)
        inds <- as.numeric(unlist(regmatches(names, gregexpr("[[:digit:]]+", names))))
        colnames(data)[inds] <- names
      }
    }

    if (pairwise) {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = T, missing = "ml")
    } else {
      fit <- lavaan::cfa(modfile$model, data, std.lv = T, orthogonal = T)
    }
  }


  sts <- lavaan::parameterestimates(fit, level = interval)
  if (fit.measures) {
    modfile$fit.measures <- lavaan::fitmeasures(fit)
  }

  return(list(omhmean = sts$est[sts$label == "omega_h"], omtmean = sts$est[sts$label == "omega_t"],
              omhlow = sts$ci.lower[sts$label == "omega_h"], omhup = sts$ci.upper[sts$label == "omega_h"],
              omtlow = sts$ci.lower[sts$label == "omega_t"], omtup = sts$ci.upper[sts$label == "omega_t"],
              modfile = modfile))
}
