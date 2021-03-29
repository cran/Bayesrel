# gives freq omega, and loadings and errors
#

omegaFreqData <- function(
  data,
  interval,
  omega.int.analytic,
  pairwise,
  n.boot = 1e3,
  callback = function(){},
  parametric = FALSE){

  p <- ncol(data)
  n <- nrow(data)
  file <- lavOneFile(data)
  colnames(data) <- file$names

  lam_names <- paste("l", 1:p, sep = "")
  err_names <- paste("e", 1:p, sep = "")
  model <- paste0("f1 =~ ")
  loadings <- paste(paste(lam_names, "*", file$names, sep = ""),
                       collapse = " + ")
  errors <- paste(paste(file$names, " ~~ ", err_names, "*",
                           file$names, sep = ""), collapse = "\n")
  sum_loads <- paste("loading :=", paste(lam_names, collapse = " + "),
                      "\n")
  sum_errs <- paste("error :=", paste(err_names, collapse = " + "),
                    "\n")
  omega <- "omega := (loading^2) / ((loading^2) + error) \n"
  mod <- paste(model, loadings, "\n", errors,
                 "\n", sum_loads, sum_errs, omega)

  if (pairwise) {
    fit <- fitmodel_mis(mod, data)
  } else {
    fit <- fitmodel(mod, data)
  }

  if (is.null(fit)) {
    return(list(omega = NA, fit.object = NULL))
  } else {
    params <- lavaan::parameterestimates(fit, level = interval)
    omega <- params$est[params$lhs=="omega"]
    if (omega.int.analytic) {
      om_low <- params$ci.lower[params$lhs=="omega"]
      om_up <- params$ci.upper[params$lhs=="omega"]
      om_obj <- NA
    } else { # omega cfa with bootstrapping:
      if (parametric) {
        if (pairwise) {
          cc <- cov(data, use = "pairwise.complete.obs")
        } else {
          cc <- cov(data)
        }
        om_obj <- numeric(n.boot)
        for (i in 1:n.boot){
          boot_data <- MASS::mvrnorm(n, colMeans(data, na.rm = TRUE), cc)
          fit <- fitmodel(mod, boot_data)
          callback()
          if (!is.null(fit)) {
            params <- lavaan::parameterestimates(fit, level = interval)
            om_obj[i] <- params$est[params$lhs=="omega"]
          } else {
            om_obj[i] <- NA
          }
        }

        if (sum(!is.na(om_obj)) > 1) {
          om_low <- quantile(om_obj, prob = (1 - interval) / 2, na.rm = TRUE)
          om_up <- quantile(om_obj, prob = interval + (1 - interval) / 2, na.rm = TRUE)
        } else {
          om_low <- NA
          om_up <- NA
          om_obj <- NA
        }

      } else { # bootstrap non parametric

        om_obj <- numeric(n.boot)
        for (i in 1:n.boot){
          boot_data <- as.matrix(data[sample.int(n, size = n, replace = TRUE), ])
          if (pairwise) {
            fit <- fitmodel_mis(mod, boot_data)
          } else {
            fit <- fitmodel(mod, boot_data)
          }
          callback()
          if (!is.null(fit)) {
            params <- lavaan::parameterestimates(fit, level = interval)
            om_obj[i] <- params$est[params$lhs=="omega"]
          } else {
            om_obj[i] <- NA
          }
        }

        if (sum(!is.na(om_obj)) > 1) {
          om_low <- quantile(om_obj, prob = (1 - interval) / 2, na.rm = TRUE)
          om_up <- quantile(om_obj, prob = interval + (1 - interval) / 2, na.rm = TRUE)
        } else {
          om_low <- NA
          om_up <- NA
          om_obj <- NA
        }
      }
    }

    fit_tmp <- lavaan::fitMeasures(fit)
    indic <- c(fit_tmp["chisq"], fit_tmp["df"], fit_tmp["pvalue"],
               fit_tmp["rmsea"], fit_tmp["rmsea.ci.lower"], fit_tmp["rmsea.ci.upper"],
               fit_tmp["srmr"])
  }
  return(list(omega = omega, omega_lower = om_low, omega_upper = om_up, indices = indic, fit.object = fit,
              omega_boot = om_obj))
}


fitmodel <- function(mod, data) {
  out <- tryCatch(
    {
      lavaan::cfa(mod, data, std.lv = T)
    },
    error = function(cond) {
      return(NULL)
    },
    warning = function(cond) {
      return(NULL)
    },
    finally = {}
  )
  return(out)
}

fitmodel_mis <- function(mod, data) {
  out <- tryCatch(
    {
      lavaan::cfa(mod, data, std.lv = T, missing = "ML")
    },
    error = function(cond) {
      return(NULL)
    },
    warning = function(cond) {
      return(NULL)
    },
    finally = {}
  )
  return(out)
}
