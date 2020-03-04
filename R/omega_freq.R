# gives freq omega, and loadings and errors
#

omegaFreqData <- function(data, pairwise){
  p <- ncol(data)
  file <- lavOneFile(data)
  colnames(data) <- file$names

  lam_names <- paste("l", 1:p, sep = "")
  err_names <- paste("e", 1:p, sep = "")
  model <- paste0("f1 =~ ")
  loadings <- paste(paste(lam_names, "*", file$names, sep = ""),
                       collapse = " + ")
  factors <- "f1 ~~ 1*f1\n"
  errors <- paste(paste(file$names, " ~~ ", err_names, "*",
                           file$names, sep = ""), collapse = "\n")
  sum_loads <- paste("loading :=", paste(lam_names, collapse = " + "),
                      "\n")
  sum_errs <- paste("error :=", paste(err_names, collapse = " + "),
                    "\n")
  omega <- "omega := (loading^2) / ((loading^2) + error) \n"
  mod <- paste(model, loadings, "\n", factors, errors,
                 "\n", sum_loads, sum_errs, omega)

  if (pairwise) {
    fit <- try(lavaan::cfa(mod, data, std.lv = T, missing = "ML"), silent = TRUE)
  } else {
    fit <- try(lavaan::cfa(mod, data, std.lv = T), silent = TRUE)
  }
  params <- try(lavaan::parameterestimates(fit), silent = TRUE)
  if ("try-error" %in% class(params)) {
    load <- resid <- omega <- om_low <- om_up <- fit_tmp <- indic <- NA
  } else {
    omega <- params$est[params$lhs=="omega"]
    om_low <- params$ci.lower[params$lhs=="omega"]
    om_up <- params$ci.upper[params$lhs=="omega"]


    fit_tmp <- lavaan::fitMeasures(fit)
    indic <- c(fit_tmp["chisq"], fit_tmp["df"], fit_tmp["pvalue"],
               fit_tmp["rmsea"], fit_tmp["rmsea.ci.lower"], fit_tmp["rmsea.ci.upper"],
               fit_tmp["srmr"])
  }

  return(list(omega = omega, omega_lower = om_low, omega_upper = om_up, indices = indic))
}

