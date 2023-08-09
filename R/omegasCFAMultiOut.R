


omegasCFAMultiOut <- function(data, n.factors, interval, fiml, model, model.type, fit.measures) {

  out <- list()
  om_out <- omegaMultiFreq(data, n.factors, interval, fiml, model, model.type, fit.measures)
  out$omega_t$est <- om_out$omtmean
  out$omega_t$conf <- c(om_out$omtlow, om_out$omtup)

  out$omega_h$est <- om_out$omhmean
  out$omega_h$conf <- c(om_out$omhlow, om_out$omhup)

  out$loadings$specific <- om_out$lambda
  out$residuals$specific <- om_out$theta
  out$loadings$general <- om_out$gloads
  out$residuals$general <- om_out$psi

  out$model <- om_out$modfile
  if (fit.measures) {
    out$fit.measures <- om_out$fit.measures
  }

  return(out)
}
