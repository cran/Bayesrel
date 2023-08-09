

bomegasMultiOut <- function(data, n.factors, n.iter, n.burnin, thin, n.chains,
                            interval, model, impute, prior.params,
                            param.out, callback, pbtick, model.type) {


  out <- list()
  om_out <- omegaMultiBayes(data, n.factors, n.iter, n.burnin,
                         n.chains, thin, model, impute, prior.params, param.out, callback, pbtick, model.type)
  out$omega_t$chains <- om_out$omt
  out$omega_t$mean <- mean(om_out$omt)
  out$omega_t$cred <- coda::HPDinterval(coda::mcmc(c(om_out$omt)),
                                        prob = interval)
  if (model.type != "correlated") {
    out$omega_h$chains <- om_out$omh
    out$omega_h$mean <- mean(om_out$omh)
    out$omega_h$cred <- coda::HPDinterval(coda::mcmc(c(om_out$omh)),
                                          prob = interval)
  }

  out$implCovs <- apply(om_out$impl_covs, c(3, 4), as.vector)
  out$imputed_data <- om_out$imputed_values
  out$model <- list()
  out$model$file <- om_out$modfile
  if (param.out) {
    out$model$lambda <- om_out$lambda
    out$model$beta <- om_out$beta
    out$model$theta <- om_out$theta
    out$model$psi <- om_out$psi
    out$model$phi <- om_out$phi
  }
  return(out)
}
