

bomegasMultiOut <- function(data, n.factors, n.iter, n.burnin, thin, n.chains,
                            interval, model, pairwise, callback) {


  out <- list()
  om_out <- omegaMultiB(data, n.factors, n.iter, n.burnin, n.chains, thin, model, pairwise, callback)
  out$omega_t$chains <- om_out$omt
  out$omega_t$mean <- mean(om_out$omt)
  out$omega_t$cred <- coda::HPDinterval(coda::mcmc(c(om_out$omt)), prob = interval)

  out$omega_h$chains <- om_out$omh
  out$omega_h$mean <- mean(om_out$omh)
  out$omega_h$cred <- coda::HPDinterval(coda::mcmc(c(om_out$omh)), prob = interval)

  out$implCovs <- apply(om_out$impl_covs, c(3, 4), as.vector)
  out$imputed_data <- om_out$imputed_values

  out$model <- om_out$modfile

  return(out)
}
