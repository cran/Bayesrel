# this function calls on other functions in order to return the sampled estimates
# and the credible intervals together with the posterior distribution objects
# to be passed on for forther analysis

gibbsFun <- function(data, n.iter, n.burnin, estimates, interval, item.dropped, pairwise){
  p <- ncol(data)
  res <- list()
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates){
    C <- covSamp(data, n.iter, n.burnin, pairwise)
    if (item.dropped) {
      Ctmp <- array(0, c(p, n.iter - n.burnin, p - 1, p - 1))
      for (i in 1:p){
        Ctmp[i, , , ] <- C[, -i, -i]
      }
    }
  } else {
    C = NULL
  }
  res$covsamp <- C

  if ("alpha" %in% estimates){
    res$samp$Bayes_alpha <- coda::mcmc(apply(C, MARGIN = 1, applyalpha))
    int <- coda::HPDinterval(res$samp$Bayes_alpha, prob = interval)
    res$cred$low$Bayes_alpha <- int[1]
    res$cred$up$Bayes_alpha <- int[2]
    res$est$Bayes_alpha<- median(res$samp$Bayes_alpha)
    if (item.dropped){
      res$ifitem$samp$alpha <- coda::mcmc(apply(Ctmp, c(2, 1), applyalpha))
      res$ifitem$est$alpha <- apply(res$ifitem$samp$alpha, 2, median)
      res$ifitem$cred$alpha <- coda::HPDinterval(res$ifitem$samp$alpha, prob = interval)
    }
  }

  if ("lambda2" %in% estimates){
    res$samp$Bayes_lambda2 <- coda::mcmc(apply(C, MARGIN = 1, applylambda2))
    int <- coda::HPDinterval(res$samp$Bayes_lambda2, prob = interval)
    res$cred$low$Bayes_lambda2 <- int[1]
    res$cred$up$Bayes_lambda2 <- int[2]
    res$est$Bayes_lambda2<- median(res$samp$Bayes_lambda2)
    if (item.dropped){
      res$ifitem$samp$lambda2 <- coda::mcmc(apply(Ctmp, c(2, 1), applylambda2))
      res$ifitem$est$lambda2 <- apply(res$ifitem$samp$lambda2, 2, median)
      res$ifitem$cred$lambda2 <- coda::HPDinterval(res$ifitem$samp$lambda2, prob = interval)

    }
  }

  if ("lambda4" %in% estimates){
    res$samp$Bayes_lambda4 <- coda::mcmc(apply(C, MARGIN = 1, applylambda4))
    int <- coda::HPDinterval(res$samp$Bayes_lambda4, prob = interval)
    res$cred$low$Bayes_lambda4 <- int[1]
    res$cred$up$Bayes_lambda4 <- int[2]
    res$est$Bayes_lambda4<- median(res$samp$Bayes_lambda4)
    if (item.dropped){
      res$ifitem$samp$lambda4 <- coda::mcmc(apply(Ctmp, c(2, 1), applylambda4))
      res$ifitem$est$lambda4 <- apply(res$ifitem$samp$lambda4, 2, median)
      res$ifitem$cred$lambda4 <- coda::HPDinterval(res$ifitem$samp$lambda4, prob = interval)

    }
  }

  if ("lambda6" %in% estimates){
    res$samp$Bayes_lambda6 <- coda::mcmc(apply(C, MARGIN = 1, applylambda6))
    int <- coda::HPDinterval(res$samp$Bayes_lambda6, prob = interval)
    res$cred$low$Bayes_lambda6 <- int[1]
    res$cred$up$Bayes_lambda6 <- int[2]
    res$est$Bayes_lambda6<- median(res$samp$Bayes_lambda6)
    if (item.dropped){
      res$ifitem$samp$lambda6 <- coda::mcmc(apply(Ctmp, c(2, 1), applylambda6))
      res$ifitem$est$lambda6 <- apply(res$ifitem$samp$lambda6, 2, median)
      res$ifitem$cred$lambda6 <- coda::HPDinterval(res$ifitem$samp$lambda6, prob = interval)

    }
  }

  if ("glb" %in% estimates){
    res$samp$Bayes_glb <- coda::mcmc(glbOnArray(C))
    if (sum(is.na(res$samp$Bayes_glb) > 0)) {
      int <- c(NA, NA)
    } else {
      int <- coda::HPDinterval(res$samp$Bayes_glb, prob = interval)
    }
    res$cred$low$Bayes_glb <- int[1]
    res$cred$up$Bayes_glb <- int[2]
    res$est$Bayes_glb<- median(res$samp$Bayes_glb)
    if (item.dropped){
      res$ifitem$samp$glb <- coda::mcmc(apply(Ctmp, c(2, 1), glbOnArray))
      res$ifitem$est$glb <- apply(res$ifitem$samp$glb, 2, median)
      res$ifitem$cred$glb <- coda::HPDinterval(res$ifitem$samp$glb, prob = interval)

    }
  }

  # special case omega -----------------------------------------------------------------
  if ("omega" %in% estimates){
    om_samp <- omegaSampler(data, n.iter, n.burnin, pairwise)
    res$samp$Bayes_omega <- coda::mcmc(om_samp$omega)
    res$loadings <- apply(om_samp$lambda, 2, mean)
    res$resid_var <- apply(om_samp$psi, 2, mean)
    # res$loadings <- coda::mcmc(om_samp$lambda)
    # res$resid_var <- coda::mcmc(om_samp$psi)
    int <- coda::HPDinterval(res$samp$Bayes_omega, prob = interval)
    res$cred$low$Bayes_omega <- int[1]
    res$cred$up$Bayes_omega<- int[2]
    res$est$Bayes_omega <- mean(res$samp$Bayes_omega)

    if (item.dropped){
      om_samp_ifitem <- matrix(0, n.iter - n.burnin, p)
      for (i in 1:p){
        tmp <- data[-i, -i]
        om_samp_ifitem[, i] <- omegaSampler(tmp, n.iter, n.burnin, pairwise)$omega
      }
      res$ifitem$samp$omega <- coda::mcmc(om_samp_ifitem)
      res$ifitem$est$omega <- apply(om_samp_ifitem, 2, mean)
      res$ifitem$cred$omega <- coda::HPDinterval(res$ifitem$samp$omega, prob = interval)

    }
  }

  return(res)

}
