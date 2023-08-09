# create lavaan cfa one factor model file from data

lavOneFile <- function(data) {
  p <- ncol(data)
  v <- 0
  for (i in 1:p) {
    v[i] <- paste0("x", i)
  }
  v <- paste0(v, collapse = "+")
  mod <- paste0("g=~", v) # dynamic lavaan model file
  mod <- paste0(mod, "; g ~~ 1*g") # fix the factor variance to 1

  # column names specify
  names <- 0
  for (i in 1:p) {
    names[i] <- paste0("x", i)
  }
  return(list(names = names, model = mod))
}


# calculate omega from loadings and residual (error variances)

omegaBasic <- function(l, e) {
  o <- sum(l)^2 / (sum(l)^2 + sum(e))
  return(o)
}

# omegas
omegasSeco <- function (lambda, beta, theta, psi) {
  gl <- lambda %*% beta
  sl <- lambda %*% sqrt(psi)
  sum_gl_glt <- sum_X_Xt(gl)
  sum_sl_slt <- sum_X_Xt(sl)
  sum_theta  <- sum(theta)
  temp <- sum_gl_glt + sum_sl_slt
  omh <- sum_gl_glt / (temp + sum_theta)
  omt <- temp / (temp + sum_theta)
  return(c(omh, omt))
}

omegasBif <- function(sl, gl, e, psi = 1) {
  sl <- sl * psi
  ll <- cbind(gl, sl)
  glProd <- sum(tcrossprod(gl))
  slProd <- sum(tcrossprod(sl))
  omh <- glProd / (glProd + slProd + sum(e))
  omt <- (glProd + slProd) / (glProd + slProd + sum(e))
  return(c(omh, omt))
}

# multivariate normal data with matrix of means and V
genNormDataTweak <- function(n, m, Sigma){
  p <- ncol(Sigma)
  randomData <- matrix(rnorm(n*p), n, p)
  cc <- chol(Sigma)
  out <- randomData %*% cc
  out <- out + matrix(m, nrow = n, ncol = ncol(Sigma), byrow = FALSE)
  return(out)
}

genNormDataLegit <- function(n, m, Sigma){
  p <- ncol(Sigma)
  randomData <- matrix(rnorm(n*p), n, p)
  cc <- chol(Sigma)
  out <- randomData %*% cc
  out <- out + matrix(m, nrow = n, ncol = ncol(Sigma), byrow = T)
  return(out)
}

# get model implied covariance matrix
implCovMultiSeco <- function(s, b, theta, psi) {
  i <- diag(nrow(b))
  ib_inv <- solve(i-b)
  out <- s %*% ib_inv %*% psi %*% t(ib_inv) %*% t(s) + theta
  out <- make_symmetric(out)
  return(out)
}

# get model implied covariance matrix
implCovMultiBif <- function(s, b, theta, psi) {
  mat_bi <- cbind(b, s)
  out <- mat_bi %*% psi %*% t(mat_bi) + theta
  return(out)
}

# generate balanced and simple model file for lavaan
lavMultiFileSeco <- function(k, ns) {

  beta_names <- paste("gl", 1:ns, sep = "")
  lambda_names <- matrix(paste("sl", 1:k, sep = ""), k/ns, ns)
  names <- matrix(paste("x", 1:k, sep = ""), k/ns, ns)
  theta_names <- paste("e", 1:k, sep = "")

  gen <- paste(paste(beta_names, "*", "s", 1:ns, sep = ""), collapse = " + " )
  gen <- paste0("g =~ ", gen)

  mod <- NULL
  for (i in 1:ns) {
    mod <- paste0(mod, "s", i, " =~ ", paste(paste(lambda_names[, i], "*", names[, i], sep = ""),
                                             collapse = " + "), "\n")
  }
  mod <- paste0(gen, " \n", mod)

  # extra stuff for omega calc
  mod <- paste0(mod, "g ~~ 1*g\n")
  mod <- paste0(mod, paste(paste(c(names), " ~~ ", theta_names, "*",
                                 c(names), sep = ""), collapse = "\n"), "\n")

  # define extra parameters
  gl <- NULL
  sl <- NULL
  for (i in 1:ns) {
    gl <- c(gl, paste0(beta_names[i] , "*", lambda_names[, i]))
    sl <- paste0(sl, "ssum", i, " := ", paste0(lambda_names[, i], collapse = " + "), "\n")
  }
  gl <- paste0(gl, collapse = "+")
  sum_gl <- paste("g_loading :=", gl, "\n")


  sum_sl <- paste0("spec_loading := ", paste(paste0("ssum", 1:ns, "^2"), collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_h <- "omega_h := (g_loading^2) / (g_loading^2 + spec_loading + residual_var) \n"
  omega_t <- "omega_t := (g_loading^2 + spec_loading) / (g_loading^2 + spec_loading + residual_var) \n"
  mod <- paste0(mod, sl, sum_gl, sum_sl, sum_errs, omega_h, omega_t)

  return(list(names = names, model = mod))
}

# generate model file for lavaan
lavMultiFileSecoSyntax <- function(k, ns, model, colnams, model_opts) {

  beta_names <- paste("gl", 1:ns, sep = "")

  mod_out <- modelSyntaxExtract(model, colnams)
  mod <- mod_out$mod
  if (ns != mod_out$mod_n.factors) {
    warning("n.factors is unequal to specified factors in model syntax")
    ns <- mod_out$mod_n.factors
  }

  imat <- model_opts$imat
  idex <- model_opts$idex
  lambda_names <- list()
  names <- list()
  for (i in seq_len(length(mod))) {
    tmp_mod <- unlist(strsplit(mod[[i]], " "))
    idex[[i]] <- which(colnams %in% tmp_mod)
    imat[idex[[i]], i] <- TRUE
    lambda_names[[i]] <- paste0("l", i, seq_len(length(idex[[i]])))
    names[[i]] <- tmp_mod[nchar(tmp_mod) > 0]
  }

  theta_names <- paste("e", 1:k, sep = "")

  gen <- paste(paste(beta_names, "*", mod_out$fac_names, sep = ""), collapse = " + " )
  gen <- paste0("g =~ ", gen)

  modfile <- NULL
  for (i in 1:ns) {
    modfile <- paste0(modfile, mod_out$fac_names[i], " =~ ", paste(paste(lambda_names[[i]], "*", names[[i]], sep = ""),
                                                     collapse = " + "), "\n")
  }
  modfile <- paste0(gen, " \n", modfile)

  # extra stuff for omega calc
  modfile <- paste0(modfile, "g ~~ 1*g\n")
  modfile <- paste0(modfile, paste(paste(unique(unlist(names)), " ~~ ", theta_names, "*",
                                         unique(unlist(names)), sep = ""), collapse = "\n"), "\n")

  # define extra parameters
  gl <- NULL
  sl <- NULL
  for (i in 1:ns) {
    gl <- c(gl, paste0(beta_names[i], "*", lambda_names[[i]]))
    sl <- paste0(sl, "ssum", i, " := ", paste0(lambda_names[[i]], collapse = " + "), "\n")
  }
  gl <- paste0(gl, collapse = " + ")
  sum_gl <- paste("g_loading :=", gl, "\n")

  sum_sl <- paste0("spec_loading := ", paste(paste0("ssum", 1:ns, "^2"), collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_h <- "omega_h := (g_loading^2) / (g_loading^2 + spec_loading + residual_var) \n"
  omega_t <- "omega_t := (g_loading^2 + spec_loading) / (g_loading^2 + spec_loading + residual_var) \n"
  model_out <- paste0(modfile, sl, sum_gl, sum_sl, sum_errs, omega_h, omega_t)

  return(out <- list(model = model_out, factor_names = mod_out$fac_names))
}


lavMultiFileBif <- function(k, ns) {

  beta_names <- paste("gl", 1:k, sep = "")
  lambda_names <- matrix(paste("sl", 1:k, sep = ""), k/ns, ns)
  names <- matrix(paste("x", 1:k, sep = ""), k/ns, ns)
  theta_names <- paste("e", 1:k, sep = "")

  gen <- paste("x", 1:k, sep="")
  gen <- paste(paste(beta_names, "*", c(names), sep = ""),
               collapse = " + ")
  gen <- paste0("g =~ ", gen)

  mod <- NULL
  for (i in 1:ns) {
    mod <- paste0(mod, "s", i, " =~ ", paste(paste(lambda_names[, i], "*", names[, i], sep = ""),
                                             collapse = " + "), "\n")
  }
  mod <- paste0(gen, " \n", mod)

  # extra stuff for omega calc
  mod <- paste0(mod, "g ~~ 1*g\n")
  for (i in 1:ns) {
    mod <- paste0(mod, "s", i, " ~~ 1*s", i, "\n")
  }

  mod <- paste0(mod, paste(paste(c(names), " ~~ ", theta_names, "*",
                                 c(names), sep = ""), collapse = "\n"), "\n")

  sum_gl <- paste("g_loading :=", paste0(beta_names, collapse = " + "),
                  "\n")

  load_sum <- NULL
  for (i in 1:ns) {
    load_sum <- paste0(load_sum, "ssum", i, " := ", paste0(lambda_names[, i], collapse = " + "), "\n")
  }

  sum_sl <- paste0("spec_loading := ", paste(paste0("ssum", 1:ns, "^2"), collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_h <- "omega_h := (g_loading^2) / (g_loading^2 + spec_loading + residual_var) \n"
  omega_t <- "omega_t := (g_loading^2 + spec_loading) / (g_loading^2 + spec_loading + residual_var) \n"
  mod <- paste0(mod, load_sum, sum_gl, sum_sl, sum_errs, omega_h, omega_t)

  return(list(names = c(names), model = mod))
}


lavMultiFileBifSyntax <- function(k, ns, model, colnams) {

  beta_names <- paste("gl", 1:k, sep = "")

  mod_out <- modelSyntaxExtract(model, colnams)
  mod <- mod_out$mod
  if (ns != mod_out$mod_n.factors) {
    warning("n.factors is unequal to specified factors in model syntax")
    ns <- mod_out$mod_n.factors
  }

  imat <- matrix(FALSE, k, ns)
  idex <- list()
  lambda_names <- list()
  names <- list()
  for (i in seq_len(length(mod))) {
    tmp_mod <- unlist(strsplit(mod[[i]], " "))
    idex[[i]] <- which(colnams %in% tmp_mod)
    imat[idex[[i]], i] <- TRUE
    lambda_names[[i]] <- paste0("l", i, seq_len(length(idex[[i]])))
    names[[i]] <- tmp_mod[nchar(tmp_mod) > 0]
  }

  theta_names <- paste("e", 1:k, sep = "")

  gen <- paste("x", 1:k, sep="")
  gen <- paste(paste(beta_names, "*", unlist(names), sep = ""),
               collapse = " + ")
  gen <- paste0("g =~ ", gen)

  modfile <- NULL
  for (i in 1:ns) {
    modfile <- paste0(modfile, "s", i, " =~ ", paste(paste(lambda_names[[i]], "*", names[[i]], sep = ""),
                                                     collapse = " + "), "\n")
  }
  modfile <- paste0(gen, " \n", modfile)

  # extra stuff for omega calc
  modfile <- paste0(modfile, "g ~~ 1*g\n")
  for (i in 1:ns) {
    modfile <- paste0(modfile, "s", i, " ~~ 1*s", i, "\n")
  }

  modfile <- paste0(modfile, paste(paste(unlist(names), " ~~ ", theta_names, "*",
                                         unlist(names), sep = ""), collapse = "\n"), "\n")

  sum_gl <- paste("g_loading :=", paste0(beta_names, collapse = " + "),
                  "\n")

  load_sum <- NULL
  for (i in 1:ns) {
    load_sum <- paste0(load_sum, "ssum", i, " := ", paste0(lambda_names[[i]], collapse = " + "), "\n")
  }

  sum_sl <- paste0("spec_loading := ", paste(paste0("ssum", 1:ns, "^2"), collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_h <- "omega_h := (g_loading^2) / (g_loading^2 + spec_loading + residual_var) \n"
  omega_t <- "omega_t := (g_loading^2 + spec_loading) / (g_loading^2 + spec_loading + residual_var) \n"
  model_out <- paste0(modfile, load_sum, sum_gl, sum_sl, sum_errs, omega_h, omega_t)

  return(list(names = names, model = model_out))
}


indexMatrix <- function(model = NULL, k, ns = NULL, colnams = NULL) {

  if (is.null(model)) {
    tmp <- matrix(seq(1:k), ns, k/ns, byrow=T)
    imat <- matrix(FALSE, k, ns)
    idex <- as.list(as.data.frame(t(tmp)))
    for (i in 1:ns) {
      imat[tmp[i, ], i] <- TRUE
    }
    mod_out <- NULL
  } else { # extract variables from model syntax
    mod_out <- modelSyntaxExtract(model, colnams)
    mod <- mod_out$mod
    ns <- mod_out$mod_n.factors

    imat <- matrix(FALSE, k, ns)
    idex <- list()
    for (i in seq_len(length(mod))) {
      tmp_mod <- unlist(strsplit(mod[i], " "))
      idex[[i]] <- which(colnams %in% tmp_mod)
      imat[idex[[i]], i] <- TRUE
    }
  }

  return(list(idex = idex, imat = imat, mod_out = mod_out, n.factors = ns))
}


modelSyntaxExtract <- function(model, colnams) {
  # get the factor and variable names
  mod_tmp <- unlist(strsplit(model, "[\n]|[=~]|[+]|[ ]"))
  # cleanup
  all_names<- mod_tmp[nchar(mod_tmp) > 0]
  var_names <- unique(all_names[!is.na(match(all_names, colnams))])
  fac_names <- unique(all_names[is.na(match(all_names, colnams))])

  # confirm that only the group factors are specified
  if (sum(all_names %in% fac_names) > length (fac_names)) {
    stop("Please only specify the group factors and their relations to the manifest variables in the syntax.")
  }

  # check the number of specified factors by checking the frequency of "=~"
  mod_n.factors <- length(unlist(gregexpr("=~", model)))
  var_count <- length(unique(all_names)) - mod_n.factors
  if (var_count > length(var_names)) {
    stop("It seems some of the item names in the model syntax do not equal the variable names in the data. Check the spelling.")
  }

  # cleanup the model to further find the paths
  mod_separated1 <- unlist(strsplit(model, "[\n]|[+]|[ ]"))
  mod_separated_tmp <- mod_separated1[nchar(mod_separated1) > 0]
  mod_separated2 <- mod_separated_tmp[-(grep("=~", mod_separated_tmp)-1)]
  mod_separated_tmp2 <- unlist(strsplit(paste(mod_separated2, collapse = " "), "[=~]"))
  mod_separated3 <- mod_separated_tmp2[nchar(mod_separated_tmp2) > 0]

  return(list(mod = mod_separated3, mod_n.factors = mod_n.factors, var_names = var_names,
              fac_names = fac_names))
}




SRMR <- function(cdat, impl) {
  nvar <- ncol(cdat)
  e <- (nvar * (nvar + 1)) / 2
  sqrt.d <- 1 / sqrt(diag(cdat))
  D <- diag(sqrt.d, ncol = length(sqrt.d))
  R <- D %*% (cdat - impl) %*% D
  srmr <- sqrt(sum(R[lower.tri(R, diag = TRUE)]^2) / e)
  return(srmr)
}


LRblav <- function(data, cmat, basell) {
  tmpll <- dmultinorm(data, cmat)
  out <- 2 * (basell - sum(tmpll))
  return(out)
}


BRMSEA <- function(chisq, p, pD, n) {
  dChisq <- (chisq - p)
  dChisq[dChisq < 0] <- 0
  rmsea <- sqrt(dChisq / ((p - pD) * n))
  return(rmsea)
}

# borrowed that from mnormt package
dmultinorm <- function(x, varcov, mm = 0) {
  d <- ncol(varcov)
  X <- t(x - mm)
  varcov <- (varcov + t(varcov))/2
  u <- chol(varcov, pivot = FALSE)
  inv <- chol2inv(u)
  logdet <- 2 * sum(log(diag(u)))
  Q <- colSums((inv %*% X) * X)
  Q <- Q[!is.na(Q)]
  logPDF <- as.vector(Q + d * logb(2 * pi) + logdet)/(-2)
  return(logPDF)
}

implCovUni <- function(ll) {
  out <- ll[1, ] %*% t(ll[1, ]) + diag(ll[2, ])
  return(make_symmetric(out))
}


Xt_A_X <- function(A, X) {
  crossprod(X, A %*% X)
}

X_A_Xt <- function(A, X) {
  tcrossprod(X %*% A, X)
}

Xt_invA_X_0 <- function(A, x) {
  t(x) %*% solve(A) %*% x
}

Xt_invA_X_1 <- function(A, x) {
  crossprod(x, solve(A, x))
}

Xt_invA_X_2 <- function(A, x) {
  Xt_invChol_X_2(chol(A), x)
}

Xt_invChol_X_2 <- function(L, x) {
  crossprod(backsolve(L, x, transpose = TRUE))
}

diag_Xt_Y <- function(X, Y) {
  # based on https://stackoverflow.com/a/42569902/ but does diag(t(X) %*% Y) rather than diag(X %*% Y)
  .colSums(X * Y, nrow(X), ncol(Y))
}
diag_Xt_X <- function(X) diag_Xt_Y(X, X)

diag_X_Yt <- function(X, Y){
  .rowSums(X * Y, nrow(X), ncol(Y))
}

diag_X_Xt <- function(X) diag_X_Yt(X, X)

sum_X_Xt <- function(X) {
  crossprod(.colSums(X, nrow(X), ncol(X)))
}


# convert a index matrix to a lavaan syntax -> helps with simulations
# model file for seco model with crossloadings
# also works for the correlated model
lavMultiFileCross <- function(imat) {
  k <- nrow(imat)
  ns <- ncol(imat)
  names <- 0
  for(i in 1:k){
    names[i] <- paste0("x", i)
  }

  str <- list()
  mod <- NULL
  for (i in 1:ns) {
    str[[i]] <- paste0(names[imat[,i]], collapse = " + ")
    mod <- paste0(mod, "s", i, " =~ ", str[[i]], " \n ")
  }

  return(out <- list(names = names, model = mod))
}

omegaCorr <- function(ll, phi, tt) {
  return(sum(colSums(ll) * colSums(ll) %*% phi) / (sum(colSums(ll) * colSums(ll) %*% phi) + sum(tt)))
}

implCovCorr <- function(lambda, phi, theta) {
  return(lambda %*% phi %*% t(lambda) + theta)
}


# generate lav model file from no input

lavMultiFileCorr <- function(k, ns) {

  lambda_names <- matrix(paste("sl", 1:k, sep = ""), k/ns, ns)
  names <- matrix(paste("x", 1:k, sep = ""), k/ns, ns)
  theta_names <- paste("e", 1:k, sep = "")

  mod <- NULL
  ss <- paste0("s", 1:ns)
  for (i in 1:ns) {
    mod <- paste0(mod, ss[i], " =~ ", paste(paste(lambda_names[, i], "*", names[, i], sep = ""),
                                                                   collapse = " + "), "\n")
  }

  ss_comb <- combn(ss, m = 2)
  rho_names <- paste("rl", 1:ncol(ss_comb), sep = "")
  rhos <- NULL
  vars <- paste0(ss, " ~~ ", "1*", ss, "\n", collapse = "")
  for (i in 1:ncol(ss_comb)) {
    rhos <- paste0(rhos, ss_comb[1, i], " ~~ ", paste0(rho_names[i], "*", ss_comb[2, i]), "\n")
  }

  mod <- paste0(mod, "\n", rhos, "\n", vars, "\n")

  # correlation matrix for latents:
  cor_mat <- matrix(NA, ns, ns)
  diag(cor_mat) <- 1
  cor_mat[lower.tri(cor_mat)] <- rho_names
  cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]

  # extra stuff for omega calc
  mod <- paste0(mod, paste(paste(c(names), " ~~ ", theta_names, "*",
                                 c(names), sep = ""), collapse = "\n"), "\n")

  # define extra parameters
  sl <- NULL
  tp <- NULL
  for (i in 1:ns) {
    sl <- paste0(sl, "ssum", i, " := ", paste0(lambda_names[, i], collapse = " + "), "\n")
    tp <- paste0(tp, "sm", i, " := (", paste0("ssum", 1:ns, "*", cor_mat[, i], collapse = " + "),
                 ")*ssum", i, "\n")
  }

  sum_sl <- paste0("spec_loading := ", paste0("sm", 1:ns, collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_t <- "omega_t := (spec_loading) / (spec_loading + residual_var) \n"
  mod <- paste0(mod, sl, tp, sum_sl, sum_errs, omega_t)

  return(list(names = names, model = mod))
}


# generate model file for lavaan from input model file
lavMultiFileCorrSyntax <- function(k, ns, model, colnams, model_opts) {

  mod_out <- model_opts$mod_out
  mod <- mod_out$mod

  if (ns != model_opts$mod_out$mod_n.factors) {
    warning("n.factors is unequal to specified factors in model syntax")
    ns <- model_opts$mod_out$mod_n.factors
  }

  imat <- model_opts$imat
  idex <- model_opts$idex
  lambda_names <- list()
  names <- list()
  for (i in seq_len(length(mod))) {
    tmp_mod <- unlist(strsplit(mod[[i]], " "))
    idex[[i]] <- which(colnams %in% tmp_mod)
    imat[idex[[i]], i] <- TRUE
    lambda_names[[i]] <- paste0("l", i, seq_len(length(idex[[i]])))
    names[[i]] <- tmp_mod[nchar(tmp_mod) > 0]
  }

  theta_names <- paste("e", 1:k, sep = "")

  modfile <- NULL
  for (i in 1:ns) {
    modfile <- paste0(modfile, mod_out$fac_names[i], " =~ ", paste(paste(lambda_names[[i]], "*", names[[i]], sep = ""),
                                                                   collapse = " + "), "\n")
  }

  vars <- mod_out$fac_names
  vars <- paste0(vars, " ~~ ", "1*", vars, "\n", collapse = "")
  ss <- mod_out$fac_names
  ss_comb <- combn(ss, m = 2)
  rho_names <- paste("rl", 1:ncol(ss_comb), sep = "")

  rhos <- NULL
  for (i in 1:ncol(ss_comb)) {
    rhos <- paste0(rhos, ss_comb[1, i], " ~~ ", paste0(rho_names[i], "*", ss_comb[2, i]), "\n")
  }

  modfile <- paste(modfile, "\n", rhos, "\n", vars, "\n")

  # correlation matrix for latents:
  cor_mat <- matrix(NA, ns, ns)
  diag(cor_mat) <- 1
  cor_mat[lower.tri(cor_mat)] <- rho_names
  cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]

  # extra stuff for omega calc
  modfile <- paste0(modfile, paste(paste(unique(unlist(names)), " ~~ ", theta_names, "*",
                                         unique(unlist(names)), sep = ""), collapse = "\n"), "\n")

  sl <- NULL
  tp <- NULL
  for (i in 1:ns) {
    sl <- paste0(sl, "ssum", i, " := ", paste0(lambda_names[[i]], collapse = " + "), "\n")
    tp <- paste0(tp, "sm", i, " := (", paste0("ssum", 1:ns, "*", cor_mat[, i], collapse = " + "),
                 ")*ssum", i, "\n")
  }

  sum_sl <- paste0("spec_loading := ", paste0("sm", 1:ns, collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_t <- "omega_t := (spec_loading) / (spec_loading + residual_var) \n"
  model_out <- paste0(modfile, sl, tp, sum_sl, sum_errs, omega_t)

  return(out <- list(model = model_out, factor_names = mod_out$fac_names))
}


