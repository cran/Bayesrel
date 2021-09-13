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
omegasSeco <- function(lambda, beta, theta, psi) {
  gl <- lambda %*% beta
  sl <- lambda %*% sqrt(psi)
  omh <- sum(gl %*% t(gl)) / (sum(gl %*% t(gl)) + sum(sl %*% t(sl)) + sum(theta))
  omt <- (sum(gl %*% t(gl)) + sum(sl %*% t(sl))) / (sum(gl %*% t(gl)) + sum(sl %*% t(sl)) + sum(theta))

  return(c(
    omh, omt
  ))
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

# get model implied covariance matrix
implCovMulti <- function(s, b, theta, psi) {
  i <- diag(nrow(b))
  ib_inv <- solve(i-b)
  out <- s %*% ib_inv %*% psi %*% t(ib_inv) %*% t(s) + theta
  out <- make_symmetric(out)
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

# generate balanced and simple model file for lavaan with specified syntax
lavMultiFileSecoSyntax <- function(k, ns, model, colnams) {

  beta_names <- paste("gl", 1:ns, sep = "")

  mod_out <- modelSyntaxExtract(model, colnams)
  mod <- mod_out$mod
  if (ns != mod_out$mod_n.factors) {
    warning("n.factors is unequal to specified factors in model syntax")
    ns <- mod_out$mod_n.factors
  }

  imat <- matrix(FALSE, k, ns)
  if (mod_out$mod_col_names) { # if the data variable names are specified in the modelsyntax
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

  } else {
    # get the digits out:
    out <- regmatches(mod, gregexpr("[[:digit:]]+", mod))
    idex <- lapply(out, as.numeric)
    lambda_names <- list()
    names <- list()
    for (i in seq_len(length(mod))) {
      imat[idex[[i]], i] <- TRUE
      lambda_names[[i]] <- paste0("l", i, seq_len(length(idex[[i]])))
      tmp <- unlist(strsplit(mod[[i]], " "))
      names[[i]] <- tmp[nchar(tmp) > 0]
    }
  }

  theta_names <- paste("e", 1:k, sep = "")

  gen <- paste(paste(beta_names, "*", "s", 1:ns, sep = ""), collapse = " + " )
  gen <- paste0("g =~ ", gen)

  modfile <- NULL
  for (i in 1:ns) {
    modfile <- paste0(modfile, "s", i, " =~ ", paste(paste(lambda_names[[i]], "*", names[[i]], sep = ""),
                                                     collapse = " + "), "\n")
  }
  modfile <- paste0(gen, " \n", modfile)

  # extra stuff for omega calc
  modfile <- paste0(modfile, "g ~~ 1*g\n")
  modfile <- paste0(modfile, paste(paste(unlist(names), " ~~ ", theta_names, "*",
                                         unlist(names), sep = ""), collapse = "\n"), "\n")

  # define extra parameters
  gl <- NULL
  sl <- NULL
  for (i in 1:ns) {
    gl <- c(gl, paste0(beta_names[i], "*", lambda_names[[i]]))
    sl <- paste0(sl, "ssum", i, " := ", paste0(lambda_names[[i]], collapse = " + "), "\n")
  }
  gl <- paste0(gl, collapse = "+")
  sum_gl <- paste("g_loading :=", gl, "\n")


  sum_sl <- paste0("spec_loading := ", paste(paste0("ssum", 1:ns, "^2"), collapse = " + "), "\n")
  sum_errs <- paste("residual_var :=", paste(theta_names, collapse = " + "), "\n")
  omega_h <- "omega_h := (g_loading^2) / (g_loading^2 + spec_loading + residual_var) \n"
  omega_t <- "omega_t := (g_loading^2 + spec_loading) / (g_loading^2 + spec_loading + residual_var) \n"
  model_out <- paste0(modfile, sl, sum_gl, sum_sl, sum_errs, omega_h, omega_t)

  return(out <- list(names = names, model = model_out, colnames = mod_out$mod_col_names))
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
  if (mod_out$mod_col_names) { # if the data variable names are specified in the modelsyntax
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

  } else {
    # get the digits out:
    out <- regmatches(mod, gregexpr("[[:digit:]]+", mod))
    idex <- lapply(out, as.numeric)
    lambda_names <- list()
    names <- list()
    for (i in seq_len(length(mod))) {
      imat[idex[[i]], i] <- TRUE
      lambda_names[[i]] <- paste0("l", i, seq_len(length(idex[[i]])))
      tmp <- unlist(strsplit(mod[[i]], " "))
      names[[i]] <- tmp[nchar(tmp) > 0]
    }
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

  return(list(names = names, model = model_out, colnames = mod_out$mod_col_names))
}


indexMatrix <- function(model, k, ns, colnams) {
  if (model == "balanced") {
    tmp <- matrix(seq(1:k), ns, k/ns, byrow=T)
    imat <- matrix(FALSE, k, ns)
    idex <- as.list(as.data.frame(t(tmp)))
    for (i in 1:ns) {
      imat[tmp[i, ], i] <- TRUE
    }
  } else { # extract variables from model syntax
    mod_out <- modelSyntaxExtract(model, colnams)
    mod <- mod_out$mod
    if (ns != mod_out$mod_n.factors) {
      warning("n.factors is unequal to specified factors in model syntax")
    }

    imat <- matrix(FALSE, k, ns)
    if (mod_out$mod_col_names) { # if the data variable names are specified in the modelsyntax
      idex <- list()
      for (i in seq_len(length(mod))) {
        tmp_mod <- unlist(strsplit(mod[i], " "))
        idex[[i]] <- which(colnams %in% tmp_mod)
        imat[idex[[i]], i] <- TRUE
      }

    } else {
      # get the digits out:
      out <- regmatches(mod, gregexpr("[[:digit:]]+", mod))
      idex <- lapply(out, as.numeric)
      for (i in seq_len(length(mod))) {
        imat[idex[[i]], i] <- TRUE
      }
    }
  }

  return(list(idex = idex, imat = imat))
}


modelSyntaxExtract <- function(model, colnams) {
  # get the factor and variable names
  mod_tmp <- unlist(strsplit(model, "[\n]|[=~]|[+]|[ ]"))
  # cleanup
  mod_names<- mod_tmp[nchar(mod_tmp) > 0]
  # check if colnames of the data are the varibale names in the syntax
  if (is.null(colnams)) {
    mod_col_names <- FALSE
  } else {
    if (all(colnams %in% mod_names)) {
      mod_col_names <- TRUE
    } else {
      mod_col_names <- FALSE
    }
  }

  # check the number of specified factors by checking the frequency of "=~"
  mod_n.factors <- length(unlist(gregexpr("=~", model)))

  # cleanup the model to further find the paths
  mod_separated1 <- unlist(strsplit(model, "[\n]|[+]|[ ]"))
  mod_separated_tmp <- mod_separated1[nchar(mod_separated1) > 0]
  mod_separated2 <- mod_separated_tmp[-(grep("=~", mod_separated_tmp)-1)]
  mod_separated_tmp2 <- unlist(strsplit(paste(mod_separated2, collapse = " "), "[=~]"))
  mod_separated3 <- mod_separated_tmp2[nchar(mod_separated_tmp2) > 0]

  return(list(mod = mod_separated3, mod_n.factors = mod_n.factors, mod_col_names = mod_col_names,
              mod_names = mod_names))
}


