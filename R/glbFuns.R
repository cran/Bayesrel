
# adjusted code from Rcsdp package and from psych Package,
# small adjustment to easy_spd function from CSDP library:
# Revelle, W. (2018) psych: Procedures for Personality and Psychological Research,
# Northwestern University, Evanston, Illinois, USA,
# https://CRAN.R-project.org/package=psych Version = 1.8.4.
# Rcsdp Package of Hector Corrada Bravo
# CSDP Library by Brian Borchers
glbOnArrayCustom <- function(Cov, callback = function(){}, printlevel = 0) {

  d <- dim(Cov)
  if (length(d) == 2L) { # turn it into an array if it is a matrix
    d <- c(1L, d)
    dim(Cov) <- d
  }

  p <- d[2L]

  opt <- rep.int(1L, p)
  A <- vector("list", p)
  for (i in seq_len(p)) {
    b <- rep(0, p)
    b[i] <- 1
    A[[i]] <- list(diag(b), -b, b)
  }

  K <- list(type = c("s", "l", "l"), size = rep(p, 3))

  prob_info <- get.prob.info2(K, length(b))
  Lo_Bounds <- rep(0, p)

  cv <- Cov[1L, , ]
  Var <- diag(cv)
  C <- list(diag(Var) - cv, -Var, Lo_Bounds)

  # make the Rcsdp object once instead of each iteration
  prob_data <- list(
    C = blkmatrix_R2csdp2(C, prob_info),
    A = constraints_R2csdp2(A, prob_info),
    b = as.double(c(0, opt))
  )

  arg1 <- as.integer(sum(prob_info$block.sizes))
  arg2 <- as.integer(prob_info$nconstraints)
  arg3 <- as.integer(prob_info$nblocks)
  arg4 <- as.integer(c(0, prob_info$block.types))
  arg5 <- as.integer(c(0, prob_info$block.sizes))


  ret <- csdpArma(
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    prob_data$C,
    prob_data$A,
    prob_data$b,
    Cov,
    callback,
    printlevel
  )
  return(ret)


}

get.prob.info2 <- function(K, m) {
  block.types <- (K$type != "s") + 1L
  nblocks <- length(K$type)
  block.sizes <- K$size
  nconstraints <- m
  ret <- list(nblocks = nblocks, nconstraints = nconstraints, block.types = block.types,
              block.sizes = block.sizes)
  return(ret)
}

vector_R2csdp <- function(x) c(0, x)

simple_triplet_sym_matrix <- function(i, j, v, n = max(c(i, j)), check.ind = FALSE) {
  if (check.ind && any(i < j)) {
    stop("Index arguments 'i' and 'j' do not point to the lower triangle. Swapping indices")
  }
  structure(list(i = i, j = j, v = v, n = n), class = "simple_triplet_sym_matrix")
}

.simple_triplet_zero_sym_matrix <- function(n, mode = "double") {
  simple_triplet_sym_matrix(integer(), integer(), vector(mode, 0L), n)
}

as.simple_triplet_sym_matrix.matrix <- function(x, check.sym = FALSE) {
  if (prod(dim(x)) == 0L)
    return(simple_triplet_sym_matrix(integer(), integer(),
                                     c(x), n = nrow(x)))
  if (nrow(x) != ncol(x))
    stop("Argument 'x' must be a square matrix")
  if (check.sym && (sum(abs(x - t(x)) / 2) != 0))
    stop("Argument 'x' must be a symmetric matrix")
  ind <- which(x != vector(typeof(x), 1L), arr.ind = TRUE)
  if (length(ind) == 0)
    return(.simple_triplet_zero_sym_matrix(nrow(x)))
  ind <- ind[ind[, 1L] >= ind[, 2L], , drop = FALSE]
  simple_triplet_sym_matrix(ind[, 1L], ind[, 2L], x[ind], n = nrow(x))
}

blkmatrix_R2csdp2 <- function(X, prob.info) {
  do.one.block <- function(blocknum) {
    cur.block <- X[[blocknum]]
    cur.type <- prob.info$block.types[blocknum]
    cur.size <- prob.info$block.sizes[blocknum]
    if (cur.type == 1L) {
      data <- as.double(cur.block)
    }
    else {
      data <- vector_R2csdp(cur.block)
    }
    structure(list(blocksize = as.integer(cur.size), blockcategory = as.integer(cur.type),
                   data = as.double(data)), class = "csdpBlkMat")
  }
  nblocks <- prob.info$nblocks
  list(nblocks = nblocks, blocks = lapply(seq_along(X), do.one.block))
}

constraints_R2csdp2 <- function(A, prob.info) {
  nblocks <- prob.info$nblocks
  do.one.constraint <- function(constraintnum) {
    Ai <- A[[constraintnum]]
    ret <- vector("list", nblocks)
    k <- 0L
    for (j in seq_len(nblocks)) {
      blocknum <- j
      Aij <- Ai[[blocknum]]
      cur.type <- prob.info$block.types[blocknum]
      cur.size <- prob.info$block.sizes[blocknum]
      if (cur.type == 1L) {
        tmp <- as.simple_triplet_sym_matrix.matrix(Aij)
        if (length(tmp$v) == 0L)
          next
        iindices <- vector_R2csdp(tmp$i)
        jindices <- vector_R2csdp(tmp$j)
        entries <- vector_R2csdp(tmp$v)
      }
      else {
        nnz <- which(Aij != 0L)
        if (length(nnz) == 0L)
          next
        iindices <- vector_R2csdp(nnz)
        jindices <- vector_R2csdp(nnz)
        entries <- vector_R2csdp(Aij[nnz])
      }
      k <- k + 1L
      ret[[k]] <- structure(list(iindices = as.integer(iindices),
                                 jindices = as.integer(jindices), entries = as.double(entries),
                                 blocknum = as.integer(blocknum), blocksize = as.integer(cur.size),
                                 constraintnum = as.integer(constraintnum), numentries = as.integer(length(entries) -
                                                                                                      1L)), class = "csdpConstrMat")
    }
    ret[1L:k]
  }
  lapply(seq_along(A), do.one.constraint)
}


