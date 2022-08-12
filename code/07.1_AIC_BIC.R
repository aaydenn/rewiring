
calc_AIC <- function(nk, Ek, p) {
    return(2*Ek)
}

calc_BIC <- function(nk, Ek, p) {
    return(log(nk)*Ek)
}

calc_EBIC <- function(nk, Ek, p) {
    return(log(nk)*Ek + 4/2*Ek*log(p))
    # return((log(nk) + 2*log(p))*Ek)
}

#' @param mat.list.t.gau a list of matrices
IC_select <- function(mat.list.t.gau, lam1 = 0.1, lam2 = 0.1, returnJGL = F,
                        calc_IC = calc_AIC) {
  S <- lapply(mat.list.t.gau, cov) # a list of covariance matrices
  jgl <- JGL::JGL(mat.list.t.gau, lambda1 = lam1, lambda2 = lam2,
                        return.whole.theta = T) # a list (?) of partial covariances
  aic <- 0
  p <- ncol(jgl$theta[[1]]) # number of variables (?)
  for (k in 1:length(mat.list.t.gau)) {
    nk <- nrow(mat.list.t.gau[[k]]) # number of observations
    traceST <- sum(diag(S[[k]] %*% jgl$theta[[k]]))
    Ek <- sum(! jgl$theta[[k]] == 0) # numb. of non-null cells
    detT <- det(jgl$theta[[k]])
    aic <- aic + traceST*nk - log(detT)*nk + calc_IC(nk, Ek, p)/1e4
  }
  res <- c(lam1, lam2, aic)
  if (returnJGL){
    return(list(res, jgl))
  } else {
    return(res)
  }
}


#'---------------------------------------
#' Tuning parameter selection by AIC
#'
#' @param mat.list.t.gau the list of matrices, which contains the Gaussian distributed counts
#' @param lam1 numeric, Tuning parameter for the sparsity
#' @param lam2 numeric, Tuning parameter for the similarities across the conditions
#' @param returnJGL logical, whether return the JGL model results. default is F.
#' @return a list of 1) a vector of tuning parameter values and the AIC value; 2) the JGL result object
#' @export

AIC_select <- function(mat.list.t.gau, lam1 = 0.1, lam2 = 0.1, returnJGL = F){
  res.S <- lapply(mat.list.t.gau, cov)
  res.jgl <- JGL::JGL(mat.list.t.gau, lambda1 = lam1, lambda2 = lam2, return.whole.theta = T)
  aic <- 0; bic <- 0; ebic <- 0
  p <- ncol(res.jgl$theta[[1]])
  for (k in 1:length(mat.list.t.gau)){
    nk <- nrow(mat.list.t.gau[[k]])
    traceST <- sum(diag(res.S[[k]] %*% res.jgl$theta[[k]]))
    Ek <- sum(! res.jgl$theta[[k]] == 0)
    detT <- det(res.jgl$theta[[k]])
    aick <- (traceST*nk - log(detT)*nk + 2*Ek)/1e4
    aic <- aic + aick
    # bick <- (traceST*nk - log(detT)*nk + log(nk)*Ek)/1e4
    # bic <- bic + bick
    # ebick <- (traceST*nk - log(detT)*nk + log(nk)*Ek + 4/2*Ek*log(p))/1e4
    # ebic <- ebic + ebick
  }
  res <- c(lam1, lam2, aic)
  if (returnJGL){
    return(list(res,
                res.jgl))
  }else{
    return(res)
  }
}

#' Derive the Joint Graphical Lasso estimation result by AIC tuning parameter selection
#'
#' @param GauList a list of gaussian distributed matrices
#' @param l1vec a vector of candidate values for lambda1, if NULL, the values will be 1:5/20
#' @param l2vec a vector of candidate values for lambda2, if NULL, the values will be 1:5/50
#' @return aic.vec: a table of AIC values ; jgl.res: the JGL result object
#' @export
getJGLTuningParallel <- function(GauList,
                                l1vec = seq(from=0.05, to=0.25, by=0.05),
                                l2vec = seq(from=0.05, to=0.25, by=0.05)) {
  # Gaulist: samples by genes
  L1 <- rep(l1vec, each  = length(l2vec))
  L2 <- rep(l2vec, times = length(l1vec))
  jgl.res <- lapply(1:length(L1), function(i) {
    AIC_select(mat.list.t.gau = GauList, lam1 = L1[i], lam2 = L2[i],
              returnJGL = T)
  })
  return(jgl.res)
}

