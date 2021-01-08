
#' Function to decon using the CIBERSORT math
#' 
#' Performs decon using v-SVR (nu support vector regression)
#' Notes from cibersort paper:
#' Our current implementation of CIBERSORT executes v-SVR using the \code{svm} function in the R package,
#' e1071. Regression coefficients are extracted with the following R command:
#'  \code{coef <- t(model$coefs) %*% model$SV}. They use best of nu = 0.25, 0.5, 0.75}.
#' 
#' @param Y p-length expression vector or p * N expression matrix - the actual (linear-scale) data
#' @param X p * K Training matrix 
#' @param bg scalar or matrix of expected background counts per data point. 
#' @param weights The same as the weights argument used by lm
#' @return a list with one element, beta
#' @import e1071
deconSVR = function(Y, X, bg = 0, weights = NULL) {
  if (length(bg) == 1) {
    bg = matrix(bg, nrow(Y), ncol(Y))
  }
  # init results object:
  beta = matrix(NA, ncol(X), ncol(Y))
  rownames(beta) = colnames(X)
  svrfits = list()
  # run SVR one sample at a time:
  for (i in 1:ncol(Y)) {
    y = pmax(Y[, i ] - bg[, i], 0)
    use = !is.na(y)
    s = svm(y[use] ~ X[use, ] - 1, 
            nu = 0.75,
            type = "nu-regression",
            kernel = "linear",
            scale = FALSE)
    svrfits[[i]] = s
    # extract coefs:
    coef <- t(s$coefs) %*% s$SV
    # rescale and thresh at 0 ala cibersort:
    coef = pmax(coef, 0)
    beta[, i] = coef 
  }
  
  out = list(beta = beta, svrfits = svrfits)
  return(out)
}



#' Function to decon using lm() with a non-negativity constraint
#' 
#' Uses optim() to minimize squared resids (i.e. run a linear model) under the constraint that coefs are non-negative.
#' 
#' @param Y p-length expression vector or p * N expression matrix - the actual (linear-scale) data
#' @param X p * K Training matrix 
#' @param bg scalar, vector or matrix of expected background counts per data point. 
#' @param weights The same as the weights argument used by lm
#' @return a list with one element, beta
deconUsingConstrainedLM = function(Y, X, bg = 0, weights = NULL) {
  if (is.vector(Y)) {
    Y = matrix(Y, nrow = length(Y))
  }
  if (length(bg) == 1) {
    bg = matrix(bg, nrow(Y), ncol(Y))
  }
  beta = matrix(NA, ncol(X), ncol(Y), dimnames = list(colnames(X), colnames(Y)))
  for (i in 1:ncol(Y)) {
    y = Y[, i]
    b = bg[, i]
    use = !is.na(y)
    
    # loss function
    scoreSqErrLoss = function(beta) {
      sqrt(mean((y[use] - (X[use, ] %*% beta + b[use]))^2))
    }
    # init with a lm:
    init = lm(y ~ X - 1)$coefficients
    # minimze loss subject to non-negativity constraint:
    op = optim(par = init, fn = scoreSqErrLoss, lower = 0, method = "L-BFGS-B")
    beta[, i] = op$par
  }
  out = list(beta = beta)
  return(out)
}



#' Function to decon using logNormReg package to run linear mean model and log error model
#' 
#' Calls lognlm() to optimize the model. 
#' 
#' @param Y p-length expression vector or p * N expression matrix - the actual (linear-scale) data
#' @param X p * K Training matrix 
#' @param bg scalar or matrix of expected background counts per data point. 
#' @param weights The same as the weights argument used by lm
#' @return a list: beta (estimate), sigmas (covariance matrix of estimate, derived by inverting the hessian from lognlm)
#' @import logNormReg
#' @export
deconLNR = function(Y, X, bg = 0, weights = NULL, epsilon = NULL) {
  if (length(weights) == 0) {
    weights = replace(Y, TRUE, 1)
  }
  if (is.vector(Y)) {
    Y = matrix(Y, nrow = length(Y))
  }
  if (length(bg) == 1) {
    bg = matrix(bg, nrow(Y), ncol(Y))
  }
  # choose "epsilon": a very small non-zero number to make fits well-behaved
  if (length(epsilon) == 0) {
    epsilon = min(replace(Y, (Y == 0) & !is.na(Y), NA), na.rm = T)
  }
  beta = matrix(NA, ncol(X), ncol(Y), dimnames = list(colnames(X), colnames(Y)))
  sigmas = array(NA, dim = c(ncol(X), ncol(X), ncol(Y)), dimnames = list(colnames(X), colnames(X), colnames(Y)))
  for (i in 1:ncol(Y)) {
    y = Y[, i]
    b = bg[, i]
    wts = weights[, i]
    
    # remove NA data:
    use = !is.na(y)
    y = y[use]
    b = b[use]
    Xtemp = X[use,  , drop = F]
    wts = wts[use]
    
    init = rep(mean(y) / (mean(X) * ncol(X)), ncol(X)); names(init) = colnames(X)
    
    # run lognlm:
    fit = lognlm(pmax(y, epsilon) ~ b + Xtemp - 1,
                 lik = FALSE,
                 weights = wts,
                 #start = rep(1, ncol(X) + 1),
                 start = c(1, init),
                 method = "L-BFGS-B", 
                 lower = c(1, rep(0, ncol(Xtemp))), 
                 upper = c(1, rep(Inf, ncol(Xtemp))),
                 opt = "optim")
    # save estimate, excluding the intercept:
    beta[, i] = fit$coefficients[-1]
    # save vcov, excluding the intercept:
    sigmas[, , i] = solve(fit$hessian)[-1, -1]
  }
  out = list(beta = pmax(beta, 0), sigmas = sigmas)
  return(out)
}


#' decon using DWLS
deconUsingDWLS = function(Y, X, bg = 0, weights = NULL) {
  if (is.vector(Y)) {
    Y = matrix(Y, nrow = length(Y))
  }
  if (length(bg) == 1) {
    bg = matrix(bg, nrow(Y), ncol(Y))
  }
  beta = matrix(NA, ncol(X), ncol(Y), dimnames = list(colnames(X), colnames(Y)))
  for (i in 1:ncol(Y)) {
    
    #y = Y[, i]
    y = pmax(Y[, i ] - bg[, i], 0)
    b = bg[, i]
    use = !is.na(y)
    
    beta[, i] = solveDampenedWLS(S = X[use, ] * 1e-3, B = y[use])

  }
  out = list(beta = beta)
  return(out)
}
