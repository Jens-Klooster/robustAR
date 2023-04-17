#' Robust AR Statistic
#'
#' Computes the robust AR statistic introduced in Klooster and Zhelonkin (202x).
#'
#' @param Y A vector of dependent variables.
#' @param X A vector of endogenous covariates.
#' @param Z A matrix/vector of instrumental variables.
#' @param W A matrix/vector of control variables.
#' @param beta0 A hypothesized beta0 value.
#' @param type Specifies which downweighting function to use: Huber, Tukey or OLS.
#' @param weighting Specifies which weighting method to use: mcd, hat matrix or no weights.
#' @param homoskedasticity Assume homoskedasticity: yes or no.
#' @returns A numeric scalar test statistic.
#' @export
robustAR  <-  function(Y, X, Z, W, beta0,
                       type = c("Huber", "Tukey", "OLS"),
                       weighting = c("hat", "mcd", "no"),
                       homoskedasticity = "no"
){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  W <- as.matrix(W)

  WZ <- cbind(W,Z)
  n <- length(Y)
  k <- dim(Z)[2]
  p <- dim(W)[2]

  if(weighting == "hat"){
    W.hat <- W
    # check if there is a column of ones. If there is, we remove it.
    if(length(which(apply(W.hat, 2, function(a) length(unique(a))==1) == TRUE)) > 0){
      constant = length(which(apply(W.hat, 2, function(a) length(unique(a))==1) == TRUE))
      column.of.ones <- which(apply(W.hat, 2, function(a) length(unique(a))==1) == TRUE)
      W.hat <- W.hat[, -column.of.ones]
    }
    WZ_ <- cbind(W.hat, Z)
    PWZ_ = WZ_ %*% solve(t(WZ_) %*% WZ_) %*% t(WZ_)
    covariate.weights <- sqrt(1 - diag(PWZ_))
  }

  if(weighting == "mcd"){
    W.mcd <- W
    # check if there is a column of ones. If there is, we remove it.
    constant = 0
    if(length(which(apply(W.mcd, 2, function(a) length(unique(a))==1) == TRUE)) > 0){
      constant = length(which(apply(W.mcd, 2, function(a) length(unique(a))==1) == TRUE))
      column.of.ones <- which(apply(W.mcd, 2, function(a) length(unique(a))==1) == TRUE)
      W.mcd <- W.mcd[, -column.of.ones]
    }
    ZW_ <- cbind(Z,W.mcd)
    mcd.obj <- MASS::cov.mcd(ZW_)
    rob.center <- mcd.obj$center
    rob.cov <- mcd.obj$cov
    mahal.obj <- mahalanobis(ZW_, rob.center, rob.cov, seed = 1)
    covariate.weights <- rep(1, n)
    chi.value <- qchisq(0.95, p + k - constant)
    for(i in 1:n){
      if(mahal.obj[i] > chi.value){
        covariate.weights[i] <- chi.value/mahal.obj[i]
      }
    }
  }

  if(weighting == "no"){
    covariate.weights <- rep(1,n)
  }


  # this function will determine which penalization function will be used: Huber, Tukey, or OLS.
  penfun <- function(x, type = c("Huber", "Tukey", "OLS")){
    if(type == "Huber"){return(huber_function(x))}
    if(type == "Tukey"){return(tukey_function(x))}
    if(type == "OLS"){return(id_function(x))}
  }
  penfunderiv <- function(x, type = c("Huber", "Tukey", "OLS")){
    if(type == "Huber"){return(huber_function_derivative(x))}
    if(type == "Tukey"){return(tukey_function_derivative(x))}
    if(type == "OLS"){return(id_function_derivative(x))}
  }


  restricted_fun <- function(Y, X, Z, W, beta0, type = c("Huber", "Tukey", "OLS"), weights){
    dep_var <- Y - X*beta0
    if(type == "Huber"){return(MASS::rlm(dep_var ~ 0 + W, weights = weights, scale.est = "proposal 2",  maxit = 100, acc = 1e-10))}
    if(type == "Tukey"){return(MASS::rlm(dep_var ~ 0 + W, weights = weights, psi = MASS::psi.bisquare, scale.est = "proposal 2",  maxit = 100, acc = 1e-10))}
    if(type == "OLS"){return(lm(dep_var ~ 0 + W))}
  }


  #the beta0-restricted model object.
  restricted_obj <- restricted_fun(Y = Y, X = X, Z = Z, W = W, beta0 = beta0, type = type, weights = covariate.weights)
  restricted_resid <- as.matrix(residuals(restricted_obj))

  # scale estimates of the beta0-restricted model object.
  if(type == "Huber" | type == "Tukey"){
    restricted_scale <- restricted_obj$s
  }
  if(type == "OLS"){
    restricted_scale <- sigma(restricted_obj)
  }

  # standardized residuals to which the penalization function is applied
  stan_resids_restricted <- apply(restricted_resid/restricted_scale, 1, penfun, type = type)
  stan_resids_restricted_deriv <- apply(restricted_resid/restricted_scale, 1, penfunderiv, type = type) * (1/restricted_scale)

  Z.matrix <- (1/n) * (covariate.weights * stan_resids_restricted) %*% Z

  M.matrix <- matrix(rep(0), nrow = p+k, ncol = p+k)
  Q.matrix <- matrix(rep(0), nrow = p+k, ncol = p+k)


  if(homoskedasticity == "yes"){
    for(j in 1:n){
      M.matrix <- M.matrix + -covariate.weights[j]  * WZ[j,] %*%  t(WZ[j,])
      Q.matrix <- Q.matrix + covariate.weights[j]^2 * WZ[j,] %*%  t(WZ[j,])
    }

    M.matrix <- (1/n) * M.matrix * mean(stan_resids_restricted_deriv)
    Q.matrix <- (1/n) * Q.matrix * (1/(n-k-p)) * sum(stan_resids_restricted^2)
  }

  if(homoskedasticity == "no"){
    for(j in 1:n){
      M.matrix <- M.matrix + covariate.weights[j] * -stan_resids_restricted_deriv[j] * WZ[j,] %*%  t(WZ[j,])
      Q.matrix <- Q.matrix + covariate.weights[j]^2 * stan_resids_restricted[j]^2 * WZ[j,] %*% t(WZ[j,])
    }
    M.matrix <- (1/n) * M.matrix
    Q.matrix <- (1/n) * Q.matrix
  }

  V.matrix <- solve(M.matrix) %*% Q.matrix %*% solve(t(M.matrix))
  V22 <- V.matrix[(p+1):(p+k), (p+1):(p+k)]

  M22 <- M.matrix[(p+1):(p+k), (p+1):(p+k)]
  M21 <- M.matrix[(p+1):(p+k), 1:p]
  M12 <- M.matrix[1:p, (p+1):(p+k)]
  M11 <- M.matrix[1:p, 1:p]
  M22.1 <- M22 - M21 %*% solve(M11) %*% M12

  U.matrix <- M22.1 %*% V22 %*% t(M22.1)

  AR.stat <- Z.matrix %*% solve(U.matrix) %*% t(Z.matrix)
  return(n * AR.stat)
}


#' (Classical) AR Statistic
#'
#' Computes the (classical) AR statistic introduced in Anderson and Rubin (1949).
#'
#' @param Y A vector of dependent variables.
#' @param X A vector of endogenous covariates.
#' @param Z A matrix/vector of instrumental variables.
#' @param W A matrix/vector of control variables.
#' @param beta0 A hypothesized beta0 value.
#' @returns A numeric scalar test statistic.
#' @export
AR <- function(Y, X, Z, W, beta0){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  W <- as.matrix(W)
  p <- dim(W)[2]
  k <- dim(Z)[2]
  n <- dim(Y)[1]

  # partial out the control variables
  MW <- diag(1, nrow = n, ncol = n) - W %*% solve(t(W) %*% W) %*% t(W)
  Z_ <- MW %*% Z
  X_ <- MW %*% X
  Y_ <- MW %*% Y

  # construct projection and annihilator matrix for the transformed Z variable
  PZ_ <- Z_ %*% solve(t(Z_) %*% Z_) %*% t(Z_)
  MZ_ <- diag(1, nrow = n, ncol = n) - PZ_

  # calculate AR test statistic
  AR.stat <- ((t(Y_ - beta0 * X_) %*% PZ_ %*% (Y_ - beta0 * X_))/k)/((t(Y_ - beta0 * X_) %*% MZ_ %*% (Y_ - beta0 * X_))/(n - p - k))
  return(AR.stat)
}

#' Robust AR Test confidence set.
#'
#' Automatically computes a confidence set for the robust AR test.
#'
#' @param Y A vector of dependent variables.
#' @param X A vector of endogenous covariates.
#' @param Z A matrix/vector of instrumental variables.
#' @param W A matrix/vector of control variables.
#' @param betagrid A grid of beta values
#' @param alpha A scalar significance level between 0 and 1
#' @param type Specifies which downweighting function to use: Huber, Tukey or OLS.
#' @param weighting Specifies which weighting method to use: mcd, hat matrix or no weights.
#' @param homoskedasticity Assume homoskedasticity: yes or no.
#' @returns A list with a confidence set and corresponding p-values.
#' @export
robustAR.conf <-  function(Y, X, Z, W, betagrid, alpha,
                           type = c("Huber", "Tukey", "OLS"),
                           weighting = c("hat", "mcd", "no"),
                           homoskedasticity = "no"){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  W <- as.matrix(W)

  WZ <- cbind(W,Z)
  n <- length(Y)
  k <- dim(Z)[2]
  p <- dim(W)[2]

  if(weighting == "hat"){
    W.hat <- W
    # check if there is a column of ones. If there is, we remove it.
    if(length(which(apply(W.hat, 2, function(a) length(unique(a))==1) == TRUE)) > 0){
      constant = length(which(apply(W.hat, 2, function(a) length(unique(a))==1) == TRUE))
      column.of.ones <- which(apply(W.hat, 2, function(a) length(unique(a))==1) == TRUE)
      W.hat <- W.hat[, -column.of.ones]
    }
    WZ <- cbind(W.hat, Z)
    PWZ = WZ %*% solve(t(WZ) %*% WZ) %*% t(WZ)
    covariate.weights <- sqrt(1 - diag(PWZ))
  }

  if(weighting == "mcd"){
    W.mcd <- W
    # check if there is a column of ones. If there is, we remove it.
    constant = 0
    if(length(which(apply(W.mcd, 2, function(a) length(unique(a))==1) == TRUE)) > 0){
      constant = length(which(apply(W.mcd, 2, function(a) length(unique(a))==1) == TRUE))
      column.of.ones <- which(apply(W.mcd, 2, function(a) length(unique(a))==1) == TRUE)
      W.mcd <- W.mcd[, -column.of.ones]
    }
    ZW_ <- cbind(Z,W.mcd)
    mcd.obj <- MASS::cov.mcd(ZW_)
    rob.center <- mcd.obj$center
    rob.cov <- mcd.obj$cov
    mahal.obj <- mahalanobis(ZW_, rob.center, rob.cov)
    covariate.weights <- rep(1, n)
    chi.value <- qchisq(0.95, p + k - constant)
    for(i in 1:n){
      if(mahal.obj[i] > chi.value){
        covariate.weights[i] <- chi.value/mahal.obj[i]
        # covariate.weights[i] <- 0
      }
    }
  }

  if(weighting == "no"){
    covariate.weights <- rep(1,n)
  }


  # Huber penalization function and its derivative
  huber_function <- function(x){
    c = 1.345
    if(abs(x) < c){
      return(x)
    }
    else{
      return(c * sign(x))
    }
  }
  huber_function_derivative <- function(x){
    c = 1.345
    if(abs(x) < c){
      return(1)
    }
    else{
      return(0)
    }
  }

  # Tukey penalization function and its derivative
  tukey_function <- function(x){
    c=4.685
    if(abs(x) < c){
      return(x*(1 - x^2/c^2)^2)
    }
    else{
      return(0)
    }
  }
  tukey_function_derivative <- function(x){
    c = 4.685
    if(abs(x) < c){
      return(1 - 6*(x^2/c^2) + 5*(x^4/c^4))
    }
    else{
      return(0)
    }
  }

  id_function <- function(x){
    return(x)
  }
  id_function_derivative <- function(x){
    return(1)
  }

  # this function will determine which penalization function will be used: Tukey or Huber
  penfun <- function(x, type = c("Huber", "Tukey", "OLS")){
    if(type == "Huber"){return(huber_function(x))}
    if(type == "Tukey"){return(tukey_function(x))}
    if(type == "OLS"){return(id_function(x))}
  }
  penfunderiv <- function(x, type = c("Huber", "Tukey", "OLS")){
    if(type == "Huber"){return(huber_function_derivative(x))}
    if(type == "Tukey"){return(tukey_function_derivative(x))}
    if(type == "OLS"){return(id_function_derivative(x))}
  }


  conf_set_robust = rep(-Inf, length(betagrid))
  p_values =  rep(-Inf, length(betagrid))

  for(i in 1:length(betagrid)){
    beta0 <- betagrid[i]

    restricted_fun <- function(Y, X, Z, W, beta0, type = c("Huber", "Tukey", "OLS"), weights){
      dep_var <- Y - X*beta0
      if(type == "Huber"){return(MASS::rlm(dep_var ~ 0 + W, weights = weights, scale.est = "proposal 2",  maxit = 100, acc = 1e-10))}
      if(type == "Tukey"){return(MASS::rlm(dep_var ~ 0 + W, weights = weights, psi = MASS::psi.bisquare, scale.est = "proposal 2",  maxit = 100, acc = 1e-10))}
      if(type == "OLS"){return(lm(dep_var ~ 0 + W))}
    }

    #the beta0-restricted model object.
    restricted_obj <- restricted_fun(Y = Y, X = X, Z = Z, W = W, beta0 = beta0, type = type, weights = covariate.weights)

    restricted_resid <- as.matrix(residuals(restricted_obj))

    # scale estimates of the beta0-restricted model object.
    if(type == "Huber" | type == "Tukey"){
      restricted_scale <- restricted_obj$s
    }
    if(type == "OLS"){
      restricted_scale <- sigma(restricted_obj)
    }

    # standardized residuals to which the penalization function is applied
    stan_resids_restricted <- apply(restricted_resid/restricted_scale, 1, penfun, type = type)
    stan_resids_restricted_deriv <- apply(restricted_resid/restricted_scale, 1, penfunderiv, type = type) * (1/restricted_scale)

    Z.matrix <- (1/n) * (covariate.weights * stan_resids_restricted) %*% Z

    M.matrix <- matrix(rep(0), nrow = p+k, ncol = p+k)
    Q.matrix <- matrix(rep(0), nrow = p+k, ncol = p+k)


    if(homoskedasticity == "yes"){
      for(j in 1:n){
        M.matrix <- M.matrix + -covariate.weights[j]  * WZ[j,] %*%  t(WZ[j,])
        Q.matrix <- Q.matrix + covariate.weights[j]^2 * WZ[j,] %*%  t(WZ[j,])
      }

      M.matrix <- (1/n) * M.matrix * mean(stan_resids_restricted_deriv)
      Q.matrix <- (1/n) * Q.matrix * (1/(n-k-p)) * sum(stan_resids_restricted^2)
    }

    if(homoskedasticity == "no"){
      for(j in 1:n){
        M.matrix <- M.matrix + covariate.weights[j] * -stan_resids_restricted_deriv[j] * WZ[j,] %*%  t(WZ[j,])
        Q.matrix <- Q.matrix + covariate.weights[j]^2 * stan_resids_restricted[j]^2 * WZ[j,] %*% t(WZ[j,])
      }
      M.matrix <- (1/n) * M.matrix
      Q.matrix <- (1/n) * Q.matrix
    }

    V.matrix <- solve(M.matrix) %*% Q.matrix %*% solve(t(M.matrix))
    V22 <- V.matrix[(p+1):(p+k), (p+1):(p+k)]

    M22 <- M.matrix[(p+1):(p+k), (p+1):(p+k)]
    M21 <- M.matrix[(p+1):(p+k), 1:p]
    M12 <- M.matrix[1:p, (p+1):(p+k)]
    M11 <- M.matrix[1:p, 1:p]
    M22.1 <- M22 - M21 %*% solve(M11) %*% M12

    U.matrix <- M22.1 %*% V22 %*% t(M22.1)

    AR.stat <- Z.matrix %*% solve(U.matrix) %*% t(Z.matrix)
    if(n * AR.stat < qchisq(1 - alpha, k)){
      conf_set_robust[i] <- beta0
    }
    p_values[i] <- pchisq(n * AR.stat, 1, lower.tail = FALSE)
  }
  return(list("conf.set" = conf_set_robust, "pvalues" = p_values))
}


#' Huber function
#'
#' @param x A scalar
#' @returns A scalar
huber_function <- function(x){
  c = 1.345
  if(abs(x) < c){
    return(x)
  }
  else{
    return(c * sign(x))
  }
}
#' Derivative of the Huber function
#'
#' @param x A scalar
#' @returns A scalar
huber_function_derivative <- function(x){
  c = 1.345
  if(abs(x) < c){
    return(1)
  }
  else{
    return(0)
  }
}

#' Tukey's biweight function
#'
#' @param x A scalar
#' @returns A scalar
tukey_function <- function(x){
  c=4.685
  if(abs(x) < c){
    return(x*(1 - x^2/c^2)^2)
  }
  else{
    return(0)
  }
}

#' Derivative of Tukey's biweight function
#'
#' @param x A scalar
#' @returns A scalar
tukey_function_derivative <- function(x){
  c = 4.685
  if(abs(x) < c){
    return(1 - 6*(x^2/c^2) + 5*(x^4/c^4))
  }
  else{
    return(0)
  }
}

#' Identity function
#'
#' @param x A scalar
#' @returns A scalar
id_function <- function(x){
  return(x)
}

#' Derivative of Identity function
#'
#' @param x A scalar
#' @returns A scalar
id_function_derivative <- function(x){
  return(1)
}
