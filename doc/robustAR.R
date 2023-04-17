## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(robustAR)

## -----------------------------------------------------------------------------
set.seed(1)
beta <- 1
n <- 250
eps1 <- rnorm(n)
eps2 <- 0.5*eps1 + 0.5*rnorm(n)
W <- rnorm(n)
Z <- rnorm(n) 
X <- 0.15 * Z + 1 * W + eps1
Y <- beta * X + 2 * W + eps2

## -----------------------------------------------------------------------------
anova(lm(X ~ 0 + W), lm(X ~ 0 + W + Z))

## -----------------------------------------------------------------------------
robustAR(Y = Y, X = X, Z = Z, W = W, beta0 = 1, 
         type = "Tukey", 
         weighting = "hat")
AR(Y = Y, X = X, Z = Z, W = W, beta0 = 1)

## -----------------------------------------------------------------------------
Y.corrupted <- Y
X.corrupted <- X
Z.corrupted <- Z
W.corrupted <- W
Y.corrupted[1] <- 2
X.corrupted[1] <- 5
Z.corrupted[1] <- 10
W.corrupted[1] <- 2

## -----------------------------------------------------------------------------
AR(Y = Y.corrupted, X = X.corrupted, Z = Z.corrupted, W = W.corrupted, beta0 = 1)
robustAR(Y = Y.corrupted, X = X.corrupted, Z = Z.corrupted, W = W.corrupted, beta0 = 1, 
         type = "Tukey", 
         weighting = "hat")

## -----------------------------------------------------------------------------
alpha = 0.05
betagrid <- seq(-6, 4, 0.1)
AR.confset <- rep(-Inf, length(betagrid))
for(i in 1:length(betagrid)){
  beta0 <- betagrid[i]
  AR.obj <- AR(Y = Y, X = X, Z = Z, W = W, beta0 = beta0)
  if(AR.obj < qchisq(0.95, 1)){
    AR.confset[i] <- beta0
  }
}
AR.confset

## -----------------------------------------------------------------------------
robustAR.confset <- robustAR.conf(Y = Y, X = X, Z = Z, W = W, betagrid = betagrid, alpha = 0.05,
                                  type = "Tukey",
                                  weighting = "hat")
robustAR.confset$conf.set

## -----------------------------------------------------------------------------
AR.confset.corrupted <- rep(-Inf, length(betagrid))
for(i in 1:length(betagrid)){
  beta0 <- betagrid[i]
  AR.obj <- AR(Y = Y.corrupted, X = X.corrupted, Z = Z.corrupted, W = W.corrupted, beta0 = beta0)
  if(AR.obj < qchisq(0.95, 1)){
    AR.confset.corrupted[i] <- beta0
  }
}
AR.confset.corrupted

robustAR.confset.corrupted <- robustAR.conf(Y = Y.corrupted, 
                                  X = X.corrupted, 
                                  Z = Z.corrupted, 
                                  W = W.corrupted, 
                                  betagrid = betagrid, 
                                  alpha = 0.05,
                                  type = "Tukey",
                                  weighting = "hat")
robustAR.confset.corrupted$conf.set

## -----------------------------------------------------------------------------
anova(lm(X.corrupted ~ 0 + W.corrupted), lm(X.corrupted ~ 0 + W.corrupted + Z.corrupted))

## -----------------------------------------------------------------------------
AR(Y = Y, X = X, Z = Z, W = cbind(1,W), beta0 = 1)
robustAR(Y = Y, X = X, Z = Z, W = cbind(1,W), beta0 = 1, 
         type = "Tukey", 
         weighting = "hat")

