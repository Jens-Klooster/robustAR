library(survey)
library(robustbase)
library(MASS)
library(AER)
library(R.matlab)
library(pim)
library(lmtest)
library(foreign)
library(mvtnorm)
library(svMisc)

source("RobustScoreTest.R")

N = 100
Panels = 3
k = 1

total.obs = N * Panels

covmat = matrix(c(1, 0.75, 0.75, 1), nrow = 2, ncol = 2)
epsV = mvrnorm(total.obs, mu = c(0,0), Sigma = covmat)
eps = epsV[,1]
V = epsV[,2]

# distributional contamination
# epsCon = rmvt(n, sigma = covmat, df = 2)
# eps[1:25] = epsCon[,1][1:25]
# V[1:25] = epsCon[,2][1:25]


Z = matrix(rnorm(total.obs*k), nrow = total.obs, ncol = k)

a = 1
Pi = a
X = Z %*% Pi + V
Y = X*0 + eps
Y[1:N] = 1 + Y[1:N]
Y[(N+1):(2*N)] = 10 + Y[(N+1):(2*N)]
Y[(2*N + 1): (3*N)] = 20 + Y[(2*N + 1): (3*N)]
plot(Y)

Y_ <- rep(0, total.obs)
Y_[1:N] <- Y[1:N] - median(Y[1:N])
Y_[(N+1):(2*N)] <- Y[(N+1):(2*N)] - median(Y[(N+1):(2*N)])
Y_[(2*N + 1): (3*N)] <- Y[(2*N + 1): (3*N)] - median(Y[(2*N + 1): (3*N)])
plot(Y_)

summary(lm(Y_ ~ X))

X <- as.numeric(X)
RobustScoreTest(lm(Y_ ~ 1 + X),
                lm(Y_ ~ 1),
                lmrob(Y_ ~ X)$scale,
                rep(1, total.obs))
summary(lm(Y_ ~ X))
