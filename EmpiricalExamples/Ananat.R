library(robustbase)
library(readstata13)
library(ivmodel)
library(robustX)
library(riv)
source("RobustScoreTest.R")

ananat.data <- read.dta13("EmpiricalExamples/aej_maindata.dta")
# ananat.data <- ananat.data[-outliers,]
# ananat.data <- ananat.data[-c(6),]


# Here we load in the variables
for(i in 1:length(names(ananat.data))){
  name = names(ananat.data)[i]
  assign(name, ananat.data[[i]])
}

#Here we use the general IV notation
Y <- lngini_w
# Y <- povrate_w
Z <- herf
X <- dism1990
W <- cbind(1, lenper)

L <- c(lenper)

# MW <- (diag(1, length(Y)) - W %*% solve(t(W) %*% W) %*% t(W))
# 
# Z <- MW %*% Z
# X <- MW %*% X
# Y <- MW %*% Y

# outlier detection
BC <- BACON(cbind(Y, X, Z, lenper),
            alpha = 0.01)
outliers <- which(BC$subset*1 == 0)
bacon.weights <- BC$subset*1
pairs(cbind(Y, X, Z, lenper),
      col = "cornflowerblue",
      labels = c("y2", "y1", "x2", "x1"))

pairs.data <- data.frame(Y, X, Z, lenper, bacon.weights)

group <- NA
group[pairs.data$bacon.weights == 0] <- 1
group[pairs.data$bacon.weights != 0] <- 2
group

pairs(pairs.data[, 1:4],
      col = c("red", "cornflowerblue", "darkred")[group],
      pch = c(8, 1)[group],
      labels = c("y2", "y1", "x2", "x1"))
      

# a1<-lmrob.control()
# a1$k.max <- 200000
# # a1$maxit.scale <- 100000
# # a1$max.it = 50000

H_matrix <- cbind(1, Z, L)
gewichtjes <- sqrt(1 - diag(H_matrix %*% solve(t(H_matrix) %*% H_matrix) %*% t(H_matrix)))^6
# gewichtjes <- sqrt(1 - diag(H_matrix %*% solve(t(H_matrix) %*% H_matrix) %*% t(H_matrix)))
# gewichtjes <- rep(1, length(herf))
# gewichtjes[6] = 0
# gewichtjes[22] <- 0

# data_mat <- cbind(Z,L)
# covMcd.obj <- covMcd(data_mat)
# rob_mean <- covMcd.obj$center
# rob_cov <- covMcd.obj$cov
# gewichtjes <- pmin(1, 3.84/sqrt(mahalanobis(data_mat, rob_mean, rob_cov)))
# gewichtjes[6] <- 0


# we can use this to find the first stage F statistic
f.stat = summary(lm(dism1990 ~ 1 + lenper + herf))


robust.f.stat = RobustScoreTest(dism1990 ~ 1 + lenper + herf,
                               dism1990 ~ 1 + lenper,
                               lmrob(dism1990 ~ 1 + lenper)$scale,
                               gewichtjes)




# Now compute the confidence interval
beta_grid = seq(-10, 10, 0.1)
conf_set_robust = rep(-Inf, length(beta_grid))
p_values = rep(1, length(beta_grid))

control.variables <-lmrob.control()
control.variables$k.max <- 100000
control.variables$refine.tol <- 1e-10

for(i in 1:length(beta_grid)){
  dep_var = Y - beta_grid[i]*X
  
  full_model <-  dep_var ~ 1 + lenper + Z
  res_model  <-  dep_var ~ 1 + lenper
  
  tau_stat = RobustScoreTest(full_model, 
                            res_model, 
                            lmrob(dep_var ~ 1 + lenper + Z,
                                  weights = gewichtjes,
                                  control = control.variables)$scale,
                            # rlm(dep_var ~ 1 + lenper + Z,
                            #     weights = as.vector(gewichtjes))$s,
                            gewichtjes, 
                            tukey.c = 4.68)
  
  p_values[i] <- tau_stat$p.value
  if(tau_stat$W.statistic < qchisq(0.95,1)){
    conf_set_robust[i] <- beta_grid[i]
  }
}
conf_set_robust
AR.test(ivmodel(Y, X, Z, L))$ci

ivmodel(Y, X, Z, L)
riv(Y, X, Xex = W[,2], Zinst = Z)$Summary.Table

plot(p_values)
lines(rep(0.05, length(beta_grid)))

beta_grid[which(p_values == max(p_values))]



