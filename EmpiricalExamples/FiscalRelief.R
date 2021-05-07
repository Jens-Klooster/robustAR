library(robustbase)
library(robustX)
library(readstata13)
library(ivmodel) # need this for AR test
library(riv) # robust IV estimator

source("RobustScoreTest.R")

fiscal.data <- read.dta13("EmpiricalExamples/FiscalRelief.dta")
fiscal.data <- fiscal.data[-outliers,]

for(i in 1:length(names(fiscal.data))){
  name = names(fiscal.data)[i]
  assign(name, fiscal.data[[i]])
}

Y <- sachange_gov_broad_pc
X <- fmap_pc
Z <- instrument_pc

# Use Bacon Algorithm to find outliers.
BC <- BACON(cbind(Y, X, Z),
            alpha = 0.01)
outliers <- which(BC$subset*1 == 0)
bacon.weights <- BC$subset*1
pairs(cbind(Y, X, Z),
      labels = c("y2", "y1", "x2"))

pairs.data <- data.frame(Y, X, Z, bacon.weights)

group <- NA
group[pairs.data$bacon.weights == 0] <- 1
group[pairs.data$bacon.weights != 0] <- 2
group

# Makes a plot where the outliers are in red
pairs(pairs.data[, 1:3],
      col = c("red", "cornflowerblue", "darkred")[group],
      pch = c(8, 1)[group],
      labels = c("y2", "y1", "x2"))


# use this to find first stage F-statistic
summary(lm(X ~ 1 + Z))



beta_grid = seq(-2, 3, 0.01)
conf_set_robust = rep(-Inf, length(beta_grid))

a1<-lmrob.control()
a1$k.max <- 100000
a1$maxit.scale <- 100000
a1$max.it <- 500

gewichtjes <- sqrt(1 - diag(Z %*% solve(t(Z) %*% Z) %*% t(Z)))

p_values = rep(1, length(beta_grid))

for(i in 1:length(beta_grid)){
  dep_var = Y - X * beta_grid[i]
  
  full_model <- dep_var ~ 1 + Z
  res_model <- dep_var ~ 1

  rob_model <- lmrob(formula = full_model, control = a1, method = "MM")
  estScale <- rob_model$scale
  
  
  tau_stat = RobustScoreTest(full_model, 
                             res_model, 
                             estScale, 
                             gewichtjes, 
                             tukey.c = 4.68)
  
  if(tau_stat$W.statistic < qchisq(0.90,1)){
    conf_set_robust[i] <- beta_grid[i]
  }
  
  p_values[i] <- tau_stat$p.value
  if(tau_stat$W.statistic < qchisq(0.95,1)){
    conf_set_robust[i] <- beta_grid[i]
  }
}
conf_set_robust
AR.test(ivmodel(Y, X, Z), alpha = 0.1)$ci

riv(Y, X, Zinst = Z)$Summary.Table
ivmodel(Y, X, Z)


plot(p_values)
lines(rep(0.1, length(beta_grid)))

beta_grid[which(p_values == max(p_values))]

max(p_values)

