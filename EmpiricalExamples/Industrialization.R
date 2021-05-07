library(robustbase) 
library(readstata13) # need this to load in stata dta file
library(ivmodel) # need this for AR test
library(robustX) # need this for Bacon algorithm
library(riv) # robust IV estimator
source("RobustScoreTest.R")

indus.data <- read.dta13("EmpiricalExamples/Industrialization.dta")
# indus.data <- indus.data[-c(270),]
indus.data <- indus.data[-outliers,]

for(i in 1:length(names(indus.data))){
  name = names(indus.data)[i]
  assign(name, indus.data[[i]])
}

Y = fac1849_total_pc
X = edu1849_adult_yos
Z = edu1816_pri_enrol
W = cbind(1, pop1849_young, pop1849_old, area1816_qkm)

# PW = W %*% solve(t(W) %*% W) %*% t(W)
# MW = diag(1, nrow = length(Y), ncol = length(Y)) - PW
# Y_ = MW %*% Y
# X_ = MW %*% X
# Z_ = MW %*% Z


BC <- BACON(cbind(Y, X, Z, cbind(pop1849_young, pop1849_old, area1816_qkm)),
            alpha = 0.01)
outliers <- which(BC$subset*1 == 0)
bacon.weights <- BC$subset*1
pairs(cbind(Y, X, Z),
      col = "cornflowerblue",
      labels = c("y2", "y1", "x2", "x1"))

pairs.data <- data.frame(Y, X, Z, bacon.weights)

group <- NA
group[pairs.data$bacon.weights == 0] <- 1
group[pairs.data$bacon.weights != 0] <- 2
group

pairs(pairs.data[, 1:3],
      col = c("red", "cornflowerblue", "darkred")[group],
      pch = c(8, 1)[group],
      labels = c("y2", "y1", "x2"))

summary(lm(X_ ~ 1 + Z_))





a1<-lmrob.control()
a1$maxit.scale<-100000
a1$k.max <- 100000



lev_mat <- cbind(edu1816_pri_enrol, pop1849_young, pop1849_old, area1816_qkm)
H_mat <- lev_mat %*% solve(t(lev_mat) %*% lev_mat) %*% t(lev_mat)
gewichtjes <- sqrt(1 - diag(H_mat))


beta_grid = seq(-0.1, 0.3, 0.01)
conf_set_robust = rep(-Inf, length(beta_grid))

p_values = rep(1, length(beta_grid))

for(i in 1:length(beta_grid)){
   dep_var <- Y - X * beta_grid[i]
   


   full_model <- lm(dep_var ~ 1 + pop1849_young + pop1849_old + area1816_qkm + edu1816_pri_enrol)
   res_model <- lm(dep_var ~ 1 + pop1849_young + pop1849_old + area1816_qkm)
   
   rob_model <- lmrob(formula = full_model, control = a1, method = "MM")
   res_model <- lmrob(formula = res_model, control = a1, weights = gewichtjes, method = "MM")
   estScale <- rob_model$scale
   
   
   
   tau_stat = RobustScoreTest(full_model, 
                              res_model, 
                              estScale, 
                              gewichtjes, 
                              tukey.c = 4.68)
   
   if(tau_stat$W.statistic < qchisq(0.9,1)){
      conf_set_robust[i] <- beta_grid[i]
   }
   
   p_values[i] <- tau_stat$p.value
   if(tau_stat$W.statistic < qchisq(0.95,1)){
      conf_set_robust[i] <- beta_grid[i]
   }
   
}
conf_set_robust
AR.test(ivmodel(Y, X, Z, W), alpha = 0.1)$ci
riv(Y, X, Xex = W[,2:4], Zinst = Z)$Summary.Table

ivmodel(Y, X, Z, W)

plot(p_values)
lines(rep(0.1, length(beta_grid)))

beta_grid[which(p_values == max(p_values))]
