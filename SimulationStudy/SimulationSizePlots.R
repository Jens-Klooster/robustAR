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
n = 250
k = 1
#e_10 = c(1, rep(1,k-1))
e_10 = c(1)
a_list = c(0.1, 1)
rho_list = seq(0, 0.95, 0.1)



outcome_matrix = matrix(rep(0, length(a_list) * length(rho_list)), nrow = length(a_list), ncol = length(rho_list))
colnames(outcome_matrix) <- rho_list
rownames(outcome_matrix) <- a_list

outcome_matrix_AR = matrix(rep(0, length(a_list) * length(rho_list)), nrow = length(a_list), ncol = length(rho_list))
colnames(outcome_matrix_AR) <- rho_list
rownames(outcome_matrix_AR) <- a_list

outcome_matrix_W = matrix(rep(0, length(a_list) * length(rho_list)), nrow = length(a_list), ncol = length(rho_list))
colnames(outcome_matrix_W) <- rho_list
rownames(outcome_matrix_W) <- a_list


sim_amount = 10000

for(i in 1:length(a_list)){
  a = a_list[i]
  for(j in 1:length(rho_list)){
    rho = rho_list[j]
    for(N in 1:sim_amount){
      covmat = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
      epsV = mvrnorm(n, mu = c(0,0), Sigma = covmat)
      eps = epsV[,1]
      V = epsV[,2]
      # distributional contamination
      # epsCon = rmvt(n, sigma = covmat, df = 2)
      # eps[1:25] = epsCon[,1][1:25]
      # V[1:25] = epsCon[,2][1:25]
      


      Z = matrix(rnorm(n*k), nrow = n, ncol = k)
      Pi = a*e_10
      X = Z %*% Pi + V
      Y = X*0 + eps
      
      # Outliers
      # Y[1] = -4
      # Z[1,1] = 4
      

      
      z1 = Z[,1]
      
      PZ = Z %*% solve(t(Z) %*% Z) %*% t(Z)
      MZ = diag(1, nrow = n, ncol = n) - PZ

      gewichtjes <- sqrt(1 - diag(PZ))



      beta_2sls = solve(t(X) %*% PZ %*% X) %*% (t(X) %*% PZ %*% Y)
      s_ee = (1/(n-1))*(t(Y - X %*% beta_2sls) %*% (Y - X %*% beta_2sls))
      var_beta_2sls = s_ee * solve(t(X) %*% PZ %*% X)
      t_stat = (beta_2sls - 0)/sqrt(var_beta_2sls)
      # beta_2sls
      # sqrt(s_ee * solve(t(X) %*% PZ %*% X))
      # t_stat
      if(abs(t_stat) > 1.96){
        outcome_matrix[i,j] <- outcome_matrix[i,j] + 1/sim_amount
      }


      # Now the Anderson Rubin test

      AR_stat = ((t(Y) %*% PZ %*% Y)/k)/((t(Y) %*% MZ %*% Y)/(n - k))
      if(AR_stat > qchisq(0.95, k)/k ){
        outcome_matrix_AR[i,j] <- outcome_matrix_AR[i,j] + 1/sim_amount
      }

      full_model <-  lm(Y ~ 1 + z1)
      res_model  <-  lm(Y ~ 1)
      
      a1<-lmrob.control()
      a1$k.max <- 10000
      
      # rob_model <- lmrob(formula = full_model, control = a1, method = "MM")
      # estScale <- res_model$scale

      W_stat = RobustScoreTest(full_model,
                                res_model, 
                                lmrob(Y ~ 1 + z1, weights = gewichtjes)$scale,
                                gewichtjes,
                                tukey.c = 4.68)
      #tau_stat = MMTauTest(full_model, res_model, scale.est.err = estScale^2, tukeyC = 4.68, OptimMethod = "Nelder-Mead")
      if(W_stat$W.statistic > qchisq(0.95, k)){
        outcome_matrix_W[i,j] <- outcome_matrix_W[i,j] + 1/sim_amount
      }
    }
    progress(j, length(rho_list))
  }
}


pairs(cbind(Y,X, Z),
      col = "cornflowerblue",
      labels = c("y2", "y1", "x2"))

outcome_matrix
outcome_matrix_W
outcome_matrix_AR

par(mfrow=c(1,2))
plot(outcome_matrix[1,],
     ylab = "Rejection Frequency",
     xlab = expression(paste(rho)),
     xaxt = "n",
     type = "l",
     main = expression(paste(pi, " = 0.1")),
     ylim = c(0,0.15),
     lty=3)
axis(1, at = 1:length(rho_list), label=rho_list)
lines(outcome_matrix_AR[1,],
      lty=2)
lines(outcome_matrix_W[1,],
      lty=1)
# legend("topleft", 
#        legend=c("t-test", "AR-test", "Robust-AR"),
#        lty=c(1,2,3))


plot(outcome_matrix[2,], ylab = "Rejection Frequency", xlab = expression(paste(rho)), xaxt = "n", type = "l", main = expression(paste(pi, " = 1")),
     ylim = c(0,0.15),
     lty=3)
axis(1, at = 1:length(rho_list), label=rho_list)
lines(outcome_matrix_AR[2,], lty=2)
lines(outcome_matrix_W[2,], lty=1)
# legend("topleft", 
#        legend=c("t-test", "AR-test", "Robust-AR"),
#        lty=c(1,2,3))



