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
source("RobustWaldTest.R")
n <- 250
k = 1
e_10 = c(1)
a_list = c(0.1, 1)
beta_list = seq(-0.3, 0.3, 0.01)

outcome_matrix = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix) <- beta_list

outcome_matrix_AR = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_AR) <- beta_list

outcome_matrix_W = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_W) <- beta_list

outcome_matrix_LM = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_LM) <- beta_list

outcome_matrix_LM_rob = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_LM_rob) <- beta_list

beta_list_weak = seq(-3, 3, 0.1)

outcome_matrix_weak = matrix(rep(0, length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_weak) <- beta_list_weak

outcome_matrix_AR_weak = matrix(rep(0, length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_AR_weak) <- beta_list_weak

outcome_matrix_W_weak = matrix(rep(0, length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_W_weak) <- beta_list_weak


sim_amount = 1000

set.seed(9)
for(N in 1:sim_amount){
  covmat = matrix(c(1, 0.25, 0.25, 1), nrow = 2, ncol = 2)
  epsV = mvrnorm(n, mu = c(0,0), Sigma = covmat)
  # epsV = rmvt(n, covmat, df = 2)
  eps = epsV[,1]
  
  V = epsV[,2]
  # distributional contamination
  # epsCon = rmvt(n, sigma = covmat, df = 2)
  # eps[1:25] = epsCon[,1][1:25]
  # V[1:25] = epsCon[,2][1:25]
  
  
  Z = matrix(rnorm(n*k), nrow = n, ncol = k)
  
 
  
  
  # Outliers:
  Z[1,1] = 4

  
  
  
  z1 = Z[,1]
  # z2 = Z[,2]
  # z3 = Z[,3]
  
  
  PZ = Z %*% solve(t(Z) %*% Z) %*% t(Z)
  MZ = diag(1, nrow = n, ncol = n) - PZ
  

  gewichtjes <- sqrt(1 - diag(PZ))
  for(number in a_list){
    a = number
    if(a == 1){
      Pi = a*e_10
      X = Z %*% Pi + V
      for(j in 1:length(beta_list)){
        beta = beta_list[j]
        Y = X*beta + eps
        
        #outlier
        Y[1] = -4
       
        
        beta_2sls = solve(t(X) %*% PZ %*% X) %*% (t(X) %*% PZ %*% Y)
        s_ee = (1/(n-1))*(t(Y - X %*% beta_2sls) %*% (Y - X %*% beta_2sls))
        var_beta_2sls = s_ee * solve(t(X) %*% PZ %*% X)
        t_stat = (beta_2sls - 0)/sqrt(var_beta_2sls)
        if(abs(t_stat) > 1.96){
          outcome_matrix[1,j] <- outcome_matrix[1,j] + 1/sim_amount
        }
        
        
        # Now the Anderson Rubin test
        
        AR_stat = ((t(Y) %*% PZ %*% Y)/k)/((t(Y) %*% MZ %*% Y)/(n - k))
        if(AR_stat > qchisq(0.95, k)/k ){
          outcome_matrix_AR[1,j] <- outcome_matrix_AR[1,j] + 1/sim_amount
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
                                weight.vector = gewichtjes,
                                tukey.c = 4.68)
        #tau_stat = MMTauTest(full_model, res_model, scale.est.err = estScale^2, tukeyC = 4.68, OptimMethod = "Nelder-Mead")
        if(W_stat$W.statistic > qchisq(0.95, k)){
          outcome_matrix_W[1,j] <- outcome_matrix_W[1,j] + 1/sim_amount
        }
      }
    }
    if(a == 0.1){
      Pi = a*e_10
      X = Z %*% Pi + V
      for(j in 1:length(beta_list_weak)){
        beta = beta_list_weak[j]
        Y = X*beta + eps
        
        # outlier
        Y[1] = -4
        
        beta_2sls = solve(t(X) %*% PZ %*% X) %*% (t(X) %*% PZ %*% Y)
        s_ee = (1/(n-1))*(t(Y - X %*% beta_2sls) %*% (Y - X %*% beta_2sls))
        var_beta_2sls = s_ee * solve(t(X) %*% PZ %*% X)
        t_stat = (beta_2sls - 0)/sqrt(var_beta_2sls)
        if(abs(t_stat) > 1.96){
          outcome_matrix_weak[1,j] <- outcome_matrix_weak[1,j] + 1/sim_amount
        }
        
        
        # Now the Anderson Rubin test
        
        AR_stat = ((t(Y) %*% PZ %*% Y)/k)/((t(Y) %*% MZ %*% Y)/(n - k))
        if(AR_stat > qchisq(0.95, k)/k ){
          outcome_matrix_AR_weak[1,j] <- outcome_matrix_AR_weak[1,j] + 1/sim_amount
        }
        
        full_model <-  lm(Y ~ 1 + z1)
        res_model  <-  lm(Y ~ 1)
        
        a1<-lmrob.control()
        a1$k.max <- 10000
        
        rob_model <- lmrob(formula = full_model, control = a1, method = "MM")
        estScale <- res_model$scale
        
        W_stat = RobustScoreTest(full_model,
                                res_model, 
                                lmrob(Y ~ 1 + z1, weights = gewichtjes)$scale,
                                gewichtjes,
                                tukey.c = 4.68)
        #tau_stat = MMTauTest(full_model, res_model, scale.est.err = estScale^2, tukeyC = 4.68, OptimMethod = "Nelder-Mead")
        if(W_stat$W.statistic > qchisq(0.95, k)){
          outcome_matrix_W_weak[1,j] <- outcome_matrix_W_weak[1,j] + 1/sim_amount
        }
      }
    }
  }
  progress(N, sim_amount)
}

par(mfrow=c(1,2))
plot(outcome_matrix_weak[1,],
     ylab = "Rejection Frequency",
     xlab = expression(paste(beta)),
     xaxt = "n",
     type = "l",
     main = expression(paste(pi, " = 0.1")),
     ylim = c(0,0.8),
     lty=3)
axis(1, at = c(1, 11, 21, 31, 41, 51, 61), label=c("-3", "-2", "-1", "0", "1", "2", "3"))
lines(outcome_matrix_AR_weak[1,],
      lty=2)
lines(outcome_matrix_W_weak[1,],
      lty=1)
# legend("bottomright", 
#        legend=c("t-test", "AR-test", "Robust-AR"),
#        lty=c(3,2,1),
#        bty="n")


plot(outcome_matrix[1,], 
     ylab = "Rejection Frequency", 
     xlab = expression(paste(beta)), 
     xaxt = "n", 
     type = "l", 
     main = expression(paste(pi, " = 1")),
     ylim = c(0,1),
     lty=3)
axis(1, at = c(1, 11, 21, 31, 41, 51, 61), label=c("-0.3", "-0.2", "-0.1", "0", "0.1", "0.2", "0.3"))
lines(outcome_matrix_AR[1,], lty=2)
lines(outcome_matrix_W[1,], lty=1)
# legend("bottomright", 
#        legend=c("t-test", "AR-test", "Robust-AR"),
#        lty=c(3,2,1),
#        bty="n")







#########################################################################
#########################################################################
#########################################################################

# Now we do the same, but for the LM test.


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
library(ivmodel)

source("RobustWaldTest.R")
source("RobustScoreTest.R")
set.seed(9)
n = 300
k = 3
#e_10 = c(1, rep(1,k-1))
#e_10 = c(1)
#a_list = c(0.1, 1)
#beta_list = c(0, 0.1)
#beta_list = c(-0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15)
#beta_list = c(-0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
beta_list = seq(-0.3, 0.3, 0.02)
# beta_list = c(0)

# We specificy a grid for the LR critical values

lr_grid <- 0:1000
grid_simul <- 10000
grid_matrix <- matrix(rep(0, grid_simul * length(lr_grid)), nrow = grid_simul, ncol = length(lr_grid))
colnames(grid_matrix) <- lr_grid
for(g in 1:length(lr_grid)){
  for(s in 1:grid_simul){
    chi1 <- rchisq(1,1)
    chi2 <- rchisq(1,k-1)
    grid_matrix[s,g] <- 0.5*(chi2 + chi1 - lr_grid[g] +
                               sqrt((chi2 + chi1 + lr_grid[g])^2 - 4*lr_grid[g]*chi2))
  }
}

sorted_grid_matrix <- apply(grid_matrix, 2, sort)


outcome_matrix = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix) <- beta_list

outcome_matrix_AR = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_AR) <- beta_list

outcome_matrix_W = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_W) <- beta_list

outcome_matrix_LM = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_LM) <- beta_list

outcome_matrix_LM_rob = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_LM_rob) <- beta_list

outcome_matrix_LR = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_LR) <- beta_list

outcome_matrix_LR_rob = matrix(rep(0, 1 * length(beta_list)), nrow = 1, ncol = length(beta_list))
colnames(outcome_matrix_LR_rob) <- beta_list





rho = 0.25
sim_amount = 200
a = 0.75
e_10 = c(a, rep(a,k-1))



for(N in 1:sim_amount){
  covmat = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  epsV = mvrnorm(n, mu = c(0,0), Sigma = covmat)
  eps = epsV[,1]
  V = epsV[,2]
  
  # distributional contamination
  epsCon = rmvt(n, sigma = covmat, df = 2)
  eps[1:30] = epsCon[,1][1:30]
  V[1:30] = epsCon[,2][1:30]
  
  Z = matrix(rnorm(n*k), nrow = n, ncol = k)
  
  # Z[1,1] <- -7
  # Z[1,2] <- -4
  # Z[1,3] <- -4
  
  
  z1 <- Z[,1]
  z2 <- Z[,2]
  z3 <- Z[,3]
  # z4 <- Z[,4]
  # z5 <- Z[,5]
  # z6 <- Z[,6]
  # z7 <- Z[,7]
  # z8 <- Z[,8]
  # z9 <- Z[,9]
  # z10 <- Z[,10]
  
  PZ = Z %*% solve(t(Z) %*% Z) %*% t(Z)
  MZ = diag(1, nrow = n, ncol = n) - PZ
  
  
  
  Pi = a*e_10
  X = Z %*% Pi + V
  
  for(j in 1:length(beta_list)){
    beta = beta_list[j]
    Y = X*beta + eps
    
    
    
    # Y[1] = 7
    
    
    
    # beta_2sls = solve(t(X) %*% PZ %*% X) %*% (t(X) %*% PZ %*% Y)
    # s_ee = (1/(n-1))*(t(Y - X %*% beta_2sls) %*% (Y - X %*% beta_2sls))
    # var_beta_2sls = s_ee * solve(t(X) %*% PZ %*% X)
    # t_stat = (beta_2sls - 0)/sqrt(var_beta_2sls)
    # if(abs(t_stat) > 1.96){
    #   outcome_matrix[1,j] <- outcome_matrix[1,j] + 1/sim_amount
    # }
    
    
    # Now the Anderson Rubin test
    
    AR_stat = ((t(Y) %*% PZ %*% Y)/k)/((t(Y) %*% MZ %*% Y)/(n - k))
    if(AR_stat > qchisq(0.95, k)/k ){
      outcome_matrix_AR[1,j] <- outcome_matrix_AR[1,j] + 1/sim_amount
    }
    
    full_model <-  lm(Y ~ 1 + z1 + z2 + z3)
    res_model  <-  lm(Y ~ 1)
    
    a1<-lmrob.control()
    a1$k.max <- 10000
    
    rob_model <- lmrob(formula = full_model, control = a1, method = "MM")
    estScale <- res_model$scale
    
    W_weights <- sqrt(1 - diag(PZ))
    
    W_stat = RobustScoreTest(full_model,
                            res_model, 
                            lmrob(Y ~ 1 + z1 + z2 + z3)$scale,
                            W_weights,
                            tukey.c = 4.68)$W.statistic
    if(W_stat > qchisq(0.95, k)){
      outcome_matrix_W[1,j] <- outcome_matrix_W[1,j] + 1/sim_amount
    }
    
    
    
    # Now the LM test
    beta0 <- 0
    eps0 <- Y - X*beta0
    sigma_hat_ee <- (1/(n-k))*(t(eps0) %*% MZ %*% eps0)
    sigma_hat_eV <- (1/(n-k))*(t(eps0) %*% MZ %*% X)
    rho_hat <- sigma_hat_eV/sigma_hat_ee

    pi_tilde <- solve(t(Z) %*% Z) %*% t(Z) %*% (X - eps0 %*% rho_hat)
    PZpi <- Z %*% pi_tilde
    projPZpi <- PZpi %*% solve(t(PZpi) %*% PZpi) %*% t(PZpi)

    LM_stat <- (1/sigma_hat_ee)*(t(eps0) %*% projPZpi %*% eps0)
    if(LM_stat > qchisq(0.95,1)){
      outcome_matrix_LM[1,j] <- outcome_matrix_LM[1,j] + 1/sim_amount
    }
    
    
    first.pi = unname(coefficients(lmrob(X ~ 0 + Z)))
    V_hat = unname(residuals(lmrob(X ~ 0 + Z)))
    eps0_hat = unname(residuals(lmrob(eps0 ~ 0 + Z))) 
    cov_eps0.V = (1/(n-k))*(t(eps0_hat) %*% V_hat)
    cov_eps0eps0 = (1/(n-k))*t(eps0_hat) %*% eps0_hat
    rho_hat1 = cov_eps0.V/cov_eps0eps0
    
    X.transformed <- X - eps0 %*% rho_hat1
    pi_tilde_rob <- coefficients(lmrob(X.transformed ~ 0 + Z))
    
    Z.pi_rob <- Z %*% pi_tilde_rob
    model_res <- lm(eps0 ~ 1)
    model_full <- lm(eps0 ~ 1 + Z.pi_rob)
    gewichtjes <- sqrt(1 - diag(Z.pi_rob %*% solve(t(Z.pi_rob) %*% Z.pi_rob) %*% t(Z.pi_rob)))
    
    rob_lm <- RobustScoreTest(model_full,
                             model_res,
                             lmrob(eps0 ~ Z)$scale,
                             gewichtjes,
                             4.68)$W.statistic 
    
    if(rob_lm > qchisq(0.95,1)){
      outcome_matrix_LM_rob[1,j] <- outcome_matrix_LM_rob[1,j] + 1/sim_amount
    }
    
    
    sigma_VV = (1/(n-k))*(t(X) %*% MZ %*% X)
    sigma_VV.e = sigma_VV - sigma_hat_eV^2/sigma_hat_ee
    r_beta0 = (1/sigma_VV.e) * t(PZpi) %*% PZpi
    # Now the LR test
    LR_stat = 0.5*(k*AR_stat - r_beta0 + sqrt((k*AR_stat + r_beta0)^2 - 4*r_beta0 *(k*AR_stat - LM_stat)))
    if(LR_stat > sorted_grid_matrix[0.95*grid_simul, 1 + min(max(lr_grid), round(r_beta0))]){
      outcome_matrix_LR[1,j] <- outcome_matrix_LR[1,j] + 1/sim_amount
    }
    
    # LR_stat = CLR(ivmodel(Y, X, Z))$p.value
    # if(LR_stat < 0.05){
    #   outcome_matrix_LR[1,j] <- outcome_matrix_LR[1,j] + 1/sim_amount
    # }
    
    # Now we try the robust LR test
    
    rob_AR <- W_stat
    rob_LM <- rob_lm
    
    sigma_hat_VV_rob <- (1/(n-k)) * (t(V_hat) %*% V_hat)
    sigma_hat_VV.e_rob <- sigma_hat_VV_rob - cov_eps0.V^2/cov_eps0eps0
    full_mod_r <- X.transformed ~ 1 + z1 + z2 + z3
    res_mod_r <- X.transformed ~ 1
    r_beta0_rob <- RobustWaldTest(full_mod_r,
                                  res_mod_r,
                                  lmrob(full_mod_r)$scale,
                                  #as.numeric(sigma_hat_VV.e_rob),
                                  W_weights,
                                  5.51)$W.statistic
    #r_beta0_rob <- (1/sigma_hat_VV.e_rob) * (t(Z.pi_rob) %*% Z.pi_rob)
    
    rob_LR <- 0.5*(rob_AR - r_beta0_rob + sqrt((rob_AR + r_beta0_rob)^2 -4*r_beta0_rob*(rob_AR - rob_LM)))
    if(rob_LR > sorted_grid_matrix[0.95*grid_simul, 1 + min(max(lr_grid), round(r_beta0_rob))]){
      outcome_matrix_LR_rob[1,j] <- outcome_matrix_LR_rob[1,j] + 1/sim_amount
    }
  }
  progress(N, sim_amount)
}

par(mfrow=c(1,2))
plot(outcome_matrix_LR[1,], ylab = "Rejection Frequency", xlab = expression(paste(beta)), xaxt = "n", type = "l", main = expression(paste(pi, " = 0.5")),
     ylim = c(0,1),
     lty=3)
# axis(1, at = 1:length(beta_list), label=beta_list)
axis(1, at = c(1, 6, 11, 16, 21, 26, 31), label=c("-0.3", "-0.2", "-0.1", "0", "0.1", "0.2", "0.3"))
# lines(outcome_matrix_LM[1,], lty = 2)
# lines(outcome_matrix_LM_rob[1,], lty = 1)
lines(outcome_matrix_LR[1,],lty=2, col = "red")
lines(outcome_matrix_LR_rob[1,], lty = 1, col = "blue")

# legend("bottomright", 
#        legend=c("t-test", "AR-test", "Robust-AR"),
#        lty=c(3,2,1),
#        bty="n")



beta_list_weak = seq(-3, 3, 0.2)
# beta_list_weak = c(0)

outcome_matrix_weak = matrix(rep(0, length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_weak) <- beta_list_weak

outcome_matrix_AR_weak = matrix(rep(0, length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_AR_weak) <- beta_list_weak

outcome_matrix_W_weak = matrix(rep(0, length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_W_weak) <- beta_list_weak

outcome_matrix_LM_weak = matrix(rep(0, 1 * length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_LM_weak) <- beta_list_weak

outcome_matrix_LM_rob_weak = matrix(rep(0, 1 * length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_LM_rob_weak) <- beta_list_weak

outcome_matrix_LR_weak = matrix(rep(0, 1 * length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_LR_weak) <- beta_list_weak

outcome_matrix_LR_rob_weak = matrix(rep(0, 1 * length(beta_list_weak)), nrow = 1, ncol = length(beta_list_weak))
colnames(outcome_matrix_LR_rob_weak) <- beta_list_weak

a = 0.3
e_10 = c(a, rep(a,k-1))

for(N in 1:sim_amount){
  covmat = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  epsV = mvrnorm(n, mu = c(0,0), Sigma = covmat)
  eps = epsV[,1]
  V = epsV[,2]
  
  # distributional contamination
  epsCon = rmvt(n, sigma = covmat, df = 2)
  eps[1:30] = epsCon[,1][1:30]
  V[1:30] = epsCon[,2][1:30]
  
  Z = matrix(rnorm(n*k), nrow = n, ncol = k)
  
  # Z[1,1] <- -7
  # Z[1,2] <- -4
  # Z[1,3] <- -4
  
  z1 <- Z[,1]
  z2 <- Z[,2]
  z3 <- Z[,3]
  # z4 <- Z[,4]
  # z5 <- Z[,5]
  # z6 <- Z[,6]
  # z7 <- Z[,7]
  # z8 <- Z[,8]
  # z9 <- Z[,9]
  # z10 <- Z[,10]
  
  PZ = Z %*% solve(t(Z) %*% Z) %*% t(Z)
  MZ = diag(1, nrow = n, ncol = n) - PZ
  
  
  
  Pi = a*e_10
  X = Z %*% Pi + V
  
  for(j in 1:length(beta_list_weak)){
    beta = beta_list_weak[j]
    Y = X*beta + eps
    # Y[1] <- 7
    
    
    
    
    
    
    # Now the Anderson Rubin test
    
    AR_stat = ((t(Y) %*% PZ %*% Y)/k)/((t(Y) %*% MZ %*% Y)/(n - k))
    if(AR_stat > qchisq(0.95, k)/k ){
      outcome_matrix_AR_weak[1,j] <- outcome_matrix_AR_weak[1,j] + 1/sim_amount
    }
    
    full_model <-  lm(Y ~ 1 + z1 + z2 + z3)
    res_model  <-  lm(Y ~ 1)
    
    a1<-lmrob.control()
    a1$k.max <- 10000
    
    rob_model <- lmrob(formula = full_model, control = a1, method = "MM")
    estScale <- res_model$scale
    
    W_weights <- sqrt(1 - diag(PZ))
    
    W_stat = RobustScoreTest(full_model,
                            res_model, 
                            lmrob(Y ~ 1 + z1 + z2 + z3)$scale,
                            W_weights,
                            tukey.c = 5.51)$W.statistic
    if(W_stat > qchisq(0.95, k)){
      outcome_matrix_W_weak[1,j] <- outcome_matrix_W_weak[1,j] + 1/sim_amount
    }
    
    
    # Now the LM test
    beta0 <- 0
    eps0 <- Y - X*beta0
    sigma_hat_ee <- (1/(n-k))*(t(eps0) %*% MZ %*% eps0)
    sigma_hat_eV <- (1/(n-k))*(t(eps0) %*% MZ %*% X)
    rho_hat <- sigma_hat_eV/sigma_hat_ee

    pi_tilde <- solve(t(Z) %*% Z) %*% t(Z) %*% (X - eps0 %*% rho_hat)
    PZpi <- Z %*% pi_tilde
    projPZpi <- PZpi %*% solve(t(PZpi) %*% PZpi) %*% t(PZpi)

    LM_stat <- (1/sigma_hat_ee)*(t(eps0) %*% projPZpi %*% eps0)
    if(LM_stat > qchisq(0.95,1)){
      outcome_matrix_LM_weak[1,j] <- outcome_matrix_LM_weak[1,j] + 1/sim_amount
    }
    
    # Now the robust LM test
    
    first.pi = unname(coefficients(lmrob(X ~ 0 + Z)))
    V_hat = unname(residuals(lmrob(X ~ 0 + Z)))
    eps0_hat = unname(residuals(lmrob(eps0 ~ 0 + Z))) 
    cov_eps0.V = (1/(n-k))*(t(eps0_hat) %*% V_hat)
    cov_eps0eps0 = (1/(n-k))*t(eps0_hat) %*% eps0_hat
    rho_hat1 = cov_eps0.V/cov_eps0eps0
    
    X.transformed <- X - eps0 %*% rho_hat1
    pi_tilde_rob <- coefficients(lmrob(X.transformed ~ 0 + Z))
    
    Z.pi_rob <- Z %*% pi_tilde_rob
    model_res <- lm(eps0 ~ 1)
    model_full <- lm(eps0 ~ 1 + Z.pi_rob)
    
    gewichtjes <- sqrt(1 - diag(Z.pi_rob %*% solve(t(Z.pi_rob) %*% Z.pi_rob) %*% t(Z.pi_rob)))
    
    rob_lm <- RobustScoreTest(model_full,
                             model_res,
                             lmrob(eps0 ~ Z)$scale,
                             gewichtjes,
                             4.68)$W.statistic 
    
    if(rob_lm > qchisq(0.95,1)){
      outcome_matrix_LM_rob_weak[1,j] <- outcome_matrix_LM_rob_weak[1,j] + 1/sim_amount
    }
    
    
    sigma_VV = (1/(n-k))*(t(X) %*% MZ %*% X)
    sigma_VV.e = sigma_VV - sigma_hat_eV^2/sigma_hat_ee
    r_beta0 = (1/sigma_VV.e) * t(PZpi) %*% PZpi
    # Now the LR test
    LR_stat = 0.5*(k*AR_stat - r_beta0 + sqrt((k*AR_stat + r_beta0)^2 - 4*r_beta0 *(k*AR_stat - LM_stat)))
    if(LR_stat > sorted_grid_matrix[0.95*grid_simul, 1 + min(max(lr_grid), round(r_beta0))]){
      outcome_matrix_LR_weak[1,j] <- outcome_matrix_LR_weak[1,j] + 1/sim_amount
    }
    # LR_stat = CLR(ivmodel(Y, X, Z))$p.value
    # if(LR_stat < 0.05){
    #   outcome_matrix_LR_weak[1,j] <- outcome_matrix_LR_weak[1,j] + 1/sim_amount
    # }
    
    # Now LR robust
    rob_AR <- W_stat
    rob_LM <- rob_lm
    sigma_hat_VV_rob <- (1/(n-k)) * (t(V_hat) %*% V_hat)
    sigma_hat_VV.e_rob <- sigma_hat_VV_rob - cov_eps0.V^2/cov_eps0eps0
    r_weights <- sqrt(1 - diag(PZ))
    r_beta0_rob <- RobustWaldTest(X.transformed ~ 1 + z1 + z2 + z3,
                                  X.transformed ~ 1,
                                  lmrob(X.transformed ~ 1 + z1 + z2 + z3)$scale,
                                  #as.numeric(sigma_hat_VV.e_rob),
                                  r_weights,
                                  5.51)$W.statistic
    #r_beta0_rob <- (1/sigma_hat_VV.e_rob) * (t(Z.pi_rob) %*% Z.pi_rob)
    
    rob_LR <- 0.5*(rob_AR - r_beta0_rob + sqrt((rob_AR + r_beta0_rob)^2 -4*r_beta0_rob*(rob_AR - rob_LM)))
    if(rob_LR > sorted_grid_matrix[0.95*grid_simul, 1 + min(max(lr_grid), round(r_beta0_rob))]){
      outcome_matrix_LR_rob_weak[1,j] <- outcome_matrix_LR_rob_weak[1,j] + 1/sim_amount
    }
  }
  progress(N, sim_amount)
}

sizeline <- rep(0.05, length(beta_list_weak))


par(mfrow=c(1,2))
plot(outcome_matrix_LR_weak[1,],
     ylab = "Rejection Frequency",
     xlab = expression(paste(beta)),
     xaxt = "n",
     type = "l",
     main = expression(paste(pi, " = 0.3")),
     ylim = c(0, 1),
     lty=2)
# axis(1, at = 1:length(beta_list_weak), label=beta_list_weak)
axis(1, at = c(1, 6, 11, 16, 21, 26, 31), label=c("-3", "-2", "-1", "0", "1", "2", "3"))
# lines(outcome_matrix_LM_weak[1,], lty = 4)
# lines(outcome_matrix_LM_rob_weak[1,], lty=1)
# lines(outcome_matrix_LR_weak[1,],lty=2, col = "red")
# lines(outcome_matrix_W_weak[1,], lty = 4)
lines(outcome_matrix_LR_rob_weak[1,], lty = 1)
lines(sizeline, lty = 3)


# legend("bottomright", 
#        legend=c("t-test", "AR-test", "Robust-AR"),
#        lty=c(3,2,1),
#        bty="n")



plot(outcome_matrix_LR[1,], ylab = "Rejection Frequency", xlab = expression(paste(beta)), xaxt = "n", type = "l", main = expression(paste(pi, " = 0.75")),
     ylim = c(0,1),
     lty=2)
# axis(1, at = 1:length(beta_list), label=beta_list)
axis(1, at = c(1, 6, 11, 16, 21, 26, 31), label=c("-0.3", "-0.2", "-0.1", "0", "0.1", "0.2", "0.3"))
# lines(outcome_matrix_LM[1,], lty = 4)
# lines(outcome_matrix_LM_rob[1,], lty = 1)
# lines(outcome_matrix_LR[1,],lty=2, col = "red")
lines(outcome_matrix_LR_rob[1,], lty = 1)
# lines(outcome_matrix_W[1,], lty = 4)
lines(sizeline, lty = 3)




