library(R.matlab)
library(robustbase) 
library(readstata13) # need this to load in stata dta file
library(ivmodel) # need this for AR test
library(robustX) # need this for Bacon algorithm
library(riv) # robust IV estimator
source("RobustScoreTest.R")
source("RobustWaldTest.R")

card.data <- readMat("EmpiricalExamples/assignmentweakinstruments.mat")

for(i in 1:length(names(card.data))){
  name = names(card.data)[i]
  assign(name, card.data[[i]])
}

Y <- wage
X <- ed
Z <- cbind(nearc2, nearc4, nearc4a, nearc4b)
W <- cbind(1, exper, exper2/100, south, smsa, race)

MW <- (diag(1, nrow = length(Y)) - W %*% solve(t(W) %*% W) %*% t(W))

Z <- MW %*% Z

MZ <- (diag(1, nrow = length(Y)) - Z %*% solve(t(Z) %*% Z) %*% t(Z))

ZW = cbind(Z, W)
MZW = (diag(1, nrow = length(Y)) - ZW %*% solve(t(ZW) %*% ZW) %*% t(ZW))

nearc2 <- Z[,1]
nearc4 <- Z[,2]
nearc4a <- Z[,3]
nearc4b <- Z[,4]

Y <- MW %*% Y
X <- MW %*% X

n = length(Y)
k = length(Z[1,])
p = length(W[1,])

H_matrix <- cbind(W, Z)
gewichtjes <- sqrt(1 - diag(H_matrix %*% solve(t(H_matrix) %*% H_matrix) %*% t(H_matrix)))

beta_grid = seq(0.05, 0.30, 0.01)
conf_set_robust = rep(-Inf, length(beta_grid))
conf_set_lm = rep(-Inf, length(beta_grid))
conf_set_klm = rep(-Inf, length(beta_grid))
conf_set_rob_lr = rep(-Inf, length(beta_grid))

p_values = rep(1, length(beta_grid))
p_values_lm =  rep(1, length(beta_grid))
p_values_klm = rep(1, length(beta_grid))

control.variables <-lmrob.control()
control.variables$k.max <- 1000
control.variables$refine.tol <- 1e-10


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

for(i in 1:length(beta_grid)){
  dep_var = Y - beta_grid[i]*X
  
  full_model <-  dep_var ~ 1 + exper + exper2 + south + smsa + race + nearc2 + nearc4 + nearc4a + nearc4b
  res_model  <-  dep_var ~ 1 + exper + exper2 + south + smsa + race
  
  tau_stat = RobustScoreTest(full_model, 
                             res_model, 
                             lmrob(dep_var ~ 1 + exper + exper2 + south + smsa + race + nearc2 + nearc4 + nearc4a + nearc4b,
                                   weights = gewichtjes,
                                   control = control.variables)$scale, 
                             gewichtjes, 
                             tukey.c = 4.68)
  
  p_values[i] <- tau_stat$p.value
  if(tau_stat$W.statistic < qchisq(0.95,4)){
    conf_set_robust[i] <- beta_grid[i]
  }
  

  # For comparison I also check the LM test
  
  eps0 <- dep_var
  sigma_hat_ee <- (1/(n - k - p))*(t(eps0) %*% MZ %*% eps0)
  sigma_hat_eV <- (1/(n - k - p))*(t(eps0) %*% MZ %*% X)
  rho_hat <- sigma_hat_eV/sigma_hat_ee
  
  pi_tilde <- solve(t(Z) %*% Z) %*% t(Z) %*% (X - eps0 %*% rho_hat)
  PZpi <- Z %*% pi_tilde
  projPZpi <- PZpi %*% solve(t(PZpi) %*% PZpi) %*% t(PZpi)
  
  LM_stat <- (1/sigma_hat_ee)*(t(eps0) %*% projPZpi %*% eps0)
  AR_stat <- (1/sigma_hat_ee)*(t(eps0) %*% (Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% eps0)
  JLM_stat <- AR_stat - LM_stat
  
  if(LM_stat < qchisq(0.95,1)){
    conf_set_klm[i] <- beta_grid[i]
  }
  
  
  
  
  
  # Now robust LM test
  
  # dep_var <- Y - beta_grid[i]*X
  eps0 <- dep_var
  
  # first.pi = unname(coefficients(lmrob(X ~ 1 + exper + exper2 + south + smsa + race + nearc2 + nearc4 + nearc4a + nearc4b)))
  V_hat = unname(residuals(lmrob(X ~ 0 + nearc2 + nearc4 + nearc4a + nearc4b,
                                 control = control.variables)))
  eps0_hat = unname(residuals(lmrob(dep_var ~ 0 + nearc2 + nearc4 + nearc4a + nearc4b,
                                    control = control.variables))) 
  cov_eps0.V = as.numeric(cov.mcd(cbind(V_hat, eps0_hat))$cov[1,2])
  cov_eps0eps0 = as.numeric(cov.mcd(eps0_hat)$cov)
  rho_hat1 = cov_eps0.V/cov_eps0eps0
  
  X.transformed <- X - eps0 %*% rho_hat1
  pi_tilde_rob <- coefficients(lmrob(X.transformed ~ 1 + exper + exper2 + south + smsa + race + nearc2 + nearc4 + nearc4a + nearc4b,
                                     control = control.variables))
  
  Z.pi_rob <- Z %*% pi_tilde_rob[7:10]
  model_res <- lm(eps0 ~ 1 + exper + exper2 + south + smsa + race)
  model_full <- lm(eps0 ~ 1 + exper + exper2 + south + smsa + race + Z.pi_rob)
  
  ZW_mat <- cbind(W, Z.pi_rob)
  gewicht <- sqrt(1 - diag(ZW_mat%*% solve(t(ZW_mat) %*% ZW_mat) %*% t(ZW_mat)))
  # gewicht <- rep(1, length(Y))
  
  rob_lm <- RobustScoreTest(model_full,
                            model_res,
                            lmrob(eps0 ~ 1 + exper + exper2 + south + smsa + race + nearc2 + nearc4 + nearc4a + nearc4b,
                                  control = control.variables)$scale,
                            gewicht,
                            4.68) 
  
  p_values_lm[i] <- rob_lm$p.value
  
  if(rob_lm$p.value > 0.05){
    conf_set_lm[i] <- beta_grid[i]
  }
  
  
  # lets try robust LR test
  
  rob_AR <- tau_stat$W.statistic
  rob_LM <- rob_lm$W.statistic

  sigma_hat_VV_rob <- as.numeric(cov.mcd(V_hat)$cov)
  sigma_hat_VV.e_rob <- sigma_hat_VV_rob - cov_eps0.V^2/cov_eps0eps0
  full_mod_r <- X.transformed ~ 1 + nearc2 + nearc4 + nearc4a + nearc4b
  res_mod_r <- X.transformed ~ 1
  r_beta0_rob <- RobustWaldTest(full_mod_r,
                                res_mod_r,
                                lmrob(full_mod_r)$scale,
                                #as.numeric(sigma_hat_VV.e_rob),
                                rep(1, length(Y)),
                                4.68)$W.statistic
  #r_beta0_rob <- (1/sigma_hat_VV.e_rob) * (t(Z.pi_rob) %*% Z.pi_rob)

  rob_LR <- 0.5*(rob_AR - r_beta0_rob + sqrt((rob_AR + r_beta0_rob)^2 -4*r_beta0_rob*(rob_AR - rob_LM)))
  if(rob_LR < sorted_grid_matrix[0.95*grid_simul, 1 + min(max(lr_grid), round(r_beta0_rob))]){
    conf_set_rob_lr[i] <- beta_grid[i]
    # outcome_matrix_LR_rob[1,j] <- outcome_matrix_LR_rob[1,j] + 1/sim_amount
  }
}

conf_set_robust
conf_set_lm
conf_set_klm
conf_set_rob_lr
ivmodel(Y, X, Z, W)
AR.test(ivmodel(Y, X, Z, W), alpha = 0.05)$ci
# LIML(ivmodel(Y, X, Z, W, beta0 = 0.17))

plot(p_values)
plot(p_values_lm)
lines(rep(0.05, length(beta_grid)))

beta_grid[which(p_values == max(p_values))]