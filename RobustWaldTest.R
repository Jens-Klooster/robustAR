library(robustbase)

RobustWaldTest <- function(full.model,
                            restricted.model,
                            scale.est,
                            weight.vector,
                            tukey.c = 4.68){
  # This block obtains the syntax of the full model and the restricted model
  full.model <- as.formula(full.model)
  restricted.model <- as.formula(restricted.model)
  weight.vector <- as.matrix(weight.vector)
  variables.full.model <-  all.vars(full.model)
  variables.restricted.model <- all.vars(restricted.model)
  
  # p is the amount of variables in the full model (it includes the dependent variable in the count)
  p <- length(variables.full.model)
  if(attr(terms(full.model), "intercept") == 1){    # here we check if the full model has an intercept
    p <- p + 1
  }
  
  # q is the amount of variables in the restricted model (it includes the dependent variable in the count)
  q <- length(variables.restricted.model)
  if(attr(terms(restricted.model), "intercept") == 1){
    q <- q + 1
  }
  
  # n is the amount of observations 
  n <-  length(get(all.vars(full.model)[1]))
  
  # here we collect all the data
  data.full.model <- apply(t(all.vars(full.model)), 2, get)
  colnames(data.full.model) <- all.vars(full.model)
  
  # here we construct the independent variable matrix of the full model
  X.full <- data.full.model
  X.full <- X.full[,-1]
  if(attr(terms(full.model), "intercept") == 1){
    X.full <- cbind(1, X.full)
  }
  
  # here we construct the independent variable matrix of the restricted model
  data.restricted.model <- apply(t(all.vars(restricted.model)), 2, get)
  colnames(data.restricted.model) <- all.vars(data.restricted.model)
  
  X.restricted <- data.restricted.model
  X.restricted <- data.restricted.model[,-1]
  if(attr(terms(restricted.model), "intercept") == 1){
    X.restricted <- cbind(1, X.restricted)
  }
  if(q == 1 & attr(terms(restricted.model), "intercept") == 0){
    X.restricted <-  rep(0, n)
  }
  
  # here we get the data of the dependent variable
  y <- data.full.model[,1]
  
  # we increase the max iteration from 200 to 10000 for our robust regression estimates
  #weights <- weight.vector
  # control.variables <-lmrob.control()
  # control.variables$k.max <- 10000
  #control.variables$setting <- "KS2011"
  
  # we estimate the restricted model using robust MM estimators 
  
  mm.estimator.object <- lmrob(formula = full.model, 
                               method = "MM")
  full.model.estimate <- unname(coefficients(mm.estimator.object))
  # full.model.estimate[q:(p-1)] <- rep(0, p-q)
  # full.model.estimate2 <- unname(coefficients(mm.estimator.object))
  
  
  # We now calculate the residuals in the full model
  
  residuals.full.model <- residuals(mm.estimator.object)/scale.est
  
  # if(length(full.model.estimate) > 1){
  #   residuals.full.model <- (y - X.full %*% full.model.estimate)/scale.est
  # }
  # if(length(full.model.estimate) == 1){
  #   residuals.full.model <- (y - full.model.estimate * X.full)/scale.est
  # }
  
  # Here we calculate the first part (the "sum" part) of our W test statistic 
  # which has the form W = 1/n * biweight.sum * invU * biweight.sum
  # We start by defining the biweight function and its derivative
  tukey_biweight_function <- function(x,c){
    if(abs(x) < c){
      return(x*(1 - x^2/c^2)^2)
    }
    else{
      return(0)
    }
  }
  tukey_biweight_derivative <- function(x,c){
    if(abs(x) < c){
      return(1 - 6*(x^2/c^2) + 5*(x^4/c^4))
    }
    else{
      return(0)
    }
  }
  tukey.residuals.full = apply(as.array(residuals.full.model), c(1), FUN = tukey_biweight_function, c = tukey.c)
  
  # now we have to construct the U matrix
  
  # to calculate the U matrix, we first need to calculate the M and Q matrix
  squared.tukey.residuals.full = weight.vector^2 *  tukey.residuals.full^2
  derivative.tukey.residuals.full = weight.vector * apply(as.array(residuals.full.model), c(1), FUN = tukey_biweight_derivative, c = tukey.c)
  
  Q.matrix.init <- matrix(rep(0, (p-1)^2), p-1, p-1)
  M.matrix.init <- matrix(rep(0, (p-1)^2), p-1, p-1)
  for(number in 1:n){
    M.matrix.init <- M.matrix.init + derivative.tukey.residuals.full[number] * X.full[number,] %*% t(X.full[number,])
    Q.matrix.init <- Q.matrix.init + squared.tukey.residuals.full[number] * X.full[number,] %*% t(X.full[number,]) 
  }
  Q.matrix <- (1/n) * Q.matrix.init
  M.matrix <- (1/n) * (1/scale.est) * M.matrix.init
  
  
  # with the Q and M matrix, we can now calculate the U matrix
  V.matrix <- solve(M.matrix) %*% Q.matrix %*% t(solve(M.matrix))
  
  W.stat =  n * t(full.model.estimate[q:(p-1)]) %*% solve(V.matrix[q:(p-1), q:(p-1)]) %*% full.model.estimate[q:(p-1)]

  pval = 1 - pchisq(W.stat, (p-q))
  return(list(W.statistic = W.stat,
              p.value = pval))
}

# library(robustbase)
# library(survey)
# n=1000
# count = 0
# sim_amount = 200
# for(i in 1:sim_amount){
#   z1 = rnorm(n)
#   z2 = rnorm(n)
#   z3 = rnorm(n)
#   z4 = rnorm(n)
#   eps = rnorm(n, 0, 3)
#   # z1[1] = 100
#   # z2[1] = 100
#   
#   
#   #gewichten = sqrt(1 - diag(z1 %*% solve(t(z1) %*% z1) %*% t(z1)))
#   # mcd.obj <- covMcd(cbind(z1,z2))
#   # mean <- mcd.obj$center
#   # cov <- mcd.obj$cov
#   # gewichten = pmin(1, 3.84/sqrt(mahalanobis(cbind(z1,z2), mean, cov)))
#   # gewichten = c(1, rep(1, n-1))
#   y = 1 + 2*z1 + 0*z2 + 0*z3 + eps
#   z2[1] = -10000
#   y[1] = 10000
#   
#   full_model = y ~ 1 + z1 + z2 + z3
#   restricted_model = y ~ 1 + z1
#   scale_est = lmrob(full_model)$scale
#   Z <- cbind(z1, z2, z3)
#   weights = sqrt(1 - diag(Z %*% solve(t(Z) %*% Z) %*% t(Z)))
#   tukey_c = 4.68
#   
#   # full.model = full_model
#   # restricted.model = restricted_model
#   # scale.est = scale_est
#   # tukey.c = tukey_c
#   # weight.vector <- sqrt(1 - diag(Z %*% solve(t(Z) %*% Z) %*% t(Z)))
#   
#   rob_wald <- RobustWaldTest2(full_model, restricted_model, scale_est, weights)
#   rob_wald
#   
#   if(rob_wald < 0.05){
#     count = count + 1
#   }
# }
# count/sim_amount


