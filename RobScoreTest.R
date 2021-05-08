# In this code we will implement the robust Wald test introduced by Markatou and Hettmansperger


# You can use this as an example
# n=100
# w1 = rnorm(n)
# w2 = rnorm(n)
# z1 = rnorm(n)
# z2 = rnorm(n)
# 
# X = 1 + 1*w1 + 2*w2 + 0*z1 + 0*z2 + rnorm(n)
# Y = 1*X + 2*w1 + 1*w2 + rnorm(n)
# 
# Y = Y
# X = X
# Z = cbind(z1, z2)
# W = cbind(w1, w2)
# beta0 =0
# weight.vector = rep(1,n)
# tukey.c = 4.68

RobScoreTest(Y,
             X,
             cbind(z1, z2),
             cbind(w1, w2),
             0,
             rep(1,n))


RobScoreTest <- function(Y,
                            X,
                            Z,
                            W,
                            beta0,
                            weight.vector,
                            tukey.c = 4.68){
  library(robustbase)
  
  # This block obtains the syntax of the full model and the restricted model
  dep_var <<- Y - beta0 * X
  full.model <- dep_var ~ 1 + W + Z
  restricted.model <- dep_var ~ 1 + W
  scale.est <- lmrob(full.model)$scale
  
  # full.model <- as.formula(full.model)
  # restricted.model <- as.formula(restricted.model)
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
  n <- length(get(all.vars(full.model)[1]))
  
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
  control.variables <-lmrob.control()
  control.variables$k.max <- 10000
  #control.variables$setting <- "KS2011"
  
  gewichten <<- weight.vector
  # we estimate the restricted model using robust MM estimators 
  if(q > 1){
    mm.estimator.object <- lmrob(formula = restricted.model, 
                                 method = "MM",  
                                 weights = gewichten,
                                 control = control.variables)
    restricted.model.estimate <- unname(coefficients(mm.estimator.object))
    restricted.model.estimate  <- c(restricted.model.estimate, rep(0,p-q))
  }
  if(q == 1){
    restricted.model.estimate <- c(rep(0,p-q))
  }
  
  # We now calculate the residuals in the restricted model
  if(length(restricted.model.estimate) > 1){
    residuals.restricted.model <- (y - X.full %*% restricted.model.estimate)/scale.est
  }
  if(length(restricted.model.estimate) == 1){
    residuals.restricted.model <- (y - restricted.model.estimate * X.full)/scale.est
  }
  
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
      return((1 - x^2/c^2)^2 - 2*(x^2/c^2)*(1 - x^2/c^2))
    }
    else{
      return(0)
    }
  }
  tukey.residuals.restricted = apply(as.array(residuals.restricted.model), c(1), FUN = tukey_biweight_function, c = tukey.c)
  if(p - q == 1 & q > 1){
    biweight.sum = sum(weight.vector * tukey.residuals.restricted * X.full[, q:(p-1)])
  }
  if(p - q == 1 & q == 1){
    biweight.sum = sum(weight.vector * tukey.residuals.restricted * X.full)
  }
  if(p - q > 1){
    biweight.sum =  apply(as.matrix(as.numeric(weight.vector) * tukey.residuals.restricted * X.full[, q:(p-1)]), 2, sum)
  }
  
  # now we have to construct the U matrix
  # to calculate the U matrix, we first need to calculate the M and Q matrix
  squared.tukey.residuals.restricted = weight.vector^2 * tukey.residuals.restricted^2
  derivative.tukey.residuals.restricted = weight.vector * apply(as.array(residuals.restricted.model), c(1), FUN = tukey_biweight_derivative, c = tukey.c)
  
  # Q.matrix.init <- matrix(rep(0, (p-1)^2), p-1, p-1)
  # M.matrix.init <- matrix(rep(0, (p-1)^2), p-1, p-1)
  # for(number in 1:n){
  #   M.matrix.init <- M.matrix.init + derivative.tukey.residuals.restricted[number] * X.full[number,] %*% t(X.full[number,])
  #   Q.matrix.init <- Q.matrix.init + squared.tukey.residuals.restricted[number] * X.full[number,] %*% t(X.full[number,])
  # }
  # Q.matrix <- (1/n) * Q.matrix.init
  # M.matrix <- (1/n) * (1/scale.est) * M.matrix.init
  
  Data.Prep.M.Matrix <- cbind(derivative.tukey.residuals.restricted, X.full)
  Data.Prep.Q.Matrix <- cbind(squared.tukey.residuals.restricted, X.full)
  M.Q.function <- function(x){
    return(x[1] * x[2:length(x)] %*% t(x[2:length(x)]))
  }
  Q.matrix1 <- (1/n) * apply(Data.Prep.Q.Matrix, 1, M.Q.function)
  M.matrix1 <- (1/n) * (1/scale.est) * apply(Data.Prep.M.Matrix, 1, M.Q.function)
  
  Q.matrix <- matrix(apply(Q.matrix1, 1, sum), nrow = p-1, ncol = p-1)
  M.matrix <- matrix(apply(M.matrix1, 1, sum), nrow = p-1, ncol = p-1)
  
  if(q > 1){
    invM11 <- solve(M.matrix[1:(q-1),1:(q-1)])
    Q22 <- Q.matrix[q:(p-1), q:(p-1)]
    Q11 <- Q.matrix[1:(q-1), 1:(q-1)]
    M12 <- M.matrix[1:(q-1), q:(p-1)]
    M21 <- M.matrix[q:(p-1), 1:(q-1)]
    Q12 <- Q.matrix[1:(q-1), q:(p-1)]
    Q21 <- Q.matrix[q:(p-1), 1:(q-1)]
    
    U.matrix <- Q22 - (M21 %*% invM11 %*% Q12) - (Q21 %*% invM11 %*% M12) + (M21 %*% invM11 %*% Q11 %*% invM11 %*% M12)
    W.stat <- (1/n) * t(biweight.sum) %*% solve(U.matrix) %*% biweight.sum
  }
  
  if(q == 1){
    Q22 <- Q.matrix[q:(p-1), q:(p-1)]
    U.matrix <- Q22
    W.stat <- (1/n) * t(biweight.sum) %*% solve(U.matrix) %*% biweight.sum
  }
  
  output <- list(W.statistic = W.stat,
                 p.value = pchisq(W.stat, p - q, lower.tail = FALSE))
  return(output)
}