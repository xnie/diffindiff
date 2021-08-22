rm(list = ls())
library(Matrix)
library(foreach)
library(magrittr)
library(grf)
library(glmnet)
library(rlearner)
library(diffindiff)
library(EQL)
library(rpart)
source("simulations/utils.R")

start_time <- Sys.time()

args = commandArgs(trailingOnly = TRUE)
setup = as.character(args[1])
n = as.numeric(args[2])
p = as.numeric(args[3])
NREP = as.numeric(args[4])

eta = 0.1

if (setup == 'A'){
  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    tau = X[,4] + 0.5*X[,5]
    TAU = 0
    f = 1/(1+exp(X[,3])) + 4*X[,5]^2
    h = 1/(1+exp(X[,4])) + 3*X[,6]^2
    probabilities = matrix(0,4,n)
    P_T1 = 0.4
    P_S1 =  0.6
    P_T1_S1 = P_T1*P_S1
    P_T1_S0 = P_T1*(1-P_S1)
    P_T0_S1 = (1-P_T1)*P_S1
    P_T0_S0 = (1-P_T1)*(1-P_S1)
    Ti = rbinom(n,1,P_T1)
    Si = rbinom(n,1,P_S1)
    b = pmax(X[,1] + X[,2], 0) + 4* X[,6]^2
    Y = b + Ti*f + Si*h + Ti*Si*tau + rnorm(n)
    list(X=X, TAU=TAU, tau=tau, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  }
} else if (setup == 'B'){
  get.params = function(){
    X = matrix(rnorm(n * p), n, p)
    tau = 0.5*(X[,1] + X[,2] + X[,3] )
    TAU = 2
    f = 5 * (sin(pi*X[,1]*X[,2]) + 2*(X[,3]-0.5)^2 )
    h = 5 * (sin(pi*X[,1]*X[,2]) +2*X[,5]^2)
    probabilities = matrix(0,4,n)
    P_T1 = 0.5
    P_S1 = 0.5
    P_T1_S1 = P_T1*P_S1
    P_T1_S0 = P_T1*(1-P_S1)
    P_T0_S1 = (1-P_T1)*P_S1
    P_T0_S0 = (1-P_T1)*(1-P_S1)
    Ti = rbinom(n,1,P_T1)
    Si = rbinom(n,1,P_S1)
    b = 0
    Y = b + Ti*f + Si*h + Ti*Si*tau + rnorm(n)
    return(list(X=X, TAU=TAU, tau=tau, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0))
  }
} else if (setup == 'C'){
  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    TAU = 1
    f = 0
    h = 0
    probabilities = matrix(0,4,n)
    locations = matrix(0,4,n)
    for (i in 1:n) {
      T1_S1 = 0.5 + 0.5*(1-6*eta)*sin(X[i,1]*1.5)
      T1_S0 = (1 - T1_S1)/3
      T0_S1 = (1-T1_S1)/3
      T0_S0 = 1 - T1_S1 - T1_S0 - T0_S1
      probabilities[,i] = c(T1_S1, T1_S0, T0_S1, T0_S0)
      locations[,i] = rmultinom(1, 1, c(T1_S1, T1_S0, T0_S1, T0_S0))
    }
    P_T1_S1 = probabilities[1,]
    P_T1_S0 = probabilities[2,]
    P_T0_S1 = probabilities[3,]
    P_T0_S0 = probabilities[4,]
    Ti = locations[1,] + locations[2,]
    Si = locations[1,] + locations[3,]
    b = 2*sin(X[,1]*1.5)
    Y = b + Ti*f + Si*h + Ti*Si*TAU + rnorm(n)
    list(X=X, TAU=TAU, tau=TAU, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  } else if (setup == 'D') {
  get.params = function(){
     X = matrix(rnorm(n * p), n, p)
    tau = 3*X[,1] + 2*X[,4]
    TAU=2
    f = 2 * X[,5]
    h = 0
    probabilities = matrix(0,4,n)
    P_T1 = pmin(pmax(eta, (1 / (1+exp(0.5*X[,2])))), 1-eta)
    P_S1 = pmin(pmax(eta, (1-1 /(1+exp(0.5*X[,3])))), 1-eta)
    P_T1_S1 = P_T1*P_S1
    P_T1_S0 = P_T1*(1-P_S1)
    P_T0_S1 = (1-P_T1)*P_S1
    P_T0_S0 = (1-P_T1)*(1-P_S1)
    Ti = rbinom(n,1,P_T1)
    Si = rbinom(n,1,P_S1)
    b = pmax(X[,1] + X[,2] + X[,4] + X[,6], 0)
    Y = b + Ti*f + Si*h + Ti*Si*tau + rnorm(n)
    list(X=X, TAU=TAU, tau=tau, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  }
} else {
  stop("bad setup")
}



results = lapply(1:NREP, function(i){
  params.train = get.params()
  params.test = get.params()

  # alg == "OLS"
  input_data = as.data.frame(cbind(params.train$Ti*params.train$Si*cbind(params.train$X,1),params.train$X, params.train$Ti*cbind(params.train$X,1), params.train$Si*cbind(params.train$X,1), params.train$Y))
  colnames(input_data)[dim(input_data)[2]] =  "Y"
  tau_fit = lm(Y~., input_data)
  tau_hat_coef = tau_fit$coefficients[2:(p+2)]
  tau_pred_ols = as.matrix(cbind(params.test$X, 1)) %*% as.matrix(tau_hat_coef, p+1, 1)

  # alg == "T"
  S1T1_fit = grf::regression_forest(X = params.train$X[params.train$Si==1 & params.train$Ti==1,], Y=params.train$Y[params.train$Si==1 & params.train$Ti==1], tune.parameters = "all")
  S1T0_fit = grf::regression_forest(X = params.train$X[params.train$Si==1 & params.train$Ti==0,], Y=params.train$Y[params.train$Si==1 & params.train$Ti==0], tune.parameters = "all")
  S0T1_fit = grf::regression_forest(X = params.train$X[params.train$Si==0 & params.train$Ti==1,], Y=params.train$Y[params.train$Si==0 & params.train$Ti==1], tune.parameters = "all")
  S0T0_fit = grf::regression_forest(X = params.train$X[params.train$Si==0 & params.train$Ti==0,], Y=params.train$Y[params.train$Si==0 & params.train$Ti==0], tune.parameters = "all")
  tau_pred_t = predict(S1T1_fit, params.test$X)$predictions - predict(S1T0_fit, params.test$X)$predictions - predict(S0T1_fit, params.test$X)$predictions + predict(S0T0_fit, params.test$X)$predictions

  # alg == "RS"
  S1_fit = grf::causal_forest(X = params.train$X[params.train$Si==1,], W=params.train$Ti[params.train$Si==1], Y=params.train$Y[params.train$Si==1], tune.parameters = "all")
  S0_fit = grf::causal_forest(X = params.train$X[params.train$Si==0,], W=params.train$Ti[params.train$Si==0], Y=params.train$Y[params.train$Si==0], tune.parameters = "all")
  tau_pred_rs = predict(S1_fit, params.test$X)$predictions - predict(S0_fit, params.test$X)$predictions

  # alg == "RT"
  T1_fit = grf::causal_forest(X = params.train$X[params.train$Ti==1,], W=params.train$Si[params.train$Ti==1], Y=params.train$Y[params.train$Ti==1], tune.parameters = "all")
  T0_fit = grf::causal_forest(X = params.train$X[params.train$Ti==0,], W=params.train$Si[params.train$Ti==0], Y=params.train$Y[params.train$Ti==0], tune.parameters = "all")
  tau_pred_rt = predict(T1_fit, params.test$X)$predictions - predict(T0_fit, params.test$X)$predictions

  # alg == "DiD"
  nu_fit = grf::causal_forest(X = params.train$X, W=params.train$Ti, Y=params.train$Y, tune.parameters="all")
  m_hat = nu_fit$Y.hat
  t_hat = nu_fit$W.hat
  nu_hat = predict(nu_fit)$predictions

  sigma_fit = grf::causal_forest(X = params.train$X, W=params.train$Si, Y=params.train$Y, tune.parameters="all")
  sigma_hat = predict(sigma_fit)$predictions
  s_hat = sigma_fit$W.hat

  delta_fit = grf::regression_forest(X = params.train$X, Y = params.train$Ti * params.train$Si, tune.parameters = "all")
  delta_hat = predict(delta_fit)$predictions - s_hat * t_hat

  scaling = 1 - (delta_hat^2 / (s_hat *(1 - s_hat) * t_hat * (1 - t_hat)))
  A_hat = (params.train$Ti - t_hat - (delta_hat * (params.train$Si - s_hat)) / (s_hat * (params.train$Si - s_hat))) / scaling
  B_hat = (params.train$Si - s_hat - (delta_hat * (params.train$Ti - t_hat)) / (t_hat * (params.train$Ti - t_hat))) / scaling
  C_hat = params.train$Si * params.train$Ti - (s_hat + delta_hat / t_hat) * A_hat - (t_hat + delta_hat / s_hat) * B_hat - (s_hat * t_hat + delta_hat)

  S_hat = params.train$Y - m_hat - A_hat*nu_hat - B_hat*sigma_hat
  tau_fit = grf::causal_forest(X=params.train$X, W = C_hat, Y=S_hat, tune.parameters = "all")
  tau_pred_did = predict(tau_fit, params.test$X)$predictions

  c(mean((params.test$tau - tau_pred_did)^2),
    mean((params.test$tau - tau_pred_rt)^2),
    mean((params.test$tau - tau_pred_rs)^2),
    mean((params.test$tau - tau_pred_t)^2),
    mean((params.test$tau - tau_pred_ols)^2))
})

end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

fnm = paste("results/results_non_const/output", setup, n, p, NREP, "full.csv", sep="-")
results = t(matrix(unlist(results, use.names=FALSE), 5, NREP))
write.csv(results, file=fnm)

}
