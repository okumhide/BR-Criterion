BRC.fit {BRC}	

R Documentation
This function implements simultaneously bias reduction and model selection using the R package stats. The loss function consist of the (-1) times likelihood function and discontinuous penalty functions. For the minimization, the approximate penalties are employed iteratively. The function also implements the best subset selection by AIC and BIC.
Description
This function implements simultaneously bias reduction and model selection using the R package stats. The loss function consist of the (-1) times likelihood function and discontinuous penalty functions. For the minimization, the approximate penalties are employed iteratively. The function also implements the best subset selection by AIC and BIC.

Usage
BRC.fit(
  X,
  y,
  family = c("binomial", "poisson", "gaussian"),
  penalty = c("AIC", "BIC", "BRC"),
  K_max = 50,
  k_scale = 2
)

Arguments
X	
The n*p design matrix without an intercept. n is sample size. p is the number of predictor variables. Note that continuous predictors standardized are used.

y	
The n*1 response vector. binary data for family="binomial", count data for family="poisson", quantitative data for family="gaussian".

family	
The type of models to fit the data. family="binomial" and family="poisson" and family='gaussian' are provided corresponding to the above type of responses.

penalty	
The type of penalty. Selection by The Akaike information criterion for penalty = "AIC", The Bayes information criterion for penalty = "BIC", and the proposed bias reduction criterion for penalty = "BRC".

K_max	
The maximum number of iterations in the BR method. the default is 50.

k_scale	
The parameter increase rate that approximates a discontinuous penalty function with a continuous function. the default is 2.

Value
Returns an object with

theta	
The estimator includeing the intercept.

theta_cnt	
The vector of indexes selected by BR with regularization step.

iter_cnt	
number of iterations.

Author(s)
Hidenori Okumura

References
Maritines, J. M. (2002). Minimization of discontinuous cost functions by smoothing. Acta Appl. Math., 71, 245–260. Okumura, H. (2021). Bias reduction and model selection in misspecified models. Commun. Stat. Theory Methods. doi: 10.1080/03610926.2021.1959613.

Examples
Run examples

library(BRC)
n = 200; p = 13; rho = 0.5
corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
corrmat[,4] = 0.25
corrmat[4, ] = 0.25
corrmat[4,4] = 1
corrmat[,5] = 0
corrmat[5, ] = 0
corrmat[5,5] = 1
cholmat = chol(corrmat)
x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
x = x%*%cholmat
x = apply(x[,1:5], 2, scale)
set.seed(1)
b = c(1,-1,2,2,2)
eta = x[,1:5]%*%b
prob = 1-exp(-exp(eta))
y = rbinom(n, 1, prob)
Model <- "binomial"
result <- BRC.fit(x, y, family = Model, penalty ="BRC")
cat("BRC",result$theta,"\n")
cat("BRC",result$theta_ix,"\n")
cat("BRC",result$iter_cnt,"\n")
