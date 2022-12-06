#  The bias reduced criterion (BRC) for Bias reduction and model selection in misspecified generalized linear models
#'
#' This function implements simultaneously bias reduction and model selection using the R package \pkg{stats}.
#' The loss function consist of the (-1) times likelihood function and discontinuous penalty functions.
#' For the minimization, the approximate penalties are employed iteratively. The function also implements the best subset selection by AIC and BIC.
#' @export
#' @importFrom stats optim
#' @param X The n*p design matrix without an intercept. n is sample size. p is the number of
#' predictor variables. Note that continuous predictors standardized are used.
#' @param y The n*1 response vector. binary data for \code{family="binomial"}, count data for \code{family="poisson"},
#' quantitative data for \code{family="gaussian"}.
#' @param family The type of models to fit the data. \code{family="binomial"} and \code{family="poisson"}
#' and \code{family='gaussian'} are provided corresponding to the above type of responses.
#' @param penalty The type of penalty. Selection by The Akaike information criterion for \code{penalty = "AIC"},
#' The Bayes information criterion for \code{penalty = "BIC"}, and the proposed bias reduction criterion for \code{penalty = "BRC"}.
#' @param K_max The maximum number of iterations in the BR method. the default is 50.
#' @param k_scale The parameter increase rate that approximates a discontinuous penalty function with a continuous function.
#' the default is 2.
#' @return Returns an object with \item{theta}{The estimator includeing the intercept.} \item{theta_cnt}{ The vector of indexes selected by
#' BR with regularization step.} \item{iter_cnt}{number of iterations.}
#' @author Hidenori Okumura
#' @references
#' Maritines, J. M. (2002). Minimization of discontinuous cost functions by smoothing. \emph{Acta Appl. Math.}, \bold{71}, 245--260.
#' Okumura, H. (2021). Bias reduction and model selection in misspecified models. \emph{Commun. Stat. Theory Methods.} doi: 10.1080/03610926.2021.1959613.
#' @examples
#' library(BRC)
#' n = 200; p = 13; rho = 0.5
#' corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
#' corrmat[,4] = 0.25
#' corrmat[4, ] = 0.25
#' corrmat[4,4] = 1
#' corrmat[,5] = 0
#' corrmat[5, ] = 0
#' corrmat[5,5] = 1
#' cholmat = chol(corrmat)
#' x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
#' x = x%*%cholmat
#' x = apply(x[,1:5], 2, scale)
#' set.seed(1)
#' b = c(1,-1,2,2,2)
#' eta = x[,1:5]%*%b
#' prob = 1-exp(-exp(eta))
#' y = rbinom(n, 1, prob)
#' Model <- "binomial"
#' result <- BRC.fit(x, y, family = Model, penalty ="BRC")
#' cat("BRC",result$theta,"\n")
#' cat("BRC",result$theta_ix,"\n")
#' cat("BRC",result$iter_cnt,"\n")



BRC.fit <- function(X, y, family = c("binomial", "poisson", "gaussian"), penalty = c("AIC", "BIC", "BRC"), K_max = 50, k_scale = 2){
  if (is.null(X) || is.null(y))
    stop("The data is missing!")
  if ( nrow(X) != length(y))
    stop("The data length is wrong!")
  family = match.arg(family)
  penalty = match.arg(penalty)
  if (mode(K_max) != "numeric"){stop("K_max must be numeric!")}
  if (mode(k_scale) != "numeric"){stop("k_scale must be numeric!")}


  n <- length(y)
  p <- length(X[1, ])
  D <- p + 1
  accel_rate <- 1
  eps <- 1e-4
  EPS <- 1e-6
  half_log_n <- log(n) * 0.5
  X <- cbind(rep(1, n), X)
  X_ <- X
  theta_ <- numeric(length(X[1, ]))

  if (penalty == "AIC") {
    M_k <- function(theta_) {
      -likelihood(X_, y, theta_, family) + penalty1(theta_, k)
    }
    M <- function(theta_){
      -likelihood(X_, y, theta_, family)
    }
  } else if(penalty == "BIC") {
    M_k <- function(theta_){
      -likelihood(X_, y, theta_, family) + half_log_n * penalty1(theta_, k)
    }
    M <- function(theta_){
      -likelihood(X_, y, theta_, family)
    }
  } else {
    M_k <- function(theta_){
      -likelihood(X_, y, theta_, family) + penalty0_k(X_, y, theta_, k, family) + half_log_n * penalty1(theta_, k)
    }
    M <- function(theta_){
      -likelihood(X_, y, theta_, family) + penalty0(X_, y, theta_, family)
    }
  }

  theta_ <- numeric(D)

  k <- 1
  theta_ix <- c(1:D)
  X_ <- X
  for(d in 1:K_max){
    #    old_theta_ix <- theta_ix
    #    old_nb <- length(old_theta_ix)
    old_theta_ <- theta_
    resoptim <- optim(theta_, M_k, method = "BFGS")
    theta_ <- resoptim$par
    #    u <- theta_ix
    nb <- length(theta_ix)
    u <- numeric(nb)
    u[abs(theta_) >= eps] <- 1
    u[1] <- 1
    theta_ix <- theta_ix[theta_ix*u>0]
    ix <- c(1:nb)
    ix <- ix[ix * u > 0]
    #    if(nb == old_nb){
    k <- k * accel_rate
    if(norm(old_theta_- theta_, "2") < EPS){
      theta_ <- theta_[ix]
      resoptim <- optim(theta_, M ,method = "BFGS")
      theta_ <- resoptim$par
      break;
    }
    #   }
    #    else{k <- k/accel_rate}

    theta_ <- theta_[ix]
    X_ <-  X_[,ix]

    k <- k*k_scale
  }
  theta <- numeric(D)
  theta[theta_ix] <- theta_

  return(list(theta = theta, theta_ix = theta_ix, iter_cnt = d))

}


H <- function(t, k){
  u <- k *t^2
  u / (1 + u)
}

H_ <- function(t,k){
  u <- exp(-k*t*t)
  2/(1+u)-1
}



likelihood <- function(X_, y, theta_, family = c("binomial", "poisson", "gaussian")){
  a <- X_ %*% theta_
  if (family == "binomial") {
    lik <- sum(y * a - log(1 + exp(a)))
  } else if (family == "poisson") {
    lik <-  sum(y * a - exp(a))
  } else {
    lik <-  sum(y * a - a^2 / 2)
  }
  return(lik)
}


penalty0_k <- function(X_, y, theta_, k, family = c("binomial", "poisson", "gaussian")){
  n <- length(y)
  a <- X_ %*% theta_
  h <- c(1, H(theta_[-1], k))
  T <- h %*% t(h)
  diag(T) <- 1
  AN <- BN <- matrix(numeric(n*n),n,n)
  if (family == "binomial") {
    eps_1 <- 1 / (2 * n)
    eps_2 <- 1 - eps_1
    p <- 1 / (1 + exp(-a))
    p[p<eps_1] <- eps_1
    p[p>eps_2] <- eps_2
    diag(AN) <- (y - p)^2
    diag(BN) <- p * (1 - p)
  } else if (family == "poisson") {
    Y_max <- max(y)
    p <- exp(a)
    p[p>Y_max] <- Y_max
    diag(AN) <- (y - p)^2
    diag(BN) <- p
  } else {
    diag(AN) <- (y - a)^2
    diag(BN) <- rep(1, n)
  }

  B <- t(X_) %*% BN %*% X_ * T
  B <- solve(B)
  A <- t(X_) %*% AN %*% X_ * T

  pen0 <- 0.5 * sum(diag(B %*% A) * h)
  return(pen0)
}


penalty0 <- function(X_, y, theta_, family = c("binomial", "poisson", "gaussian")){
  n <- length(y)
  a <- X_ %*% theta_
  AN <- BN <- matrix(numeric(n*n),n,n)
  if (family == "binomial") {
    eps_1 <- 1 / (2 * n)
    eps_2 <- 1 - eps_1
    p <- 1 / (1 + exp(-a))
    p[p<eps_1] <- eps_1
    p[p>eps_2] <- eps_2
    diag(AN) <- (y - p)^2
    diag(BN) <- p * (1 - p)
  } else if (family == "poisson") {
    Y_max <- 2 * max(y)
    p <- exp(a)
    p[p>Y_max] <- Y_max
    diag(AN) <- (y - p)^2
    diag(BN) <-  p
  } else {
    diag(AN) <- (y - a)^2
    diag(BN) <- rep(1, n)
  }

  B <- t(X_) %*% BN %*% X_
  B <- solve(B)
  A <- t(X_) %*% AN %*% X_

  pen0 <- 0.5 * sum(diag(B %*% A))
  return(pen0)
}


penalty1 <- function(theta_, k){
  #  h <- c(1, H(theta_[-1], k))
  h_ <- H(theta_[-1], k)
  pen1 <- sum(h_)
  return(pen1)
}
