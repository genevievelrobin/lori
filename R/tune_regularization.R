#' Computes the threshold $\lambda_{\text{QUT}}$ with parametric bootstrap when covariates are available.
#' If you don't have any covariates, use the function \code{lambda_QUT} which will be significantly faster.
#' @param Y A matrix of counts (contingency table).
#' @param projection A projection function on the space orthogonal to covariates. By default centers by rows and columns
#' @param quantile A number between \code{0} and \code{1}. The quantile of the distribution of $\lambda_{\text{QUT}}$ to take.
#' @param n An integer. The number of parametric bootstrap samples to draw.
#' @return the value of $\lambda_{\text{QUT}}$ to use in LoRI.
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' lambda = lambda_QUT_covariates(Y)
lambda_QUT_covariates = function(Y, projection = default_projection, q = 0.95, n = 1e4){
  m1 = nrow(Y)
  m2 = ncol(Y)
  W = Y
  W[W < 1] = 1e-6
  null_estimator = admm_algorithm(W, cov = T, lambda = 1e5, projection = projection, max_it = 1e4)
  X_0 = null_estimator$X - null_estimator$Theta
  lambdas = rep(0, n)
  proc = Sys.time()
  for(i in 1:n)
  {
    print(i)
    Y_simul = matrix(rpois(n = m1 * m2, exp(c(X_0))), nrow = m1)
    Y_simul[Y_simul <= 0] = 1e-6
    null_estimator_simul = admm_algorithm(Y_simul, cov = TRUE, lambda=1e5, projection = projection, max_it = 1e4)
    X_0_simul = null_estimator_simul$X - null_estimator_simul$Theta
    lambdas[i] = (1 / (m1 * m2)) * svd::propack.svd(projection(Y_simul - exp(X_0_simul)), neig = 1, opts = list(maxiter = 1e5))$d
  }
  Sys.time() - proc
  return(quantile(lambdas, q)[[1]])
}

#' Computes the threshold $\lambda_{\text{QUT}}$ with parametric bootstrap when  NO covariates are available.
#' If you don't have any covariates, use this function instead of \code{lambda_QUT_covariates} which will be significantly slower.
#' @param Y A matrix of counts (contingency table).
#' @param quantile A number between \code{0} and \code{1}. The quantile of the distribution of $\lambda_{\text{QUT}}$ to take.
#' @param n An integer. The number of parametric bootstrap samples to draw.
#' @return the value of $\lambda_{\text{QUT}}$ to use in LoRI.
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' lambda = lambda_QUT(Y)

lambda_QUT = function(Y, q = 0.95, n = 1e4){
  m1 = nrow(Y)
  m2 = ncol(Y)
  W = Y
  W[W < 1] = 1e-6
  null_estimator = estimate_null_model(W)
  X_0 = null_estimator$X
  lambdas = rep(0,n)
  for(i in 1:n)
  {
    Y_simul = matrix(rpois(n = m1 * m2, exp(c(X_0))), nrow = m1)
    Y_simul[Y_simul <= 0] = 1e-6
    null_estimator_simul = estimate_null_model(Y_simul)
    lambdas[i] = (1 / (m1 * m2)) * svd::propack.svd(default_projection(Y_simul - exp(null_estimator_simul$X)),
                                             neig = 1, opts = list(maxiter = 1e5))$d
  }
  return(quantile(lambdas, q)[[1]])
}


#' Selects the parameter $\lambda_{\text{CV}}$ with cross-validation.
#'
#' @param Y A matrix of counts (contingency table).
#' @param cov A boolean \code{TRUE} if covariate matrices are provided, \code{FALSE} otherwise. Default is \code{FALSE}.
#' @param projection A projection function on the space orthogonal to covariates. By default centers by rows and columns
#' @param X_init A matrix. Initial Poisson parameter matrix, same size as Y.
#' @param Theta_init A matrix. Initial interaction matrix, same size as X_init.
#' @param gamma_init A matrix. Initial dual variable, same size as X_init.
#' @param tau A number (augmented Lagrangian parameter).
#' @param epsilon A number. Convergence tolerance of ADMM, by default \code{1e-6}.
#' @param tol A number. Convergence tolerance of gradient descent in ADMM iteration in \code{X}. By default\code{1e-12}.
#' @param max_it An integer. Maximum allowed number of iterations.
#' @param upper upper bound on the values of \code{X}.
#' @param lower lower bound on the values of \code{X}.
#' @param K An integer. The number of folds of the cross-validation.
#' @return the value of $\lambda_{\text{QUT}}$ to use in LoRI.
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' lambda = lambda_QUT(Y)
lambda_cv = function(Y, cov = FALSE, projection = default_projection, gamma_init = NULL, X_init = NULL,
                     Theta_init = NULL, tau = 0.1, epsilon = 1e-6,
                     tol = 1e-12, max_it = 5 * 1e5, upper = -log(1e-6), lower = log(1e-6), K = 10)
{
  null_estimator=estimate_null_model(Y)
  m1 = nrow(Y)
  m2 = ncol(Y)
  mu = null_estimator$mu
  mu = matrix(mu, nrow = m1, ncol = m2)
  alpha = null_estimator$alpha
  alpha = matrix(rep(alpha,m2), nrow = m1,ncol = m2)
  beta = null_estimator$beta
  beta = matrix(rep(beta, m1), nrow = m1, ncol = m2, byrow = TRUE)
  Y_na_rm = Y
  Y_na_rm[is.na(Y)] = 0
  n = sum(!is.na(Y))
  lambda_max = (1 / (m1 * m2)) * max(svd(projection(Y_na_rm - exp(mu + alpha + beta)))$d)
  lambda_grid = exp(seq(log(lambda_max), log(lambda_max / 1e3), length.out = 10))
  Pi = c(Y / sum(Y))
  error = rep(0, K)
  for(i in 1:K)
  {
    N = floor(8 * (n / 10))
    R = rep(0 , m1 * m2)
    R[which(!is.na(Y))[sample(1:n, N)]] = 1
    R = matrix(R, nrow = m1)
    Y_sample = Y
    Y_sample[R == 0] = NA
    estimator_list = list()
    estimator_list[[1]] = admm_algorithm(Y_sample, cov = cov, lambda = NULL, projection, gamma_init, X_init, Theta_init, tau, epsilon, tol, max_it, upper, lower)
    indices_to_predict = 1*(!is.na(Y)) - 1*(!is.na(Y_sample))
    error[1] = error[1] + norm(exp(estimator_list[[1]]$X)[indices_to_predict > 0]-Y[indices_to_predict>0], type="2")
    for(k in 2:(length(lambda_grid) - 1)){
      estimator_list[[k]] = admm_algorithm(Y_sample, cov = cov, lambda = lambda_grid[k], projection, gamma_init, X_init, Theta_init, tau, epsilon, tol, max_it, upper, lower)
      error[k] = error[k] + norm(exp(estimator_list[[k]]$X)[indices_to_predict > 0] - Y[indices_to_predict > 0], type="2")
    }
    
  }
  return(lambda_grid[which(error == min(error))])
}
