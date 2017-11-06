#' Computes the value of the objective function in X to be optimized at each iteration of ADMM.
#'
#' @param X A matrix of Poisson parameters.
#' @param Y A matrix of counts (same size as X).
#' @param projection A projection function - identifiability constraint - by default centers by rows and columns
#' @param Theta A matrix interaction (same size as X).
#' @param gamma A matrix (dual variable, same size as X).
#' @param tau A number (augmented Lagrangian parameter)
#' @return The value of the augmented Lagrangian (with fixed \code{Theta} and \code{gamma}) taken at \code{X}.
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Theta = default_projection(matrix(rnorm(rep(0, 15)), 5))
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' gamma = matrix(rnorm(rep(0, 15)), 5)
#' tau = 1e-3
#' objective_function(X, Y,  projection = default_projection, Theta, gamma, tau)
objective_function = function(X, Y,  projection = default_projection, Theta, gamma, tau) {
  m1 = nrow(Y)
  m2 = ncol(Y)
  Omega = 1-1*is.na(Y)
  n = sum(Omega > 0)
  X = matrix(X, nrow = m1, ncol = m2)
  X_projected = projection(X)
  Y_rm_na = Y
  Y_rm_na[is.na(Y_rm_na)] = 0
  return(-1 / (n) * sum(Y_rm_na * X - Omega * exp(X)) + tr(t(gamma) %*% X_projected) +
           (tau / 2) * norm(X_projected - Theta, type="F")^2)
}

#' Estimates the GAMMIT model under the constraint Theta = 0.
#'
#' @param Y A matrix of counts (contingency table).
#' @return the null estimator under constraint Theta = 0
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' res = estimate_null_model(Y)
estimate_null_model = function(Y){
  Y[Y == 0] = 1e-6
  m1 = nrow(Y)
  m2 = ncol(Y)
  Omega = 1-1*is.na(Y)
  mu = (1 / m1) * sum(log(rowSums(Y, na.rm = TRUE))) + (1 / m2) * sum(log(colSums(Y, na.rm = TRUE))) -
    log(sum(Y, na.rm = TRUE))
  alpha = log(rowSums(Y, na.rm = TRUE)) - (1 / m1) * sum(log(rowSums(Y, na.rm = TRUE)))
  beta = log(colSums(Y, na.rm = TRUE)) - (1 / m2) * sum(log(colSums(Y, na.rm = TRUE)))
  mu_matrix = matrix(mu, nrow = m1, ncol = m2)
  alpha_matrix = matrix(rep(alpha, m2),nrow = m1, ncol = m2)
  beta_matrix = matrix(rep(beta, m1), nrow = m1, ncol = m2, byrow = TRUE)
  return(structure(list(mu = mu, alpha = alpha, beta = beta, X = mu_matrix + alpha_matrix + beta_matrix)))
}

#' Computes deviance of Poisson parameter matrix with count data matrix.
#'
#' @param X A Poisson parameter matrix.
#' @param Y A matrix of counts (contingency table, same size as X).
#' @return the deviance between \code{X} and \code{Y}
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' deviance(Y, X)
deviance = function(Y, X){
  y = as.matrix(Y)
  x = as.matrix(X)
  m = y * log(y / exp(x))
  m[is.na(m)] = 0
  2 * sum(m) + 2 * sum(exp(x) - y, na.rm = TRUE)
}


#' Computes the gradient of the objective function in X for gradient descent step in ADMM.
#'
#' @param X A matrix of Poisson parameters.
#' @param Y A matrix of counts (same size as X).
#' @param projection A projection function - identifiability constraint - by default centers by rows and columns
#' @param Theta A matrix interaction (same size as X).
#' @param gamma A matrix (dual variable, same size as X).
#' @param tau A number (augmented Lagrangian parameter)
#' @return The gradient of the augmented Lagrangian (with fixed \code{Theta} and \code{gamma}) taken at \code{X}.
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Theta = default_projection(matrix(rnorm(rep(0, 15)), 5))
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' gamma = matrix(rnorm(rep(0, 15)), 5)
#' tau = 1e-3
#' gradient(X, Y,  projection = default_projection, Theta, gamma, tau)

gradient = function(X, Y, projection = default_projection, Theta, gamma, tau) {
  m1 = nrow(Y)
  m2 = ncol(Y)
  Omega = 1-1*is.na(Y)
  n = sum(!is.na(Y))
  X = matrix(X, nrow = m1, ncol = m2)
  X_projected = projection(X)
  Y_rm_na = Y
  Y_rm_na[is.na(Y_rm_na)] = 0
  return((-1 / (n)) * (Y_rm_na - Omega * exp(X)) + gamma + tau * (X_projected - Theta))
}

#' Default projection on covariate space: center by row and column.
#'
#' @param X A matrix.
#' @return The projection of X
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' default_projection(X)

default_projection = function(M){
  M_projected = as.matrix(sweep(M, 2, colMeans(M)))
  M_projected = as.matrix(sweep(M_projected, 1, rowMeans(M_projected)))
  return(M_projected)
}

#' Performs the Alternating Descent Method of Multipliers (ADMM) algorithm to estimate the GAMMIT parameters.
#'
#' @param Y A matrix of counts (same size as X).
#' @param cov A boolean \code{TRUE} if covariate matrices are provided, \code{FALSE} otherwise. Default is \code{FALSE}.
#' @param lambda A number, the regularization parameter.
#' @param projection A projection function - identifiability constraint - by default centers by rows and columns.
#' @param X_init A matrix - initial Poisson parameter matrix, same size as Y.
#' @param Theta_init A matrix - initial interaction matrix, same size as X_init.
#' @param gamma_init A matrix - initial dual variable, same size as X_init.
#' @param tau A number (augmented Lagrangian parameter).
#' @param epsilon A number - convergence tolerance of ADMM, by default \code{1e-6}.
#' @param tol A number - convergence of gradient descent in ADMM iteration in \code{X}. By default\code{1e-12}.
#' @param max_it An integer. Maximum allowed number of iterations.
#' @param upper upper bound on the values of \code{X}.
#' @param lower lower bound on the values of \code{X}.
#' @return The value of the augmented Lagrangian (with fixed \code{Theta} and \code{gamma}) taken at \code{X}.
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y = matrix(rpois(length(c(X)), exp(c(X))), 5)
#' res = admm_algorithm(Y)

admm_algorithm = function(Y, cov = FALSE, lambda = NULL, projection = default_projection, gamma_init = NULL, X_init = NULL,
                          Theta_init = NULL, tau = 0.1, epsilon = 1e-6,
                          tol = 1e-12, max_it = 5 * 1e5, upper = -log(1e-6), lower = log(1e-6))
{
  Y = as.matrix(Y)
  m1 = nrow(Y)
  m2 = ncol(Y)
  Omega = 1-1*is.na(Y)
  n = sum(!is.na(Y))
  if(is.null(X_init)){
    X = matrix(0, m1, m2)
    X[which(Omega > 0)] = log(Y[which(Omega > 0)] + 1e-5)
    X[which(Omega <= 0)] = mean(X[which(Omega > 0)])
  } else{
    X = X_init
  }
  X_projected = projection(X)
  if(is.null(lambda)){
    if(cov == FALSE){
      lambda = lambda_QUT(Y, quantile = 0.95, n = 1e3)
    } else{
      lambda = lambda_QUT_covariates(Y, projection = projection, quantile = 0.95, n = 1e3)
    }

  }
  if(is.null(gamma_init)){
    gamma = matrix(0, nrow=m1, ncol=m2)
    gamma = projection(gamma)
  } else{
    gamma = gamma_init
  }
  if(is.null(Theta_init)){
    estim = estimate_null_model(Y)
    mu = estim$mu
    mu = matrix(mu, nrow=m1, ncol=m2)
    alpha = estim$alpha
    alpha = matrix(rep(alpha,m2),nrow=m1,ncol=m2)
    beta = estim$beta
    beta = matrix(rep(beta, m1), nrow = m1, ncol = m2, byrow = TRUE)
    Y_na_rm = Y
    Y_na_rm[is.na(Y_na_rm)] = 0
    Theta = (1 / (n)) * (Y_na_rm - Omega * exp(mu + alpha + beta))
    Theta=projection(Theta)
  } else{
    Theta=Theta_init
  }
  error = 1
  d = svd(Theta)$d
  objective = list(-(1 / (n)) * sum(Y_na_rm * X - Omega * exp(X)) + lambda * sum(d)
                   + tr(t(gamma) %*% (X_projected - Theta)) + (tau / 2) * norm(matrix(X_projected - Theta),
                                                                               type = "F")^2)
  count = 1
  while(error > epsilon && count < max_it){
    X_tmp = X
    X =  optim(par = X_tmp, method="L-BFGS-B", fn = objective_function, gr = gradient,
               Y, projection, Theta, gamma, tau, lower = matrix(lower, nrow = m1, ncol = m2),
               upper = matrix(upper, nrow = m1, ncol = m2), control = list(pgtol = tol,maxit = 1e5))$par
    X_projected = projection(X)
    Theta_tmp = Theta
    X_projected_svd = svd(X_projected + (1 / tau) * gamma)
    d = X_projected_svd$d[X_projected_svd$d > lambda / tau]
    u = X_projected_svd$u[, X_projected_svd$d > lambda / tau]
    v = X_projected_svd$v[, X_projected_svd$d > lambda / tau]
    if(length(d) == 1){
      Theta = u %*% ((d - lambda / tau) * t(v))
    } else{
      Theta = u %*% diag(d - lambda / tau) %*% t(v)
    }
    gamma = gamma + tau * (X_projected - Theta)
    d = svd(Theta)$d
    count = count + 1
    objective[[count]] = -(1 / (m1 * m2)) * sum(Y_na_rm * X - Omega * exp(X)) + lambda * sum(d) +
      tr(t(gamma) %*% (X_projected - Theta)) + (tau / 2) * norm(X_projected - Theta, type = "F")^2
    residual_1 = X_projected - Theta
    residual_2 = Theta_tmp - Theta
    residual_1 = norm(residual_1, type="F")
    residual_2 = norm(residual_2, type="F")
    if(residual_1 > 10 * residual_2)
    {
      tau = 2 * tau
    } else if (residual_2 > 10 * residual_1)
    {
      tau = tau / 2
    }
    error = (residual_1 + residual_2) / 2
    error_1 = norm(X - X_tmp)
    error_2 = objective[[count]] - objective[[count - 1]]
    if(count %% 1e3 == 0){
      print(paste("iteration:", count))
      print(paste("objective step:", error_2))
    }
  }
  X = X - projection(X) + Theta
  mu = mean(X)
  alpha = rowMeans(X - mu)
  beta = colMeans(X - mu)
  rank = sum(svd(Theta)$d > 5e-06)
  deviance = deviance(Y, X)
  return(structure(list(X = X, Theta = Theta, gamma = gamma, mu = mu, alpha = alpha, beta = beta,
                        objective = unlist(objective), deviance = deviance,
                        iter = count, rank = rank, convergence = count < max_it)))
}

