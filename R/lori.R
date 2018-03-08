#' Main function to be used to fit the LORI model
#' @param Y a matrix of counts (m1 x m2).
#' @param R a matrix of row covariates (m1 x K1)
#' @param C a matrix of column covariates (m2 x K2)
#' @param lambda a positive number, the regularization parameter
#' @param projection a projection function, by default centers by rows and columns
#' @param Theta_init an initial value for the matrix of interactions (same size as Y).
#' @param gamma_init an initial value for the dual matrix (dual variable, same size as Y).
#' @param X_init an initial value for the parameter matrix (dual variable, same size as Y).
#' @param tau a positive number (augmented Lagrangian parameter)
#' @param epsilon a positive number, the convergence criterion
#' @param tol a positive number, the convergence criterion for optimization step inside iterations
#' @param max_it a positive number, the maximum number of iterations
#' @param upper a number, the upper bound on the entries of matrix X
#' @param lower a number, the lower bound on the entries of matrix X
#' @return The value of the augmented Lagrangian (with fixed \code{Theta} and \code{gamma}) taken at \code{X}.
#' @export
#' @examples
#' X = matrix(rnorm(rep(0, 15)), 5)
#' Y <- matrix(rpois(length(c(X)), exp(c(X))), 5)
#' res_lori <- lori(Y)
lori = function(Y, R = NULL, C = NULL, lambda = NULL, projection = default_projection, gamma_init = NULL, X_init = NULL,
                  Theta_init = NULL, tau = 0.1, epsilon = 1e-6,
                  tol = 1e-12, max_it = 5 * 1e5, upper = -log(1e-6), lower = log(1e-6)){
  Y = as.matrix(Y)
  m1 = nrow(Y)
  m2 = ncol(Y)
  n = sum(!is.na(Y))
  cov = FALSE
  if(!is.null(R) || !is.null(C)){
    cov = TRUE
    R=qr(as.matrix(R))
    R=qr.Q(R)
    C=qr(as.matrix(C))
    C=qr.Q(C)
    projection = function(X){
      return(X - R %*% t(R) %*% X - X %*% C %*% t(C) + R %*% t(R) %*% X %*% C %*% t(C))
    }
  }
  res = admm_algorithm(Y, cov, lambda, projection, gamma_init, X_init, Theta_init, tau, epsilon, tol, max_it, upper, lower)
  X = res$X
  colnames(X) <- colnames(Y)
  rownames(X) <- rownames(Y)
  if(is.null(R)){
    R = (1 / sqrt(m1)) * matrix(1, m1, 1)
  }
  if(is.null(C)){
    C = (1 / sqrt(m2)) * matrix(1, 1, m2)
  }
  if(length(C) == ncol(X)) C = t(C)
  mu = solve(t(R) %*% R) %*% t(R) %*% X %*% C %*% solve(t(C) %*% C)
  R_mu_C = R %*% mu %*% t(C)
  beta = solve(t(R) %*% R) %*% t(R) %*% (X - R_mu_C)
  R_beta = R %*% beta
  alpha = (X - R_mu_C) %*% C %*% solve(t(C) %*% C)
  alpha_C = alpha %*% t(C)
  Theta <- res$Theta
  colnames(Theta) <- colnames(Y)
  rownames(Theta) <- rownames(Y)
  colnames(R_mu_C) <- colnames(Y)
  rownames(R_mu_C) <- rownames(Y)
  colnames(R_beta) <- colnames(Y)
  rownames(R_beta) <- rownames(Y)
  colnames(alpha_C) <- colnames(Y)
  rownames(alpha_C) <- rownames(Y)
  return(structure(list(X = X, Theta = Theta, d = res$d, gamma = res$gamma, mu = mu, alpha = alpha, beta = beta,
                        R_mu_C = R_mu_C, R_beta = R_beta, alpha_C = alpha_C,
                        objective = res$objective,
                        iter = res$count, rank = res$rank, convergence = res$iter < max_it)))


}
