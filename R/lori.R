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
                        objective = res$objective, deviance = res$deviance,
                        iter = res$count, rank = res$rank, convergence = res$iter < max_it)))


}
