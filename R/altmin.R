altmin <-
  # internal function to do alternating minimization
  function(Y,
           lambda1,
           lambda2,
           cov = NULL,
           rank.max = 2,
           thresh = 1e-5,
           maxit = 100,
           trace.it = T,
           intercept = F,
           reff = T,
           ceff = T,
           parallel=F)
  {
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(cov)
    rank.max <- min(n - 1, p - 1, rank.max)
    Omega <- 1 * (!is.na(Y))
    m <- sum(Omega)
    mu <- 0
    theta <- matrix(0, nrow = n, ncol = p)
    alpha <- rep(0, n)
    beta <- rep(0, p)
    if (!is.null(cov)){
      epsilon <- rep(0, ncol(cov))
    } else{
      epsilon <- 0
    }
    error = 1
    objective <- NULL
    count = 1
    Y2 <- Y
    Y2[is.na(Y)] <- 0
    X <- mu + matrix(cov %*% epsilon, nrow = n) + theta
    X <- sweep(X, 1, alpha, "+")
    X <- sweep(X, 2, beta, "+")
    Ytalp <- c(Y[!is.na(Y)])
    alpha.step <- 1
    theta.step <- 1
    while (error > thresh && count < maxit) {
      mu.tmp <- mu
      alpha.tmp <- alpha
      beta.tmp <- beta
      epsilon.tmp <- epsilon
      d.tmp <- svd(theta, nu=rank.max, nv=rank.max)$d
      if(!intercept) mu <- 0 else mu <- log(sum(Y, na.rm=T)/sum(Omega*exp(X-mu.tmp)))
      if(!reff) alpha <- rep(0,n)
      if(!ceff) beta <- rep(0,p)
      grad <- grad(Y, cov, mu, alpha, beta, epsilon, theta)
      if(!reff) grad[1:n] <- 0
      if(!ceff) grad[(n+1):(n+p)] <- 0
      flag <- TRUE
      alpha.step <- 1
      Xtmp <- mu + theta + matrix(cov%*% epsilon, nrow=n)
      Xtmp <- sweep(Xtmp, 1, alpha, "+")
      Xtmp <- sweep(Xtmp, 2, beta, "+")
      ref.obj <-
        sum((-Y * (Xtmp) + exp(Xtmp)),na.rm = T) + lambda1 * sum(d.tmp) + sum(lambda2 * c(abs(alpha), abs(beta),abs(epsilon)))
      while(flag){
        #print(step)
        maineff <- sign(c(alpha.tmp, beta.tmp, epsilon.tmp)-alpha.step*grad)*pmax(abs(c(alpha.tmp, beta.tmp, epsilon.tmp)-alpha.step*grad)-alpha.step*lambda2, 0)
        alpha <- maineff[1:n]
        beta <- maineff[(n+1):(n+p)]
        epsilon <- maineff[(n+p+1):(n+p+q)]
        Xtmp <- mu + theta + matrix(cov%*% epsilon, nrow=n)
        Xtmp <- sweep(Xtmp, 1, alpha, "+")
        Xtmp <- sweep(Xtmp, 2, beta, "+")
        diff <-
          sum((-Y * (Xtmp) + exp(Xtmp)), na.rm = T)  + lambda1 * sum(d.tmp) + sum(lambda2 * c(abs(alpha), abs(beta),abs(epsilon))) - ref.obj
        flag <- diff > thresh * abs(ref.obj)
        alpha.step <- 0.5*alpha.step
      }
      X <- Xtmp
      grad_theta <- (-Y2 + exp(X))
      grad_theta <- sweep(grad_theta, 2, colMeans(grad_theta))
      grad_theta <- sweep(grad_theta, 1, rowMeans(grad_theta))
      theta.step <- 1
      flag <- TRUE
      ref.obj <-
        sum((-Y * (X) + exp(X)), na.rm = T) + lambda1 * sum(d.tmp) + sum(lambda2 * c(abs(alpha), abs(beta), abs(epsilon)))
      Xtmp <- mu + matrix(cov%*% epsilon, nrow=n)
      Xtmp <- sweep(Xtmp, 1, alpha, "+")
      Xtmp <- sweep(Xtmp, 2, beta, "+")
      while (flag) {
        #print(step)
        svd_theta <- svd(theta - theta.step * grad_theta, nu=rank.max, nv=rank.max)
        d <- svd_theta$d[1:rank.max]
        u <- svd_theta$u[, 1:rank.max]
        v <- svd_theta$v[, 1:rank.max]
        d <- d[d > (lambda1 * theta.step)] - lambda1 * theta.step
        if (length(d) == 0) {
          theta2 <- matrix(rep(0, n * p), nrow = n)
        } else if (length(d) == 1) {
          u <- svd_theta$u[, 1:length(d)]
          v <- svd_theta$v[, 1:length(d)]
          theta2 <- d * u %*% t(v)
        } else {
          u <- svd_theta$u[, 1:length(d)]
          v <- svd_theta$v[, 1:length(d)]
          theta2 <- u %*% diag(d) %*% t(v)
        }
        diff <-
          sum((-Y* (Xtmp + theta2) + exp(Xtmp+ theta2)), na.rm = T) + lambda1 * sum(d) + sum(lambda2 * c(abs(alpha), abs(beta),abs(epsilon)))- ref.obj
        flag <- diff > thresh * abs(ref.obj)
        theta.step <- 0.5 * theta.step
      }
      theta <- theta2
      X <- Xtmp + theta
      objective = c(objective,-sum((Y * X - exp(X)), na.rm = T)+ lambda1 *
                      sum(d) + sum(lambda2 * c(abs(alpha), abs(beta),abs(epsilon))))
      if(count==1){
        error <- 1
      } else{
        error = abs(objective[count] - objective[count - 1]) / abs(objective[count])
      }
      count = count + 1
      epsilon <- as.vector(epsilon)
      names(epsilon) <- colnames(cov)
      rownames(X) <- rownames(Y)
      colnames(X) <- colnames(Y)
      rownames(theta) <- rownames(Y)
      colnames(theta) <- colnames(Y)
    }
    return(structure(
      list(
        X = X,
        mu = mu,
        alpha = alpha,
        beta = beta,
        epsilon = epsilon,
        theta = theta,
        objective = unlist(objective),
        iter = count,
        convergence = count < maxit
      )
    ))
  }
