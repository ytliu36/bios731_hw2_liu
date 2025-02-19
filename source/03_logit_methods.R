newton_logit = function(X, y, beta0 = c(0,0), tol = 1e-6, max_iter = 100) {
  
  beta_cur = beta0
  beta_history = gradient_vec = matrix(NA, nrow = max_iter, 
                                       ncol = length(beta0))
  for (iter in 1:max_iter) {
    
    # store results
    beta_history[iter,] = beta_cur
    
    # Compute the gradient and hessian
    gradient = as.numeric(t(y-exp(X%*%beta_cur)/(1+exp(X%*%beta_cur)))%*%X)
    
    hessian_ls = as.list(rep(NA, length.out = length(y)))
    for(i in 1:length(y)){
      hessian_ls[[i]] = as.numeric(exp(X[i,]%*%beta_cur)/(1+exp(X[i,]%*%beta_cur))^2) * tcrossprod(X[i,], X[i,])
    }
    
    hessian <- -1 * Reduce("+", hessian_ls)
    
    gradient_vec[iter,] = gradient
    
    # Check stopping criterion
    if(sqrt(sum(gradient^2)) < tol){
      message("Converged in", iter, "iterations.\n")
      break
    }
    
    # Update the solution
    beta_cur = beta_cur - solve(hessian) %*% gradient
  }
  hessian_ls = as.list(rep(NA, length.out = length(y)))
  for(i in 1:length(y)){
    hessian_ls[[i]] = as.numeric(exp(X[i,]%*%beta_cur)/(1+exp(X[i,]%*%beta_cur))^2) * tcrossprod(X[i,], X[i,])
  }
  
  var <- solve(Reduce("+", hessian_ls))
  
  return(list(solution = beta_cur, 
              beta_history = beta_history,
              gradient = gradient_vec,
              var = var,
              converged = (iter < max_iter),
              niter = iter))
}


mm_logit = function(X, y, beta0 = c(0,0), tol = 1e-6, max_iter = 200) {
  
  beta_cur = beta0
  beta_history = gradient_vec = matrix(NA, nrow = max_iter, 
                                       ncol = length(beta0))
  inner_iter <- 0
  for (iter in 1:max_iter) {
    
    # store results
    beta_history[iter,] = beta_cur
    eta <- X %*% beta_cur
    mu <- exp(eta) / (1 + exp(eta))
    
    # set a value for beta_0 because X_i1 = 1 all the time
    beta_0 = log(sum(y)/sum((exp(eta)* exp(-2 *beta_cur[1]))/(1 + exp(eta))))/2
    
    # Newton method to calculate beta_1
    inner_diff <- Inf
    beta_1 = beta_cur[2]
    
    while (inner_diff > tol) {
      exp_term = exp(-2 * X[, 2] * beta_cur[2])
      denom = sum((mu * X[, 2] * exp_term) * exp(2 * X[, 2] * beta_1))
      num = sum(y * X[, 2])
      
      # Compute gradient and Hessian for Newton's update
      f_beta_1 = -denom + num
      f_prime_beta_1 = -2 * sum((mu * X[, 2]^2 * exp_term) * exp(2 * X[, 2] * beta_1))
      
      # Newton update
      beta_1_new = beta_1 - f_beta_1 / f_prime_beta_1
      
      inner_diff = abs(beta_1_new - beta_1)
      beta_1 = beta_1_new
      inner_iter = inner_iter + 1
    }
    
    # Check stopping criterion
    if(sqrt(sum((c(beta_0,beta_1)-beta_cur)^2)) < tol){
      message("Converged in", iter, "iterations.\n")
      break
    }
    
    # Update the solution
    beta_cur = c(beta_0,beta_1)
  }
  
  
  # Compute Fisher information matrix
  eta <- X %*% beta_cur
  mu <- exp(eta) / (1 + exp(eta))
  W <- diag(as.vector(mu * (1 - mu)), nrow = length(y))
  fisher_info <- t(X) %*% W %*% X
  var_cov_matrix <- solve(fisher_info)
  
  return(list(solution = beta_cur, 
              beta_history = beta_history,
              converged = (iter < max_iter),
              niter = iter,
              inner_iter = inner_iter,
              var = var_cov_matrix))
  
}


