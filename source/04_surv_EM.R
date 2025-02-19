EM_exp <- function(times, status, tol = 1e-6, max_iter = 1000) {
  lambda <- mean(times[status == 1])  
  n <- length(times)
  iter <- 0
  diff <- Inf
  
  while (diff > tol && iter < max_iter) {
    
    # M-step: Update lambda
    lambda_new <- sum(times+(1-status)*lambda)/n
    
    # Check convergence
    diff <- abs(lambda_new - lambda)
    iter <- iter + 1
    
    lambda = lambda_new
  }
  
  list(lambda = lambda, iterations = iter, converged = (iter < max_iter))
}