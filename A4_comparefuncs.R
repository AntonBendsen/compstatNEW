
gradient_descent <- function(theta_init, fn, gr, tol = 1e-8, max_iter = 1000, 
                             gamma_0 = 1.0, c = 0.1, d = 0.5) {
  theta <- theta_init

  for (i in seq_len(max_iter)) {
    grad <- gr(theta)
    grad_norm_sq <- sum(grad^2)

    if (sqrt(grad_norm_sq) <= tol) {
      message("Converged successfully after ", i, " iterations.")
      return(list(theta = theta, gamma = gamma))
    }

    gamma <- gamma_0
    f_curr <- fn(theta)

    while (fn(theta - gamma * grad) > f_curr - c * gamma * grad_norm_sq) {
      gamma <- d * gamma

      if (gamma < 1e-16) {
        warning("Line search failed to find a sufficient decrease step.")
        break
      }
    }

    theta <- theta - gamma * grad
  }

  warning("Maximum number of iterations reached without full convergence.")
  return(list(theta = theta, gamma = gamma))
}

damped_newton <- function(theta_init, fn, gr, hess, tol = 1e-8, max_iter = 100, 
                          c = 0.1, d = 0.5) {
  theta <- theta_init

  for (i in seq_len(max_iter)) {

    grad <- gr(theta)
    grad_norm <- sqrt(sum(grad^2))

    if (grad_norm <= tol) {
      message("Converged successfully after ", i, " iterations.")
      return(list(theta = theta, gamma = gamma))
    }

    H_mat <- hess(theta)
    v <- solve(H_mat, -grad)

    if (sum(grad * v) >= 0) {
      warning("Newton direction is not a descent direction. Falling back to negative gradient.")
      v <- -grad
    }

    gamma <- 1.0
    f_curr <- fn(theta)
    directional_derivative <- sum(grad * v) 

    while (fn(theta + gamma * v) > f_curr + c * gamma * directional_derivative) {
      gamma <- d * gamma

      if (gamma < 1e-16) {
        warning("Line search failed to find a sufficient decrease step.")
        break
      }
    }

    theta <- theta + gamma * v
  }

  warning("Maximum number of iterations reached without full convergence.")
  return(list(theta = theta, gamma = gamma))
}



# test_fct <- function(theta) sum((theta - c(2, -1))^2)
# test_fct_prime <- function(theta) 2 * (theta - c(2, -1))
# test_fct_hess <- function(theta) 2 * diag(length(theta))

# gradient_descent(c(0, 0), test_fct, test_fct_prime)
# damped_newton(c(0, 0), test_fct, test_fct_prime, test_fct_hess)
