
gradient_descent <- function(theta_init, fn, gr, tol = 1e-5, max_iter = 1000,
                            gamma_0 = 2.0, c = 0.25, d = 0.75, cb = NULL) {
  theta <- theta_init
  gamma <- gamma_0 # Initialize gamma to avoid errors on early convergence
  if (!is.null(cb)) cb() # Record initial theta at time zero.

  for (i in seq_len(max_iter)) {
    grad <- gr(theta)
    grad_norm <- sqrt(sum(grad^2))

    if (grad_norm <= tol) {
      message("Converged successfully after ", i, " iterations.")
      return(list(theta = theta, gamma = gamma))
    }

    gamma <- gamma_0
    f_curr <- fn(theta)

    # Backtracking line search (Armijo condition)
    while (fn(theta - gamma * grad) > f_curr - c * gamma * grad_norm^2) {
      gamma <- d * gamma

      if (gamma < 1e-16) {
        warning("Line search failed to find a sufficient decrease step.")
        break
      }
    }

    theta <- theta - gamma * grad
    if (!is.null(cb)) cb() # Record theta after the completed update.
  }

  warning("Maximum number of iterations reached without full convergence.")
  return(list(theta = theta, gamma = gamma))
}

damped_newton <- function(theta_init, fn, gr, hess, tol = 1e-5, max_iter = 100,
                          gamma_0 = 1.0, c = 0.25, d = 0.5, cb = NULL) {
  theta <- theta_init
  gamma <- gamma_0
  if (!is.null(cb)) cb() # Record initial theta at time zero.

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

    gamma <- gamma_0 # Restart backtracking at the initial step size each iteration.
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
    if (!is.null(cb)) cb() # Record theta after the completed update.
  }

  warning("Maximum number of iterations reached without full convergence.")
  return(list(theta = theta, gamma = gamma))
}

stochastic_gradient_descent <- function(theta_init, x, y, batch_size = 200,
                                         gamma_0 = 0.5, K = 2500, a = 1,
                                         max_epochs = 1000,
                                         cb = NULL) {
  theta <- theta_init
  n <- length(y)
  step <- 0
  if (!is.null(cb)) cb() # Record initial theta before shuffling or updates.

  for (epoch in seq_len(max_epochs)) {
    shuffled_indices <- sample.int(n)

    for (start in seq.int(1, n, by = batch_size)) {
      step <- step + 1
      learning_rate <- gamma_0 * K / (K + step^a)
      batch_indices <- shuffled_indices[
        start:min(start + batch_size - 1L, n)
      ]
      batch_grad <- numeric(length(theta))

      for (i in batch_indices) {
        batch_grad <- batch_grad + grad_i(theta, x[i], y[i])
      }
      batch_grad <- batch_grad / length(batch_indices)

      theta <- theta - learning_rate * batch_grad
    }
    if (!is.null(cb)) cb()
  }
  list(theta = theta, epochs = max_epochs)
}

dose_response <- function(x, theta) {
  alpha <- theta[1]
  beta  <- theta[2]
  gamma <- theta[3]
  delta <- theta[4]

  p <- plogis(-beta * (log(x) - alpha))

  gamma + (delta - gamma) * p
}

set.seed(123)

theta_true <- c(
  alpha = 2.5,
  beta  = 1.4,
  gamma = 0.5,
  delta = 4.0
)

# 3. Observation-level loss: 1/2 * residual^2
loss_i <- function(theta, x, y) {
  residual <- y - dose_response(x, theta)
  0.5 * residual^2
}

# 4. Gradient of one observation's loss
grad_i <- function(theta, x, y) {
  alpha <- theta[1]
  beta  <- theta[2]
  gamma <- theta[3]
  delta <- theta[4]

  lx <- log(x)
  p <- plogis(-beta * (lx - alpha)) # as plogis is the logistic function so plogis(z) = 1 / (1 + exp(-z))
  q <- p * (1 - p)

  fit <- gamma + (delta - gamma) * p
  residual <- y - fit

  df <- c(
    (delta - gamma) * beta * q,
    -(delta - gamma) * (lx - alpha) * q,
    1 - p,
    p
  )

  -residual * df
}

full_loss <- function(theta, x, y) {
  mean(vapply(
    seq_along(y),
    function(i) loss_i(theta, x[i], y[i]),
    numeric(1)
  ))
}

full_gradient <- function(theta, x, y) {
  g <- numeric(4)

  for (i in seq_along(y)) {
    g <- g + grad_i(theta, x[i], y[i])
  }

  g / length(y)
}

hess_i <- function(theta, x, y) {
  alpha <- theta[1]
  beta  <- theta[2]
  gamma <- theta[3]
  delta <- theta[4]

  u <- log(x) - alpha
  p <- plogis(-beta * u)
  q <- p * (1 - p)
  q_prime <- q * (1 - 2 * p)

  residual <- y - dose_response(x, theta)
  df <- c(
    (delta - gamma) * beta * q,
    -(delta - gamma) * u * q,
    1 - p,
    p
  )

  hess_f <- matrix(0, nrow = 4, ncol = 4)
  hess_f[1, 1] <- (delta - gamma) * q_prime * beta^2
  hess_f[1, 2] <- hess_f[2, 1] <-
    (delta - gamma) * (q - beta * u * q_prime)
  hess_f[2, 2] <- (delta - gamma) * u^2 * q_prime
  hess_f[1, 3] <- hess_f[3, 1] <- -beta * q
  hess_f[1, 4] <- hess_f[4, 1] <- beta * q
  hess_f[2, 3] <- hess_f[3, 2] <- u * q
  hess_f[2, 4] <- hess_f[4, 2] <- -u * q

  tcrossprod(df) - residual * hess_f
}

full_hessian <- function(theta, x, y) {
  H <- matrix(0, nrow = 4, ncol = 4)

  for (i in seq_along(y)) {
    H <- H + hess_i(theta, x[i], y[i])
  }

  H / length(y)
}


# 6. Gradient Descent
library(CSwR)
library(tidyverse)

compare_optimizers <- function(theta_init, N, sigma, ...) {
  # 15 grid points on the log-dose scale
  log_grid <- seq(1, 5, length.out = 15)
  x <- exp(sample(log_grid, N, replace = TRUE))
  y <- dose_response(x, theta_true) +
    rnorm(N, mean = 0, sd = sigma)

  fn <- function(theta) full_loss(theta, x, y)
  gr <- function(theta) full_gradient(theta, x, y)
  hess <- function(theta) full_hessian(theta, x, y)

  theta_ref <- damped_newton(
    theta_true, fn, gr, hess, tol = 1e-12, max_iter = 1000
  )$theta
  reference_loss <- fn(theta_ref)

  gd_tracer <- tracer(c("theta"), Delta = 0)
  gradient_descent(theta_init, fn, gr, cb = gd_tracer$tracer)
  dn_tracer <- tracer(c("theta"), Delta = 0)
  damped_newton(theta_init, fn, gr, hess, cb = dn_tracer$tracer)
  sgd_tracer <- tracer(c("theta"), Delta = 0)
  stochastic_gradient_descent(
    theta_init, x, y, cb = sgd_tracer$tracer, ...
  )

  add_loss_error <- function(summary_data) {
    summary_data %>%
      rowwise() %>%
      mutate(loss_error = abs(
        fn(c_across(starts_with("theta."))) - reference_loss
      )) %>%
      ungroup()
  }

  gd_summary <- add_loss_error(summary(gd_tracer))
  dn_summary <- add_loss_error(summary(dn_tracer))
  sgd_summary <- add_loss_error(summary(sgd_tracer))

  plot_data <- rbind(
    transform(gd_summary, algorithm = "Gradient descent"),
    transform(dn_summary, algorithm = "Damped Newton"),
    transform(sgd_summary, algorithm = "Stochastic gradient descent")
  )

  convergence_plot <- ggplot(
    plot_data,
    aes(x = .time, y = loss_error, color = algorithm, group = algorithm)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_y_log10() +
    scale_color_manual(values = c(
      "Gradient descent" = "steelblue",
      "Damped Newton" = "firebrick",
      "Stochastic gradient descent" = "darkgreen"
    )) +
    labs(
      x = "Time",
      y = "|H(theta) - H(theta_ref)|",
      color = "Algorithm"
    ) +
    theme_minimal()

  list(
    gd_summary = gd_summary,
    dn_summary = dn_summary,
    sgd_summary = sgd_summary,
    convergence_plot = convergence_plot,
    theta_ref = theta_ref,
    reference_loss = reference_loss
  )
}

comparison_2000 <- compare_optimizers(
  c(2, 1, 3, 5), N = 2000, sigma = 0.25, batch_size = 50
)

comparison_2000$convergence_plot

comparison_10000 <- compare_optimizers(
  c(2, 1, 3, 5), N = 10000, sigma = 0.25, batch_size = 50
)

comparison_10000$convergence_plot