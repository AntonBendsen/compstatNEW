
gradient_descent <- function(theta_init, fn, gr, tol = 1e-5, max_iter = 1000,
                             gamma_0 = 1.0, c = 0.05, d = 0.5, cb = NULL) {
  theta <- theta_init
  gamma <- gamma_0
  if (!is.null(cb)) cb()

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
                          gamma_0 = 1.0, c = 0.05, d = 0.5, cb = NULL) {
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

sgd <- function(theta_init, x, y, grad,
                gamma_0 = 0.5, K = 2500, a = 1,
                batch_size = 200, max_epochs = 500,
                cb = NULL) {
  theta <- theta_init
  n <- length(y)
  step <- 0
  if (!is.null(cb)) cb() # Record initial theta

  for (epoch in seq_len(max_epochs)) {
    shuffled_indices <- sample.int(n)

    for (start in seq.int(1, n, by = batch_size)) {
      step <- step + 1
      learning_rate <- gamma_0 * K / (K + step^a)
      batch_indices <- shuffled_indices[
        start:min(start + batch_size - 1L, n)
      ]
      batch_grad <- grad(theta, x[batch_indices], y[batch_indices])

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

  p <- plogis(-beta * log(x) + alpha)
  gamma + (delta - gamma) * p
}

theta_true <- c(
  alpha = 2.5,
  beta  = 1.4,
  gamma = 0.5,
  delta = 4.0
)

N <- 5000
sd = 0.25
log_grid <- seq(1, 5, length.out = 15)
set.seed(123)
x <- exp(sample(log_grid, N, replace = TRUE))
set.seed(123)
y <- dose_response(x, theta_true) +
  rnorm(N, mean = 0, sd = sd)

loss_func <- function(theta, x, y) {
  residual <- y - dose_response(x, theta)
  sum(residual^2) / (2 * length(y))
}

gradient_func <- function(theta, x, y) {
  alpha <- theta[[1]]
  beta  <- theta[[2]]
  gamma <- theta[[3]]
  delta <- theta[[4]]

  lx <- log(x)
  p <- plogis(alpha - beta * lx)
  q <- p * (1 - p)

  residual <- y - (gamma + (delta - gamma) * p)
  rq <- residual * q

  -c(
    (delta - gamma) * sum(rq),
    -(delta - gamma) * sum(rq * lx),
    sum(residual * (1 - p)),
    sum(residual * p)
  ) / length(y)
}

hess_full <- function(theta, x, y) {
  alpha <- theta[[1]]
  beta  <- theta[[2]]
  gamma <- theta[[3]]
  delta <- theta[[4]]

  lx <- log(x)
  p <- plogis(alpha - beta * lx)
  q <- p * (1 - p)
  q_prime <- q * (1 - 2 * p)

  residual <- y - (gamma + (delta - gamma) * p)
  # Each row contains the dose-response gradient for one observation.
  J <- cbind(
    alpha = (delta - gamma) * q,
    beta = -(delta - gamma) * lx * q,
    gamma = 1 - p,
    delta = p
  )

  # Sum of residual_i * Hessian(f_i), required for the exact loss Hessian.
  curvature <- residual * (delta - gamma) * q_prime
  rq <- residual * q
  correction <- matrix(0, nrow = 4, ncol = 4)
  correction[1, 1] <- sum(curvature)
  correction[1, 2] <- correction[2, 1] <- -sum(lx * curvature)
  correction[2, 2] <- sum(lx^2 * curvature)
  correction[1, 3] <- correction[3, 1] <- -sum(rq)
  correction[1, 4] <- correction[4, 1] <- sum(rq)
  correction[2, 3] <- correction[3, 2] <- sum(lx * rq)
  correction[2, 4] <- correction[4, 2] <- -sum(lx * rq)

  # Hessian of sum(residual^2) / (2 * n).
  (crossprod(J) - correction) / length(y)
}

fn <- function(theta) loss_func(theta, x, y)
gr <- function(theta) gradient_func(theta, x, y)
hess <- function(theta) hess_full(theta, x, y)

# gradient_descent(c(2.5, 1.5, 0.3, 4), fn, gr)
# damped_newton(c(3, 2, 0, 2), fn, gr, hess)
# stochastic_gradient_descent(c(2.5, 1.5, 0.3, 4), x, y, gradient_func)

# 6. Gradient Descent
library(CSwR)
library(tidyverse)

compare_optimizers <- function(theta_init, N, sigma, bs_vec = 200, ...) {
  if (!is.numeric(bs_vec) || length(bs_vec) == 0L ||
      any(!is.finite(bs_vec)) || any(bs_vec < 1 | bs_vec != floor(bs_vec))) {
    stop("bs_vec must be a non-empty vector of positive whole numbers.")
  }
  bs_vec <- unique(bs_vec)
  sgd_labels <- paste0("SGD (batch size = ", bs_vec, ")")

  # Keep the tracer's parameter columns consistent from the initial state.
  names(theta_init) <- c("alpha", "beta", "gamma", "delta")
  # 15 grid points on the log-dose scale
  log_grid <- seq(1, 5, length.out = 15)
  x <- exp(sample(log_grid, N, replace = TRUE))
  y <- dose_response(x, theta_true) +
    rnorm(N, mean = 0, sd = sigma)

  fn <- function(theta) loss_func(theta, x, y)
  gr <- function(theta) gradient_func(theta, x, y)
  hess <- function(theta) hess_full(theta, x, y)

  theta_ref <- damped_newton(
    theta_true, fn, gr, hess, tol = 1e-12, max_iter = 1000
  )$theta
  reference_loss <- fn(theta_ref)

  gd_tracer <- tracer(c("theta"), Delta = 0)
  gradient_descent(theta_init, fn, gr, cb = gd_tracer$tracer)
  dn_tracer <- tracer(c("theta"), Delta = 0)
  damped_newton(theta_init, fn, gr, hess, cb = dn_tracer$tracer)

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
  sgd_summary <- bind_rows(lapply(bs_vec, function(batch_size) {
    sgd_tracer <- tracer(c("theta"), Delta = 0)
    stochastic_gradient_descent(
      theta_init, x, y, gradient_func, batch_size = batch_size,
      cb = sgd_tracer$tracer, ...
    )
    add_loss_error(summary(sgd_tracer)) %>%
      mutate(
        batch_size = batch_size,
        algorithm = paste0("SGD (batch size = ", batch_size, ")")
      )
  }))

  plot_data <- bind_rows(
    transform(gd_summary, algorithm = "Gradient descent"),
    transform(dn_summary, algorithm = "Damped Newton"),
    sgd_summary
  )
  algorithm_levels <- c("Gradient descent", "Damped Newton", sgd_labels)
  plot_data$algorithm <- factor(plot_data$algorithm, levels = algorithm_levels)
  algorithm_colors <- c(
    "Gradient descent" = "black",
    "Damped Newton" = "firebrick",
    setNames(hcl.colors(length(bs_vec), palette = "Dark 3"), sgd_labels)
  )

  convergence_plot <- ggplot(
    plot_data,
    aes(x = .time, y = loss_error, color = algorithm, group = algorithm)
  ) +
    geom_line(linewidth = 1) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = algorithm_colors) +
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

set.seed(123)
comparison_1000 <- compare_optimizers(
  c(2, 1, 3, 5),
  N = 1000,
  sigma = 0.25,
  bs_vec = c(1, 50, 200)
)

comparison_1000$convergence_plot

set.seed(123)
comparison_10000 <- compare_optimizers(
  c(2, 1, 3, 5),
  N = 10000,
  sigma = 0.25,
  bs_vec = c(1, 50, 200)
)

comparison_10000$convergence_plot

set.seed(123)
comparison_100000 <- compare_optimizers(
  c(2, 1, 3, 5),
  N = 100000,
  sigma = 0.25,
  bs_vec = c(1, 50, 200)
)

comparison_100000$convergence_plot
