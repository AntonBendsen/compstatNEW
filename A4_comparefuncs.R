# Refresh definitions without running the examples:
# source("A4_load_functions.R", chdir = TRUE)

gradient_descent <- function(theta_init, fn, gr, tol = 1e-4, max_iter = 1000,
                             gamma_0 = 1.0, c = 0.05, d = 0.5, cb = NULL) {
  theta <- theta_init
  gamma <- gamma_0
  if (!is.null(cb)) cb()

  for (i in seq_len(max_iter)) {
    grad <- gr(theta)
    grad_norm <- sqrt(sum(grad^2))

    if (grad_norm <= tol) {
      message("Converged successfully after ", i - 1, " iterations.")
      return(list(theta = theta, gamma = gamma))
    }

    gamma <- gamma_0
    f_curr <- fn(theta)

    # Backtracking line search (Armijo condition)
    while (fn(theta - gamma * grad) > f_curr - c * gamma * grad_norm^2) {
      gamma <- d * gamma

      if (gamma < 1e-16) {
        warning("Line search failed to find a sufficient decrease step.")
        return(list(theta = theta, gamma = gamma))
      }
    }

    theta <- theta - gamma * grad
    if (!is.null(cb)) cb() # Record theta after the completed update.
  }

  warning("Maximum number of iterations reached without full convergence.")
  return(list(theta = theta, gamma = gamma))
}

damped_newton <- function(theta_init, fn, gr, hess, tol = 1e-4, max_iter = 100,
                          gamma_0 = 1.0, c = 0.05, d = 0.5, cb = NULL) {
  theta <- theta_init
  gamma <- gamma_0
  if (!is.null(cb)) cb() # Record initial theta at time zero.

  for (i in seq_len(max_iter)) {

    grad <- gr(theta)
    grad_norm <- sqrt(sum(grad^2))

    if (grad_norm <= tol) {
      message("Converged successfully after ", i - 1, " iterations.")
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
        return(list(theta = theta, gamma = gamma))
      }
    }

    theta <- theta + gamma * v
    if (!is.null(cb)) cb() # Record theta after the completed update.
  }

  warning("Maximum number of iterations reached without full convergence.")
  return(list(theta = theta, gamma = gamma))
}

sgd <- function(theta_init, x, y, grad,
                gamma_0 = 0.5, K = 100, a = 1,
                batch_size = 100, max_epochs = 500,
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
  list(theta = theta)
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
    sgd(
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
      y = "|L(theta) - L(theta_ref)|",
      color = "Algorithm"
    ) +
    theme_minimal()

  list(
    gd_summary = gd_summary,
    dn_summary = dn_summary,
    sgd_summary = sgd_summary,
    convergence_plot = convergence_plot,
    data = data.frame(x = x, y = y),
    theta_ref = theta_ref,
    reference_loss = reference_loss
  )
}

compare_optimizers_contour <- function(theta_init, N, sigma, bs_vec = 200,
                                       dims = c(1, 2), ...) {
  if (!is.numeric(theta_init) || length(theta_init) != 4L ||
      any(!is.finite(theta_init))) {
    stop("theta_init must contain four finite numeric parameters.")
  }
  if (!is.numeric(bs_vec) || length(bs_vec) == 0L ||
      any(!is.finite(bs_vec)) || any(bs_vec < 1 | bs_vec != floor(bs_vec))) {
    stop("bs_vec must be a non-empty vector of positive whole numbers.")
  }
  if (!is.numeric(dims) || length(dims) != 2L ||
      any(!is.finite(dims)) || any(dims != floor(dims)) ||
      any(dims < 1 | dims > 4) || anyDuplicated(dims)) {
    stop("dims must contain two distinct parameter indices between 1 and 4.")
  }

  parameter_names <- c("alpha", "beta", "gamma", "delta")
  names(theta_init) <- parameter_names
  bs_vec <- unique(bs_vec)
  axis_names <- parameter_names[dims]
  fixed_dims <- setdiff(seq_len(4), dims)
  sgd_labels <- paste0("SGD (batch size = ", bs_vec, ")")

  # The contours and every optimizer use these same observations.
  log_grid <- seq(1, 5, length.out = 15)
  x <- exp(sample(log_grid, N, replace = TRUE))
  y <- dose_response(x, theta_true) + rnorm(N, mean = 0, sd = sigma)
  fn <- function(theta) loss_func(theta, x, y)
  gr <- function(theta) gradient_func(theta, x, y)
  hess <- function(theta) hess_full(theta, x, y)

  theta_ref <- damped_newton(
    theta_true, fn, gr, hess, tol = 1e-12, max_iter = 1000
  )$theta
  names(theta_ref) <- parameter_names

  gd_tracer <- CSwR::tracer("theta", Delta = 0)
  gradient_descent(theta_init, fn, gr, cb = gd_tracer$tracer)
  gd_summary <- summary(gd_tracer)

  dn_tracer <- CSwR::tracer("theta", Delta = 0)
  damped_newton(theta_init, fn, gr, hess, cb = dn_tracer$tracer)
  dn_summary <- summary(dn_tracer)

  sgd_summary <- dplyr::bind_rows(lapply(bs_vec, function(batch_size) {
    sgd_tracer <- CSwR::tracer("theta", Delta = 0)
    sgd(
      theta_init, x, y, gradient_func, batch_size = batch_size,
      cb = sgd_tracer$tracer, ...
    )
    result <- summary(sgd_tracer)
    result$batch_size <- batch_size
    result$algorithm <- paste0("SGD (batch size = ", batch_size, ")")
    result
  }))

  path_data <- dplyr::bind_rows(
    transform(gd_summary, algorithm = "Gradient descent"),
    transform(dn_summary, algorithm = "Damped Newton"),
    sgd_summary
  )
  names(path_data)[match(paste0("theta.", parameter_names), names(path_data))] <-
    parameter_names
  path_data$algorithm <- factor(
    path_data$algorithm,
    levels = c("Gradient descent", "Damped Newton", sgd_labels)
  )

  padded_range <- function(values) {
    limits <- range(values, finite = TRUE)
    padding <- 0.08 * max(diff(limits), 0.1 * max(1, abs(limits)))
    limits + c(-padding, padding)
  }
  x_range <- padded_range(c(path_data[[axis_names[1]]], theta_ref[dims[1]]))
  y_range <- padded_range(c(path_data[[axis_names[2]]], theta_ref[dims[2]]))
  contour_grid <- expand.grid(
    seq(x_range[1], x_range[2], length.out = 80),
    seq(y_range[1], y_range[2], length.out = 80)
  )
  names(contour_grid) <- axis_names

  # Hold the two unplotted parameters at the reference fit.
  # fn() captures the numeric x and y generated above.
  contour_grid$loss <- vapply(seq_len(nrow(contour_grid)), function(i) {
    theta <- theta_ref
    theta[dims] <- c(contour_grid[i, 1], contour_grid[i, 2])
    fn(theta)
  }, numeric(1))

  markers <- as.data.frame(rbind(theta_init[dims], theta_ref[dims]))
  names(markers) <- axis_names
  markers$point <- c("Initial", "Reference fit")
  algorithm_colors <- c(
    "Gradient descent" = "black",
    "Damped Newton" = "firebrick",
    setNames(hcl.colors(length(bs_vec), palette = "Dark 3"), sgd_labels)
  )
  fixed_values <- paste(sprintf(
    "%s = %.3g", parameter_names[fixed_dims], theta_ref[fixed_dims]
  ), collapse = ", ")

  contour_plot <- ggplot2::ggplot() +
    ggplot2::geom_contour(
      data = contour_grid,
      ggplot2::aes(x = .data[[axis_names[1]]], y = .data[[axis_names[2]]], z = loss),
      bins = 20, color = "grey50", linewidth = 0.7
    ) +
    ggplot2::geom_path(
      data = path_data,
      ggplot2::aes(x = .data[[axis_names[1]]], y = .data[[axis_names[2]]],
                   color = algorithm, group = algorithm),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = markers,
      ggplot2::aes(x = .data[[axis_names[1]]], y = .data[[axis_names[2]]],
                   shape = point),
      size = 3, stroke = 1.2, color = "black"
    ) +
    ggplot2::scale_color_manual(values = algorithm_colors) +
    ggplot2::scale_shape_manual(
      values = c("Initial" = 4, "Reference fit" = 8), name = NULL
    ) +
    ggplot2::labs(
      title = "Optimizer paths on full-loss contours",
      subtitle = paste0("Contours fix ", fixed_values, " at the reference fit."),
      caption = paste(
        "Paths project optimization over all four parameters.",
        "SGD paths connect epoch endpoints."
      ),
      x = axis_names[1], y = axis_names[2], color = "Algorithm"
    ) +
    ggplot2::theme_minimal()

  list(
    gd_summary = gd_summary,
    dn_summary = dn_summary,
    sgd_summary = sgd_summary,
    contour_plot = contour_plot,
    contour_grid = contour_grid,
    path_data = path_data,
    data = data.frame(x = x, y = y),
    dims = dims,
    theta_slice = theta_ref,
    theta_ref = theta_ref,
    reference_loss = fn(theta_ref)
  )
}

# set.seed(123)
# comparison_1000 <- compare_optimizers(
#   c(2, 1, 3, 5),
#   N = 1000,
#   sigma = 0.25,
#   bs_vec = c(1, 50, 200)
# )

# comparison_1000$convergence_plot

# set.seed(123)
# comparison_10000 <- compare_optimizers(
#   c(2, 1, 3, 5),
#   N = 10000,
#   sigma = 0.25,
#   bs_vec = c(1, 50, 200),
# )

# comparison_10000$convergence_plot

# set.seed(123)
# comparison_100000 <- compare_optimizers(
#   c(2, 1, 3, 5),
#   N = 100000,
#   sigma = 0.25,
#   bs_vec = c(1, 50, 200)
# )

# comparison_100000$convergence_plot


# #NEWTON AND GD GET STUCK:

# set.seed(123)
# comparison_10000_p2 <- compare_optimizers(
#   c(1, 4, 0, 2),
#   N = 10000,
#   sigma = 0.25,
#   bs_vec = c(1, 50, 200),
# )

# comparison_1000_p2$convergence_plot

paths_10000 <- compare_optimizers_contour(
  c(2, 1, 3, 5),
  N = 10000,
  sigma = 0.25,
  bs_vec = c(1, 50, 200),
  dims = c(1, 2)
)

paths_10000$contour_plot

paths_10000_p2 <- compare_optimizers_contour(
  c(1, 4, 0, 2),
  N = 10000,
  sigma = 0.25,
  bs_vec = c(1, 50, 200),
  dims = c(1, 2)
)

paths_10000_p2$contour_plot
