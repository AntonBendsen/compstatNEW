gradient_descent <- function(theta_init, fn, gr, tol = 1e-4, max_iter = 600,
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
                gamma_0 = 0.5, K = 100, a = 0.9,
                batch_size = 50, max_epochs = 500,
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

sgd_update_traced <- function(theta_init, x, y, grad,
                gamma_0 = 0.5, K = 100, a = 0.9,
                batch_size = 50, max_epochs = 500,
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
      if (!is.null(cb)) cb()
    }
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

set.seed(123)

N <- 10000
sigma <- 0.25

log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))

y <- dose_response(x, theta_true) +
  rnorm(N, mean = 0, sd = sigma)

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


make_trace_loss_error <- function(x, y, theta_ref, block_size = 10000L) {
  # All observations at a dose share the same prediction. Within-dose
  # variation cancels in L(theta) - L(theta_ref), leaving just group means.
  doses <- sort(unique(x))
  dose_group <- match(x, doses)
  weights <- tabulate(dose_group, nbins = length(doses)) / (2 * length(y))
  mean_y <- vapply(seq_along(doses), function(i) mean(y[dose_group == i]), numeric(1))
  log_doses <- log(doses)
  reference_prediction <- dose_response(doses, theta_ref)
  reference_residual <- mean_y - reference_prediction
  theta_columns <- paste0("theta.", c("alpha", "beta", "gamma", "delta"))

  function(summary_data) {
    loss_error <- numeric(nrow(summary_data))
    if (nrow(summary_data) > 0L) {
      for (start in seq.int(1L, nrow(summary_data), by = block_size)) {
        rows <- seq.int(start, min(start + block_size - 1L, nrow(summary_data)))
        theta <- as.matrix(summary_data[rows, theta_columns, drop = FALSE])
        prediction <- theta[, 3] + (theta[, 4] - theta[, 3]) *
          plogis(theta[, 1] - outer(theta[, 2], log_doses))
        change <- sweep(prediction, 2, reference_prediction, "-")
        # (res_ref - change)^2 - res_ref^2 avoids subtracting two full losses.
        difference <- change * sweep(change, 2, 2 * reference_residual, "-")
        loss_error[rows] <- abs(drop(difference %*% weights))
      }
    }
    summary_data$loss_error <- loss_error
    summary_data
  }
}

converged_trace_endpoints <- function(plot_data, fits, gradient, tolerances) {
  converged <- vapply(seq_along(fits), function(i) {
    norm <- sqrt(sum(gradient(fits[[i]]$theta)^2))
    is.finite(norm) && norm <= tolerances[i]
  }, logical(1))
  # Select the last point visible on the log axes for converged algorithms.
  plot_data %>%
    filter(as.character(algorithm) %in% names(fits)[converged],
           is.finite(.time), .time > 0, is.finite(loss_error), loss_error > 0) %>%
    group_by(algorithm) %>%
    slice_tail(n = 1) %>%
    ungroup()
}

validate_comparison_epochs <- function(mepoch_vec, variable_vec) {
  if (is.null(mepoch_vec)) return(invisible(NULL))
  if (!is.numeric(mepoch_vec) || length(mepoch_vec) != length(variable_vec) ||
      any(!is.finite(mepoch_vec)) || any(mepoch_vec <= 0) ||
      any(mepoch_vec != floor(mepoch_vec))) {
    stop("mepoch_vec must contain one positive whole-number epoch limit per variable_vec entry.")
  }
  if (anyDuplicated(variable_vec)) {
    stop("variable_vec must have unique values when mepoch_vec is supplied.")
  }
  invisible(NULL)
}

make_optimizer_comparison_plots <- function(plot_data, converged_endpoints,
    curve_data, observations, algorithm_colors, algorithm_levels, curve_levels,
    best_so_far) {
  # Keep optimizer callbacks out of the plotting environment.
  error_column <- if (best_so_far) "best_loss_error" else "loss_error"
  x <- observations$x
  y <- observations$y
  convergence_plot <- ggplot(
    plot_data,
    aes(x = .time, y = .data[[error_column]], color = algorithm, group = algorithm)
  ) +
    geom_step(
      data = plot_data[!as.character(plot_data$algorithm) %in%
                         c("Gradient descent", "Damped Newton"), ],
      direction = "hv", linewidth = 1
    ) +
    geom_step(
      data = plot_data[as.character(plot_data$algorithm) %in%
                        c("Gradient descent", "Damped Newton"), ],
      direction = "hv", linewidth = 1
    ) +
    geom_point(
      data = converged_endpoints, aes(x = .time, y = .data[[error_column]]),
      inherit.aes = FALSE, color = "black", shape = 4,
      size = 3, stroke = 1.1, show.legend = FALSE
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = algorithm_colors) +
    labs(
      x = "Time (s)",
      y = if (best_so_far) "Best recorded loss error" else "|L(theta) - L(theta_ref)|",
      color = "Algorithm"
    ) +
    theme_minimal()
  curve_plot <- ggplot(
    curve_data, aes(x = x, y = response, color = algorithm, linetype = algorithm)
  ) +
    geom_point(
      data = data.frame(x = x, y = y), aes(x = x, y = y),
      inherit.aes = FALSE, color = "grey40", alpha = 0.25, size = 1.5
    ) +
    geom_line(linewidth = 1.5) +
    scale_x_log10(breaks = scales::breaks_log(n = 5),
                  labels = scales::label_number()) +
    scale_color_manual(values = c(algorithm_colors, "Data generating curve" = "black")) +
    scale_linetype_manual(values = setNames(
      c(rep("solid", length(algorithm_levels)), "longdash"), curve_levels
    )) +
    labs(x = "Dose (log scale)", y = "Response", color = "Curve",
         linetype = "Curve") +
    theme_minimal(base_size = 13) +
    theme(axis.title = element_text(color = "black"))
  list(convergence_plot = convergence_plot, curve_plot = curve_plot)
}


compare_optimizers <- function(theta_init, N, sigma, variable_vec = 200,
                               variable = "batch_size", mepoch_vec = NULL,
                               best_so_far = FALSE, ...) {
  variable <- match.arg(variable, c("batch_size", "gamma_0"))
  if (!is.numeric(variable_vec) || length(variable_vec) == 0L ||
      any(!is.finite(variable_vec)) || any(variable_vec <= 0)) {
    stop("variable_vec must be a non-empty vector of positive finite numbers.")
  }
  if (variable == "batch_size" && any(variable_vec != floor(variable_vec))) {
    stop("variable_vec must contain whole numbers when varying batch_size.")
  }
  validate_comparison_epochs(mepoch_vec, variable_vec)
  variable_vec <- unique(variable_vec)
  variable_label <- if (variable == "batch_size") "batch size" else "gamma_0"
  sgd_labels <- paste0("SGD (", variable_label, " = ", variable_vec, ")")

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
  gd_fit <- gradient_descent(theta_init, fn, gr, cb = gd_tracer$tracer)
  dn_tracer <- tracer(c("theta"), Delta = 0)
  dn_fit <- damped_newton(theta_init, fn, gr, hess, cb = dn_tracer$tracer)

  add_loss_error <- make_trace_loss_error(x, y, theta_ref)

  gd_summary <- add_loss_error(summary(gd_tracer))
  dn_summary <- add_loss_error(summary(dn_tracer))
  sgd_results <- lapply(seq_along(variable_vec), function(i) {
    value <- variable_vec[i]
    sgd_tracer <- CSwR::tracer("theta", Delta = 0)
    sgd_args <- list(
      theta_init = theta_init, x = x, y = y, grad = gradient_func,
      cb = sgd_tracer$tracer, ...
    )
    sgd_args[[variable]] <- value
    # A per-run limit overrides max_epochs supplied through ... .
    if (!is.null(mepoch_vec)) sgd_args$max_epochs <- mepoch_vec[i]
    fit <- do.call(sgd_update_traced, sgd_args)
    epoch_limit <- if (is.null(sgd_args$max_epochs)) {
      formals(sgd_update_traced)$max_epochs
    } else sgd_args$max_epochs
    batch_size <- if (is.null(sgd_args$batch_size)) {
      formals(sgd_update_traced)$batch_size
    } else sgd_args$batch_size
    updates_per_epoch <- ceiling(N / batch_size)
    result <- summary(sgd_tracer)
    expected_points <- epoch_limit * updates_per_epoch + 1L
    if (nrow(result) != expected_points) {
      stop("SGD recorded ", nrow(result), " points; expected ", expected_points,
           ". The callback must run after every update inside the batch loop.")
    }
    # Every callback point is evaluated and plotted, without thinning.
    result <- add_loss_error(result)
    result$.update <- seq_len(nrow(result)) - 1L
    result$.epoch <- result$.update / updates_per_epoch
    result[[variable]] <- value
    result$max_epochs <- epoch_limit
    result$algorithm <- paste0("SGD (", variable_label, " = ", value, ")")
    result$epochs_completed <- epoch_limit
    result$updates <- nrow(result) - 1L
    list(theta = fit$theta, summary = result)
  })
  sgd_summary <- bind_rows(lapply(sgd_results, function(result) result$summary))

  plot_data <- bind_rows(
    transform(gd_summary, algorithm = "Gradient descent"),
    transform(dn_summary, algorithm = "Damped Newton"),
    sgd_summary
  )
  algorithm_levels <- c("Gradient descent", "Damped Newton", sgd_labels)
  plot_data$algorithm <- factor(plot_data$algorithm, levels = algorithm_levels)
  # Keep the observed error and theta together for selecting fitted curves.
  # The envelope is a separate column, so a repeated minimum never borrows
  # parameters from a later, worse checkpoint.
  plot_data <- plot_data %>%
    arrange(algorithm, .time) %>%
    group_by(algorithm) %>%
    mutate(best_loss_error = cummin(ifelse(is.finite(loss_error), loss_error, Inf))) %>%
    ungroup()
  error_column <- if (best_so_far) "best_loss_error" else "loss_error"
  algorithm_colors <- c(
    "Gradient descent" = "#e0cc11",
    "Damped Newton" = "firebrick",
    setNames(hcl.colors(length(variable_vec), palette = "Dark 3"), sgd_labels)
  )

  converged_endpoints <- converged_trace_endpoints(
    plot_data,
    list("Gradient descent" = gd_fit, "Damped Newton" = dn_fit), gr,
    c(formals(gradient_descent)$tol, formals(damped_newton)$tol)
  )


  final_thetas <- c(
    list("Gradient descent" = gd_fit$theta, "Damped Newton" = dn_fit$theta),
    setNames(lapply(sgd_results, function(result) result$theta), sgd_labels),
    list("Data generating curve" = theta_true)
  )
  x_grid <- exp(seq(min(log_grid), max(log_grid), length.out = 300))
  curve_data <- bind_rows(lapply(names(final_thetas), function(algorithm) {
    data.frame(
      x = x_grid,
      response = dose_response(x_grid, final_thetas[[algorithm]]),
      algorithm = algorithm
    )
  }))
  curve_levels <- c(algorithm_levels, "Data generating curve")
  curve_data$algorithm <- factor(curve_data$algorithm, levels = curve_levels)
  plots <- make_optimizer_comparison_plots(
    plot_data, converged_endpoints, curve_data, data.frame(x = x, y = y),
    algorithm_colors, algorithm_levels, curve_levels, best_so_far
  )

  list(
    gd_summary = gd_summary,
    dn_summary = dn_summary,
    sgd_summary = sgd_summary,
    sgd_settings = unique(sgd_summary[c("algorithm", variable, "max_epochs",
                                        "epochs_completed", "updates")]),
    trace_mode = "full",
    best_so_far = best_so_far,
    convergence_plot = plots$convergence_plot,
    converged_endpoints = converged_endpoints,
    curve_plot = plots$curve_plot,
    data = data.frame(x = x, y = y),
    theta_ref = theta_ref,
    reference_loss = reference_loss
  )
}

compare_optimizers_ut <- function(theta_init, N, sigma, variable_vec = 200,
                               variable = "batch_size", mepoch_vec = NULL,
                               best_so_far = FALSE, ...) {

  compare_optimizers(
    theta_init, N, sigma, variable_vec, variable, mepoch_vec,
    best_so_far = best_so_far, ...
  )
}

best_curves_before_time <- function(comparison, time_limit = 1e-3) {
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      !is.finite(time_limit) || time_limit <= 0) {
    stop("time_limit must be one positive finite number.")
  }

  # The convergence data contains every update and its already computed error.
  trace_data <- comparison$convergence_plot$data
  theta_columns <- paste0("theta.", c("alpha", "beta", "gamma", "delta"))
  selected_traces <- trace_data %>%
    filter(is.finite(.time), .time < time_limit, is.finite(loss_error)) %>%
    group_by(algorithm) %>%
    arrange(loss_error, .time, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(algorithm, .time, loss_error, all_of(theta_columns))

  missing_algorithms <- setdiff(
    as.character(unique(trace_data$algorithm)),
    as.character(selected_traces$algorithm)
  )
  if (length(missing_algorithms) > 0L) {
    stop("No finite trace before time_limit for: ",
         paste(missing_algorithms, collapse = ", "))
  }

  # Reuse the original grid, generating curve, observations, and styling.
  curve_plot <- comparison$curve_plot
  curve_data <- curve_plot$data
  for (i in seq_len(nrow(selected_traces))) {
    curve_rows <- as.character(curve_data$algorithm) ==
      as.character(selected_traces$algorithm[i])
    theta <- as.numeric(unlist(selected_traces[i, theta_columns],
                               use.names = FALSE))
    curve_data$response[curve_rows] <- dose_response(curve_data$x[curve_rows], theta)
  }
  curve_plot$data <- curve_data

  list(curve_plot = curve_plot, selected_traces = selected_traces)
}

compare_optimizers_contour <- function(theta_init, N, sigma, variable_vec = 200,
                                       dims = c(1, 2), variable = "batch_size", ...) {
  if (!is.numeric(theta_init) || length(theta_init) != 4L ||
      any(!is.finite(theta_init))) {
    stop("theta_init must contain four finite numeric parameters.")
  }
  variable <- match.arg(variable, c("batch_size", "gamma_0"))
  if (!is.numeric(variable_vec) || length(variable_vec) == 0L ||
      any(!is.finite(variable_vec)) || any(variable_vec <= 0)) {
    stop("variable_vec must be a non-empty vector of positive finite numbers.")
  }
  if (variable == "batch_size" && any(variable_vec != floor(variable_vec))) {
    stop("variable_vec must contain whole numbers when varying batch_size.")
  }
  if (!is.numeric(dims) || length(dims) != 2L ||
      any(!is.finite(dims)) || any(dims != floor(dims)) ||
      any(dims < 1 | dims > 4) || anyDuplicated(dims)) {
    stop("dims must contain two distinct parameter indices between 1 and 4.")
  }

  parameter_names <- c("alpha", "beta", "gamma", "delta")
  names(theta_init) <- parameter_names
  variable_vec <- unique(variable_vec)
  axis_names <- parameter_names[dims]
  fixed_dims <- setdiff(seq_len(4), dims)
  variable_label <- if (variable == "batch_size") "batch size" else "gamma_0"
  sgd_labels <- paste0("SGD (", variable_label, " = ", variable_vec, ")")

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

  sgd_summary <- dplyr::bind_rows(lapply(variable_vec, function(value) {
    sgd_tracer <- CSwR::tracer("theta", Delta = 0)
    sgd_args <- list(
      theta_init = theta_init, x = x, y = y, grad = gradient_func,
      cb = sgd_tracer$tracer, ...
    )
    sgd_args[[variable]] <- value
    do.call(sgd, sgd_args)
    result <- summary(sgd_tracer)
    result[[variable]] <- value
    result$algorithm <- paste0("SGD (", variable_label, " = ", value, ")")
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
    setNames(hcl.colors(length(variable_vec), palette = "Dark 3"), sgd_labels)
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
