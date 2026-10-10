

library(tidyverse)
library(dplyr)
library(CSwR)

rnorm(100)


test <- rnorm(100)
test



# 1. Model
dose_response <- function(x, theta) {
  alpha <- theta[1]
  beta  <- theta[2]
  gamma <- theta[3]
  delta <- theta[4]

  p <- plogis(-beta * (log(x) - alpha))

  gamma + (delta - gamma) * p
}


# 2. Simulate data
set.seed(123)

theta_true <- c(
  alpha = 2.5,
  beta  = 1.4,
  gamma = 0.5,
  delta = 4.0
)

N <- 200
sigma <- 0.25

# 15 grid points on the log-dose scale
log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))

y <- dose_response(x, theta_true) +
  rnorm(N, mean = 0, sd = sigma)


# 3. Observation-level loss: 1/2 * residual^2
loss_i <- function(theta, x, y) {
  residual <- y - dose_response(x, theta)
  residual^2
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


# 5. Full empirical loss and gradient (for all observations)
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


### We now implement the optimization algorithms: gradient descent, stochastic gradient descent, and mini-batch stochastic gradient descent.

#First we look at SGD with gamma and delta fixed at their true values. This is a simpler problem, and we can visualize the optimization path in the alpha-beta plane.




# Explore SGD with gamma and delta fixed at their true values.
sgd_two_parameters <- function(x, y, theta_ab, gamma_fixed, delta_fixed,
                               lr0, max_epochs = 200, a = 1,
                               K = 2500, theta_ab_true, seed = 1) {
  set.seed(seed)
  N <- length(y)
  theta <- c(
    alpha = theta_ab[1],
    beta = theta_ab[2],
    gamma = gamma_fixed,
    delta = delta_fixed
  )
  max_steps <- max_epochs * N
  trajectory <- matrix(NA_real_, nrow = max_steps + 1, ncol = 2)
  trajectory[1, ] <- theta[1:2]
  parameter_mse <- numeric(max_epochs)
  objective <- numeric(max_epochs)
  elapsed_time <- numeric(max_epochs)
  start_time <- proc.time()[["elapsed"]]
  step <- 0

  for (epoch in seq_len(max_epochs)) {
    for (j in seq_len(N)) {
      index <- sample.int(N, 1)
      step <- step + 1
      learning_rate <- lr0 / (1 + step^a / K)
      gradient <- grad_i(theta, x[index], y[index])[1:2]
      theta[1:2] <- theta[1:2] - learning_rate * gradient
      trajectory[step + 1, ] <- theta[1:2]
    }

    parameter_mse[epoch] <- mean((theta[1:2] - theta_ab_true)^2)
    objective[epoch] <- full_loss(theta, x, y)
    elapsed_time[epoch] <- proc.time()[["elapsed"]] - start_time
  }

  list(
    theta = theta[1:2],
    trajectory = trajectory[seq_len(step + 1), , drop = FALSE],
    parameter_mse = parameter_mse,
    elapsed_time = elapsed_time,
    objective = objective,
    epochs = max_epochs,
    final_parameter_mse = parameter_mse[max_epochs],
    final_loss = objective[max_epochs]
  )
}

theta_ab_start <- c(alpha = 2.0, beta = 1.0)
learning_rates <- c(0.001, 0.005, 0.01, 0.05)

sgd_experiments <- lapply(learning_rates, function(rate) {
  timing <- system.time(result <- sgd_two_parameters(
    x = x,
    y = y,
    theta_ab = theta_ab_start,
    gamma_fixed = theta_true["gamma"],
    delta_fixed = theta_true["delta"],
    theta_ab_true = theta_true[c("alpha", "beta")],
    lr0 = rate,
    max_epochs = 200,
    a = 1,
    K = 2500,
    seed = 1
  ))
  result$runtime_seconds <- unname(timing[["elapsed"]])
  result
})
names(sgd_experiments) <- paste0("lr0 = ", learning_rates)

sgd_experiment_summary <- bind_rows(lapply(
  seq_along(sgd_experiments),
  function(i) {
    result <- sgd_experiments[[i]]
    data.frame(
      learning_rate = learning_rates[i],
      epochs = result$epochs,
      runtime_seconds = result$runtime_seconds,
      final_loss = result$final_loss,
      final_parameter_mse = result$final_parameter_mse
    )
  }
))
print(sgd_experiment_summary)



# Compare the decay parameters while holding the initial learning rate fixed.
# Each heatmap cell reports the median final parameter MSE across seeds.
a_grid <- c(0.55, 0.7, 0.85, 1.0)
K_grid <- c(100, 500, 2500, 10000, 50000)
heatmap_lr0 <- 0.01
heatmap_epochs <- 200
heatmap_seeds <- 1:3

decay_grid <- expand.grid(
  a = a_grid,
  K = K_grid,
  seed = heatmap_seeds
)

decay_results <- bind_rows(lapply(seq_len(nrow(decay_grid)), function(i) {
  run <- sgd_two_parameters(
    x = x,
    y = y,
    theta_ab = theta_ab_start,
    gamma_fixed = theta_true["gamma"],
    delta_fixed = theta_true["delta"],
    theta_ab_true = theta_true[c("alpha", "beta")],
    lr0 = heatmap_lr0,
    max_epochs = heatmap_epochs,
    a = decay_grid$a[i],
    K = decay_grid$K[i],
    seed = decay_grid$seed[i]
  )

  data.frame(
    a = decay_grid$a[i],
    K = decay_grid$K[i],
    seed = decay_grid$seed[i],
    final_parameter_mse = run$final_parameter_mse,
    runtime_seconds = tail(run$elapsed_time, 1)
  )
}))

decay_heatmap_data <- decay_results %>%
  group_by(a, K) %>%
  summarise(
    median_parameter_mse = median(final_parameter_mse),
    median_runtime_seconds = median(runtime_seconds),
    .groups = "drop"
  )

decay_heatmap <- ggplot(
  decay_heatmap_data,
  aes(x = factor(K, levels = K_grid), y = factor(a, levels = rev(a_grid)),
      fill = median_parameter_mse)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = formatC(median_parameter_mse, format = "e", digits = 1)),
    size = 3
  ) +
  scale_fill_gradient(
    trans = "log10",
    low = "lightblue",
    high = "darkblue"
  ) +
  labs(
    title = "SGD parameter MSE by decay parameters",
    subtitle = paste0(
      "Median over ", length(heatmap_seeds),
      " seeds; lr0 = ", heatmap_lr0,
      ", epochs = ", heatmap_epochs
    ),
    x = "K",
    y = "Decay exponent a",
    fill = "Median parameter MSE"
  ) +
  theme_minimal()

print(decay_heatmap)

runtime_heatmap <- ggplot(
  decay_heatmap_data,
  aes(
    x = factor(K, levels = K_grid),
    y = factor(a, levels = rev(a_grid)),
    fill = median_runtime_seconds
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = formatC(median_runtime_seconds, format = "f", digits = 3)),
    size = 3
  ) +
  scale_fill_gradient(
    low = "lightyellow",
    high = "darkred"
  ) +
  labs(
    title = "SGD runtime by decay parameters",
    subtitle = paste0(
      "Median over ", length(heatmap_seeds),
      " seeds; lr0 = ", heatmap_lr0,
      ", epochs = ", heatmap_epochs
    ),
    x = "K",
    y = "Decay exponent a",
    fill = "Median runtime (s)"
  ) +
  theme_minimal()

print(runtime_heatmap)

# We see that setting a=1 and K=2500 is a good tradeoff between speed and accuracy. We will use these decay parameters in the following experiments.





sgd_parameter_mse_time_data <- bind_rows(lapply(
  seq_along(sgd_experiments),
  function(i) {
    result <- sgd_experiments[[i]]
    data.frame(
      elapsed_time = result$elapsed_time,
      parameter_mse = result$parameter_mse,
      learning_rate = names(sgd_experiments)[i]
    )
  }
)) %>%
  filter(elapsed_time > 0, parameter_mse > 0)

sgd_parameter_mse_time_plot <- ggplot(
  sgd_parameter_mse_time_data,
  aes(x = elapsed_time, y = parameter_mse, color = learning_rate)
) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Parameter MSE versus elapsed time",
    x = "Elapsed time (seconds)",
    y = "MSE to the true alpha and beta",
    color = "Learning rate"
  ) +
  theme_minimal()

print(sgd_parameter_mse_time_plot)

sgd_path_data <- bind_rows(lapply(
  seq_along(sgd_experiments),
  function(i) {
    result <- sgd_experiments[[i]]
    trajectory <- result$trajectory
    keep <- unique(round(seq(
      1,
      nrow(trajectory),
      length.out = min(nrow(trajectory), 1500)
    )))
    data.frame(
      alpha = trajectory[keep, 1],
      beta = trajectory[keep, 2],
      iteration = keep - 1,
      learning_rate = names(sgd_experiments)[i]
    )
  }
))
sgd_path_data$learning_rate <- factor(
  sgd_path_data$learning_rate,
  levels = names(sgd_experiments)
)
sgd_path_endpoints <- sgd_path_data %>%
  group_by(learning_rate) %>%
  slice_max(iteration, n = 1, with_ties = FALSE) %>%
  ungroup()

sgd_alpha_range <- range(sgd_path_data$alpha, theta_true["alpha"])
sgd_beta_range <- range(sgd_path_data$beta, theta_true["beta"])
sgd_alpha_padding <- 0.08 * diff(sgd_alpha_range)
sgd_beta_padding <- 0.08 * diff(sgd_beta_range)
sgd_alpha_range <- sgd_alpha_range +
  c(-sgd_alpha_padding, sgd_alpha_padding)
sgd_beta_range <- sgd_beta_range +
  c(-sgd_beta_padding, sgd_beta_padding)

sgd_contour_grid <- expand.grid(
  alpha = seq(sgd_alpha_range[1], sgd_alpha_range[2], length.out = 80),
  beta = seq(sgd_beta_range[1], sgd_beta_range[2], length.out = 80)
)
sgd_contour_grid$loss <- vapply(
  seq_len(nrow(sgd_contour_grid)),
  function(i) {
    theta_slice <- c(
      alpha = sgd_contour_grid$alpha[i],
      beta = sgd_contour_grid$beta[i],
      gamma = theta_true["gamma"],
      delta = theta_true["delta"]
    )
    full_loss(theta_slice, x, y)
  },
  numeric(1)
)

sgd_two_parameter_path_plot <- ggplot() +
  geom_contour(
    data = sgd_contour_grid,
    aes(x = alpha, y = beta, z = loss),
    color = "grey70"
  ) +
  geom_path(
    data = sgd_path_data,
    aes(
      x = alpha, y = beta,
      color = learning_rate,
      group = learning_rate
    ),
    linewidth = 0.7,
    alpha = 0.85
  ) +
  geom_point(
    data = sgd_path_endpoints,
    aes(x = alpha, y = beta, color = learning_rate, shape = learning_rate),
    size = 2.5
  ) +
  geom_point(
    data = data.frame(
      alpha = theta_ab_start[1],
      beta = theta_ab_start[2]
    ),
    aes(x = alpha, y = beta),
    shape = 4,
    size = 3,
    color = "black"
  ) +
  geom_point(
    data = data.frame(
      alpha = theta_true["alpha"],
      beta = theta_true["beta"]
    ),
    aes(x = alpha, y = beta),
    shape = 8,
    size = 3,
    color = "black"
  ) +
  coord_cartesian(
    xlim = sgd_alpha_range,
    ylim = sgd_beta_range
  ) +
  labs(
    title = "SGD paths with gamma and delta fixed",
    subtitle = "Comparing learning rates; cross = start, star = true parameters",
    x = expression(alpha),
    y = expression(beta),
    color = "Learning rate",
    shape = "Final point"
  ) +
  theme_minimal()

sgd_parameter_mse_data <- bind_rows(lapply(
  seq_along(sgd_experiments),
  function(i) {
    data.frame(
      epoch = seq_along(sgd_experiments[[i]]$parameter_mse),
      parameter_mse = sgd_experiments[[i]]$parameter_mse,
      learning_rate = names(sgd_experiments)[i]
    )
  }
))

sgd_parameter_mse_plot <- ggplot(
  sgd_parameter_mse_data,
  aes(x = epoch, y = parameter_mse, color = learning_rate)
) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(
    title = "MSE to the true parameters during SGD",
    x = "Epoch",
    y = expression(frac((hat(alpha)-alpha[true])^2 +
                          (hat(beta)-beta[true])^2, 2)),
    color = "Learning rate"
  ) +
  theme_minimal()

print(sgd_two_parameter_path_plot)
print(sgd_parameter_mse_plot)


sgd_experiment_summary

#taking initial learning rate of 0.005 seems to result in the best convergence compared to time and accuracy. We will use this learning rate in the following experiments.

# Compare SGD trajectories across five random seeds with fixed hyperparameters.
trajectory_seeds <- 1:5
trajectory_lr0 <- 0.005
trajectory_a <- 1
trajectory_K <- 2500
trajectory_epochs <- 200

seed_trajectories <- lapply(trajectory_seeds, function(seed_value) {
  sgd_two_parameters(
    x = x,
    y = y,
    theta_ab = theta_ab_start,
    gamma_fixed = theta_true["gamma"],
    delta_fixed = theta_true["delta"],
    theta_ab_true = theta_true[c("alpha", "beta")],
    lr0 = trajectory_lr0,
    max_epochs = trajectory_epochs,
    a = trajectory_a,
    K = trajectory_K,
    seed = seed_value
  )
})
names(seed_trajectories) <- paste("Seed", trajectory_seeds)

seed_trajectory_data <- bind_rows(lapply(
  seq_along(seed_trajectories),
  function(i) {
    trajectory <- seed_trajectories[[i]]$trajectory
    keep <- unique(round(seq(
      1,
      nrow(trajectory),
      length.out = min(nrow(trajectory), 1500)
    )))
    data.frame(
      alpha = trajectory[keep, 1],
      beta = trajectory[keep, 2],
      iteration = keep - 1,
      seed = names(seed_trajectories)[i]
    )
  }
))

seed_trajectory_endpoints <- bind_rows(lapply(
  seq_along(seed_trajectories),
  function(i) {
    data.frame(
      alpha = seed_trajectories[[i]]$theta["alpha"],
      beta = seed_trajectories[[i]]$theta["beta"],
      seed = names(seed_trajectories)[i]
    )
  }
))

seed_alpha_range <- range(seed_trajectory_data$alpha, theta_true["alpha"])
seed_beta_range <- range(seed_trajectory_data$beta, theta_true["beta"])
seed_alpha_padding <- 0.08 * diff(seed_alpha_range)
seed_beta_padding <- 0.08 * diff(seed_beta_range)
seed_alpha_range <- seed_alpha_range +
  c(-seed_alpha_padding, seed_alpha_padding)
seed_beta_range <- seed_beta_range +
  c(-seed_beta_padding, seed_beta_padding)

seed_contour_grid <- expand.grid(
  alpha = seq(seed_alpha_range[1], seed_alpha_range[2], length.out = 80),
  beta = seq(seed_beta_range[1], seed_beta_range[2], length.out = 80)
)
seed_contour_grid$loss <- vapply(
  seq_len(nrow(seed_contour_grid)),
  function(i) {
    theta_slice <- c(
      alpha = seed_contour_grid$alpha[i],
      beta = seed_contour_grid$beta[i],
      gamma = theta_true["gamma"],
      delta = theta_true["delta"]
    )
    full_loss(theta_slice, x, y)
  },
  numeric(1)
)

sgd_seed_trajectory_plot <- ggplot() +
  geom_contour(
    data = seed_contour_grid,
    aes(x = alpha, y = beta, z = loss),
    color = "grey70"
  ) +
  geom_path(
    data = seed_trajectory_data,
    aes(x = alpha, y = beta, color = seed, group = seed),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    data = seed_trajectory_endpoints,
    aes(x = alpha, y = beta, color = seed),
    size = 2
  ) +
  geom_point(
    data = data.frame(
      alpha = theta_ab_start["alpha"],
      beta = theta_ab_start["beta"]
    ),
    aes(x = alpha, y = beta),
    shape = 4,
    size = 3,
    color = "black"
  ) +
  geom_point(
    data = data.frame(
      alpha = theta_true["alpha"],
      beta = theta_true["beta"]
    ),
    aes(x = alpha, y = beta),
    shape = 8,
    size = 3,
    color = "black"
  ) +
  coord_cartesian(xlim = seed_alpha_range, ylim = seed_beta_range) +
  labs(
    title = "Two-parameter SGD trajectories across five seeds",
    subtitle = paste0(
      "lr0 = ", trajectory_lr0,
      ", a = ", trajectory_a,
      ", K = ", trajectory_K,
      "; cross = start, star = true parameters"
    ),
    x = expression(alpha),
    y = expression(beta),
    color = "Random seed"
  ) +
  theme_minimal()

print(sgd_seed_trajectory_plot)












































# 6. Gradient Descent
fn <- function(theta) full_loss(theta, x, y)
gr <- function(theta) full_gradient(theta, x, y)


gradient_descent <- function(theta_init, fn, gr, tol = 1e-8, max_iter = 1000,
                             gamma_0 = 1.0, c = 0.1, d = 0.5,
                             cb = NULL) {
  theta <- theta_init
  gamma <- gamma_0
  history <- numeric(max_iter + 1)
  history[1] <- fn(theta)
  trajectory <- matrix(NA_real_, nrow = max_iter + 1, ncol = length(theta))
  trajectory[1, ] <- theta
  time_history <- numeric(max_iter + 1)
  start_time <- proc.time()[["elapsed"]]
  gradient_norm_history <- numeric(0)
  gradient_time_history <- numeric(0)
  iterations <- 0

  for (i in seq_len(max_iter)) {
    grad <- gr(theta)
    grad_norm_sq <- sum(grad^2)
    gradient_norm_history <- c(gradient_norm_history, sqrt(grad_norm_sq))
    gradient_time_history <- c(
      gradient_time_history,
      proc.time()[["elapsed"]] - start_time
    )

    if (!is.null(cb)) {
      cb()
    }

    if (sqrt(grad_norm_sq) <= tol) {
      message("Converged successfully after ", i, " iterations.")
      break
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
    iterations <- i
    history[i + 1] <- fn(theta)
    trajectory[i + 1, ] <- theta
    time_history[i + 1] <- proc.time()[["elapsed"]] - start_time
  }

  if (iterations == max_iter) {
    warning("Maximum number of iterations reached without full convergence.")
  }

  list(
    theta = theta,
    gamma = gamma,
    history = history[seq_len(iterations + 1)],
    trajectory = trajectory[seq_len(iterations + 1), , drop = FALSE],
    time_history = time_history[seq_len(iterations + 1)],
    gradient_norm_history = gradient_norm_history,
    gradient_time_history = gradient_time_history
  )
}


# 7. Stochastic Gradient Descent
sgd <- function(x, y, theta, lr0 = 0.001,
                max_epochs = 200, a = 0.6, K = 100,
                tol = 1e-4, seed = 1, cb = NULL) {
  set.seed(seed)
  N <- length(y)
  history <- numeric(max_epochs)
  step <- 0
  trajectory <- matrix(
    NA_real_,
    nrow = max_epochs * N + 1,
    ncol = length(theta)
  )
  trajectory[1, ] <- theta
  time_history <- numeric(max_epochs * N + 1)
  start_time <- proc.time()[["elapsed"]]
  gradient_norm_history <- numeric(max_epochs * N)
  gradient_time_history <- numeric(max_epochs * N)
  convergence_gradient_history <- numeric(max_epochs)
  convergence_time_history <- numeric(max_epochs)
  epochs_run <- 0
  converged <- FALSE
  full_gradient_norm <- Inf

  for (epoch in seq_len(max_epochs)) {
    # N random draws WITH replacement per epoch
    for (j in seq_len(N)) {
      i <- sample.int(N, 1)
      step <- step + 1

      # Robbins-Monro style decreasing learning rate
      lr <- lr0 / (1 + step / K)^a

      g <- grad_i(theta, x[i], y[i])
      grad_norm_sq <- sum(g^2)
      gradient_norm_history[step] <- sqrt(grad_norm_sq)

      if (!is.null(cb)) {
        cb()
      }
      gradient_time_history[step] <- proc.time()[["elapsed"]] - start_time

      theta <- theta - lr * g

      if (any(!is.finite(theta))) {
        stop("SGD diverged. Try a smaller lr0.")
      }

      trajectory[step + 1, ] <- theta
      time_history[step + 1] <- proc.time()[["elapsed"]] - start_time
    }

    history[epoch] <- full_loss(theta, x, y)
    epochs_run <- epoch
    full_gradient_norm <- sqrt(sum(full_gradient(theta, x, y)^2))
    convergence_gradient_history[epoch] <- full_gradient_norm
    convergence_time_history[epoch] <- proc.time()[["elapsed"]] - start_time

    if (full_gradient_norm < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "SGD did not reach the gradient-norm tolerance before max_epochs."
    )
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history[seq_len(epochs_run)],
    trajectory = trajectory[seq_len(step + 1), , drop = FALSE],
    time_history = time_history[seq_len(step + 1)],
    gradient_norm_history = gradient_norm_history[seq_len(step)],
    gradient_time_history = gradient_time_history[seq_len(step)],
    convergence_gradient_history =
      convergence_gradient_history[seq_len(epochs_run)],
    convergence_time_history = convergence_time_history[seq_len(epochs_run)],
    epochs = epochs_run,
    full_gradient_norm = full_gradient_norm,
    converged = converged
  )
}


# 8. Mini-batch SGD
mini_batch_sgd <- function(x, y, theta, batch_size = 32,
                           lr0 = 0.001, max_epochs = 200,
                           a = 0.6, K = 100, tol = 1e-4, seed = 1,
                           cb = NULL) {
  set.seed(seed)
  N <- length(y)
  history <- numeric(max_epochs)
  step <- 0
  max_steps <- max_epochs * ceiling(N / batch_size)
  trajectory <- matrix(NA_real_, nrow = max_steps + 1, ncol = length(theta))
  trajectory[1, ] <- theta
  time_history <- numeric(max_steps + 1)
  start_time <- proc.time()[["elapsed"]]
  gradient_norm_history <- numeric(max_steps)
  gradient_time_history <- numeric(max_steps)
  convergence_gradient_history <- numeric(max_epochs)
  convergence_time_history <- numeric(max_epochs)
  epochs_run <- 0
  converged <- FALSE
  full_gradient_norm <- Inf

  for (epoch in seq_len(max_epochs)) {
    # Shuffle: every observation used once per epoch
    indices <- sample.int(N)
    batches <- split(
      indices,
      ceiling(seq_along(indices) / batch_size)
    )

    for (batch in batches) {
      step <- step + 1
      lr <- lr0 / (1 + step / K)^a

      g <- Reduce(
        `+`,
        lapply(batch, function(i) grad_i(theta, x[i], y[i])
        )
      ) / length(batch)
      grad_norm_sq <- sum(g^2)
      gradient_norm_history[step] <- sqrt(grad_norm_sq)

      if (!is.null(cb)) {
        cb()
      }
      gradient_time_history[step] <- proc.time()[["elapsed"]] - start_time

      theta <- theta - lr * g

      if (any(!is.finite(theta))) {
        stop("Mini-batch SGD diverged. Try a smaller lr0.")
      }

      trajectory[step + 1, ] <- theta
      time_history[step + 1] <- proc.time()[["elapsed"]] - start_time
    }

    history[epoch] <- full_loss(theta, x, y)
    epochs_run <- epoch
    full_gradient_norm <- sqrt(sum(full_gradient(theta, x, y)^2))
    convergence_gradient_history[epoch] <- full_gradient_norm
    convergence_time_history[epoch] <- proc.time()[["elapsed"]] - start_time

    if (full_gradient_norm < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "Mini-batch SGD did not reach the gradient-norm tolerance before max_epochs."
    )
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history[seq_len(epochs_run)],
    trajectory = trajectory[seq_len(step + 1), , drop = FALSE],
    time_history = time_history[seq_len(step + 1)],
    gradient_norm_history = gradient_norm_history[seq_len(step)],
    gradient_time_history = gradient_time_history[seq_len(step)],
    convergence_gradient_history =
      convergence_gradient_history[seq_len(epochs_run)],
    convergence_time_history = convergence_time_history[seq_len(epochs_run)],
    epochs = epochs_run,
    full_gradient_norm = full_gradient_norm,
    converged = converged
  )
}





theta_start <- c(
  alpha = 2.0,
  beta  = 1.0,
  gamma = 0.0,
  delta = 3.0
)

# Trace the parameter vector and squared gradient norm at each update.
gd_tracer <- tracer(c("theta", "grad_norm_sq"), Delta = 0)
sgd_tracer <- tracer(c("theta", "grad_norm_sq"), Delta = 0)
mb_tracer <- tracer(c("theta", "grad_norm_sq"), Delta = 0)

# Gradient Descent
result_gd <- gradient_descent(
  theta_init = theta_start,
  fn = fn,
  gr = gr,
  cb = gd_tracer$tracer
)

# SGD
result_sgd <- sgd(
  x, y, theta_start,
  lr0 = 0.001,
  max_epochs = 2000,
  a = 0.6,
  K = 10000,
  seed = 1,
  cb = sgd_tracer$tracer
)

# Mini-batch SGD
result_mb <- mini_batch_sgd(
  x, y, theta_start,
  batch_size = 32,
  lr0 = 0.001,
  max_epochs = 2000,
  a = 0.6,
  K = 10000,
  seed = 1,
  cb = mb_tracer$tracer
)

#summary(gd_tracer)
#summary(sgd_tracer)
#summary(mb_tracer)

# Sammenlign parameterestimater
results <- rbind(
  True = theta_true,
  GD = result_gd$theta,
  SGD = result_sgd$theta,
  MiniBatch = result_mb$theta
)

print(round(results, 4))
print(data.frame(
  method = c("SGD", "Mini-batch SGD"),
  epochs = c(result_sgd$epochs, result_mb$epochs),
  full_gradient_norm = c(
    result_sgd$full_gradient_norm,
    result_mb$full_gradient_norm
  ),
  converged = c(result_sgd$converged, result_mb$converged)
))



### --- Plots ---















# Compare Euclidean parameter distance to the true parameters.
distance_history <- function(result, method) {
  data.frame(
    iteration = seq_len(nrow(result$trajectory) - 1),
    time = result$time_history[-1],
    distance = sqrt(rowSums(
      (result$trajectory[-1, , drop = FALSE] -
         matrix(theta_true, nrow = nrow(result$trajectory) - 1,
                ncol = length(theta_true), byrow = TRUE))^2
    )),
    method = method
  )
}

distance_data <- rbind(
  distance_history(result_gd, "GD"),
  distance_history(result_sgd, "SGD"),
  distance_history(result_mb, "Mini-batch SGD")
)
distance_data <- distance_data[
  distance_data$distance > 0,
]

distance_by_iteration_plot <- ggplot(
  distance_data,
  aes(x = iteration, y = distance, color = method)
) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Parameter distance to truth by iteration",
    x = "Iteration (parameter update)",
    y = "Euclidean distance to true parameters",
    color = "Method"
  ) +
  theme_minimal()

gradient_norm_data <- bind_rows(
  data.frame(
    time = result_gd$gradient_time_history,
    gradient_norm = result_gd$gradient_norm_history,
    method = "GD"
  ),
  data.frame(
    time = result_sgd$convergence_time_history,
    gradient_norm = result_sgd$convergence_gradient_history,
    method = "SGD"
  ),
  data.frame(
    time = result_mb$convergence_time_history,
    gradient_norm = result_mb$convergence_gradient_history,
    method = "Mini-batch SGD"
  )
) %>%
  filter(time > 0, gradient_norm > 0)

gradient_norm_by_time_plot <- ggplot(
  gradient_norm_data,
  aes(x = time, y = gradient_norm, color = method)
) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Full-data gradient norm over elapsed time",
    x = "Elapsed time (seconds)",
    y = "Norm of full-data gradient",
    color = "Method"
  ) +
  theme_minimal()

print(distance_by_iteration_plot)
print(gradient_norm_by_time_plot)






# Plot the optimization paths projected onto the alpha-beta plane.
make_path_data <- function(result, method, max_points = 1e1000) {
  trajectory <- result$trajectory
  keep <- unique(round(seq(
    from = 1,
    to = nrow(trajectory),
    length.out = min(nrow(trajectory), max_points)
  )))

  data.frame(
    alpha = trajectory[keep, 1],
    beta = trajectory[keep, 2],
    iteration = keep - 1,
    method = method
  )
}

trajectory_plot_data <- bind_rows(
  make_path_data(result_gd, "GD"),
  make_path_data(result_sgd, "SGD"),
  make_path_data(result_mb, "Mini-batch SGD")
)

alpha_range <- range(trajectory_plot_data$alpha, theta_true["alpha"])
beta_range <- range(trajectory_plot_data$beta, theta_true["beta"])
alpha_padding <- 0.08 * diff(alpha_range)
beta_padding <- 0.08 * diff(beta_range)
alpha_range <- alpha_range + c(-alpha_padding, alpha_padding)
beta_range <- beta_range + c(-beta_padding, beta_padding)

contour_grid <- expand.grid(
  alpha = seq(alpha_range[1], alpha_range[2], length.out = 80),
  beta = seq(beta_range[1], beta_range[2], length.out = 80)
)
contour_grid$loss <- vapply(
  seq_len(nrow(contour_grid)),
  function(i) {
    theta_slice <- c(
      alpha = contour_grid$alpha[i],
      beta = contour_grid$beta[i],
      gamma = theta_true["gamma"],
      delta = theta_true["delta"]
    )
    mean((y - dose_response(x, theta_slice))^2)
  },
  numeric(1)
)

path_plot <- ggplot() +
  geom_contour(
    data = contour_grid,
    aes(x = alpha, y = beta, z = loss),
    color = "grey70"
  ) +
  geom_path(
    data = trajectory_plot_data,
    aes(x = alpha, y = beta, color = method, group = method),
    linewidth = 0.7,
    alpha = 0.85
  ) +
  geom_point(
    data = data.frame(alpha = theta_start[1], beta = theta_start[2]),
    aes(x = alpha, y = beta),
    shape = 4,
    size = 3,
    stroke = 1.2,
    color = "black"
  ) +
  geom_point(
    data = data.frame(alpha = theta_true[1], beta = theta_true[2]),
    aes(x = alpha, y = beta),
    shape = 8,
    size = 3,
    stroke = 1.2,
    color = "black"
  ) +
  annotate(
    "text",
    x = theta_start[1],
    y = theta_start[2],
    label = "Start",
    hjust = -0.15,
    vjust = 1.2
  ) +
  annotate(
    "text",
    x = theta_true[1],
    y = theta_true[2],
    label = "True",
    hjust = -0.15,
    vjust = -0.5
  ) +
  coord_cartesian(xlim = alpha_range, ylim = beta_range) +
  labs(
    title = "Optimization paths toward the true parameters",
    subtitle = "Paths projected onto alpha-beta; contours fix gamma and delta at their true values",
    x = expression(alpha),
    y = expression(beta),
    color = "Method"
  ) +
  theme_minimal()

print(path_plot)



