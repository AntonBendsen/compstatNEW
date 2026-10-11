## This document is the "storyline" where things are presented in order. 

library(tidyverse)
library(dplyr)
library(CSwR)

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

N <- 5000
sigma <- 0.25

# 15 grid points on the log-dose scale
log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))

y <- dose_response(x, theta_true) +
  rnorm(N, mean = 0, sd = sigma)


#data generating curve:

dose_grid <- exp(seq(min(log(x)), max(log(x)), length.out = 300))
data_generating_curve <- data.frame(
  dose = dose_grid,
  response = dose_response(dose_grid, theta_true)
)

simulated_data_plot <- ggplot(data.frame(dose = x, response = y),
                              aes(x = dose, y = response)) +
  geom_point(alpha = 0.3, color = "grey35") +
  geom_line(
    data = data_generating_curve,
    aes(x = dose, y = response),
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_x_log10() +
  labs(
    title = "Simulated dose-response data",
    subtitle = "Dashed line is the data-generating curve",
    x = "Dose (log scale)",
    y = "Response"
  ) +
  theme_minimal()

print(simulated_data_plot)



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

hess_i <- function(theta, x, y) {
  alpha <- theta[1]
  beta <- theta[2]
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
  hess_fit <- matrix(0, nrow = 4, ncol = 4)
  hess_fit[1, 1] <- (delta - gamma) * q_prime * beta^2
  hess_fit[1, 2] <- hess_fit[2, 1] <-
    (delta - gamma) * (q - beta * u * q_prime)
  hess_fit[2, 2] <- (delta - gamma) * u^2 * q_prime
  hess_fit[1, 3] <- hess_fit[3, 1] <- -beta * q
  hess_fit[1, 4] <- hess_fit[4, 1] <- beta * q
  hess_fit[2, 3] <- hess_fit[3, 2] <- u * q
  hess_fit[2, 4] <- hess_fit[4, 2] <- -u * q

  tcrossprod(df) - residual * hess_fit
}

full_hessian <- function(theta, x, y) {
  hessian <- matrix(0, nrow = 4, ncol = 4)

  for (i in seq_along(y)) {
    hessian <- hessian + hess_i(theta, x[i], y[i])
  }

  hessian / length(y)
}

#
#We use newton to find the "sample theta values"


damped_newton <- function(theta_init, fn, gr, hess, tol = 1e-3, max_iter = 100, 
                          c = 0.1, d = 0.5, cb = NULL) {
  theta <- theta_init

  for (i in seq_len(max_iter)) {

    grad <- gr(theta)
    grad_norm <- sqrt(sum(grad^2))

    if (!is.null(cb)) {
      cb() # Record the current state for later inspection.
    }

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

fn <- function(theta) full_loss(theta, x, y)
gr <- function(theta) full_gradient(theta, x, y)
hess <- function(theta) full_hessian(theta, x, y)

expand_sample_theta <- function(theta_ab) {
  c(
    alpha = unname(theta_ab[1]),
    beta = unname(theta_ab[2]),
    gamma = unname(theta_true["gamma"]),
    delta = unname(theta_true["delta"])
  )
}

fn_sample_ab <- function(theta_ab) {
  full_loss(expand_sample_theta(theta_ab), x, y)
}
gr_sample_ab <- function(theta_ab) {
  full_gradient(expand_sample_theta(theta_ab), x, y)[1:2]
}
hess_sample_ab <- function(theta_ab) {
  full_hessian(expand_sample_theta(theta_ab), x, y)[1:2, 1:2]
}

set.seed(123)
newton_sample_fit <- damped_newton(
  theta_init = c(alpha = 2.0, beta = 1.0),
  fn = fn_sample_ab,
  gr = gr_sample_ab,
  hess = hess_sample_ab,
  tol = 1e-12
)
theta_ab_sample <- newton_sample_fit$theta
theta_sample <- expand_sample_theta(theta_ab_sample)
damped_newton_alpha <- theta_sample[1]
damped_newton_beta <- theta_sample[2]
damped_newton_gamma <- theta_sample[3]
damped_newton_delta <- theta_sample[4]


#

### We now implement the optimization algorithms: gradient descent, stochastic gradient descent, and mini-batch stochastic gradient descent.

#First we look at SGD with gamma and delta fixed at their true values. This is a simpler problem, and we can visualize the optimization path in the alpha-beta plane.


# Explore SGD with gamma and delta fixed at their true values.
sgd_two_parameters <- function(x, y, theta_ab, gamma_fixed, delta_fixed,
                               lr0, max_epochs = 200, a = 1,
                               K = 2500, theta_ab_target, seed = 1) {
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

    parameter_mse[epoch] <- mean((theta[1:2] - theta_ab_target)^2)
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
    theta_ab_target = theta_ab_sample,
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


### ---- Heatmap ----

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
    theta_ab_target = theta_ab_sample,
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


### --- Learning rates ---

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
    y = "MSE to Newton sample alpha and beta",
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
      alpha = theta_ab_sample["alpha"],
      beta = theta_ab_sample["beta"]
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
    subtitle = "Comparing learning rates; cross = start, star = Newton sample estimate",
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
    title = "MSE to Newton sample parameters during SGD",
    x = "Epoch",
    y = expression(frac((hat(alpha)-alpha[sample])^2 +
                (hat(beta)-beta[sample])^2, 2)),
    color = "Learning rate"
  ) +
  theme_minimal()

print(sgd_two_parameter_path_plot)
print(sgd_parameter_mse_plot)



### --- Different seeds ----


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
    theta_ab_target = theta_ab_sample,
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
      alpha = theta_ab_sample["alpha"],
      beta = theta_ab_sample["beta"]
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
      "; cross = start, star = Newton sample estimate"
    ),
    x = expression(alpha),
    y = expression(beta),
    color = "Random seed"
  ) +
  theme_minimal()

print(sgd_seed_trajectory_plot)


### ---- Differnet starting points ---


starting_values <- data.frame(
  alpha = c(1.5, 2.0, 2.0, 3.0, 3.5,0),
  beta = c(0.5, 1.0, 2.0, 1.0, 2.0,0)
)
starting_values$start <- paste0(
  "(", starting_values$alpha, ", ", starting_values$beta, ")"
)

start_value_trajectories <- lapply(seq_len(nrow(starting_values)), function(i) {
  sgd_two_parameters(
    x = x,
    y = y,
    theta_ab = c(
      alpha = starting_values$alpha[i],
      beta = starting_values$beta[i]
    ),
    gamma_fixed = theta_true["gamma"],
    delta_fixed = theta_true["delta"],
    theta_ab_target = theta_ab_sample,
    lr0 = trajectory_lr0,
    max_epochs = trajectory_epochs,
    a = trajectory_a,
    K = trajectory_K,
    seed = trajectory_seeds[1]
  )
})
names(start_value_trajectories) <- starting_values$start

start_value_trajectory_data <- bind_rows(lapply(
  seq_along(start_value_trajectories),
  function(i) {
    trajectory <- start_value_trajectories[[i]]$trajectory
    keep <- unique(round(seq(
      1,
      nrow(trajectory),
      length.out = min(nrow(trajectory), 1500)
    )))
    data.frame(
      alpha = trajectory[keep, 1],
      beta = trajectory[keep, 2],
      iteration = keep - 1,
      start = names(start_value_trajectories)[i]
    )
  }
))

start_value_trajectory_endpoints <- bind_rows(lapply(
  seq_along(start_value_trajectories),
  function(i) {
    final_theta <- start_value_trajectories[[i]]$theta
    data.frame(
      alpha = unname(final_theta["alpha"]),
      beta = unname(final_theta["beta"]),
      start = names(start_value_trajectories)[i]
    )
  }
))

start_value_alpha_range <- range(
  start_value_trajectory_data$alpha,
  theta_true["alpha"]
)
start_value_beta_range <- range(
  start_value_trajectory_data$beta,
  theta_true["beta"]
)
start_value_alpha_padding <- 0.08 * diff(start_value_alpha_range)
start_value_beta_padding <- 0.08 * diff(start_value_beta_range)
start_value_alpha_range <- start_value_alpha_range +
  c(-start_value_alpha_padding, start_value_alpha_padding)
start_value_beta_range <- start_value_beta_range +
  c(-start_value_beta_padding, start_value_beta_padding)

start_value_contour_grid <- expand.grid(
  alpha = seq(
    start_value_alpha_range[1],
    start_value_alpha_range[2],
    length.out = 80
  ),
  beta = seq(
    start_value_beta_range[1],
    start_value_beta_range[2],
    length.out = 80
  )
)
start_value_contour_grid$loss <- vapply(
  seq_len(nrow(start_value_contour_grid)),
  function(i) {
    theta_slice <- c(
      alpha = start_value_contour_grid$alpha[i],
      beta = start_value_contour_grid$beta[i],
      gamma = theta_true["gamma"],
      delta = theta_true["delta"]
    )
    full_loss(theta_slice, x, y)
  },
  numeric(1)
)

sgd_start_value_trajectory_plot <- ggplot() +
  geom_contour(
    data = start_value_contour_grid,
    aes(x = alpha, y = beta, z = loss),
    color = "grey70"
  ) +
  geom_path(
    data = start_value_trajectory_data,
    aes(x = alpha, y = beta, color = start, group = start),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    data = starting_values,
    aes(x = alpha, y = beta, color = start),
    shape = 4,
    size = 3
  ) +
  geom_point(
    data = start_value_trajectory_endpoints,
    aes(x = alpha, y = beta, color = start),
    size = 2
  ) +
  geom_point(
    data = data.frame(
      alpha = unname(theta_ab_sample["alpha"]),
      beta = unname(theta_ab_sample["beta"])
    ),
    aes(x = alpha, y = beta),
    shape = 8,
    size = 3,
    color = "black"
  ) +
  coord_cartesian(
    xlim = start_value_alpha_range,
    ylim = start_value_beta_range
  ) +
  labs(
    title = "Two-parameter SGD trajectories from five starting values",
    subtitle = paste0(
      "Same seed and hyperparameters: lr0 = ", trajectory_lr0,
      ", a = ", trajectory_a,
      ", K = ", trajectory_K,
      "; crosses = starts, circles = final points, star = Newton sample estimate"
    ),
    x = expression(alpha),
    y = expression(beta),
    color = "Starting (alpha, beta)"
  ) +
  theme_minimal()

print(sgd_start_value_trajectory_plot)




#### --- Mini batches ---



# Compare ordinary SGD with mini-batch SGD for the two-parameter model.
comparison_lr0 <- 0.005
comparison_a <- 1
comparison_K <- 2500
comparison_epochs <- 200
mini_batch_sizes <- c(10, 50, 250, 1000, 2500)

gradient_ab <- function(theta_ab, x_batch, y_batch) {
  alpha <- theta_ab[1]
  beta <- theta_ab[2]
  log_x <- log(x_batch)
  p <- plogis(-beta * (log_x - alpha))
  q <- p * (1 - p)
  fit <- theta_true["gamma"] +
    (theta_true["delta"] - theta_true["gamma"]) * p
  residual <- y_batch - fit
  d_alpha <- (theta_true["delta"] - theta_true["gamma"]) * beta * q
  d_beta <- -(theta_true["delta"] - theta_true["gamma"]) *
    (log_x - alpha) * q

  c(
    alpha = mean(-residual * d_alpha),
    beta = mean(-residual * d_beta)
  )
}

loss_ab <- function(theta_ab) {
  theta <- c(
    alpha = theta_ab[1],
    beta = theta_ab[2],
    gamma = theta_true["gamma"],
    delta = theta_true["delta"]
  )
  mean(0.5 * (y - dose_response(x, theta))^2)
}

mini_batch_sgd_two_parameters <- function(x, y, theta_ab, batch_size,
                                          lr0, max_epochs, a, K, seed = 1) {
  set.seed(seed)
  n <- length(y)
  batches_per_epoch <- ceiling(n / batch_size)
  max_updates <- max_epochs * batches_per_epoch
  theta <- theta_ab
  trajectory <- matrix(NA_real_, nrow = max_updates + 1, ncol = 2)
  colnames(trajectory) <- c("alpha", "beta")
  trajectory[1, ] <- theta
  elapsed_time <- numeric(max_epochs)
  gradient_norm <- numeric(max_epochs)
  objective <- numeric(max_epochs)
  start_time <- proc.time()[["elapsed"]]
  update <- 0

  for (epoch in seq_len(max_epochs)) {
    indices <- sample.int(n)

    for (first in seq.int(1, n, by = batch_size)) {
      batch_indices <- indices[first:min(first + batch_size - 1L, n)]
      update <- update + 1
      learning_rate <- lr0 / (1 + update^a / K)
      batch_gradient <- gradient_ab(
        theta,
        x[batch_indices],
        y[batch_indices]
      )
      theta <- theta - learning_rate * batch_gradient
      trajectory[update + 1, ] <- theta
    }

    gradient_norm[epoch] <- sqrt(sum(gradient_ab(theta, x, y)^2))
    objective[epoch] <- loss_ab(theta)
    elapsed_time[epoch] <- proc.time()[["elapsed"]] - start_time
  }

  list(
    theta = theta,
    trajectory = trajectory[seq_len(update + 1), , drop = FALSE],
    elapsed_time = elapsed_time,
    gradient_norm = gradient_norm,
    objective = objective,
    epochs = max_epochs,
    final_loss = objective[max_epochs],
    final_gradient_norm = gradient_norm[max_epochs]
  )
}

mini_batch_results <- lapply(mini_batch_sizes, function(batch_size) {
  mini_batch_sgd_two_parameters(
    x = x,
    y = y,
    theta_ab = theta_ab_start,
    batch_size = batch_size,
    lr0 = comparison_lr0,
    max_epochs = comparison_epochs,
    a = comparison_a,
    K = comparison_K,
    seed = 1
  )
})
names(mini_batch_results) <- paste("Batch", mini_batch_sizes)

sgd_reference <- sgd_experiments[[which(learning_rates == comparison_lr0)]]
sgd_expected_rows <- comparison_epochs * length(y) + 1L
sgd_history_is_current <-
  nrow(sgd_reference$trajectory) == sgd_expected_rows &&
  length(sgd_reference$elapsed_time) == comparison_epochs &&
  length(sgd_reference$objective) == comparison_epochs

if (!sgd_history_is_current) {
  sgd_timing <- system.time(sgd_reference <- sgd_two_parameters(
    x = x,
    y = y,
    theta_ab = theta_ab_start,
    gamma_fixed = theta_true["gamma"],
    delta_fixed = theta_true["delta"],
    theta_ab_target = theta_ab_sample,
    lr0 = comparison_lr0,
    max_epochs = comparison_epochs,
    a = comparison_a,
    K = comparison_K,
    seed = 1
  ))
  sgd_reference$runtime_seconds <- unname(sgd_timing[["elapsed"]])
}

sgd_epoch_rows <- seq_len(comparison_epochs) * length(y) + 1L
sgd_gradient_norm <- vapply(sgd_epoch_rows, function(row) {
  sqrt(sum(gradient_ab(sgd_reference$trajectory[row, ], x, y)^2))
}, numeric(1))

mini_batch_summary <- bind_rows(lapply(seq_along(mini_batch_results), function(i) {
  result <- mini_batch_results[[i]]
  data.frame(
    method = names(mini_batch_results)[i],
    batch_size = mini_batch_sizes[i],
    epochs = result$epochs,
    elapsed_seconds = tail(result$elapsed_time, 1),
    final_loss = result$final_loss,
    final_gradient_norm = result$final_gradient_norm
  )
}))
sgd_summary <- data.frame(
  method = "Ordinary SGD",
  batch_size = 1L,
  epochs = sgd_reference$epochs,
  elapsed_seconds = sgd_reference$runtime_seconds,
  final_loss = sgd_reference$final_loss,
  final_gradient_norm = tail(sgd_gradient_norm, 1)
)
print(bind_rows(sgd_summary, mini_batch_summary))

path_keep <- function(trajectory, max_points = 1500) {
  unique(round(seq(
    from = 1,
    to = nrow(trajectory),
    length.out = min(nrow(trajectory), max_points)
  )))
}

comparison_paths <- bind_rows(
  lapply(seq_along(mini_batch_results), function(i) {
    result <- mini_batch_results[[i]]
    keep <- path_keep(result$trajectory)
    data.frame(
      alpha = result$trajectory[keep, 1],
      beta = result$trajectory[keep, 2],
      method = names(mini_batch_results)[i]
    )
  }),
  list({
    keep <- path_keep(sgd_reference$trajectory)
    data.frame(
      alpha = sgd_reference$trajectory[keep, 1],
      beta = sgd_reference$trajectory[keep, 2],
      method = "Ordinary SGD"
    )
  })
)
comparison_paths$method <- factor(
  comparison_paths$method,
  levels = c("Ordinary SGD", names(mini_batch_results))
)

theta_ab_newton <- unname(theta_ab_sample)
contour_alpha <- range(comparison_paths$alpha, theta_ab_newton[1])
contour_beta <- range(comparison_paths$beta, theta_ab_newton[2])
contour_grid <- expand.grid(
  alpha = seq(contour_alpha[1], contour_alpha[2], length.out = 80),
  beta = seq(contour_beta[1], contour_beta[2], length.out = 80)
)
contour_grid$loss <- vapply(seq_len(nrow(contour_grid)), function(i) {
  loss_ab(c(contour_grid$alpha[i], contour_grid$beta[i]))
}, numeric(1))

mini_batch_contour_plot <- ggplot() +
  geom_contour(
    data = contour_grid,
    aes(x = alpha, y = beta, z = loss),
    color = "grey75"
  ) +
  geom_path(
    data = comparison_paths,
    aes(x = alpha, y = beta, color = method, group = method),
    linewidth = 0.7,
    alpha = 0.9
  ) +
  geom_point(
    data = data.frame(alpha = theta_ab_newton[1], beta = theta_ab_newton[2]),
    aes(x = alpha, y = beta),
    shape = 8,
    size = 3,
    color = "black"
  ) +
  facet_wrap(~method, ncol = 3) +
  labs(
    title = "SGD and mini-batch SGD paths",
    subtitle = paste0(
      "a = ", comparison_a,
      ", K = ", comparison_K,
      ", lr0 = ", comparison_lr0,
      "; star = constrained Newton estimate"
    ),
    x = expression(alpha),
    y = expression(beta),
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

initial_parameter_error <- sqrt(sum((theta_ab_start - theta_ab_newton)^2))

parameter_error_comparison <- bind_rows(
  data.frame(
    elapsed_time = c(0, sgd_reference$elapsed_time[seq_len(comparison_epochs)]),
    parameter_error = c(
      initial_parameter_error,
      sqrt(rowSums(
        (sgd_reference$trajectory[sgd_epoch_rows, 1:2, drop = FALSE] -
            matrix(theta_ab_newton, nrow = comparison_epochs, ncol = 2,
                  byrow = TRUE))^2
      ))
    ),
    epoch = 0:comparison_epochs,
    method = "Ordinary SGD"
  ),
  bind_rows(lapply(seq_along(mini_batch_results), function(i) {
    result <- mini_batch_results[[i]]
    updates_per_epoch <- ceiling(length(y) / mini_batch_sizes[i])
    epoch_rows <- seq_len(comparison_epochs) * updates_per_epoch + 1L
    data.frame(
      elapsed_time = c(0, result$elapsed_time),
      parameter_error = c(
        initial_parameter_error,
        sqrt(rowSums(
          (result$trajectory[epoch_rows, 1:2, drop = FALSE] -
               matrix(theta_ab_newton, nrow = comparison_epochs, ncol = 2,
                    byrow = TRUE))^2
        ))
      ),
      epoch = 0:comparison_epochs,
      method = names(mini_batch_results)[i]
    )
  }))
) %>%
  filter(parameter_error > 0, is.finite(parameter_error))

mini_batch_parameter_error_plot <- ggplot(
  parameter_error_comparison,
  aes(x = elapsed_time, y = parameter_error, color = method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1)) +
  scale_y_log10() +
  labs(
    title = "Distance to constrained Newton estimate over time",
    subtitle = paste0(
      "a = ", comparison_a,
      ", K = ", comparison_K,
      ", lr0 = ", comparison_lr0
    ),
    x = "Elapsed time (seconds)",
    y = "Euclidean distance to Newton estimate",
    color = "Method"
  ) +
  theme_minimal()

early_parameter_error_plot <- ggplot(
  filter(parameter_error_comparison, epoch <= 20),
  aes(x = epoch, y = parameter_error, color = method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  labs(
    title = "Distance to constrained Newton estimate: first 20 epochs",
    x = "Epoch",
    y = "Euclidean distance to Newton estimate",
    color = "Method"
  ) +
  theme_minimal()

loss_comparison <- bind_rows(
  data.frame(
    elapsed_time = c(0, sgd_reference$elapsed_time[seq_len(comparison_epochs)]),
    loss = c(
      loss_ab(theta_ab_start),
      sgd_reference$objective[seq_len(comparison_epochs)]
    ),
    method = "Ordinary SGD"
  ),
  bind_rows(lapply(seq_along(mini_batch_results), function(i) {
    result <- mini_batch_results[[i]]
    data.frame(
      elapsed_time = c(0, result$elapsed_time),
      loss = c(loss_ab(theta_ab_start), result$objective),
      method = names(mini_batch_results)[i]
    )
  }))
)

mini_batch_loss_plot <- ggplot(
  loss_comparison,
  aes(x = elapsed_time, y = loss, color = method)
) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(
    title = "Full-data loss over time",
    subtitle = paste0(
      "a = ", comparison_a,
      ", K = ", comparison_K,
      ", lr0 = ", comparison_lr0
    ),
    x = "Elapsed time (seconds)",
    y = "Mean squared-error objective",
    color = "Method"
  ) +
  theme_minimal()

print(mini_batch_contour_plot)
print(mini_batch_parameter_error_plot)
print(early_parameter_error_plot)
print(mini_batch_loss_plot)



# --- Comparing with Newton and GD


# Compare gradient descent, damped Newton, and SGD while optimizing only
# alpha and beta. Gamma and delta stay fixed at their true values.
theta_ab_start <- c(alpha = 2.0, beta = 1.0)
theta_ab_true <- theta_true[c("alpha", "beta")]
gamma_fixed <- theta_true["gamma"]
delta_fixed <- theta_true["delta"]

expand_theta_ab <- function(theta_ab) {
  c(
    alpha = theta_ab[1],
    beta = theta_ab[2],
    gamma = gamma_fixed,
    delta = delta_fixed
  )
}

fn_ab <- function(theta_ab) {
  full_loss(expand_theta_ab(theta_ab), x, y)
}

gr_ab <- function(theta_ab) {
  full_gradient(expand_theta_ab(theta_ab), x, y)[1:2]
}

hess_ab <- function(theta_ab) {
  full_hessian(expand_theta_ab(theta_ab), x, y)[1:2, 1:2]
}

gradient_descent_ab <- function(theta_init, tol = 1e-6, max_iter = 1000,
                                gamma_0 = 1, c = 0.1, d = 0.5) {
  theta <- theta_init
  trajectory <- matrix(NA_real_, nrow = max_iter + 1, ncol = 2)
  trajectory[1, ] <- theta
  time_history <- numeric(max_iter + 1)
  gradient_norm_history <- numeric(max_iter)
  start_time <- proc.time()[["elapsed"]]
  iterations <- 0

  for (i in seq_len(max_iter)) {
    gradient <- gr_ab(theta)
    gradient_norm_sq <- sum(gradient^2)
    gradient_norm_history[i] <- sqrt(gradient_norm_sq)

    if (sqrt(gradient_norm_sq) <= tol) {
      break
    }

    step_size <- gamma_0
    current_loss <- fn_ab(theta)
    while (fn_ab(theta - step_size * gradient) >
           current_loss - c * step_size * gradient_norm_sq) {
      step_size <- d * step_size
      if (step_size < 1e-16) {
        warning("Gradient descent line search failed.")
        break
      }
    }

    theta <- theta - step_size * gradient
    iterations <- i
    trajectory[i + 1, ] <- theta
    time_history[i + 1] <- proc.time()[["elapsed"]] - start_time
  }

  list(
    theta = theta,
    loss = fn_ab(theta),
    trajectory = trajectory[seq_len(iterations + 1), , drop = FALSE],
    time_history = time_history[seq_len(iterations + 1)],
    gradient_norm_history = gradient_norm_history[seq_len(iterations)],
    gradient_time_history = time_history[seq_len(iterations)]
  )
}

damped_newton_ab <- function(theta_init, tol = 1e-6, max_iter = 100,
                             c = 0.1, d = 0.5) {
  theta <- theta_init
  trajectory <- matrix(NA_real_, nrow = max_iter + 1, ncol = 2)
  trajectory[1, ] <- theta
  time_history <- numeric(max_iter + 1)
  gradient_norm_history <- numeric(max_iter)
  start_time <- proc.time()[["elapsed"]]
  iterations <- 0

  for (i in seq_len(max_iter)) {
    gradient <- gr_ab(theta)
    gradient_norm <- sqrt(sum(gradient^2))
    gradient_norm_history[i] <- gradient_norm

    if (gradient_norm <= tol) {
      break
    }

    direction <- tryCatch(
      solve(hess_ab(theta), -gradient),
      error = function(error) -gradient
    )
    if (sum(gradient * direction) >= 0) {
      direction <- -gradient
    }

    directional_derivative <- sum(gradient * direction)
    step_size <- 1
    current_loss <- fn_ab(theta)
    while (fn_ab(theta + step_size * direction) >
           current_loss + c * step_size * directional_derivative) {
      step_size <- d * step_size
      if (step_size < 1e-16) {
        warning("Newton line search failed.")
        break
      }
    }

    theta <- theta + step_size * direction
    iterations <- i
    trajectory[i + 1, ] <- theta
    time_history[i + 1] <- proc.time()[["elapsed"]] - start_time
  }

  list(
    theta = theta,
    loss = fn_ab(theta),
    trajectory = trajectory[seq_len(iterations + 1), , drop = FALSE],
    time_history = time_history[seq_len(iterations + 1)],
    gradient_norm_history = gradient_norm_history[seq_len(iterations)],
    gradient_time_history = time_history[seq_len(iterations)]
  )
}

result_gd_ab <- gradient_descent_ab(theta_ab_start)
result_newton_ab <- damped_newton_ab(theta_ab_start)
theta_ab_newton <- unname(result_newton_ab$theta)
result_sgd_ab <- sgd_two_parameters(
  x = x,
  y = y,
  theta_ab = theta_ab_start,
  gamma_fixed = gamma_fixed,
  delta_fixed = delta_fixed,
  theta_ab_target = theta_ab_newton,
  lr0 = 0.005,
  max_epochs = 200,
  a = 1,
  K = 2500,
  seed = 1
)

comparison_summary <- data.frame(
  method = c("Gradient descent", "Damped Newton", "SGD"),
  alpha = c(result_gd_ab$theta[1], result_newton_ab$theta[1],
            result_sgd_ab$theta[1]),
  beta = c(result_gd_ab$theta[2], result_newton_ab$theta[2],
           result_sgd_ab$theta[2]),
  loss = c(result_gd_ab$loss, result_newton_ab$loss,
           result_sgd_ab$final_loss),
  parameter_mse = c(
    mean((result_gd_ab$theta - theta_ab_newton)^2),
    mean((result_newton_ab$theta - theta_ab_newton)^2),
    result_sgd_ab$final_parameter_mse
  ),
  elapsed_seconds = c(
    tail(result_gd_ab$time_history, 1),
    tail(result_newton_ab$time_history, 1),
    tail(result_sgd_ab$elapsed_time, 1)
  )
)
print(comparison_summary)

sgd_epoch_rows <- seq_len(result_sgd_ab$epochs) * length(y) + 1L

comparison_parameter_mse_data <- bind_rows(
  data.frame(
    time = result_gd_ab$time_history,
    parameter_mse = rowMeans((result_gd_ab$trajectory -
      matrix(theta_ab_newton, nrow = nrow(result_gd_ab$trajectory), ncol = 2,
             byrow = TRUE))^2),
    method = "Gradient descent"
  ),
  data.frame(
    time = result_newton_ab$time_history,
    parameter_mse = rowMeans((result_newton_ab$trajectory -
      matrix(theta_ab_newton, nrow = nrow(result_newton_ab$trajectory), ncol = 2,
             byrow = TRUE))^2),
    method = "Damped Newton"
  ),
  data.frame(
    time = c(0, result_sgd_ab$elapsed_time),
    parameter_mse = c(
      mean((theta_ab_start - theta_ab_newton)^2),
      rowMeans((result_sgd_ab$trajectory[sgd_epoch_rows, , drop = FALSE] -
        matrix(theta_ab_newton, nrow = result_sgd_ab$epochs, ncol = 2,
               byrow = TRUE))^2)
    ),
    method = "SGD"
  )
) %>%
  filter(time >= 0, parameter_mse > 0, is.finite(parameter_mse))

comparison_parameter_mse_plot <- ggplot(
  comparison_parameter_mse_data,
  aes(x = time, y = parameter_mse, color = method)
) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1)) +
  scale_y_log10() +
  labs(
    title = "Parameter MSE to the Newton estimate",
    x = "Elapsed time (seconds)",
    y = "Mean squared parameter error",
    color = "Method"
  ) +
  theme_minimal()

comparison_path_data <- bind_rows(
  data.frame(result_gd_ab$trajectory, method = "Gradient descent"),
  data.frame(result_newton_ab$trajectory, method = "Damped Newton"),
  data.frame(
    result_sgd_ab$trajectory[unique(round(seq(
      1, nrow(result_sgd_ab$trajectory),
      length.out = min(1500, nrow(result_sgd_ab$trajectory))
    ))), , drop = FALSE],
    method = "SGD"
  )
)
names(comparison_path_data)[1:2] <- c("alpha", "beta")

comparison_alpha_range <- range(comparison_path_data$alpha, theta_ab_newton[1])
comparison_beta_range <- range(comparison_path_data$beta, theta_ab_newton[2])
comparison_contour <- expand.grid(
  alpha = seq(comparison_alpha_range[1], comparison_alpha_range[2], length.out = 80),
  beta = seq(comparison_beta_range[1], comparison_beta_range[2], length.out = 80)
)
comparison_contour$loss <- vapply(seq_len(nrow(comparison_contour)), function(i) {
  fn_ab(c(comparison_contour$alpha[i], comparison_contour$beta[i]))
}, numeric(1))

comparison_path_plot <- ggplot() +
  geom_contour(
    data = comparison_contour,
    aes(x = alpha, y = beta, z = loss),
    color = "grey70"
  ) +
  geom_path(
    data = comparison_path_data,
    aes(x = alpha, y = beta, color = method, group = method),
    linewidth = 0.8
  ) +
  geom_point(
    data = data.frame(alpha = theta_ab_start[1], beta = theta_ab_start[2]),
    aes(x = alpha, y = beta),
    shape = 4,
    size = 3,
    color = "black"
  ) +
  geom_point(
    data = data.frame(alpha = theta_ab_newton[1], beta = theta_ab_newton[2]),
    aes(x = alpha, y = beta),
    shape = 8,
    size = 3,
    color = "black"
  ) +
  labs(
    title = "Optimization paths in the alpha-beta plane",
    subtitle = "Gamma and delta are fixed at their true values; star = Newton estimate",
    x = expression(alpha),
    y = expression(beta),
    color = "Method"
  ) +
  theme_minimal()

print(comparison_parameter_mse_plot)
print(comparison_path_plot)


