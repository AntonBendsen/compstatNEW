library(tidyverse)

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

N <- 2000
sigma <- 0.25

# 15 grid points on the log-dose scale
log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))

y <- dose_response(x, theta_true) +
  rnorm(N, mean = 0, sd = sigma)


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


# 6. Gradient Descent
fn <- function(theta) full_loss(theta, x, y)
gr <- function(theta) full_gradient(theta, x, y)


gradient_descent <- function(theta_init, fn, gr, tol = 1e-8, max_iter = 1000,
                             gamma_0 = 1.0, c = 0.1, d = 0.5) {
  theta <- theta_init
  gamma <- gamma_0
  history <- numeric(max_iter + 1)
  history[1] <- fn(theta)
  trajectory <- matrix(NA_real_, nrow = max_iter + 1, ncol = length(theta))
  trajectory[1, ] <- theta
  time_history <- numeric(max_iter + 1)
  start_time <- proc.time()[["elapsed"]]
  iterations <- 0

  for (i in seq_len(max_iter)) {
    grad <- gr(theta)
    grad_norm_sq <- sum(grad^2)

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
    time_history = time_history[seq_len(iterations + 1)]
  )
}

# 7. Stochastic Gradient Descent
sgd <- function(x, y, theta, lr0 = 0.001,
                max_epochs = 200, a = 0.6, K = 100,
                seed = 1) {
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

  for (epoch in seq_len(max_epochs)) {
    # N random draws WITH replacement per epoch
    for (j in seq_len(N)) {
      i <- sample.int(N, 1)
      step <- step + 1

      # Robbins-Monro style decreasing learning rate
      lr <- lr0 / (1 + step / K)^a

      g <- grad_i(theta, x[i], y[i])
      theta <- theta - lr * g

      if (any(!is.finite(theta))) {
        stop("SGD diverged. Try a smaller lr0.")
      }

      trajectory[step + 1, ] <- theta
      time_history[step + 1] <- proc.time()[["elapsed"]] - start_time
    }

    history[epoch] <- full_loss(theta, x, y)
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history,
    trajectory = trajectory[seq_len(step + 1), , drop = FALSE],
    time_history = time_history[seq_len(step + 1)]
  )
}


# 8. Mini-batch SGD
mini_batch_sgd <- function(x, y, theta, batch_size = 32,
                           lr0 = 0.001, max_epochs = 200,
                           a = 0.6, K = 100, seed = 1) {
  set.seed(seed)
  N <- length(y)
  history <- numeric(max_epochs)
  step <- 0
  max_steps <- max_epochs * ceiling(N / batch_size)
  trajectory <- matrix(NA_real_, nrow = max_steps + 1, ncol = length(theta))
  trajectory[1, ] <- theta
  time_history <- numeric(max_steps + 1)
  start_time <- proc.time()[["elapsed"]]

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

      theta <- theta - lr * g

      if (any(!is.finite(theta))) {
        stop("Mini-batch SGD diverged. Try a smaller lr0.")
      }

      trajectory[step + 1, ] <- theta
      time_history[step + 1] <- proc.time()[["elapsed"]] - start_time
    }

    history[epoch] <- full_loss(theta, x, y)
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history,
    trajectory = trajectory[seq_len(step + 1), , drop = FALSE],
    time_history = time_history[seq_len(step + 1)]
  )
}





theta_start <- c(
  alpha = 2.0,
  beta  = 1.0,
  gamma = 0.0,
  delta = 3.0
)

# Gradient descent
result_gd <- gradient_descent(
  theta_init = theta_start,
  fn = fn,
  gr = gr
)

# SGD
result_sgd <- sgd(
  x, y, theta_start,
  lr0 = 0.001,
  max_epochs = 200,
  a = 0.6,
  K = 10000,
  seed = 1
)

# Mini-batch SGD
result_mb <- mini_batch_sgd(
  x, y, theta_start,
  batch_size = 32,
  lr0 = 0.001,
  max_epochs = 200,
  a = 0.6,
  K = 10000,
  seed = 1
)

# Sammenlign parameterestimater
results <- rbind(
  True = theta_true,
  GD = result_gd$theta,
  SGD = result_sgd$theta,
  MiniBatch = result_mb$theta
)

print(round(results, 4))

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

distance_by_time_plot <- ggplot(
  distance_data[distance_data$time > 0, ],
  aes(x = time, y = distance, color = method)
) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Parameter distance to truth over elapsed time",
    x = "Elapsed time (seconds)",
    y = "Euclidean distance to true parameters",
    color = "Method"
  ) +
  theme_minimal()

print(distance_by_iteration_plot)
print(distance_by_time_plot)
