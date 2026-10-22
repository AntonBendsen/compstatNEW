library(here)
library(Rcpp)
library(bench)
library(gridExtra)
library(tidyverse)
library(CSwR)


sourceCpp(here::here("minibatch.cpp"))

# Data and helpers needed to run this script on its own.
theta_true <- c(alpha = 2.5, beta = 1.4, gamma = 0.5, delta = 4.0)
set.seed(123)
N <- 5000
sigma <- 0.25
log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))
y <- dose_response(x, theta_true) + rnorm(N, sd = sigma)

full_loss <- function(theta, x, y) {
  mean(vapply(seq_along(y), function(i) loss_i(theta, x[i], y[i]), numeric(1)))
}

full_gradient <- function(theta, x, y) {
  gradient <- numeric(length(theta))
  for (i in seq_along(y)) {
    gradient <- gradient + grad_i(theta, x[i], y[i])
  }
  gradient / length(y)
}

mini_batch_sgd <- function(x, y, theta, batch_size = 32, lr0 = 0.001,
                           max_epochs = 200, a = 0.6, K = 100,
                           tol = 1e-4, seed = 1, cb = NULL) {
  set.seed(seed)
  N <- length(y)
  max_steps <- max_epochs * ceiling(N / batch_size)
  trajectory <- matrix(NA_real_, nrow = max_steps + 1, ncol = length(theta))
  trajectory[1, ] <- theta
  time_history <- numeric(max_steps + 1)
  gradient_norm_history <- numeric(max_steps)
  gradient_time_history <- numeric(max_steps)
  convergence_gradient_history <- numeric(max_epochs)
  convergence_time_history <- numeric(max_epochs)
  history <- numeric(max_epochs)
  start_time <- proc.time()[["elapsed"]]
  step <- 0
  converged <- FALSE

  for (epoch in seq_len(max_epochs)) {
    indices <- sample.int(N)
    batches <- split(indices, ceiling(seq_along(indices) / batch_size))

    for (batch in batches) {
      step <- step + 1
      gradient <- Reduce(
        `+`, lapply(batch, function(i) grad_i(theta, x[i], y[i]))
      ) / length(batch)
      gradient_norm_history[step] <- sqrt(sum(gradient^2))
      grad_norm_sq <- sum(gradient^2)
      gradient_time_history[step] <- proc.time()[["elapsed"]] - start_time
      if (!is.null(cb)) cb()

      theta <- theta - lr0 / (1 + step / K)^a * gradient
      if (any(!is.finite(theta))) {
        stop("Mini-batch SGD diverged. Try a smaller lr0.")
      }
      trajectory[step + 1, ] <- theta
      time_history[step + 1] <- proc.time()[["elapsed"]] - start_time
    }

    history[epoch] <- full_loss(theta, x, y)
    convergence_gradient_history[epoch] <- sqrt(sum(full_gradient(theta, x, y)^2))
    convergence_time_history[epoch] <- proc.time()[["elapsed"]] - start_time
    if (convergence_gradient_history[epoch] < tol) {
      converged <- TRUE
      break
    }
  }

  epochs_run <- epoch
  if (!converged) {
    warning("Mini-batch SGD did not reach the gradient-norm tolerance before max_epochs.")
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history[seq_len(epochs_run)],
    trajectory = trajectory[seq_len(step + 1), , drop = FALSE],
    time_history = time_history[seq_len(step + 1)],
    gradient_norm_history = gradient_norm_history[seq_len(step)],
    gradient_time_history = gradient_time_history[seq_len(step)],
    convergence_gradient_history = convergence_gradient_history[seq_len(epochs_run)],
    convergence_time_history = convergence_time_history[seq_len(epochs_run)],
    epochs = epochs_run,
    full_gradient_norm = convergence_gradient_history[epochs_run],
    converged = converged
  )
}

make_path_data <- function(result, method) {
  data.frame(
    alpha = result$trajectory[, 1],
    beta = result$trajectory[, 2],
    iteration = seq_len(nrow(result$trajectory)) - 1L,
    method = method
  )
}

# 8.2 Mini-batch SGD (Rcpp)
mini_batch_sgd_rcpp <- function(x, y, theta, batch_size = 32,
                           lr0 = 0.001, max_epochs = 200,
                           a = 0.6, K = 100, tol = 1e-4, seed = 1,
                           cb = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
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
    
    epoch_start_time <- proc.time()[["elapsed"]] - start_time
    res <- sgd_batch(x, y, theta, batches, step, lr0, K, a, trajectory)
    first_step <- step + 1L
    theta <- res$theta
    trajectory <- res$trajectory
    step <- res$step
    steps <- seq.int(first_step, step)
    gradient_norm_history[steps] <- res$gradient_norm
    gradient_time_history[steps] <- epoch_start_time + res$gradient_time
    time_history[steps + 1L] <- epoch_start_time + res$elapsed_time

    if (!is.null(cb)) {
      grad_norm_sq <- tail(res$gradient_norm, 1)^2
      cb()
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


# Mini-batch SGD
result_mb_rcpp <- mini_batch_sgd_rcpp(
  x, y, theta_start,
  batch_size = 32,
  lr0 = 0.001,
  max_epochs = 1000,
  a = 0.6,
  K = 10000,
  #seed = 1,
  cb = mb_tracer$tracer
)


result_mb <- mini_batch_sgd(
  x, y, theta_start,
  batch_size = 32,
  lr0 = 0.001,
  max_epochs = 1000,
  a = 0.6,
  K = 10000,
  #seed = 2,
  cb = mb_tracer$tracer
)

## 34% tidsbesparelse! Det er sgu meget godt.
1 - max(result_mb_rcpp$time_history, na.rm=TRUE)/max(result_mb$time_history, na.rm=TRUE)


## Compare paths:
trajectory_plot_data <- bind_rows(
  make_path_data(result_mb, "Mini-batch SGD"),
  make_path_data(result_mb_rcpp, "Mini-batch SGD (Rcpp)")
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


resid <- data.frame(
  "alpha_resid" = result_mb$trajectory[,1] -  result_mb_rcpp$trajectory[,1] ,
  "beta_resid" = result_mb$trajectory[,2] -  result_mb_rcpp$trajectory[,2] ,
  "gamma_resid" = result_mb$trajectory[,3] -  result_mb_rcpp$trajectory[,3] ,
  "delta_resid" = result_mb$trajectory[,4] -  result_mb_rcpp$trajectory[,4]
) |>
  mutate(
    step = row_number()
  )

## Nedenstående viser, at der IKKE er tale om, at der er "samme" path igen og igen - tror bare fejlen er så lille (som man kan se nedenfor) at man ikke kan se forskellen i paths!

alpha_resid_plot <- ggplot(resid, aes(x = step, y = alpha_resid)) + geom_line(color = "red", alpha = 0.6) + 
  geom_hline(yintercept=0, linetype = "dashed") + ylab("Alpha, residual")
beta_resid_plot <- ggplot(resid, aes(x = step, y = beta_resid)) + geom_line(color = "blue", alpha = 0.6) + 
  geom_hline(yintercept=0, linetype = "dashed") + ylab("Beta, residual")
gamma_resid_plot <- ggplot(resid, aes(x= step, y = gamma_resid)) + geom_line(color = "green", alpha = 0.6) + 
  geom_hline(yintercept=0, linetype = "dashed")  + ylab("Gamma, residual")
delta_resid_plot <- ggplot(resid, aes(x = step, y = delta_resid)) + geom_line(color = "yellow", alpha = 0.6) + 
  geom_hline(yintercept=0, linetype = "dashed") + ylab("Delta, residual")

resid_plot <- grid.arrange(alpha_resid_plot, beta_resid_plot, gamma_resid_plot, delta_resid_plot, nrow = 2)

bench <- bench::mark(
  mini_batch = mini_batch_sgd(x, y, theta_start,
                              batch_size = 32,
                              lr0 = 0.001,
                              max_epochs = 2000,
                              a = 0.6,
                              K = 10000,
                              #seed = 2,
                              cb = mb_tracer$tracer) ,
  mini_batch_rcpp = mini_batch_sgd_rcpp(x, y, theta_start,
                                batch_size = 32,
                                lr0 = 0.001,
                                max_epochs = 2000,
                                a = 0.6,
                                K = 10000,
                                #seed = 2,
                                cb = mb_tracer$tracer) ,
  check = FALSE
)

print(bench)