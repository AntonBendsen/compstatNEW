source("A4_comparefuncs.R")

N <- 10000
sigma <- 1e-5

log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))

y <- dose_response(x, theta_true) +
rnorm(N, mean = 0, sd = sigma)

theta_true <- c(2.5, 1.4, 0.5, 4.0)
theta_guess1 <- c(1, 2, 2, 4)
theta_guess2 <- c(2.3, 1.7, 0.1, 3.2)

y <- dose_response(x, theta_true) +
  rnorm(N, mean = 0, sd = sigma)

theta_hat1 <- sgd(theta_guess1, x, y, gradient_func,
                 gamma_0 = 0.5, K = 100, a = 0.9,
                 batch_size = 100, max_epochs = 100)$theta

theta_hat2 <- sgd(theta_guess2, x, y, gradient_func,
                 gamma_0 = 0.5, K = 100, a = 0.9,
                 batch_size = 100, max_epochs = 100)$theta

fn <- function(theta) loss_func(x, y, theta)
gr <- function(theta) gradient_func(x, y, theta)
hess <- function(theta) hess_full(theta, x, y)

sqrt(sum((theta_true - theta_guess)^2))
sqrt(sum((theta_true - theta_hat1)^2))
sqrt(sum((theta_true - theta_hat2)^2))
