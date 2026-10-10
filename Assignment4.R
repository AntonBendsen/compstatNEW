library(tidyverse)

rnorm(100)


test <- rnorm(100)
test



theta_true <- c(alpha = 2.5, beta = 1.4,
                gamma = 0.5, delta = 4)

# Log-logistisk model
f <- function(x, theta) {
  alpha <- theta[1]
  beta  <- theta[2]
  gamma <- theta[3]
  delta <- theta[4]

  p <- plogis(-beta * (log(x) - alpha))
  gamma + (delta - gamma) * p
}

# Simulér data på et grid
set.seed(123)

N <- 2000
sigma <- 0.25

log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))

y <- f(x, theta_true) + rnorm(N, sd = sigma)

dat <- data.frame(x = x, y = y)

# Tabsfunktion: mean squared error
loss <- function(theta, x, y) {
  mean((y - f(x, theta))^2)
}

loss(theta_true, x, y)


runif(67)
