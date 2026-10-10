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

log_grid <- seq(1, 5, length.out = 15) ## 15 different doses
x <- exp(sample(log_grid, N, replace = TRUE))

y <- f(x, theta_true) + rnorm(N, sd = sigma)

dat <- data.frame(x = x, y = y)

# Tabsfunktion: mean squared error
loss <- function(theta, x, y) {
  mean((y - f(x, theta))^2)
}

loss(theta_true, x, y)


#Calculate the gradient of the loss function

# Beregn modelværdier og Jacobimatrix
model_jacobian <- function(x, theta) {
  alpha <- theta[1]
  beta  <- theta[2]
  gamma <- theta[3]
  delta <- theta[4]

  lx <- log(x)
  p <- plogis(-beta * (lx - alpha))
  q <- p * (1 - p)

  fit <- gamma + (delta - gamma) * p

  J <- cbind(
    alpha = (delta - gamma) * beta * q,
    beta  = -(delta - gamma) * (lx - alpha) * q,
    gamma = 1 - p,
    delta = p
  ) #analytic Jacobian maybe compare to the numerical one.

  list(fit = fit, J = J)
}

# Gradient af den gennemsnitlige kvadrerede fejl
grad_loss <- function(theta, x, y) {
  obj <- model_jacobian(x, theta)
  residual <- y - obj$fit

  as.vector(-2 * crossprod(obj$J, residual) / length(y))
}



sgd <- function(x, y, theta, lr = 0.01,
                epochs = 100, decay = 0.01) {
  n <- length(y)
  history <- numeric(epochs)

  for (epoch in seq_len(epochs)) {
    idx <- sample.int(n, size = n, replace = TRUE)

    for (i in idx) {
      g <- grad_loss(theta, x[i], y[i])
      eta <- lr / (1 + decay * (epoch - 1))
      theta <- theta - eta * g
    }

    history[epoch] <- loss(theta, x, y)
  }

  list(theta = theta, loss = loss(theta, x, y),
       history = history)
}

mini_batch_sgd <- function(x, y, theta, lr = 0.01,
                           batch_size = 32,
                           epochs = 100, decay = 0.01) {
  n <- length(y)
  history <- numeric(epochs)

  for (epoch in seq_len(epochs)) {
    idx <- sample.int(n)
    batches <- split(idx, ceiling(seq_along(idx) / batch_size))

    for (batch in batches) {
      g <- grad_loss(theta, x[batch], y[batch])
      eta <- lr / (1 + decay * (epoch - 1))
      theta <- theta - eta * g
    }

    history[epoch] <- loss(theta, x, y)
  }

  list(theta = theta, loss = loss(theta, x, y),
       history = history)
}

# Test SGD
# At a fixed parameter, the average single-observation gradient equals
# the gradient of the full empirical loss.
theta_test <- c(alpha = 2.2, beta = 1.2, gamma = 0.6, delta = 3.8)
individual_gradients <- vapply(
  seq_len(N),
  function(i) grad_loss(theta_test, x[i], y[i]),
  numeric(length(theta_test))
)
average_individual_gradient <- rowMeans(individual_gradients)
full_gradient <- grad_loss(theta_test, x, y)

stopifnot(isTRUE(all.equal(
  unname(average_individual_gradient),
  unname(full_gradient),
  tolerance = 1e-10
)))

# Check that SGD runs and reduces the empirical loss.
set.seed(456)
initial_loss <- loss(theta_test, x, y)
sgd_test <- sgd(x, y, theta_test, lr = 0.001, epochs = 30)

stopifnot(
  length(sgd_test$theta) == length(theta_test),
  all(is.finite(sgd_test$theta)),
  length(sgd_test$history) == 30,
  sgd_test$loss < initial_loss
)

cat("SGD test passed. Initial loss:", initial_loss,
    "Final loss:", sgd_test$loss, "\n")








#NEW: 


# ============================================================
# B: Log-Logistic Dose-Response Curves
# Computational Statistics
# ============================================================

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
  p <- plogis(-beta * (lx - alpha))
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


# 5. Full empirical loss and gradient
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
gd <- function(x, y, theta, lr = 0.001,
               max_epochs = 200, tol = 1e-8) {
  N <- length(y)
  history <- numeric(max_epochs)

  for (epoch in seq_len(max_epochs)) {
    old_theta <- theta

    g <- full_gradient(theta, x, y)
    theta <- theta - lr * g

    history[epoch] <- full_loss(theta, x, y)

    if (any(!is.finite(theta))) {
      stop("GD diverged. Try a smaller learning rate.")
    }

    if (sqrt(sum((theta - old_theta)^2)) <
        tol * (1 + sqrt(sum(old_theta^2)))) {
      history <- history[seq_len(epoch)]
      break
    }
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history
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
    }

    history[epoch] <- full_loss(theta, x, y)
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history
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
    }

    history[epoch] <- full_loss(theta, x, y)
  }

  list(
    theta = theta,
    loss = full_loss(theta, x, y),
    history = history
  )
}









theta_start <- c(
  alpha = 2.0,
  beta  = 1.0,
  gamma = 0.0,
  delta = 3.0
)

# Gradient descent
result_gd <- gd(
  x, y, theta_start,
  lr = 0.001,
  max_epochs = 200
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
