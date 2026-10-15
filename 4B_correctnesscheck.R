source("A4_comparefuncs.R")

N <- 10000
sigma <- 0

set.seed(123) # Reproduce both the sampled doses and SGD batch shuffling.
log_grid <- seq(1, 5, length.out = 15)
x <- exp(sample(log_grid, N, replace = TRUE))

theta_true <- c(2.5, 1.4, 0.5, 4.0)
theta_guess1 <- c(1.9, 1, 2, 2)
theta_guess2 <- c(2.3, 1.7, 1, 5)

y <- dose_response(x, theta_true) +
  rnorm(N, mean = 0, sd = sigma)

theta_hat <- sgd(theta_true, x, y, gradient_func,
                 gamma_0 = 0.5, K = 100, a = 0.9,
                 batch_size = 100, max_epochs = 100)$theta

theta_hat1 <- sgd(theta_guess1, x, y, gradient_func,
                 gamma_0 = 0.5, K = 100, a = 0.9,
                 batch_size = 100, max_epochs = 100)$theta

theta_hat2 <- sgd(theta_guess2, x, y, gradient_func,
                 gamma_0 = 0.5, K = 100, a = 0.9,
                 batch_size = 100, max_epochs = 100)$theta

# Squared Euclidean errors before and after SGD (no square root).
initial_thetas <- list(True = theta_true, `Guess 1` = theta_guess1,
                      `Guess 2` = theta_guess2)
fitted_thetas <- list(True = theta_hat, `Guess 1` = theta_hat1,
                     `Guess 2` = theta_hat2)
correctness_results <- data.frame(
  start = names(initial_thetas),
  initial_squared_error = vapply(initial_thetas, function(theta) {
    sum((theta_true - theta)^2)
  }, numeric(1)),
  fitted_squared_error = vapply(fitted_thetas, function(theta) {
    sum((theta_true - theta)^2)
  }, numeric(1)),
  row.names = NULL
)
print(correctness_results, row.names = FALSE)

# Generate the slide table and curves from these same estimates.
dir.create("Pictures", showWarnings = FALSE)
theta_labels <- vapply(initial_thetas, function(theta) {
  paste0("$(", paste(theta, collapse = ", "), ")$")
}, character(1))
writeLines(c(
  "| Start | $\\theta_0$ | Before | After SGD |",
  "|:--|:--|--:|--:|",
  sprintf("| %s | %s | %.4f | %.4f |", correctness_results$start,
          theta_labels, correctness_results$initial_squared_error,
          correctness_results$fitted_squared_error)
), "Pictures/sgd-correctness-results.md")

dose_grid <- exp(seq(1, 5, length.out = 300))
correctness_curves <- do.call(rbind, lapply(names(initial_thetas), function(start) {
  curves <- list(True = theta_true, Guessed = initial_thetas[[start]],
                 Fitted = fitted_thetas[[start]])
  do.call(rbind, lapply(names(curves), function(curve) {
    data.frame(start = start, curve = curve, dose = dose_grid,
               response = dose_response(dose_grid, curves[[curve]]))
  }))
}))
correctness_curves$start <- factor(correctness_curves$start,
                                    levels = names(initial_thetas))
correctness_curves$curve <- factor(correctness_curves$curve,
                                    levels = c("True", "Guessed", "Fitted"))
correctness_plot <- ggplot2::ggplot(
  correctness_curves,
  ggplot2::aes(dose, response, color = curve, linetype = curve)
) +
  ggplot2::geom_line(
    data = correctness_curves[correctness_curves$curve != "True", ],
    linewidth = 0.9
  ) +
  # Draw the true curve last so it stays visible over the fitted curve.
  ggplot2::geom_line(
    data = correctness_curves[correctness_curves$curve == "True", ],
    linewidth = 0.9
  ) +
  ggplot2::facet_wrap(~start, ncol = 1) +
  ggplot2::scale_x_log10(breaks = c(3, 10, 30, 100)) +
  ggplot2::scale_color_manual(values = c(True = "black", Guessed = "#D55E00",
                                        Fitted = "#56B4E9")) +
  ggplot2::scale_linetype_manual(values = c(True = "longdash", Guessed = "dotted",
                                           Fitted = "solid")) +
  ggplot2::labs(x = "Dose (log scale)", y = "Response", color = NULL,
                linetype = NULL) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "bottom",
                 strip.text = ggplot2::element_text(face = "bold"))
ggplot2::ggsave("Pictures/sgd-correctness.png", correctness_plot,
                width = 6.5, height = 6.5, dpi = 160, bg = "white")
