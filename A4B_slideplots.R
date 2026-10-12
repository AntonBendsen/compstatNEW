# Plot used by the simulated-data slide. Source this file, then call:
# make_slide_data_plot(N = 1000, sigma = 0.25)
# This file does not run the optimizer comparisons.

make_slide_data_plot <- function(N = 1000, sigma = 0.25, seed = 123) {
  # Same model and true parameters as A4_comparefuncs.R.
  dose_response <- function(x, theta) {
    alpha <- theta[1]
    beta <- theta[2]
    gamma <- theta[3]
    delta <- theta[4]
    gamma + (delta - gamma) * plogis(alpha - beta * log(x))
  }
  theta_true <- c(alpha = 2.5, beta = 1.4, gamma = 0.5, delta = 4.0)

  set.seed(seed)
  log_grid <- seq(1, 5, length.out = 15)
  x <- exp(sample(log_grid, N, replace = TRUE))

  y <- dose_response(x, theta_true) +
    rnorm(N, mean = 0, sd = sigma)

  dose_grid <- exp(seq(1, 5, length.out = 300))
  generating_curve <- data.frame(
    dose = dose_grid,
    response = dose_response(dose_grid, theta_true)
  )

  ggplot2::ggplot(data.frame(dose = x, response = y),
                  ggplot2::aes(x = dose, y = response)) +
    ggplot2::geom_point(color = "grey40", alpha = 0.3, size = 1.5) +
    ggplot2::geom_line(
      data = generating_curve, color = "black",
      linetype = "longdash", linewidth = 1
    ) +
    ggplot2::scale_x_log10(
      breaks = scales::breaks_log(n = 5), labels = scales::label_number()
    ) +
    ggplot2::labs(
      title = "Simulated dose-response data",
      subtitle = "Dashed line is the data-generating curve",
      x = "Dose (log scale)", y = "Response"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, color = "black"),
      plot.subtitle = ggplot2::element_text(size = 13, color = "black"),
      axis.title = ggplot2::element_text(color = "black")
    )
}
