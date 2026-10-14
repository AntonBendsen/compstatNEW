source("A4_comparefuncs.R")


set.seed(123)
comparison_1000_ut <- compare_optimizers_ut(
  c(3, 2, 2, 5),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(200, 500, 500)
)

set.seed(123)
comparison_1000_p2_ut <- compare_optimizers_ut(
  c(1, 4, 0, 2),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(200, 500, 500)
)

# Use each algorithm's best recorded theta strictly before 5 milliseconds.
comparison_1000_ut_early <- best_curves_before_time(
  comparison_1000_ut, time_limit = 0.005
)
comparison_1000_p2_ut_early <- best_curves_before_time(
  comparison_1000_p2_ut, time_limit = 0.005
)

update_traced_comparisons <- list(comparison_1000_ut, comparison_1000_p2_ut)
update_traced_early <- list(comparison_1000_ut_early, comparison_1000_p2_ut_early)
start_labels <- c("Good start: (3, 2, 2, 5)", "Poor start: (1, 4, 0, 2)")

update_traced_convergence_plots <- lapply(seq_along(update_traced_comparisons), function(i) {
  update_traced_comparisons[[i]]$convergence_plot +
    ggplot2::labs(title = start_labels[i], x = "Time (s)") +
    ggplot2::guides(color = "none")
})
update_traced_curve_plots <- lapply(seq_along(update_traced_early), function(i) {
  update_traced_early[[i]]$curve_plot +
    ggplot2::labs(title = "Best fit before 5 ms") +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1),
      linetype = ggplot2::guide_legend(nrow = 1)
    )
})

# Display this object after sourcing to see the same 2 x 2 grid as the slide.
update_traced_performance_grid <- patchwork::wrap_plots(
  c(update_traced_convergence_plots, update_traced_curve_plots),
  ncol = 2, nrow = 2, byrow = TRUE, guides = "collect"
) &
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 13, hjust = 0.5),
    axis.title = ggplot2::element_text(size = 10),
    axis.text = ggplot2::element_text(size = 9),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(4, 6, 2, 2)
  )
