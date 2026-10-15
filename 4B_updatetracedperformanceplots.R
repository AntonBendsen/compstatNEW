source("A4_comparefuncs.R")


set.seed(123)
comparison_1000_ut <- compare_optimizers_ut(
  c(3, 2, 2, 5),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(100, 500, 500)
)

comparison_1000_ut$convergence_plot

set.seed(123)
comparison_1000_p2_ut <- compare_optimizers_ut(
  c(1, 4, 0, 2),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(100, 500, 500)
)

comparison_1000_p2_ut$convergence_plot


set.seed(123)
comparison_lr1 <- compare_optimizers_ut(
  c(1, 4, 0, 2),
  N = 10000,
  sigma = 0.25,
  variable = "gamma_0",
  variable_vec = c(0.5, 2, 3),
  mepoch_vec = c(500, 500, 500)
)

#comparison_lr1$convergence_plot

set.seed(123)
comparison_lr2 <- compare_optimizers_ut(
  c(3, 0, 2, 7),
  N = 10000,
  sigma = 0.25,
  variable = "gamma_0",
  variable_vec = c(0.5, 2, 3),
  mepoch_vec = c(500, 500, 500)
)

#comparison_lr2$convergence_plot



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

# Learning-rate slide: convergence above the best curves strictly before 5 ms.
comparison_lr1_early <- best_curves_before_time(comparison_lr1, time_limit = 0.005)
comparison_lr2_early <- best_curves_before_time(comparison_lr2, time_limit = 0.005)
learning_rate_comparisons <- list(comparison_lr1, comparison_lr2)
learning_rate_early <- list(comparison_lr1_early, comparison_lr2_early)
learning_rate_start_labels <- c("Start: (1, 4, 0, 2)", "Start: (3, 0, 2, 7)")

learning_rate_convergence_plots <- lapply(seq_along(learning_rate_comparisons), function(i) {
  learning_rate_comparisons[[i]]$convergence_plot +
    ggplot2::labs(title = learning_rate_start_labels[i], x = "Time (s)") +
    ggplot2::geom_vline(
      xintercept = 0.005, color = "grey35", linetype = "dashed",
      linewidth = 0.6, show.legend = FALSE
    ) +
    ggplot2::guides(color = "none")
})
learning_rate_curve_plots <- lapply(seq_along(learning_rate_early), function(i) {
  learning_rate_early[[i]]$curve_plot +
    ggplot2::labs(title = "Best fit before 5 ms") +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1),
      linetype = ggplot2::guide_legend(nrow = 1)
    )
})

learning_rate_performance_grid <- patchwork::wrap_plots(
  c(learning_rate_convergence_plots, learning_rate_curve_plots),
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

# Rendering the slides reads this image without rerunning the comparisons.
dir.create("Pictures", showWarnings = FALSE)
ggplot2::ggsave(
  "Pictures/performance-grid-learning-rates.png", learning_rate_performance_grid,
  width = 13.2, height = 5.6, dpi = 120, bg = "white"
)
