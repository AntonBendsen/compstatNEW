source("A4_comparefuncs.R")

# Include N = 100,000 for the complete performance slides.
run_large_comparisons <- TRUE

# Every SGD callback point is evaluated and plotted; there is no thinning.
# Fitted curves use the best update strictly before the time limit.

cat("Computing N = 1,000 (good start)...\n", file = stderr())
set.seed(123)
comparison_1000 <- compare_optimizers(
  c(3, 2, 2, 5),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(100, 500, 500)
)

#comparison_1000$convergence_plot
#comparison_1000$curve_plot

cat("Computing N = 10,000 (good start)...\n", file = stderr())
set.seed(123)
comparison_10000 <- compare_optimizers(
  c(3, 2, 2, 5),
  N = 10000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(100, 500, 500)
)

#comparison_10000$convergence_plot
# comparison_10000$curve_plot

if (run_large_comparisons) {
  cat("Computing N = 100,000 (good start)...\n", file = stderr())
  set.seed(123)
  comparison_100000 <- compare_optimizers(
    c(3, 2, 2, 5),
    N = 100000,
    sigma = 0.25,
    variable_vec = c(1, 50, 250),
    mepoch_vec = c(100, 500, 500)
  )
}

# comparison_100000$convergence_plot
# comparison_100000$curve_plot

# #NEWTON AND GD GET STUCK:

cat("Computing N = 1,000 (poor start)...\n", file = stderr())
set.seed(123)
comparison_1000_p2 <- compare_optimizers(
  c(1, 4, 0, 2),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(100, 500, 500)
)

# comparison_1000_p2$convergence_plot
# comparison_1000_p2$curve_plot

cat("Computing N = 10,000 (poor start)...\n", file = stderr())
set.seed(123)
comparison_10000_p2 <- compare_optimizers(
  c(1, 4, 0, 2),
  N = 10000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250),
  mepoch_vec = c(100, 500, 500)
)

# comparison_10000_p2$convergence_plot
# comparison_10000_p2$curve_plot

if (run_large_comparisons) {
  cat("Computing N = 100,000 (poor start)...\n", file = stderr())
  set.seed(123)
  comparison_100000_p2 <- compare_optimizers(
    c(1, 4, 0, 2),
    N = 100000,
    sigma = 0.25,
    variable_vec = c(1, 50, 250),
    mepoch_vec = c(100, 500, 500)
  )
}

# comparison_100000_p2$convergence_plot
# comparison_100000_p2$curve_plot


# Compatibility names: the main comparisons are now fully update-traced.
comparison_1000_ut <- comparison_1000
comparison_1000_p2_ut <- comparison_1000_p2

# Best loss error across EVERY early update, strictly before the time limit.
comparison_10000_early <- best_curves_before_time(
  comparison_10000, time_limit = 1e-2
)
comparison_10000_early$curve_plot

# Poor starting guess: separate time limits for the three sample sizes.
comparison_1000_p2_early <- best_curves_before_time(
  comparison_1000_p2, time_limit = 0.001
)
comparison_10000_p2_early <- best_curves_before_time(
  comparison_10000_p2, time_limit = 0.01
)
if (run_large_comparisons) {
  comparison_100000_p2_early <- best_curves_before_time(
    comparison_100000_p2, time_limit = 0.1
  )
}

#CONTOUR PLOTS

# paths_10000 <- compare_optimizers_contour(
#   c(3, 2, 2, 5),
#   N = 1000,
#   sigma = 0.25,
#   variable_vec = c(1, 50, 250),
#   dims = c(1, 2)
# )

# paths_10000$contour_plot

# paths_10000_p2 <- compare_optimizers_contour(
#   c(1, 4, 0, 2),
#   N = 10000,
#   sigma = 0.25,
#   variable_vec = c(1, 50, 250),
#   dims = c(1, 2)
# )

# paths_10000_p2$contour_plot


# Regenerate the static performance images used by Assignment4B_slides.qmd.
# Slide rendering reads these PNGs and does not source this script.
make_performance_grid <- function(comparisons) {
  sample_sizes <- c("1,000", "10,000", "100,000")
  time_limits <- c(0.001, 0.01, 0.1)
  early_comparisons <- Map(best_curves_before_time, comparisons, time_limits)
  convergence_plots <- lapply(seq_along(comparisons), function(i) {
    comparisons[[i]]$convergence_plot +
      ggplot2::labs(title = paste("N =", sample_sizes[i]), x = "Time (s)") +
      ggplot2::geom_vline(
        xintercept = time_limits[i], color = "grey35", linetype = "dashed",
        linewidth = 0.6, show.legend = FALSE
      ) +
      ggplot2::guides(color = "none")
  })
  curve_plots <- lapply(seq_along(early_comparisons), function(i) {
    early_comparisons[[i]]$curve_plot +
      ggplot2::labs(title = paste("Best fit before", time_limits[i], "s")) +
      ggplot2::guides(
        color = ggplot2::guide_legend(nrow = 1),
        linetype = ggplot2::guide_legend(nrow = 1)
      )
  })
  patchwork::wrap_plots(
    c(convergence_plots, curve_plots), ncol = 3, nrow = 2,
    byrow = TRUE, guides = "collect"
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
}

if (run_large_comparisons) {
  dir.create("Pictures", showWarnings = FALSE)
  performance_grids <- list(
    good = list(comparison_1000, comparison_10000, comparison_100000),
    poor = list(comparison_1000_p2, comparison_10000_p2, comparison_100000_p2)
  )
  for (guess in names(performance_grids)) {
    image_path <- file.path("Pictures", paste0("performance-grid-", guess, ".png"))
    grDevices::png(image_path, width = 13.2, height = 5.6, units = "in", res = 96)
    tryCatch({
      grid::grid.newpage()
      grid::grid.draw(patchwork::patchworkGrob(
        make_performance_grid(performance_grids[[guess]])
      ))
    }, finally = grDevices::dev.off())
    cat("Saved ", image_path, "\n", file = stderr(), sep = "")
  }
  # Keep the slide's parameter footer in sync if these settings change.
}
