source("A4_comparefuncs.R")

set.seed(123)
comparison_1000 <- compare_optimizers(
  c(3, 2, 2, 5),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

# comparison_1000$convergence_plot
# comparison_1000$curve_plot

set.seed(123)
comparison_10000 <- compare_optimizers(
  c(3, 2, 2, 5),
  N = 10000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

comparison_10000$convergence_plot
# comparison_10000$curve_plot

set.seed(123)
comparison_100000 <- compare_optimizers(
  c(3, 2, 2, 5),
  N = 100000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

# comparison_100000$convergence_plot
# comparison_100000$curve_plot

# #NEWTON AND GD GET STUCK:

set.seed(123)
comparison_1000_p2 <- compare_optimizers(
  c(1, 4, 0, 2),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

# comparison_1000_p2$convergence_plot
# comparison_1000_p2$curve_plot

set.seed(123)
comparison_10000_p2 <- compare_optimizers(
  c(1, 4, 0, 2),
  N = 10000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

# comparison_10000_p2$convergence_plot
# comparison_10000_p2$curve_plot

set.seed(123)
comparison_100000_p2 <- compare_optimizers(
  c(1, 4, 0, 2),
  N = 100000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

# comparison_100000_p2$convergence_plot
# comparison_100000_p2$curve_plot


set.seed(123)
comparison_1000_ut <- compare_optimizers_ut(
  c(3, 2, 2, 5),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

comparison_1000_ut$convergence_plot

set.seed(123)
comparison_1000_p2_ut <- compare_optimizers_ut(
  c(1, 4, 0, 2),
  N = 1000,
  sigma = 0.25,
  variable_vec = c(1, 50, 250)
)

comparison_1000_p2_ut$convergence_plot

# Best recorded loss error for each algorithm strictly before the time limit.
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
comparison_100000_p2_early <- best_curves_before_time(
  comparison_100000_p2, time_limit = 0.1
)

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
