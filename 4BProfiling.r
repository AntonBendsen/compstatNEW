library(profvis)
source("A4_comparefuncs.R")

profvis::profvis(
  for (i in seq_len(10L)) sgd(c(3,2,2,5), x, y, gradient_func,
                gamma_0 = 0.5, K = 100, a = 0.9,
                batch_size = 50, max_epochs = 500,
                cb = NULL)
)
