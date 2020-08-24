#!/usr/bin/env Rscript

# generate dataset with certain seed
set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/cellrouter",
  num_cells = 200,
  num_features = 300,
  model = "linear",
  normalise = FALSE
)

# add method specific args (if needed)
data$parameters <- list()
data$seed <- 1L

# write example dataset to file
file <- commandArgs(trailingOnly = TRUE)[[1]]
dynutils::write_h5(data, file)
