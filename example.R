#!/usr/local/bin/Rscript

# generate dataset with certain seed
set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/cellrouter",
  num_cells = 200,
  num_features = 300,
  model = "tree",
  normalise = FALSE
)
set.seed(1)
data$counts@x <- rpois(length(data$counts@x), 10)
data$expression@x <- runif(length(data$expression@x))

# add method specific args (if needed)
data$params <- list()
data$seed <- 1

# write example dataset to file
file <- commandArgs(trailingOnly = TRUE)[[1]]
dynutils::write_h5(data, file)
