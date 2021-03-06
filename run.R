#!/usr/bin/env Rscript

task <- dyncli::main()

# load software
library(dyncli, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)

cellrouter_root <- "/cellrouter"
source(paste0(cellrouter_root, "/CellRouter_Class.R"))

#####################################
###           LOAD DATA           ###
#####################################

# load data
parameters <- task$parameters
expression <- as.matrix(task$expression)
start_id <- task$priors$start_id

# TIMING: done with preproc
timings <- list(method_afterpreproc = Sys.time())

#####################################
###        INFER TRAJECTORY       ###
#####################################

# steps are taken from https://github.com/edroaldo/cellrouter/blob/master/Myeloid_Progenitors/CellRouter_Paul_Tutorial.md
cellrouter <- CellRouter(as.data.frame(Matrix::t(expression)))
cellrouter <- scaleData(cellrouter)

# do pca
num_pcs <- pmin(parameters$ndim_pca, ncol(expression) - 1)
cellrouter <- computePCA(cellrouter, num.pcs = num_pcs, seed = NULL) # seed is set by dyncli

# do safe tsne
if (parameters$max_iter == "Inf") {
  parameters$max_iter <- 100000
}

again <- TRUE
while (again) {
  tryCatch({
    cellrouter <- computeTSNE(
      cellrouter,
      num.pcs = parameters$ndim_tsne,
      max_iter = parameters$max_iter,
      perplexity = parameters$perplexity,
      seed = NULL
    )
    again <- FALSE
  }, error = function(e) {
    if (grepl("Perplexity is too large", e$message)) {
      # don't forget both <<-'s, otherwise this will result in an infinite loop
      again <<- TRUE
      cat("TSNE: Perplexity is too large. Reducing perplexity by 25%: ", parameters$perplexity, " -> ", parameters$perplexity * .75, "\n", sep = "")
      parameters$perplexity <<- parameters$perplexity * .75
    } else {
      stop(e)
    }
  })
}

# louvain clustering
cellrouter <- findClusters(
  cellrouter,
  method = "graph.clustering",
  num.pcs = min(parameters$ndim_pca_clustering, num_pcs),
  k = parameters$k_clustering
)

# do knn
cellrouter <- buildKNN(
  cellrouter,
  k = parameters$k_knn,
  column.ann = 'population',
  num.pcs = min(parameters$ndim_pca_knn, num_pcs),
  sim.type = parameters$sim_type
)

# create trajectory using start cells as source
outputdir <- dynutils::safe_tempdir("cellrouter")
on.exit(unlink(outputdir, recursive = TRUE))
filename <- file.path(outputdir, "cell_edge_weighted_network.txt")
write.table(cellrouter@graph$edges, file = filename, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

sources <- unique(cellrouter@sampTab$population[cellrouter@sampTab$sample_id %in% start_id])
targets <- setdiff(as.vector(cellrouter@sampTab$population), sources)

# this function uses global variables...
libdir <- paste0(cellrouter_root, "/CellRouter")
cellrouter <- findPaths(
  cellrouter,
  column = 'population',
  libdir,
  outputdir,
  method = parameters$distance_method_paths
)

# process trajectories
cellrouter <- processTrajectories(
  cellrouter,
  rownames(cellrouter@ndata),
  path.rank = parameters$ranks,
  num.cells = parameters$num_cells,
  neighs = parameters$neighs,
  column.ann = 'population',
  column.color = 'colors'
)

timings$method_aftermethod <- Sys.time()

#####################################
###     SAVE OUTPUT TRAJECTORY    ###
#####################################
# first get network of backbone cells
backbone_network <-
  map_df(
    cellrouter@paths$path,
    function(x) {
      order <- tail(head(stringr::str_split(x, "->")[[1]], -1), -1)
      tibble(
        from = order[-length(order)],
        to = lead(order)[-length(order)],
        directed = TRUE,
        length = 1
      )
    }
  )

# now get for every non-backbone cell the shortest backbone cell
backbone_cells <- unique(c(backbone_network$from, backbone_network$to))
nonbackbone_cells <- setdiff(rownames(expression), backbone_cells)

nonbackbone_network <-
  distances(cellrouter@graph$network, backbone_cells, nonbackbone_cells) %>%
  apply(2, which.min) %>%
  {backbone_cells[.]} %>%
  set_names(nonbackbone_cells) %>%
  enframe("from", "to") %>%
  mutate(length = 1, directed = TRUE)

# combine to cell_graph & remove duplicated edges
cell_graph <- bind_rows(
  backbone_network,
  nonbackbone_network
) %>%
  group_by(from, to) %>%
  filter(row_number() == 1) %>%
  ungroup()

to_keep <- backbone_cells

# dimred
dimred <-
  cellrouter@tsne %>%
  as.data.frame() %>%
  rownames_to_column("cell_id")

# save output
output <-
  wrap_data(
    cell_ids = unique(c(cell_graph$from, cell_graph$to))
  ) %>%
  add_dimred(
    dimred = dimred
  )%>%
  add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep
  )  %>%
  add_timings(
    timings = timings
  )

dyncli::write_output(output, task$output)
