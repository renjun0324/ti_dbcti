#!/usr/local/bin/Rscript

task <- dyncli::main()

library(dbcti, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)

source("/code/dbcti.R")

#-----------------------------------------------------------------------------
#
#                           satisfy r cmd check
#
#-----------------------------------------------------------------------------

# library(dynwrap)
# library(dyneval)
# library(dynmethods)
# library(dynplot)
# library(dyncli)
# library(dynparam)
# library(tidyverse)
# 
# data_path = "/share/data6/tmp/renjun/CellTrekResult/CellTrek/benchmark_datasets/synthetic/dyntoy/bifurcating/bifurcating_1.rds"
# dataset = readRDS(data_path)
# dataset = add_cell_waypoints(dataset)
# dataset_wrap = wrap_expression(expression = dataset$expression,
#                                counts = dataset$counts,
#                                prior_information = dataset$prior_information)
# expression = dataset_wrap$expression
# priors = dataset_wrap$prior_information
# start_id = priors$start_id
# parameters = list(normalized = TRUE,
#                   gene_cri = 1,
#                   cell_cri = 1,
#                   use_normalized_data = TRUE,
#                   gene_number = 50,
#                   ndraw = 100,
#                   expansion = 1.5,
#                   r = 5)

expression <- task$expression
parameters <- task$parameters
start_id <- task$priors$start_id

#-----------------------------------------------------------------------------
#
#                                  create
#
#-----------------------------------------------------------------------------

## TIMING - start
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

## create object
obj <- as.data.frame( t(as.matrix(expression)) )
obj <- create_object(obj, normalized = parameters$normalized) 
obj <- filter_data(obj, 
                   gene_cri = parameters$gene_cri, 
                   cell_cri = parameters$cell_cri, 
                   use_normalized_data = parameters$use_normalized_data) 
obj <- select_var_feature(obj, 
                          use_normalized_data = parameters$use_normalized_dat, 
                          n = parameters$gene_number)

## plot
obj <- tsneplot2(obj, 
                 use_normalized_data = parameters$use_normalized_dat, 
                 specified_gene = FALSE, 
                 pca = TRUE, 
                 perplexity = 10)
obj <- contour_plot(obj)

## Estimate distribution
obj <- distribution_estimation(obj, 
                               ndraw = parameters$ndraw, 
                               expansion = parameters$expansion, 
                               ... = 1,2,3) 

## Calculate possibility for points
obj <- point_possibility(obj, r = parameters$r)

## Connect cluster
obj <- connect_cluster(obj)

## Infer_trajectory
obj <- infer_trajectory2(obj)

## Calculate pseudotime
index = which(rownames(expression)==start_id)
start_state_name = obj@distribution_estimation$cluster_index[index] %>% as.character
if(start_state_name==0) start_state_name = "1"
obj <- dbcti::calculate_pseudotime(obj, start_state_name = start_state_name)

## TIMING - end
checkpoints$method_aftermethod <- as.numeric(Sys.time())

#-----------------------------------------------------------------------------
#
#                         process celltrek output
#
#-----------------------------------------------------------------------------

pseudotime = obj@pseudotime
trajectory = obj@trajectory
connection_matrix = obj@connect_cluster$cluster_connection
g <- igraph::graph_from_adjacency_matrix(connection_matrix, 
                                         mode = "undirected")

## cluster_network
g_edge <- igraph::as_data_frame(g)
g_edge$length <- 0
tmp = obj@connect_cluster$vague_ratio
dist_m = matrix(0, nr = length(tmp), nc = length(tmp)) %>% data.frame
rownames(dist_m) = colnames(dist_m) = as.character(1:length(tmp))
for(i in names(tmp)){
  a = strsplit(i, "_")[[1]][2] %>% as.numeric
  b = strsplit(i, "_")[[1]][3] %>% as.numeric
  dist_m[a,b] = dist_m[b,a] = tmp[[i]]
}
for(i in 1:nrow(g_edge)){
  g_edge$length[i] = dist_m[g_edge$from[i], g_edge$to[i]]
}
g_edge = g_edge %>% mutate( directed = FALSE) %>% tibble
 
## create lin_df
lapply(1:length(obj@pseudotime), function(i){
  t = data.frame(obj@pseudotime[[i]],
                 names(obj@pseudotime[[i]]),
                 rownames(expression)[as.numeric(names(obj@pseudotime[[i]]))],
                 obj@trajectory$result_list[[i]],
                 lineage = i)
  colnames(t) = c("pseudot", "index", "cell_id", "comp_1", "comp_2","lineage")
  return(t)
}) -> tmplist
lin_df = do.call(rbind, tmplist) %>% data.frame

## dimred
dimred = lin_df[,c("comp_1", "comp_2")] %>%
  magrittr::set_rownames(lin_df$cell_id)

## percentages
percentages <- tibble(cell_id = lin_df$cell_id, 
                      milestone_id = as.character(lin_df$lineage), 
                      percentage = 1)

## return output
output <-
  dynwrap::wrap_data(
    cell_ids = rownames(dimred)
  ) %>%
  dynwrap::add_trajectory(
    milestone_network = g_edge,
    milestone_percentages = percentages
  ) %>%
  dynwrap::add_dimred(
    dimred = as.matrix(dimred)
  ) %>%
  dynwrap::add_pseudotime(pseudotime = lin_df$pseudot %>%
                            magrittr::set_names(lin_df$cell_id) ) %>%
  dynwrap::add_timings(checkpoints)

dyncli::write_output(output, task$output)

## dimred
# g_comp <- igraph::components(g)
# n <- tail(sort(table(g_comp$membership)), 1) + 1
# ordered_pseudotime <- pseudotime[[1]]
# ordered_trajectory <- trajectory$result_list[[1]]
# for (i in pseudotime[-1]) ordered_pseudotime <- rbind(as.matrix(i),
#                                                       as.matrix(ordered_pseudotime))
# for (i in trajectory$result_list[-1]) ordered_trajectory <- rbind(as.matrix(i),
#                                                                   as.matrix(ordered_trajectory))
# ordered_pseudotime <- ordered_pseudotime[order(as.numeric(rownames(ordered_pseudotime))),
# ]
# ordered_trajectory <- ordered_trajectory[order(as.numeric(rownames(ordered_trajectory))), 
# ]
# dimred <- as.data.frame(ordered_trajectory)
# cell_id = rownames(expression)[as.numeric(names(ordered_pseudotime))]
# dimnames(dimred) <- list(cell_id, c("comp_1", "comp_2") )


## progressions
# pst.full = lin_df$pseudot
# lineages = obj@trajectory$lines_list
# progressions <- map_df(1:length(lineages), function(l) {
# 
#   path_red = lineages[[l]]$s
#   x = rownames(expression)[as.numeric(rownames(path_red))]
# 
#   index = which(lin_df$cell_id %in% x)
#   pst = lin_df$pseudot[index] # extract pseudotime
# 
#   lin = unique(lin_df[index,"lineage"])
#   means <- sapply(lin, function(clID){
#     cluster = rep(0, length(pst.full))
#     cluster[which(lin_df$lineage==clID)] = 1
#     stats::weighted.mean(pst.full, cluster)
#   })
#   names(means) <- lin
# 
#   non_ends <- means[-c(1,length(means))]
#   edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
#   from.l <- lin[edgeID.l] %>% as.character
#   to.l <- lin[edgeID.l + 1] %>% as.character
#   m.from <- means[from.l]
#   m.to <- means[to.l]
# 
#   pct <- (pst - m.from) / (m.to - m.from)
#   pct[pct < 0] <- 0
#   pct[pct > 1] <- 1
# 
#   tibble(cell_id = x, from = from.l, to = to.l, percentage = pct)
# })