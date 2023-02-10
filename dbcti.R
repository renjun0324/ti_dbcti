
#tsneplot #################################
#input object, perplexity, max_iter, use_normalized_data...
#if specified gene = true, specified gene will be used for tsne plot
#' Tsne plot for cell_trajectory object
#'
#' @param object a cell_trajectory object
#' @param perplexity perplexity value for tsne
#' @param max_iter max iteration number
#' @param specified_gene feature names if specified
#' @param use_normalized_data if or not use normalized data as input
#' @param pca if pca will be performed
#' @param check_duplicates if or not duplicates are checked at tsne step
#' @param title title of plot, default as none
#' @param file file the plot to be written to, defalut as empity
#'
#' @return values in the tsne_data slot of cell_trajectory object
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @examples
#' tsneplot(sample_data, use_normalized_data = TRUE, perplexity = 5)
tsneplot2 <- function (object, 
                       perplexity = 30, 
                       max_iter = 500, 
                       specified_gene = FALSE, 
                       use_normalized_data = TRUE, 
                       pca = TRUE, 
                       check_duplicates = FALSE, 
                       title = "", 
                       file = "", 
                       initial_dims = 50) {
  # object = obj
  # perplexity = 30
  # max_iter = 500
  # specified_gene = FALSE
  # use_normalized_data = TRUE
  # pca = TRUE
  # check_duplicates = FALSE
  # title = ""
  # file = ""
  # initial_dims = 50
  
  set.seed(12345)
  if (specified_gene == TRUE) {
    if (use_normalized_data == FALSE) 
      data_feature <- scale(object@raw_data)[object@specified_gene, 
      ]
    else data_feature <- scale(object@normalized_data)[object@specified_gene, 
    ]
  } else {
    if (use_normalized_data == FALSE) 
      data_feature <- scale(object@raw_data)[object@selected_feature$selected_gene, 
      ]
    else data_feature <- scale(object@normalized_data)[object@selected_feature$selected_gene, 
    ]
  }
  tsne <- Rtsne::Rtsne(t(data_feature), dims = 2, perplexity = perplexity, 
                       verbose = FALSE, max_iter = max_iter, pca = pca, check_duplicates = check_duplicates, 
                       initial_dims = initial_dims)
  if (nchar(file) > 0) {
    pdf(file = file)
    plot(tsne$Y, main = title)
    dev.off()
  }
  plot(tsne$Y, main = title)
  class(tsne) = "list"
  object@tsne_data <- tsne
  return(object)
}


#fit principle curve iteratively
#input connection matrix(n*n) from above step, cluster_df collumn as x, y and cluster index, and corresponding cell name, iter_n
#i.e. test<-infer_trajectory(connection_matrix = connection_matrix, tsne_df = dat, cluster_index = cluster_index)
#' Title
#'
#' @param object a cell_trajectory object
#' @param iter_n iteration number
#'
#' @return values in the trajectory slot of cell_trajectory object
#' @export
#' @importFrom princurve principal_curve
#' @examples
#' infer_trajectory(sample_data, iter_n = 50)
infer_trajectory2 <- function (object, iter_n = 50) {
  
  connection_matrix = object@connect_cluster$cluster_connection
  tsne_df = as.data.frame(object@tsne_data$Y)
  cluster_index = object@distribution_estimation$cluster_index
  if (nrow(connection_matrix) < 2) 
    stop("not enough clusters")
  n_cluster = length(which(unique(cluster_index)!=0))
  lines_list <- result_list <- list()
  for (i in 1:n_cluster) {
    cluster_name_i <- rownames(connection_matrix)[i]
    result_list[[cluster_name_i]] <- tsne_df[cluster_index == 
                                               cluster_name_i, ]
  }
  isolated_cluster <- c()
  for (i in 1:n_cluster) {
    cat(i, "\n")
    cluster_name_i <- rownames(connection_matrix)[i]
    if (sum(connection_matrix[, i]) == 0) {
      fit_data <- as.matrix(result_list[[cluster_name_i]])
      fit_i <- princurve::principal_curve(fit_data)
      resulted_i <- as.data.frame(fit_i$s)
      colnames(resulted_i) <- c("x", "y")
      result_list[[cluster_name_i]] <- resulted_i
      lines_list_name <- paste("fit", cluster_name_i, sep = "_")
      lines_list[[lines_list_name]] <- fit_i
      isolated_cluster <- append(isolated_cluster, cluster_name_i)
    }
  }
  for (k in 1:iter_n) {
    for (i in 1:n_cluster) {
      cluster_name_i <- rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }
      assign(paste("fc", cluster_name_i, "x", sep = "_"), 
             data.frame(matrix(0, nrow = sum(cluster_index == 
                                               cluster_name_i), ncol = 1)))
      assign(paste("fc", cluster_name_i, "y", sep = "_"), 
             data.frame(matrix(0, nrow = sum(cluster_index == 
                                               cluster_name_i), ncol = 1)))
      assign(paste("cluster_index", cluster_name_i, sep = "_"), 
             names(cluster_index)[cluster_index == cluster_name_i])
    }
    for (i in 1:(nrow(connection_matrix) - 1)) {
      for (j in (i + 1):nrow(connection_matrix)) {
        if (connection_matrix[i, j] == 1) {
          cluster_name_i <- rownames(connection_matrix)[i]
          cluster_name_j <- rownames(connection_matrix)[j]
          fit_data <- as.matrix(rbind(result_list[[cluster_name_i]], 
                                      result_list[[cluster_name_j]]))
          assign(paste("fit", cluster_name_i, cluster_name_j, 
                       sep = "_"), princurve::principal_curve(fit_data))
          fitted <- eval(parse(text = paste("fit", cluster_name_i, 
                                            cluster_name_j, sep = "_")))$s
          fitted_i_x <- fitted[eval(parse(text = paste("cluster_index", 
                                                       cluster_name_i, sep = "_"))), ][, 1]
          fitted_i_y <- fitted[eval(parse(text = paste("cluster_index", 
                                                       cluster_name_i, sep = "_"))), ][, 2]
          fitted_j_x <- fitted[eval(parse(text = paste("cluster_index", 
                                                       cluster_name_j, sep = "_"))), ][, 1]
          fitted_j_y <- fitted[eval(parse(text = paste("cluster_index", 
                                                       cluster_name_j, sep = "_"))), ][, 2]
          assign(paste("fc", cluster_name_i, "x", sep = "_"), 
                 cbind(eval(parse(text = paste("fc", cluster_name_i, 
                                               "x", sep = "_"))), fitted_i_x))
          assign(paste("fc", cluster_name_i, "y", sep = "_"), 
                 cbind(eval(parse(text = paste("fc", cluster_name_i, 
                                               "y", sep = "_"))), fitted_i_y))
          assign(paste("fc", cluster_name_j, "x", sep = "_"), 
                 cbind(eval(parse(text = paste("fc", cluster_name_j, 
                                               "x", sep = "_"))), fitted_j_x))
          assign(paste("fc", cluster_name_j, "y", sep = "_"), 
                 cbind(eval(parse(text = paste("fc", cluster_name_j, 
                                               "y", sep = "_"))), fitted_j_y))
        }
      }
    }
    for (i in 1:n_cluster) {
      cluster_name_i <- rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }
      x_df <- eval(parse(text = paste("fc", cluster_name_i, 
                                      "x", sep = "_")))
      y_df <- eval(parse(text = paste("fc", cluster_name_i, 
                                      "y", sep = "_")))
      assign(paste("fc", cluster_name_i, "x", sep = "_"), 
             as.data.frame(x_df[, -1], row.names = rownames(x_df)))
      assign(paste("fc", cluster_name_i, "y", sep = "_"), 
             as.data.frame(y_df[, -1], row.names = rownames(y_df)))
      data_x <- eval(parse(text = paste("fc", cluster_name_i, 
                                        "x", sep = "_")))
      data_y <- eval(parse(text = paste("fc", cluster_name_i, 
                                        "y", sep = "_")))
      assign(paste("loc", cluster_name_i, "x", sep = "_"), 
             apply(as.data.frame(data_x), 1, mean))
      assign(paste("loc", cluster_name_i, "y", sep = "_"), 
             apply(as.data.frame(data_y), 1, mean))
      assign(paste("loc", cluster_name_i, sep = "_"), data.frame(x = eval(parse(text = paste("loc", 
                                                                                             cluster_name_i, "x", sep = "_"))), y = eval(parse(text = paste("loc", 
                                                                                                                                                            cluster_name_i, "y", sep = "_")))))
      result_list[[cluster_name_i]] <- eval(parse(text = paste("loc", 
                                                               cluster_name_i, sep = "_")))
    }
  }
  for (i in 1:(nrow(connection_matrix) - 1)) {
    for (j in (i + 1):nrow(connection_matrix)) {
      cluster_name_i <- rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }
      cluster_name_j <- rownames(connection_matrix)[j]
      if (connection_matrix[i, j] == 1) {
        fitted <- eval(parse(text = paste("fit", cluster_name_i, 
                                          cluster_name_j, sep = "_")))
        lines_list_name <- paste("fit", cluster_name_i, 
                                 cluster_name_j, sep = "_")
        lines_list[[lines_list_name]] <- fitted
      }
    }
  }
  result <- list(result_list = result_list, lines_list = lines_list)
  object@trajectory <- result
  return(object)
}


#connection_matrix, trajectory, start_state_name(character), separate==FALSE, cluster_index
#' calculate pseudotime for each cell
#'
#' @param object a cell_trajectory object
#' @param start_state_name index for the state that starts
#'
#' @return values in the pseudotime slot of cell_trajectory object
#' @export
#' @import igraph
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @examples
#' calculate_pseudotime(sample_data, start_state_name = c('1','2'))
calculate_pseudotime2 <- function(object, start_state_name){
  connection_matrix = object@connect_cluster$cluster_connection
  trajectory = object@trajectory
  cluster_index = object@distribution_estimation$cluster_index
  
  #build graph
  g<-igraph::graph_from_adjacency_matrix(connection_matrix, mode = 'undirected')
  #if separated
  if (igraph::components(g)$no == 1) separated = FALSE else separated = TRUE
  #build cluster index for each cluster
  for(i in 1:nrow(connection_matrix)){
    cluster_name_i<-rownames(connection_matrix)[i]
    assign(paste('cluster_index', cluster_name_i, sep = '_'), names(cluster_index)[cluster_index==cluster_name_i])
  }
  #separated == false ###################
  if (separated==FALSE) {
    #calculate distance of each cluster i
    for(i in 1:nrow(connection_matrix)){
      d_list<-list()
      cluster_name_i<-rownames(connection_matrix)[i]
      fit_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
      part<-names(trajectory$lines_list)[fit_index]
      #for each part j of cluster i calculate the mean of j
      for (j in 1:length(part)) {
        part_index_j<-part[j]
        tra_lam<-trajectory$lines_list[[part_index_j]]$lambda
        d_part_j<-tra_lam[names(tra_lam) %in% eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))]
        d_range_j<-max(d_part_j)-min(d_part_j)
        d_list<-append(d_list, d_range_j)
      }
      assign(paste('d', cluster_name_i, sep = '_'), mean(unlist(d_list)))
    }
    
    #assign psdeotime
    psdeotime_list<-list()
    for(i in 1:nrow(connection_matrix)){
      cat(i, "\n")
      cluster_name_i<-rownames(connection_matrix)[i]
      path<-names(unlist(igraph::shortest_paths(g, cluster_name_i, start_state_name)[[1]]))
      
      
      cl_point<-eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))
      #find the start point of the cluster
      #situation where point is not in the cluster of start
      if (length(path)>=2) {
        if (as.numeric(path[1]) < as.numeric(path[2])) curve_name<-paste('fit', path[1], path[2], sep = '_') else curve_name<-paste('fit', path[2], path[1], sep = '_')
        
        ordered<-trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ][rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ]) %in% cl_point, ]
        name_first<-rownames(ordered)[1]
        name_last<-rownames(ordered)[nrow(ordered)]
        
        match_first<-match(name_first, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        match_last<-match(name_last, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', path[2], sep = '_'))), rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,])))
        if (avg_other_cl<=match_first) start_point<-name_first else start_point<-name_last
        
        #distance
        d_within_cluster<-abs(trajectory$lines_list[[curve_name]]$lambda[start_point]-trajectory$lines_list[[curve_name]]$lambda[cl_point])
        d_total<-d_within_cluster
        for (j in path[-1]) d_total <- d_total + eval(parse(text = paste('d', j, sep = '_')))
        
      }
      
      #situation where cluster is state 0
      if (length(path)==1){
        part_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
        part<-names(trajectory$lines_list)[part_index]
        path_1_df<-data.frame(matrix(nrow = length(cl_point)), row.names = cl_point)
        n=0
        for (j in part) {
          n=n+1
          ordered<-trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ][rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ]) %in% cl_point, ]
          name_first<-rownames(ordered)[1]
          name_last<-rownames(ordered)[nrow(ordered)]
          
          match_first<-match(name_first, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))
          match_last<-match(name_last, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))
          
          #other cluster
          if (cluster_name_i == unlist(stringr::str_split(j, '_'))[2]) other_cluster<-unlist(stringr::str_split(j, '_'))[3] else other_cluster<-unlist(stringr::str_split(j, '_'))[2]
          
          
          avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', other_cluster, sep = '_'))), rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,])))
          if (avg_other_cl<=match_first) {
            start_point<-name_first
            
            #all distance
            d_within_cluster<-eval(parse(text = paste('d', cluster_name_i, sep = '_')))-abs(trajectory$lines_list[[j]]$lambda[start_point]-trajectory$lines_list[[j]]$lambda[cl_point])
            
            path_1_df[,n]<-d_within_cluster
          } else {
            start_point<-name_first
            
            #all distance
            d_within_cluster<-trajectory$lines_list[[j]]$lambda[cl_point]
            
            path_1_df[,n]<-d_within_cluster
          }
        }
        
        
        #average distance
        d_total<-apply(path_1_df, 1, mean)
      }
      psdeotime_list[[cluster_name_i]]<-d_total
    }
  }
  #separated == true ####################
  if (separated == TRUE) {
    #if the start_state_name meet requirement
    g_comp<-igraph::components(g)
    print(g_comp)
    if (length(start_state_name)<=1) stop('start_state_name for separated graph should be a at least length of 2 vertex')
    
    
    start_cluster_list<-c()
    for (i in start_state_name) {
      separated_comp<-g_comp$membership[match(i, names(igraph::V(g)))]
      start_cluster_list<-append(start_cluster_list, separated_comp)
    }
    if (sum(duplicated(start_cluster_list)) >= 1) stop('multiple start state name from same graph part')
    
    
    #calculate distance of each cluster i
    for(i in 1:nrow(connection_matrix)){
      d_list<-list()
      cluster_name_i<-rownames(connection_matrix)[i]
      fit_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
      part<-names(trajectory$lines_list)[fit_index]
      #for each part j of cluster i calculate the mean of j
      for (j in 1:length(part)) {
        part_index_j<-part[j]
        tra_lam<-trajectory$lines_list[[part_index_j]]$lambda
        d_part_j<-tra_lam[names(tra_lam) %in% eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))]
        d_range_j<-max(d_part_j)-min(d_part_j)
        d_list<-append(d_list, d_range_j)
      }
      assign(paste('d', cluster_name_i, sep = '_'), mean(unlist(d_list)))
    }
    
    #isolated part situation##########################
    psdeotime_list<-list()
    isolated_cluster<-c()
    for (i in 1:nrow(connection_matrix)) {
      cluster_name_i<-rownames(connection_matrix)[i]
      if (sum(connection_matrix[, i]) == 0) {
        psdeotime_list[[cluster_name_i]]<-trajectory$lines_list[[paste('fit', cluster_name_i, sep = '_')]]$lambda
        isolated_cluster<-append(isolated_cluster, cluster_name_i)
      }
    }
    #
    #assign psdeotime
    
    for(i in 1:nrow(connection_matrix)){
      cluster_name_i<-rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }
      #build cluster and their start state relationship
      separated_comp_i<-g_comp$membership[match(cluster_name_i, names(igraph::V(g)))]
      for (k in start_state_name) {
        separated_comp_k_start<-g_comp$membership[match(k, names(igraph::V(g)))]
        if (separated_comp_k_start == separated_comp_i) break
      }
      
      path<-names(unlist(igraph::shortest_paths(g, cluster_name_i, names(separated_comp_k_start))[[1]]))
      cl_point<-eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))
      #line################
      #find the start point of the cluster
      #situation where point is not in the cluster of start
      if (length(path)>=2) {
        if (as.numeric(path[1]) < as.numeric(path[2])) curve_name<-paste('fit', path[1], path[2], sep = '_') else curve_name<-paste('fit', path[2], path[1], sep = '_')
        
        ordered<-trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ][rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ]) %in% cl_point, ]
        name_first<-rownames(ordered)[1]
        name_last<-rownames(ordered)[nrow(ordered)]
        
        match_first<-match(name_first, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        match_last<-match(name_last, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', path[2], sep = '_'))), rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,])))
        if (avg_other_cl<=match_first) start_point<-name_first else start_point<-name_last
        
        #distance
        d_within_cluster<-abs(trajectory$lines_list[[curve_name]]$lambda[start_point]-trajectory$lines_list[[curve_name]]$lambda[cl_point])
        d_total<-d_within_cluster
        for (j in path[-1]) d_total <- d_total + eval(parse(text = paste('d', j, sep = '_')))
        
      }
      
      #situation where cluster is state 0
      if (length(path)==1){
        part_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
        part<-names(trajectory$lines_list)[part_index]
        path_1_df<-data.frame(matrix(nrow = length(cl_point)), row.names = cl_point)
        n=0
        for (j in part) {
          n=n+1
          ordered<-trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ][rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ]) %in% cl_point, ]
          name_first<-rownames(ordered)[1]
          name_last<-rownames(ordered)[nrow(ordered)]
          
          match_first<-match(name_first, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))
          match_last<-match(name_last, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))
          
          #other cluster
          if (cluster_name_i == unlist(stringr::str_split(j, '_'))[2]) other_cluster<-unlist(stringr::str_split(j, '_'))[3] else other_cluster<-unlist(stringr::str_split(j, '_'))[2]
          
          
          avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', other_cluster, sep = '_'))), rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,])))
          if (avg_other_cl<=match_first) {
            start_point<-name_first
            
            #all distance
            d_within_cluster<-eval(parse(text = paste('d', cluster_name_i, sep = '_')))-abs(trajectory$lines_list[[j]]$lambda[start_point]-trajectory$lines_list[[j]]$lambda[cl_point])
            
            path_1_df[,n]<-d_within_cluster
          } else {
            start_point<-name_first
            
            #all distance
            d_within_cluster<-trajectory$lines_list[[j]]$lambda[cl_point]
            
            path_1_df[,n]<-d_within_cluster
          }
        }
        #average distance
        d_total<-apply(path_1_df, 1, mean)
      }
      psdeotime_list[[cluster_name_i]]<-d_total
    }
  }
  object@pseudotime <- psdeotime_list
  return(object)
}