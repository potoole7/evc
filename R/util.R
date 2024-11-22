`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}

# calculate adjacency matrix from Voronoi cells
#' @title Calculate adjacency matrix from Voronoi cells
#' @description Calculate adjacency matrix from Voronoi cells
#' @param pts Spatial points object
#' @param cut_vor Logical, whether to cut Voronoi cells to coastline
#' @param plot Logical, whether to plot Voronoi cells
#' @return Adjacency matrix
#' @export
calc_adj_mat <- \(pts, cut_vor = TRUE, plot = FALSE) {
  # Calculate voronoi partition for sites
  vor <- pts %>% 
    st_union() %>%
    st_voronoi(envelope = st_as_sfc(st_bbox(pts))) %>%
    # st_voronoi(envelope = st_as_sfc(st_bbox(areas))) %>%
    st_collection_extract(type = "POLYGON") %>% 
    st_sf() %>%  # convert from geometry set to simple feature collection
    identity()
  
  # order voronoi cells to match points
  vor <- vor[order(st_nearest_feature(st_centroid(vor), pts)), ]
  
  # cutoff voronoi cells from ocean, if desired (stops far away neighbours)
  # TODO: Cut via coastline
  if (cut_vor == TRUE) {
    # slightly smaller than areas bbox
    vor <- st_crop(
      vor, c("xmin" = -10, "ymin" = 51.6, "xmax" = -6, "ymax" = 55.2)
    )
  }
   
  # check that voronoi cells have been produced correctly
  if (plot == TRUE) {
    plot(st_geometry(areas)) 
    plot(vor, add = TRUE)
    plot(pts, col = "blue", pch = 16, add = TRUE)
    # test that voronoi cells and points ordered correctly
    # plot(vor[10, ], col = "blue", fill = NA, lwd = 5, add = TRUE)
    # plot(pts[10, ], col = "red", pch = 16, add = TRUE)
  }
  
  
  # calculate adjacency matrix from voronoi cells for stations
  return(spdep::nb2mat(spdep::poly2nb(vor), style = "B", zero.policy = TRUE))
}

# compute the total within-cluster sum of distances
# TODO: Create methods for the below functions to differ for PAM vs k-means
within_cluster_sum <- function(k, distance_matrix, fun = cluster::pam, ...) {
  clust_res <- fun(distance_matrix, k, ...)
  if (inherits(clust_res, "kmeans")) {
    return(clust_res$tot.withinss)
  } else if (inherits(clust_res, "pam")) {
    return(clust_res$objective[1])
  } else {
    stop("Clustering class not currently supported")
  }
}

# compute and produce scree plot
scree_plot <- \(dist_mat, k = 1:10, fun = cluster::pam, ...) {
  total_within_ss <- vapply(
    k, within_cluster_sum, dist_mat, fun = fun, ..., FUN.VALUE = numeric(1)
  )
  
  # scree plot
  plot(
    k, total_within_ss, type = "b", pch = 19 # , 
    # xlab = "Number of Clusters (k)", 
    # ylab = "Total Within-Cluster Sum of Distances",
    # main = "Scree Plot for K-medoids Clustering"
  )
  return(total_within_ss)
}

# plot clustering solution on map
plt_clust <- \(pts, clust_obj) {
  
  if (inherits(clust_obj, "kmeans")) {
    clust_element <- "cluster"  
    medoids <- NA
  } else if (inherits(clust_obj, "pam")) {
    clust_element <- "clustering"  
    medoids <- clust_obj$medoids
    if (inherits(medoids, "matrix")) medoids <- as.numeric(rownames(medoids))
  } else {
    stop("Clustering class not currently supported")
  }
  
  pts_plt <- cbind(pts, data.frame("clust" = clust_obj[[clust_element]])) %>% 
    mutate(row = row_number()) %>% 
    mutate(mediod = ifelse(row %in% medoids, TRUE, FALSE))
  
  ggplot(areas) + 
    geom_sf(colour = "black", fill = NA) + 
    geom_sf(
      data = pts_plt, 
      aes(colour = factor(clust), shape = mediod, size = as.numeric(mediod)), 
      alpha = 0.8
    ) + 
    scale_shape_discrete(breaks = c(1, 15)) + 
    scale_size_continuous(range = c(3, 6)) +  
    guides(shape = "none", size = "none") + 
    theme
}

evc_theme <- \() {
   theme_bw() +
     theme(
       legend.position = "bottom",
       plot.title = element_text(size = 16, hjust = 0.5),
       axis.text = element_text(size = 12),
       axis.title = element_text(size = 14, face = "bold"),
       legend.text = element_text(size = 12),
       strip.text = element_text(size = 13, face = "bold"),
       strip.background = element_rect(fill = NA, colour = "black"),
       plot.tag = element_text(size = 16, face = "bold"),
       panel.background = element_rect(fill = NA, colour = "black")
     )
}
