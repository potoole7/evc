`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}

#' @title Calculate adjacency matrix from Voronoi cells
#' @description Calculate adjacency matrix from Voronoi cells
#' @param pts Spatial points object
#' @param cut_vor Logical, whether to cut Voronoi cells to coastline
#' @param plot Logical, whether to plot Voronoi cells
#' @param areas Spatial polygons object
#' @return Adjacency matrix
#' @rdname calc_adj_mat
#' @export
calc_adj_mat <- \(pts, cut_vor = TRUE, plot = FALSE, areas = NULL) {
  # Calculate voronoi partition for sites
  vor <- pts |>
    sf::st_union() |>
    sf::st_voronoi(envelope = sf::st_as_sfc(sf::st_bbox(pts))) |>
    # sf::st_voronoi(envelope = sf::st_as_sfc(sf::st_bbox(areas))) |>
    sf::st_collection_extract(type = "POLYGON") |>
    sf::st_sf() |> # convert from geometry set to simple feature collection
    identity()

  # order voronoi cells to match points
  vor <- vor[order(sf::st_nearest_feature(sf::st_centroid(vor), pts)), ]

  # cutoff voronoi cells from ocean, if desired (stops far away neighbours)
  # TODO: Cut via coastline
  # TODO: Generalise more! Functionalise box below
  if (cut_vor == TRUE) {
    # slightly smaller than areas bbox
    vor <- sf::st_crop(
      vor, c("xmin" = -10, "ymin" = 51.6, "xmax" = -6, "ymax" = 55.2)
    )
  }

  # check that Voronoi cells have been produced correctly
  if (plot == TRUE && !is.null(areas)) {
    plot(sf::st_geometry(areas))
    plot(vor, add = TRUE)
    plot(pts, col = "blue", pch = 16, add = TRUE)
  }

  # calculate adjacency matrix from voronoi cells for stations
  return(spdep::nb2mat(spdep::poly2nb(vor), style = "B", zero.policy = TRUE))
}

#' @title Calculate within cluster sum of distances
#' @description Calculate within cluster sum of distances
#' @param k Number of clusters
#' @param distance_matrix Distance matrix
#' @param fun Clustering function
#' @return Total within-cluster sum of distances
#' @keywords internal
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

#' @title Scree plot
#' @description Produce scree plot for TWGSS. 
#' @param dist_mat Distance matrix
#' @param k Number of clusters
#' @param fun Clustering function
#' @param ... Additional arguments to clustering function
#' @return Total within-cluster sum of distances
#' @rdname scree_plot
#' @export
# TODO Could make this object method
# TODO Could expand this to work for AIC as well as TWGSS
scree_plot <- \(dist_mat, k = 1:10, fun = cluster::pam, ...) {
  total_within_ss <- vapply(
    k, within_cluster_sum, dist_mat, fun = fun, ..., FUN.VALUE = numeric(1)
  )

  # scree plot
  plot(
    k, total_within_ss, type = "b", pch = 19,
    xlab = "Number of Clusters",
    ylab = "Total Within-Cluster Sum of Distances",
    main = "Scree Plot"
  )
  return(total_within_ss)
}

#' @title Plot clustering solution on map
#' @description Plot clustering solution on map
#' @param pts Spatial points object
#' @param areas Spatial polygons object
#' @param clust_obj Clustering object
#' @return ggplot object
#' @rdname plt_clust_map
#' @export
# TODO: Could make this plot/ggplot method for object
plt_clust_map <- \(pts, areas, clust_obj) {

  name <- clust <- medoid <- NULL

  if (inherits(clust_obj, "kmeans")) {
    clust_element <- "cluster"
    medoid_locs <- NA
  } else if (inherits(clust_obj, "pam")) {
    clust_element <- "clustering"
    medoids <- clust_obj$medoids
    medoid_locs <- NA
    if (inherits(clust_obj$medoids, "character")) {
      medoid_locs <- medoids
    } else if (!is.null(rownames(medoids))) {
      medoid_locs <- rownames(medoids)
    }
  } else {
    stop("Clustering class not currently supported")
  }

  # reorder alphabetically
  # TODO: Look into this, required??
  # clust_names <- names(clust_obj[[clust_element]])
  # if (!is.null(clust_names)) {
  #   clust_obj[[clust_element]] <- clust_obj[[clust_element]][order(clust_names)]
  # }

  # TODO: Medoids doesn#t work as rows are named, fix!
  pts_plt <- cbind(pts, data.frame("clust" = clust_obj[[clust_element]])) |>
    dplyr::mutate(
      medoid = ifelse(name %in% medoid_locs, TRUE, FALSE), 
      medoid = factor(medoid, levels = c(FALSE, TRUE))
    )
  
  ggplot2::ggplot(areas) +
    ggplot2::geom_sf(colour = "black", fill = NA) +
    ggplot2::geom_sf(
      data = pts_plt,
      ggplot2::aes(
        colour = factor(clust), shape = medoid, size = as.numeric(medoid)
      ),
      alpha = 0.8
    ) +
    ggplot2::scale_shape_discrete(breaks = c(1, 15)) +
    ggplot2::scale_size_continuous(range = c(4, 8)) +
    ggplot2::guides(shape = "none", size = "none") +
    ggplot2::labs(colour = "Cluster") + 
    evc_theme() + 
    ggsci::scale_colour_nejm()
}

#' @title `ggpplot2` plotting theme
#' @description `ggpplot2` plotting theme
#' @return `ggpplot2` theme
#' @rdname evc_theme
#' @export
evc_theme <- \() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(size = 16, hjust = 0.5),
      axis.text        = ggplot2::element_text(size = 12),
      axis.title       = ggplot2::element_text(size = 14, face = "bold"),
      legend.text      = ggplot2::element_text(size = 12),
      strip.text       = ggplot2::element_text(size = 13, face = "bold"),
      strip.background = ggplot2::element_rect(fill = NA, colour = "black"),
      plot.tag         = ggplot2::element_text(size = 16, face = "bold"),
      panel.background = ggplot2::element_rect(fill = NA, colour = "black")
    )
}
