# Function to calculate and cluster on the Jensen-Shannon Divergence for
# the conditional extremes model
#' @title Cluster based on Jensen-Shannon divergence for conditional extremes
#' model
#' @description Function to calculate and cluster on the Jensen-Shannon
#' Divergence for the conditional extremes model.
#' @param dependence List of `mexDependence` objects for each location.
#' @param k Number of clusters to fit, set to NULL to produce scree plot and
#' return distance matrix.
#' @param dist_mat Distance matrix, set to NULL to calculate.
#' @param cluster_mem Optional known cluster memberships for each location,
#' which can be used to evaluate the quality of the clustering.
#' @param dat_max_mult Multiplier for the maximum multiple of the largest
#' threshold across all locations.
#' @param n_dat Number of data points to use in the Jensen-Shannon divergence.
#' @param scree_k Vector of k values to use in scree plot, if `k` is NULL.
#' @return List containing the clustering results and, if `cluster_mem` is
#' provided, the adjusted Rand index.
#' @rdname js_clust
#' @export
# TODO: Document how globally largest threshold is used as data
# TODO: Document how pam is also used here
js_clust <- \(
  dependence,
  k = NULL,
  dist_mat = NULL,
  cluster_mem = NULL,
  dat_max_mult = 2,
  n_dat = 10,
  scree_k = 1:5
) {

  # TODO: Add stopifnot clause for class of dependence

  # pull parameter values for each location
  if (is.null(dist_mat)) {
    params <- lapply(dependence, pull_params)

    # pull Laplace threshold values for each location
    # TODO: Problem, last threshold different to others, investigate
    thresh <- lapply(dependence, pull_thresh_trans)

    # take maximum Laplace thresholds; want to geenra
    thresh_max <- lapply(dplyr::bind_rows(thresh), max)

    # list of locs containing vars -> list of vars, each containing all locs
    params <- purrr::transpose(params)

    dist_mats <- lapply(seq_along(params), \(i) {
      proxy::dist(
        params[[i]],
        method = js_div,
        thresh_max = thresh_max[[i]],
        data_max = dat_max_mult * thresh_max[[i]],
        n_dat = n_dat
      )
    })

    # sum distance matrices over different variables together
    dist_mat <- do.call(`+`, dist_mats)
  }

  if (is.null(k)) {
    return(list(
      "dist_mat"        = dist_mat,
      "total_within_ss" = scree_plot(dist_mat, k = scree_k)
    ))
  }

  # cluster for rain and wind speed using PAM
  # TODO: functionalise this, exact same as in kl_sim_eval!
  pam_js_clust <- cluster::pam(dist_mat, k = k)
  ret <- list("pam" = pam_js_clust)
  # evaluate quality
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      pam_js_clust$clustering,
      cluster_mem
    )
    ret <- c(ret, list("adj_rand" = adj_rand))
  }

  return(ret)

  pam_js_clust <- cluster::pam(dist_mat, k = k)
  ret <- pam_js_clust
  # evaluate quality
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      pam_js_clust$clustering,
      cluster_mem
    )
    ret <- list("pam" = ret, "adj_rand" = adj_rand)
  }
  return(ret)
}


#' @title Pull parameters for conditional extremes model
#' @description Function to pull the parameters (a, b, m and s) for the
#' conditional extremes model.
#' @param dep List of `mexDependence` objects for each location.
#' @return List of named vectors containing the parameters for each location.
#' @keywords internal
pull_params <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  # loop through conditioning variables for single location
  return(lapply(dep, \(x) {
    # pull parameter vals
    ret <- as.vector(x$dependence$coefficients)
    names(ret) <- rownames(x$dependence$coefficients)
    # remove c and d if 0 (i.e. Laplace margins, rather than Gumbel)
    if (ret["c"] == 0 && ret["d"] == 0) {
      ret <- ret[!names(ret) %in% c("c", "d")]
    }
    return(ret)
  }))
}

#' @title Pull thresholds for conditional extremes model
#' @description Function to pull the thresholds for the conditional extremes
#' model.
#' @param dep List of `mexDependence` objects for each location.
#' @return List of named vectors containing the thresholds for each location.
#' @keywords internal
pull_thresh_trans <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  stopifnot(all(vapply(dep, class, character(1)) == "mex"))
  # return quantile of transformed data (already calculated in texmex)
  return(lapply(dep, \(x) x$dependence$dth))
}

#' @title Calculate KL divergence between two Gaussian distributions
#' @description Function to calculate the Kullback-Leibler divergence between
#' two Gaussian distributions, using the closed form solution for the case of
#' two univariate Gaussians provided at
#' \href{https://statproofbook.github.io/P/norm-kl.html}{test}.
#' @param mu1 Mean of the first Gaussian distribution.
#' @param mu2 Mean of the second Gaussian distribution.
#' @param var1 Variance of the first Gaussian distribution.
#' @param var2 Variance of the second Gaussian distribution.
#' @return The Kullback-Leibler divergence between the two Gaussian
#' distributions.
#' @keywords internal
# TODO: Have mus and variances as lists perhaps?
kl_gauss <- \(mu1, mu2, var1, var2) {
  return(1 / 2 * (
    (mu2 - mu1)^2 / var2 + (var1 / var2) - log(var1 / var2) - 1
  ))
}

#' @title Calculate Jensen-Shannon divergence between two Gaussian distributions
#' @description Function to calculate the Jensen-Shannon divergence between two
#' Gaussian distributions, using the Kullback-Leibler divergence between them.
#' This metric is symmetric, as required for clustering.
#' @param mu1 Mean of the first Gaussian distribution.
#' @param mu2 Mean of the second Gaussian distribution.
#' @param var1 Variance of the first Gaussian distribution.
#' @param var2 Variance of the second Gaussian distribution.
#' @return The Jensen-Shannon divergence between the two Gaussian
#' distributions.
#' @keywords internal
js_gauss <- \(mu1, mu2, var1, var2) {
  # calculate mean, variance for mixture distribution M = (P + Q)/2
  # Sum of normals as in https://online.stat.psu.edu/stat414/lesson/26/26.1
  mu_m <- (mu1 + mu2) / 2
  var_m <- (var1 + var2) / 4

  # calculate JS(P||Q) = ((KL(P||M)) + KL(Q||M))/2
  # TODO: Check that this is right, shouldn't it be var2 for second kl_gauss?
  # TODO: Unit test this is symmetric, previously wasn't (?)
  return(
    (kl_gauss(mu1, mu_m, var1, var_m) + kl_gauss(mu2, mu_m, var2, var_m)) / 2
  )
}

#' @title Calculate Jensen-Shannon divergence for each data point
#' @description Function to calculate Jensen-Shannon divergence for each data
#' point.
#' @param params_x Parameters for the first Gaussian distribution.
#' @param params_y Parameters for the second Gaussian distribution.
#' @param thresh_max Maximum threshold value across all locations.
#' @param data_max Maximum data value.
#' @param n_dat Number of data points to use in the Jensen-Shannon divergence.
#' return The Jensen-Shannon divergence for each data point.
#' @keywords internal
js_div <- \(
  params_x,
  params_y,
  thresh_max,
  data_max = 2 * thresh_max,
  n_dat = 10
) {

  # test that input vectors have correct conditional extremes parameters
  stopifnot(is.vector(params_x) && is.vector(params_y))
  stopifnot(
    all(c(names(params_x), names(params_y)) == rep(c("a", "b", "m", "s"), 2))
  )

  # create data sequence from specified arguments
  data <- seq(thresh_max, data_max, length = n_dat)

  # funs to calculate mu and sd for normal dist as in 5.2 of Heff & Tawn '04
  mu_fun <- \(x, data) {
    return(x[["a"]] * data + x[["m"]] * (data ^ x[["b"]]))
  }
  var_fun <- \(x, data) {
    # take square now to avoid in kl_single formula
    return((x[["s"]] * (data ^ x[["b"]]))^2)
  }

  # calculate mu and sigma for each data point
  mus <- lapply(list(params_x, params_y), mu_fun, data = data)
  vars <- lapply(list(params_x, params_y), var_fun, data = data)

  # Calculate Jensen-Shannon divergence for each data point
  # TODO: How best to summarise across all data points? Sum? Average?
  return(sum(mapply(
    js_gauss,
    mu1 = mus[[1]], mu2 = mus[[2]], var1 = vars[[1]], var2 = vars[[2]]
  )))
}
