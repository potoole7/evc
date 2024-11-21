# Function to calculate and cluster on the Jensen-Shannon Divergence for 
# the conditional extremes model
js_clust <- \(
  dependence,
  nclust = 3,
  cluster_mem = NULL,
  dat_max_mult = 2,
  n_dat = 10
) {
  
  # TODO: Add stopifnot clause for class of dependence
   
  # pull parameter values for each location
  params <- lapply(dependence, pull_params)
  
  # pull Laplace threshold values for each location
  # TODO: Problem, last threshold different to others, investigate
  thresh <- lapply(dependence, pull_thresh_trans)
  
  # take maximum Laplace thresholds; want to geenra
  # thresh_max <- apply(bind_rows(thresh), 2, max)
  thresh_max <- lapply(bind_rows(thresh), max)
  
  # list of locs -> list of len 2 of variables, each containing all locs
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
  # scree plots look to suggest 3 clusters for both
  # lapply(dist_mats, scree_plot) 
  # lapply(dist_mats, scree_plot, fun = kmeans) 
  
  # cluster for rain and wind speed using both k-means and PAM
  # js_clust <- lapply(dist_mats, \(x) {
  #   lapply(c(pam, kmeans), \(fun) {
  #     fun(x, 3)
  #   })
  # })
  # cluster for rain and wind speed using PAM
  pam_js_clust <- cluster::pam(dist_mat, k = nclust)
  ret <- list("pam" = pam_js_clust)
  # evaluate quality
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      pam_js_clust$clustering, 
      cluster_mem
    )
    # print(paste0("adjusted Rand index" = adj_rand))
    ret <- c(ret, list("adj_rand" = adj_rand))
  }
  
  return(ret)
}


# pull a, b, mu, sigma for a given location
# dep - List of mex objects (i.e. CE models for each var) for single location
pull_params <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  # stopifnot(all(vapply(dep, class, character(1)) == "mex"))
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

# pull thresholds
# dep - List of mex objects (i.e. CE models for each var) for single location
pull_thresh_trans <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  stopifnot(all(vapply(dep, class, character(1)) == "mex"))
  # return quantile of transformed data (already calculated in texmex)
  return(lapply(dep, \(x) x$dependence$dth))
}

# function to calculate KL divergence for single data point
# derivation at https://statproofbook.github.io/P/norm-kl.html
kl_gauss <- \(mu1, mu2, var1, var2) {
  return(1 / 2 * (
    (mu2 - mu1)^2 / var2 + (var1 / var2) - log(var1 / var2) - 1
  ))
}

# function to calculate Jensen-Shannon divergence metric from KL divergence
# This metric is symmetric, as required for clustering 
# https://tinyurl.com/526rwy9f
js_gauss <- \(mu1, mu2, var1, var2) {
  # calculate mean, variance for mixture distribution M = (P + Q)/2
  # Sum of normals as in https://online.stat.psu.edu/stat414/lesson/26/26.1
  # (N(u1, v1) + N(u2, v2)) / 2 = N((u1 + u2) / 2, (v1 + v2) / 2^2)
  mu_m <- (mu1 + mu2) / 2 
  var_m <- (var1 + var2) / 4
  
  # calculate JS(P||Q) = ((KL(P||M)) + KL(Q||M))/2
  # TODO: Check that this is right, shouldn't it be var2 for second kl_gauss?
  # TODO: Unit test this is symmetric, previously wasn't (?)
  return(
    # (kl_gauss(mu1, mu_m, var1, var_m) + kl_gauss(mu2, mu_m, var1, var_m)) / 2
    (kl_gauss(mu1, mu_m, var1, var_m) + kl_gauss(mu2, mu_m, var2, var_m)) / 2
  )
}

# Function to calculate Jensen-Shannon divergence for each data point
js_div <- \(params_x, params_y, thresh_max, data_max = 2 * thresh_max, n_dat) {
  
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


