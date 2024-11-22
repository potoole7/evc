# TODO: Generate beyond wind speed and rain
#' @title PAM clustering based on empirical KL divergence
#' @description Perform clustering based on empirical KL divergence between
#' bivariate extremes for two locations.
#' @param data Data matrix with rows representing locations and columns
#' @param kl_prob Extremal quantile to use.
#' @param k Number of clusters.
#' @param cluster_mem Optional true cluster membership to calculate adjusted
#' Rand index.
#' @param ... Additional arguments to pass to `proxy::dist`.
#' @return List with clustering solution and adjusted Rand index if cluster_mem
#' is provided.
#' @rdname kl_sim_eval
#' @export
# Function to calculate KL divergence for 2 cluster simulation dataset
kl_sim_eval <- \(
  data, 
  kl_prob,  
  k = 2,   
  cluster_mem = NULL, 
  ...
) {
  
  # prob must be valid probability
  stopifnot(0 <= kl_prob && kl_prob <= 1)
  
  # KL divergence between areas using Vignotto 2021 method
  kl_mat <- proxy::dist(
    data, method = emp_kl_div, print = FALSE, prob = kl_prob, ...
  )
  
  # clustering solution
  pam_kl_clust <- cluster::pam(kl_mat, k = k)
  ret <- list("pam" = pam_kl_clust)
  # evaluate quality
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      pam_kl_clust$clustering, 
      cluster_mem
    )  
    ret <- c(ret, list("adj_rand" = adj_rand))
  }
  return(ret)
}

#' @title Empirical KL divergence between two locations
#' @description Functions to calculate empirical KL divergence between two 
#' locations using the Vignotto 2021 method.
#' @param x Vector of data representing location 1
#' @param y Vector of data representing location 2
#' @param prob Extremal quantile to use.
#' @return KL divergence between x and y.
#' @keywords internal
# function to calculate KL divergence between any two locations
emp_kl_div <- \(x, y, prob = 0.9) {
  
  # split x and y in half (rain vs wind speed)
  # remove NAs from x and y from padding matrix
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n <- length(x)
  m <- length(y) # need for m???

  # TODO: Make more general for multivariate (not bivariate) case
  df_lst <- list(
  # x_df <- data.frame(
    "x" = data.frame(
      "rain"       = x[1:(n / 2)],
      "wind_speed" = x[((n / 2) + 1):n]
    ),
    "y" = data.frame(
    "rain"       = y[1:(m / 2)],
    "wind_speed" = y[((m / 2) + 1):m]
    )
  )
  
  # convert to Pareto scale
  df_lst_par <- df_lst
  df_lst_par <- lapply(df_lst, \(z) {
    data.frame(apply(z, 2, pareto_trans))
  })
  
  # partition into 3 subsets
  df_part <- lapply(df_lst_par, \(z) {
    partition_max(list(z[, 1], z[, 2]), prob = prob)
  })
  
  # fail if there are any 0s
  stopifnot(
    "No set can contain zeros, reduce `prob`" = all(unlist(df_part) != 0)
  )
  
  # calculate proportions of partitions
  df_part_prop <- lapply(df_part, \(z) {
    denom <- sum(unlist(z)) # denominator is # extreme obs
    lapply(z, `/`, denom) # find prop of extreme obs for each disjoint set
  })
  
  # calculate proportions of partitions
  x_part <- df_part_prop[[1]]
  y_part <- df_part_prop[[2]]
  if (print) {
    print(paste0("both for x: ", round(x_part$both, 3)))
    print(paste0("both for y: ", round(y_part$both, 3)))
  }
  sum_vals <- vector(length = length(x_part))
  for (i in seq_along(sum_vals)) { 
    sum_vals[i] <- (x_part[[i]] - y_part[[i]]) * 
                      log(x_part[[i]] / y_part[[i]])
  }
  if (any(is.nan(sum_vals) | is.infinite(sum_vals), na.rm = TRUE)) return(NA)
  return((1  / 2) * sum(sum_vals))
}


#' @title Pareto transformation
#' @description Function to convert a vector to Pareto scale.
#' @param x Vector to convert.
#' @return Vector in Pareto scale.
#' @keywords internal
pareto_trans <- \(x) {
  stopifnot(is.vector(x))
  
  # order and sort data
  x_ord <- order(x)
  x_sort <- x[x_ord]
  
  # calculate ECDF
  n <- length(x)
  ecdf_vals <- (seq_len(n)) / (n + 1)
  # convert back to original order
  ecdf_vals_x_ord <- numeric(n)
  ecdf_vals_x_ord[x_ord] <- ecdf_vals
  
  # pareto transform 
  return(1 / (1 - ecdf_vals_x_ord))
}

#' @title Bivariate risk function
#' @description Function to calculate bivariate risk function from Vignotto
#' 2021.
#' @param x List of vectors representing different locations.
#' @param fun Function to calculate "risk", defaults to max as easier to
#' partition. Note that `sum` has not been implemented yet.
#' @return Vector of risk values.
#' @keywords internal
risk_fun <- \(x, fun = max) {
  # test list
  stopifnot(is.list(x)) 
  # test equal length
  stopifnot(length(unique(vapply(x, length, numeric(1)))) == 1) 
  
  # TODO: Improve/speedup (could have x as a dataframe!)
  risk_vec <- vector(length = length(x[[1]]))
  for (i in seq_along(x[[1]])) {  
    risk_vec[[i]] <- fun(
      # pull ith entry in each vector in the list x
      vapply(x, `[[`, i, FUN.VALUE = numeric(1))
    )
  }
  return(risk_vec)
}

#' @title Partition into 3 subsets
#' @description Function to partition bivariate extremes into into 3 subsets
#' based on the risk function from Vignotto 2021.
#' @param x List of vectors representing different locations.
#' @param prob Extremal quantile to use.
#' @return List with counts of points in each subset. 
#' @keywords internal
# TODO: Only works for max risk fun, unsure how to do for sum...
partition_max <- \(x, prob, plot = FALSE) {
  
  # x must be a list of length 2 (rain and wind speed)
  stopifnot(is.list(x) && length(x) == 2) 
  
  # convert to dataframe, easier to subset (must follow this order!)
  df <- data.frame(
    rain       = x[[1]], 
    wind_speed = x[[2]], # TODO: Make x a dataframe already
    # TODO: Fix this, wrong somehow!!!
    R          = risk_fun(x, max)
  ) 
  
  # calculate quantile of risk function
  qu <- stats::quantile(df$R, prob = prob)[[1]]
  
  # partition into 3 subsets
  df <- df |> 
    dplyr::mutate(extreme = dplyr::case_when(
      rain > qu & wind_speed > qu   ~ "both",
      rain > qu                     ~ "rain",
      wind_speed > qu               ~ "wind_speed",
      TRUE                          ~ NA
    ))
  
  # return list with points falling into each category
  return(list(
    "rain"       = sum(df$extreme == "rain", na.rm = TRUE),
    "wind_speed" = sum(df$extreme == "wind_speed", na.rm = TRUE),
    "both"       = sum(df$extreme == "both", na.rm = TRUE)
  ))
}

#' @title Calculate proportion of each subset
#' @description Function to calculate proportion of each subset based on the
#' number of points in each subset.
#' @param x List with counts of points in each subset.
#' @return List with proportions of points in each subset.
#' @keywords internal
calc_prop <- \(x) {
  stopifnot(is.list(x) & names(x) == c("rain", "wind_speed", "both"))
  denom <- sum(unlist(x))
  ret <- lapply(x, \(y) y / denom)
  names(ret) <- names(x)
  return(ret)
}
