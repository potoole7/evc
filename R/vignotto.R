# function to calculate KL divergence between any two locations
emp_kl_div <- \(x, y, prob = 0.9) {
  # stopifnot(length(x) == length(y))
  
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
  # sum_vals[is.nan(sum_vals) | is.infinite(sum_vals)] <- 0
  if (any(is.nan(sum_vals) | is.infinite(sum_vals), na.rm = TRUE)) return(NA)
  return((1  / 2) * sum(sum_vals))
}


# Function to calculate KL divergence for 2 cluster simulation dataset
kl_sim_eval <- \(
  data_mix, # Mixture data from copulas
  kl_prob,  # Extremal quantile
  k = 2,    # number of clusters
  cluster_mem = NULL, # known cluster membership, for ARI
  ...
) {
  
  # prob must be valid probability
  stopifnot(0 <= kl_prob && kl_prob <= 1)
  
  # KL divergence between areas using Vignotto 2021 method
  kl_mat <- proxy::dist(
    data_mix, method = emp_kl_div, print = FALSE, prob = kl_prob, ...
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
    # print(paste0("adjusted Rand index" = adj_rand))
    ret <- c(ret, list("adj_rand" = adj_rand))
  }
  return(ret)
}


# function to convert to Pareto scale
# divide by N + 1, not N (so do i / n + 1)
# x: Vector representing e.g. rain
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

# function to compute bivariate risk function
# x: list of vectors representing different 
# fun: fun used to calculate "risk", defaults to max as easier to partition
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

# Partition into 3 subsets, calculate empirical proportion in each
# TODO: Only works for max risk fun, unsure how to do for sum...
# TODO: Fix, obviously wrong since it doesn't identify times where both 
# are high!
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
  qu <- quantile(df$R, prob = prob)[[1]]
  
  # partition into 3 subsets
  df <- df %>% 
    mutate(extreme = case_when(
      rain > qu & wind_speed > qu   ~ "both",
      rain > qu                     ~ "rain",
      wind_speed > qu               ~ "wind_speed",
      TRUE                          ~ NA
    ))
  
  if (plot) {
    df %>% 
      ggplot() +
      geom_point(aes(x = rain, y = wind_speed, colour = extreme)) + 
      geom_vline(xintercept = qu, linetype = "dashed") + 
      geom_hline(yintercept = qu, linetype = "dashed") + 
      ggsci::scale_colour_nejm() + 
      labs(y = "wind speed", title = paste0(
        "Hazard subsets, q_u = ", round(qu, 3), " (prob = ", prob, ")"
      )) + 
      theme
  }
  
  # return list with points falling into each category
  return(list(
      "rain"       = sum(df$extreme == "rain", na.rm = TRUE),
      "wind_speed" = sum(df$extreme == "wind_speed", na.rm = TRUE),
      "both"       = sum(df$extreme == "both", na.rm = TRUE)
    )
  )
}

# Function to calculate proportion 
calc_prop <- \(x) {
  stopifnot(is.list(x) & names(x) == c("rain", "wind_speed", "both"))
  # denom <- length(unlist(x))
  denom <- sum(unlist(x))
  ret <- lapply(x, \(y) y / denom)
  names(ret) <- names(x)
  return(ret)
}

