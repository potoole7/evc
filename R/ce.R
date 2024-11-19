#### Functions for fitting Conditional Extremes model ####

# TODO: Document functions
# TODO: Extend to more than 2 variables
# TODO: Write tests for functions

# Function to fit CE model
fit_ce <- \(
  data_mix, 
  vars = c("rain", "wind_speed"), 
  marg_prob = 0.9,
  cond_prob = 0.9,
  f = list(excess ~ name, ~ 1), # keep shape constant for now
  split_data = TRUE
) {
  
  # convert to dataframe
  # TODO: Tidy up, quite verbose
  if (split_data) {
    data_df <- bind_rows(lapply(seq_along(data_mix), \(i) {
      n <- length(data_mix[[i]])
      data.frame(
        "rain"       = data_mix[[i]][1:(n / 2)],
        "wind_speed" = data_mix[[i]][((n / 2) + 1):n]
      ) %>% 
        mutate(name = paste0("location_", i))
    }))
  } else {
    data_df <- bind_rows(lapply(seq_along(data_mix), \(i) {
      as.data.frame(data_mix[[i]]) %>% 
        mutate(name = paste0("location_", i))
    }))
    names(data_df)[1:2] <- c("rain", "wind_speed")
  }
  
  # First, calculate threshold (90th quantile across all locations)
  thresh <- apply(data_df[, 1:2], 2, quantile, marg_prob)
  
  # for each variable, calculate excess over threshold
  data_thresh <- lapply(vars, \(x) {
    data_df %>% 
      dplyr::select(matches(x), name) %>% 
      mutate(
        thresh = thresh[x],
        excess = !!sym(x) - thresh
      ) %>% 
      filter(excess > 0)
  })
  
  # Now fit evgam model for each marginal
  # TODO: Is just fitting different models by loc appropriate? (Yes I think!)
  evgam_fit <- lapply(data_thresh, \(x) {
    fit_evgam(
      data = x, 
      pred_data = x,
      # model scale and shape for each location
      # f = list(excess ~ name, ~ name) 
      f = f 
    )
  })
  
  # Join scale and shape estimates into data
  # TODO: Functionalise to work with > 2 variables
  data_gpd <- distinct(data_df, name) %>% 
    bind_cols(
      rename(
        distinct(evgam_fit[[1]]$predictions), 
        scale_rain = scale, 
        shape_rain = shape),
      rename(
        distinct(evgam_fit[[2]]$predictions), 
        scale_ws = scale, 
        shape_ws = shape),
    )
  
  # Now convert marginals to migpd (i.e. texmex format)
  marginal <- gen_marg_migpd(data_gpd, data_df)
  names(marginal) <- paste0(data_gpd$name, " - ", data_gpd$county)
  
  # Calculate dependence from marginals
  return(fit_texmex_dep(
    marginal, 
    mex_dep_args = list(dqu = cond_prob), 
    fit_no_keef = TRUE 
  ))
}

# function to fit varying threshold using quantile regression
# https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
# TODO: Vary tau and see how that effects results (QQ plots, etc)
# quantile_thresh <- function(data, response, tau = .95) {
quantile_thresh <- function(data, response, tau = .95, jitter = TRUE) {
  fm_ald <- list(
    # response ~ s(lon, lat, k = 50), # location
    formula(paste(response, " ~ s(lon, lat, k = 50)")), # location
    ~ s(lon, lat, k = 40)                               # logscale
  )
  
  # jitter, if specified, to remove 0s when calculating quantiles
  if (jitter == TRUE) {
    data <- data %>%
      mutate(across(all_of(response), ~ . + rnorm(n(), 0, 1e-6)))
  }
  
  # fit the quantile regression model at tau'th percentile
  ald_fit <- evgam::evgam(
    fm_ald, data, family = "ald", ald.args = list(tau = tau)
  )
  
  # add threshold to data and filter
  data %>% 
    mutate(
      thresh = evgam:::predict.evgam(ald_fit)$location, 
      # excess = wind_speed - thresh
      excess = !!sym(response) - thresh
    ) %>% 
    filter(excess > 0) %>% 
    return()
}

# Fit evgam model with spline for spatial location, create predictions
fit_evgam <- \(
  data, 
  pred_data,
  # formula used in evgam, fitting to both scale and shape parameters
  f = list(
    excess ~ s(lon, lat, k = 40), # increase smoothing on scale parameter
    ~ s(lon, lat, k = 40) # shape parameter
  ) 
) {
  # model formula
  # f <- list(
  #   excess ~ s(lon, lat, k = 40), # increase smoothing on scale parameter
  #   ~ s(lon, lat, k = 40) # shape parameter
  # ) 
  # fit evgam model
  m <- evgam::evgam(f, data = data, family = "gpd")

  # create predictions
  predictions <- evgam:::predict.evgam(m, pred_data, type = "response")
  
  # return model fit and predictions
  return(list(
    "m"           = m,
    "predictions" = predictions
  ))
}

# create marginal `migpd` objects from `evgam` objects for each site
# - data_gpd: scale and shape parameters for each location
# - data: Data for each location
gen_marg_migpd <- \(data_gpd, data, mqu = 0.95) {
  # Create "dummy" migpd object to fill in with evgam values
  dat_mat <- data %>% 
    filter(name == data$name[[1]]) %>% 
    dplyr::select(rain, wind_speed) %>% 
    as.matrix()
  names(dat_mat) <- c("rain", "wind_speed")
  
  temp <- texmex::migpd(dat_mat, mqu = mqu, penalty = "none")
  # m <- evm(y = rain, data = data, qu = 0.95, penalty = "none", famuly = "gpd")
  # m1 <- update(m, phi = ~lon + lat)
  
  marginal <- lapply(seq_len(nrow(data_gpd)), \(i) {
    # initialise
    # browser()
    spec_marg <- temp
    # replace data 
    spec_marg$data <- data %>% 
      filter(name == data_gpd$name[i]) %>% 
      dplyr::select(rain, wind_speed) %>% 
      as.matrix()
    names(spec_marg$data) <- c("rain", "wind_speed")
    # replace thresholds
    # spec_marg$models$rain$threshold <- thresh_rain
    # spec_marg$models$wind_speed$threshold <- thresh_wind
    spec_marg$models$rain$threshold <- data_gpd$thresh_rain[[i]]
    spec_marg$models$wind_speed$threshold <- data_gpd$thresh_wind[[i]]
    # replace coefficients
    spec_marg$models$rain$coefficients[1:2] <- c(
      # data_gpd$scale_rain[i], 
      log(data_gpd$scale_rain[i]),
      data_gpd$shape_rain[i]
    )
    spec_marg$models$wind_speed$coefficients[1:2] <- c(
      # data_gpd$scale_ws[i], 
      log(data_gpd$scale_ws[i]), 
      data_gpd$shape_ws[i]
    )
    
    return(spec_marg)
  })
  return(marginal)
}

# Fit CE dependence model for each site
fit_texmex_dep <- \(
  marginal, 
  vars = c("rain", "wind_speed"), 
  mex_dep_args = list(
    start = c(0.01, 0.01), 
    dqu = 0.7,
    fixed_b = FALSE,
    PlotLikDo = FALSE
  ),
  fit_no_keef = FALSE
) {
  dependence <- lapply(seq_along(marginal), \(i) {
    # fit for rain and wind speed
    ret <- lapply(vars, \(col) {
      mod <- do.call(
        texmex::mexDependence, 
        args = c(list(x = marginal[[i]], which = col), mex_dep_args)
      )
      
      # if model didn't optimise with Keef 2013 constrains, return NA
      ll <- mod$dependence$loglik
      if (is.na(ll) || abs(mod$dependence$loglik) > 1e9) {
        message("Model not fitting properly under Keef constraints")
        if (fit_no_keef) {
          mod <- do.call(
            texmex::mexDependence, 
            args = c(
              list(x = marginal[[i]], which = col, constrain = FALSE), 
              mex_dep_args
            )
          )
        } else {
          return(NA)
        }
      }
      return(mod)
    })
    names(ret) <- vars
    return(ret)
  })
  names(dependence) <- names(marginal)
  return(dependence)
}

