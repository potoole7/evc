#### Functions for fitting Conditional Extremes model ####

# TODO: Extend to more than 2 variables
# TODO: Write tests for functions
# TODO: Allow use of base `texmex` with `fit_ce`

#' @title Fit Conditional Extremes model
#' @description Fit Conditional Extremes model using marginal fits form `evgam`.
#' @param data List of dataframes at each location containing data for each
#' variable.
#' @param vars Variable names for each location.
#' @param marg_prob Can be a list of arguments to `quantile_thresh` if a
#' variable quantile is desired, or a numeric value for simple quantile
#' thresholding.
#' @param cond_prob Conditional quantile for dependence modelling.
#' @param f Formula for `evgam` model.
#' @param split_data if `data` has variables stacked, unstack.
#' @param output_all Logical argument for whether to return quantiles, `evgam`
#' and marginal fits.
#' @return Object of type `mexDependence` for each location.
#' @rdname fit_ce
#' @export
fit_ce <- \(
  data,
  vars = c("rain", "wind_speed"),
  marg_prob = list(
    f        = list("response ~ name", "~ name"), # must be as character
    tau      = .95,
    jitter   = TRUE
  ),
  cond_prob   = 0.9,
  f           = list(excess ~ name, ~ 1), # keep shape constant for now
  split_data  = TRUE,
  output_all = FALSE
) {

  # initialise to remove `devtools::check()` note
  name <- excess <- shape <- n <- NULL

  # if marg_prob used as args to `quantile_thresh`, check args correct
  if (is.list(marg_prob)) {
    stopifnot(all(names(marg_prob) %in% names(formals(quantile_thresh))))
    stopifnot(is.list(marg_prob$f))
    if (!all(vapply(marg_prob$f, is.character, logical(1)))) {
      stop(paste(
        "f should be a list of characters where 'response' is replaced by",
        "each specified 'vars'"
      ))
    }
  }

  # TODO: Should really be done outside of this function, shouldn't be here!
  # convert to dataframe
  # TODO: Tidy up, quite verbose
  # TODO: Need to extend to multiple variables, using length of vars
  if (!is.data.frame(data) && split_data) {
    data_df <- dplyr::bind_rows(lapply(seq_along(data), \(i) {
      n <- length(data[[i]])
      ret <- data.frame(
        "rain"       = data[[i]][1:(n / 2)],
        "wind_speed" = data[[i]][((n / 2) + 1):n]
      ) |>
        dplyr::mutate(name = paste0("location_", i))
    }))
  } else if (!is.data.frame(data)) {
    data_df <- dplyr::bind_rows(lapply(seq_along(data), \(i) {
      as.data.frame(data[[i]]) |>
        dplyr::mutate(name = paste0("location_", i))
    }))
    names(data_df)[1:2] <- c("rain", "wind_speed") # TODO: Change to var_j
  } else {
    data_df <- data
  }

  # add row number for labelling

  # First, calculate threshold (90th quantile across all locations)
  if (!is.list(marg_prob) && is.numeric(marg_prob)) {
    thresh <- apply(data_df[, c(vars)], 2, stats::quantile, marg_prob)
    # for each variable, calculate excess over threshold
    data_thresh <- lapply(vars, \(x) {
      data_df |>
        # dplyr::select(dplyr::matches(x), name) |>
        # remote other responses, will be joined together after
        dplyr::select(-dplyr::matches(vars[vars != x])) |>
        dplyr::mutate(
          thresh = thresh[x],
          excess = !!rlang::sym(x) - thresh
        ) |>
        dplyr::filter(excess > 0)
    })
  # If thresh is a list, assume it is arguments to quantile_thresh
  } else if (is.list(marg_prob)) {
    # TODO: Redo with do.call
    data_thresh <- lapply(vars, \(x) {
      # Change formula to include response in question
      marg_prob$f <- lapply(marg_prob$f, \(f_spec) {
        stats::formula(stringr::str_replace_all(f_spec, "response", x))
      })
      # Run `quantile_thresh` for each response with specified args
      do.call(
        quantile_thresh,
        # data args
        args = c(list(data = data_df, response = x),
        marg_prob
        )
      )
    })
  } else {
    stop("marg_prob must be numeric or arguments to `quantile_thresh`")
  }

  # Now fit evgam model for each marginal
  # TODO: Allow use of base `texmex` (i.e. migpd instead of evgam) as well!
  evgam_fit <- lapply(data_thresh, \(x) {
    fit_evgam(
      data = x,
      pred_data = data_df,
      f = f # formula to use in evgam::evgam, specified arg above
    )
  })

  # Join scale and shape estimates into data
  # TODO: Functionalise to work with > 2 variables and for different names
  # pull variables specified as predictors in f
  preds <- unique(as.vector(unlist(lapply(
    evgam_fit, \(x) x$m$predictor.names
  ))))
  # TODO: Recursively perform this action
  # TODO: Change pipe function
  data_gpd <- data_df |>
    # add scale and shape parameters for each location
    dplyr::bind_cols(
      dplyr::rename(
        evgam_fit[[1]]$predictions,
        scale_rain = scale,
        shape_rain = shape),
      dplyr::rename(
        evgam_fit[[2]]$predictions,
        scale_ws = scale,
        shape_ws = shape),
    ) |>
    # add thresholds and number of exceedances for each predictor combo
    dplyr::left_join(
      data_thresh[[1]] |>
        dplyr::mutate(thresh = round(thresh, 3)) |>
        dplyr::count(dplyr::across(dplyr::all_of(preds)), thresh) |>
        dplyr::rename(n_rain = n, thresh_rain = thresh),
      by = preds
    ) |>
    dplyr::left_join(
      data_thresh[[2]] |>
        dplyr::mutate(thresh = round(thresh, 3)) |>
        dplyr::count(dplyr::across(dplyr::all_of(preds)), thresh) |>
        dplyr::rename(n_wind = n, thresh_wind = thresh),
      by = preds
    ) |>
    # fill in NAs (indicating no exceedances) with 0
    dplyr::mutate(
      dplyr::across(dplyr::contains("n_"), \(x) ifelse(is.na(x), 0, x))
    ) |>
    # must have unique rows to loop through in `gen_marg_migpd`
    # TODO: Still forces us to have name column, change in gen_marg_migpd?
    dplyr::distinct(name, dplyr::across(c(
      !!preds,
      dplyr::contains("scale"),
      dplyr::contains("shape"),
      dplyr::contains("thresh"),
      dplyr::contains("n_")))
    )

  # Now convert marginals to migpd (i.e. texmex format)
  # TODO: Working, but summary/coef doesn't show new thresholds (needed?)
  marginal <- gen_marg_migpd(data_gpd, data_df)
  names(marginal) <- paste0(data_gpd$name, " - ", data_gpd$county)

  # Calculate dependence from marginals (default output object)
  ret <- fit_texmex_dep(
    marginal,
    mex_dep_args = list(dqu = cond_prob),
    fit_no_keef = TRUE
  )
  # output more than just dependence object, if desired
  if (output_all) {
    ret <- list(
      "thresh"   = data_thresh,
      "evgam"    = evgam_fit,
      "marginal" = marginal,
      "dependence" = ret
    )
  }
  return(ret)
}

#' @title `evgam` varying threshold
#' @description Fit varying threshold using quantile regression, as in
#' https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R.
#' @param data Dataframe for one location which we wish to
#' threshold.
#' @param response Name of variable to threshold.
#' @param f Formula for `evgam` model.
#' @param tau Quantile to threshold at (see \link[evgam]{evgam} for details)
#' @param jitter Add jitter to data to remove 0s.
#' @param thresh return thresholded (i.e. filtered) data if TRUE.
#' @return Dataframe with `thresh` and `excess` columns, optionally thresholded.
#' @rdname quantile_thresh
#' @export
# TODO: Vary tau and see how that effects results (QQ plots, etc)
quantile_thresh <- function(
  data,
  response,
  f = list(
    # response ~ s(lon, lat, k = 50), # location
    stats::formula(paste(response, " ~ s(lon, lat, k = 50)")), # location
    ~ s(lon, lat, k = 40)                               # logscale
  ),
  tau = .95,
  jitter = TRUE,
  thresh = TRUE
) {

  excess <- NULL

  # jitter, if specified, to remove 0s when calculating quantiles
  if (jitter == TRUE) {
    data <- data |>
      dplyr::mutate(
        dplyr::across(dplyr::all_of(response), ~ . + stats::rnorm(n(), 0, 1e-6))
      )
  }

  # fit the quantile regression model at tau'th percentile
  ald_fit <- evgam::evgam(
    f, data, family = "ald", ald.args = list(tau = tau)
  )

  # add threshold to data
  data_thresh <- data |>
    dplyr::mutate(
      thresh = stats::predict(ald_fit)$location,
      excess = !!rlang::sym(response) - thresh
    )
  # threshold if desired
  if (thresh == TRUE) {
    data_thresh <- dplyr::filter(data_thresh, excess > 0)
  }

  return(data_thresh)
}

#' @title Fit `evgam` model
#' @description Fit and generate predictions from `evgam` model
#' @param data Dataframe for one location.
#' @param pred_data Dataframe for one location to predict on.
#' @param f Formula for `evgam` model.
#' @return List with model `m` and predictions `predictions`.
#' @rdname fit_evgam
#' @export
fit_evgam <- \(
  data,
  pred_data,
  # formula used in evgam, fitting to both scale and shape parameters
  f = list(
    excess ~ s(lon, lat, k = 40), # increase smoothing on scale parameter
    ~ s(lon, lat, k = 40) # shape parameter
  )
) {
  # fit evgam model
  m <- evgam::evgam(f, data = data, family = "gpd")

  # create predictions
  predictions <- stats::predict(m, pred_data, type = "response")

  # return model fit and predictions
  return(list(
    "m"           = m,
    "predictions" = predictions
  ))
}

#' @title Generate marginal `migpd` objects
#' @description Generate marginal `migpd` objects from `evgam` objects for each
#' site.
#' @param data_gpd Scale and shape parameters for each location.
#' @param data Data for each location.
#' @param mqu Marginal quantile for thresholding.
#' @return List of `migpd` objects for each location.
#' @rdname gen_marg_migpd
#' @export
gen_marg_migpd <- \(data_gpd, data, mqu = 0.95) {

  name <- rain <- wind_speed <- NULL

  # Create "dummy" migpd object to fill in with evgam values
  # TODO: Extend to > 2 variables
  dat_mat <- data |>
    dplyr::filter(name == data$name[[1]]) |>
    dplyr::select(rain, wind_speed) |>
    as.matrix()
  names(dat_mat) <- c("rain", "wind_speed")

  temp <- texmex::migpd(dat_mat, mqu = mqu, penalty = "none")

  marginal <- lapply(seq_len(nrow(data_gpd)), \(i) {
    # initialise
    spec_marg <- temp
    # replace data
    # TODO: Extend to work for more variables than just rain and wind_speed
    spec_marg$data <- data |>
      dplyr::filter(name == data_gpd$name[i]) |>
      dplyr::select(rain, wind_speed) |>
      as.matrix()
    names(spec_marg$data) <- c("rain", "wind_speed")
    # replace thresholds
    # TODO: Fix, thresholds not included in above
    spec_marg$models$rain$threshold <- data_gpd$thresh_rain[[i]]
    spec_marg$models$wind_speed$threshold <- data_gpd$thresh_wind[[i]]
    # replace coefficients
    spec_marg$models$rain$coefficients[1:2] <- c(
      log(data_gpd$scale_rain[i]),
      data_gpd$shape_rain[i]
    )
    spec_marg$models$wind_speed$coefficients[1:2] <- c(
      log(data_gpd$scale_ws[i]),
      data_gpd$shape_ws[i]
    )

    return(spec_marg)
  })
  return(marginal)
}

#' @title Fit `texmex` dependence model
#' @description Fit `texmex` dependence model for each site.
#' @param marginal List of `migpd` objects for each site.
#' @param vars Variables to fit dependence on.
#' @param mex_dep_args Arguments to pass to \link[texmex]{mexDependence}.
#' @param fit_no_keef If model doesn't fit under Keef constraints, fit without
#' (see \link[texmex]{mexDependence} for details).
#' @return List of `mexDependence` objects for each site.
#' @rdname fit_texmex_dep
#' @export
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
