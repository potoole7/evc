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
#' @param ncores Number of cores to use for parallel computation, Default: 1.
#' @param fit_no_keef If model doesn't fit under Keef constraints, fit without
#' @param output_all Logical argument for whether to return quantiles, `evgam`
#' and marginal fits.
#' @return Object of type `mexDependence` for each location.
#' @rdname fit_ce
#' @importFrom rlang .data :=
#' @export
fit_ce <- \(
  data,
  vars = c("rain", "wind_speed"),
  marg_prob = list(
    f      = list("response ~ name", "~ name"), # must be as character
    tau    = .95,
    jitter = TRUE
  ),
  cond_prob   = 0.9,
  f           = list(excess ~ name, ~ 1), # keep shape constant for now
  ncores      = 1,
  fit_no_keef = FALSE,
  output_all  = FALSE
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
  
  # Parallel setup
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }
  
  # number of variables
  nvars <- length(vars)

  # convert to data frame, if required
  if (!is.data.frame(data) && is.list(data)) {
    data_df <- dplyr::bind_rows(lapply(seq_along(data), \(i) {
      ret <- as.data.frame(data[[i]])
      # Add name column if not in list already
      if (!"name" %in% names(ret)) {
        ret <- ret |>
          dplyr::mutate(name = paste0("location_", i))
      }
    }))
    names(data_df)[seq_len(nvars)] <- vars
  } else {
    data_df <- data
  }

  # must have names column, and all of vars must be columns in data_df
  stopifnot("Must have a `name` column" = "name" %in% names(data_df))
  stopifnot("All of `vars` must be in data" = all(vars %in% names(data_df)))

  # First, calculate threshold (90th quantile across all locations)
  if (!is.list(marg_prob) && is.numeric(marg_prob)) {
    thresh <- apply(data_df[, c(vars)], 2, stats::quantile, marg_prob)
    # for each variable, calculate excess over threshold
    data_thresh <- lapply(vars, \(x) {
      data_df |>
        # remote other responses, will be joined together after
        dplyr::select(-dplyr::matches(vars[vars != x])) |>
        dplyr::mutate(
          thresh = thresh[x],
          excess = !!rlang::sym(x) - thresh
        ) |>
        dplyr::filter(excess > 0)
    })
  # If thresh is a list, assume it is arguments to quantile_thresh
  # TODO: May be easier to just copy each vars column as response in data_df
  # Would allow for simpler formula specification
  } else if (is.list(marg_prob)) {
    data_thresh <- lapply(vars, \(x) {
      # Change formula to include response in question
      marg_prob$f <- lapply(marg_prob$f, \(f_spec) {
        stats::formula(stringr::str_replace_all(f_spec, "response", x))
      })
      # Run `quantile_thresh` for each response with specified args
      do.call(
        quantile_thresh,
        args = c(list(data = data_df, response = x), # data args
        marg_prob
        )
      )
    })
  } else {
    stop("marg_prob must be numeric or arguments to `quantile_thresh`")
  }

  # If f NULL, fit ordinary marginal models with `texmex::migpd`
  if (is.null(f)) {
    # calculate marginal fits for all locations
    marginal <- data_df |> 
      dplyr::group_split(name) |>
      loop_fun(\(x) {
        # pull marginal thresholds
        mth <- vapply(data_thresh, \(y) {
          y |>
            dplyr::filter(name == x$name[[1]]) |> # need thresh for correct loc
            dplyr::slice(1) |>
            dplyr::pull(thresh)
        }, numeric(1))
        texmex::migpd(as.matrix(x[, vars]), mth = mth)
      })
  # Now fit evgam model for each marginal
  } else {
    evgam_fit <- loop_fun(data_thresh, \(x) {
      fit_evgam(
        data = x,
        pred_data = data_df,
        f = f # formula to use in evgam::evgam, specified arg above
      )
    })
  
    # Join scale and shape estimates into data
    # pull variables specified as predictors in f
    preds <- unique(as.vector(unlist(lapply(
      evgam_fit, \(x) x$m$predictor.names
    ))))
  
    # add predictions of scale and shape parameters for each variable
    data_df_wide <- data_df |>
      dplyr::bind_cols(
        # for each variable, take predictions for scale + shape, rename
        lapply(seq_along(evgam_fit), \(i) {
          evgam_fit[[i]]$predictions |>
            dplyr::rename(
              !!paste0("scale_", vars[[i]]) := scale,
              !!paste0("shape_", vars[[i]]) := shape
            )
        })
      )
  
    # add thresholds and number of exceedances for each predictor combination
    # TODO: Replace for loop somehow?
    data_df_wide_join <- data_df_wide
    for (i in seq_len(nvars)) {
      data_df_wide_join <- data_df_wide_join |>
        dplyr::left_join(
          data_thresh[[i]] |>
            dplyr::mutate(thresh = round(thresh, 3)) |>
            dplyr::count(dplyr::across(dplyr::all_of(preds)), thresh) |>
            dplyr::rename(
              !!paste0("n_", vars[[i]]) := n,
              !!paste0("thresh_", vars[[i]]) := thresh
            ),
          by = preds # predictors supplied to evgam formula
        )
    }
  
    data_gpd <- data_df_wide_join |>
      # fill in NAs (indicating no exceedances) with 0
      dplyr::mutate(
        dplyr::across(dplyr::starts_with("n_"), ~ifelse(is.na(.), 0, .))
      ) |>
      # must have unique rows to loop through in `gen_marg_migpd`
      dplyr::distinct(name, .keep_all = TRUE)
  
    # Now convert marginals to migpd (i.e. texmex format)
    marginal <- gen_marg_migpd(data_gpd, data_df, vars, loop_fun = loop_fun)  
  }
  names(marginal) <- unique(data_df$name)

  # Calculate dependence from marginals (default output object)
  # TODO: Replace with our own conditional extremes implementation
  ret <- fit_texmex_dep(
    marginal     = marginal,
    vars         = vars,
    mex_dep_args = list(dqu = cond_prob),
    fit_no_keef  = fit_no_keef, 
    loop_fun     = loop_fun
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
quantile_thresh <- function(
  data,
  response,
  f = list(
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
      dplyr::mutate(dplyr::across(
        dplyr::all_of(response), ~ . + stats::rnorm(n(), 0, 1e-6)
      ))
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
#' @param vars Variable names for each location.
#' @param loop_fun Function to loop through each location, Default: `lapply`.
#' @return List of `migpd` objects for each location.
#' @rdname gen_marg_migpd
#' @export
gen_marg_migpd <- \(data_gpd, data, vars, loop_fun = lapply) {

  name <- NULL

  # Create "dummy" migpd object to fill in with evgam values
  dat_mat <- data |>
    dplyr::filter(name == data$name[[1]]) |>
    dplyr::select(dplyr::all_of(vars)) |>
    as.matrix()

  # create dummy `migpd` object, can replace key values with those in data_gpd
  temp <- texmex::migpd(dat_mat, mqu = 0.9, penalty = "none")

  # loop through each location
  # jmarginal <- lapply(seq_len(nrow(data_gpd)), \(i) {
  marginal <- loop_fun(seq_len(nrow(data_gpd)), \(i) {
    # initialise
    spec_marg <- temp
    # replace data
    spec_marg$data <- data |>
      dplyr::filter(name == data_gpd$name[i]) |>
      dplyr::select(dplyr::all_of(vars)) |>
      as.matrix()

    # replace thresholds and coefficients for each variable (for each location)
    spec_marg$models <- lapply(seq_along(spec_marg$models), \(j) {
      # replace thresholds for each variable
      spec_mod <- spec_marg$models[[j]]
      spec_mod$threshold <- data_gpd[[
        paste0("thresh_", vars[j])
      ]][[i]] # thresh may not be fixed quantile, so take row specific value

      # replace coefficients for each variable
      spec_mod$coefficients[1:2] <- c(
        log(data_gpd[[paste0("scale_", vars[j])]][[i]]),
        data_gpd[[paste0("shape_", vars[j])]][[i]]
      )
      class(spec_mod) <- "evmOpt"
      return(spec_mod)
    })
    names(spec_marg$models) <- vars

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
#' @param loop_fun Function to loop through each location, Default: `lapply`.
#' (see \link[texmex]{mexDependence} for details).
#' @return List of `mexDependence` objects for each site.
#' @rdname fit_texmex_dep
#' @export
fit_texmex_dep <- \(
  marginal,
  vars         = c("rain", "wind_speed"),
  mex_dep_args = list(
    start     = c(0.01, 0.01),
    dqu       = 0.7,
    fixed_b   = FALSE,
    PlotLikDo = FALSE
  ),
  fit_no_keef  = FALSE, 
  loop_fun     = lapply
) {

  # dependence <- lapply(seq_along(marginal), \(i) {
  dependence <- loop_fun(seq_along(marginal), \(i) {
    # fit for rain and wind speed
    ret <- lapply(vars, \(col) {
      mod <- do.call(
        texmex::mexDependence,
        args = c(list(x = marginal[[i]], which = col), mex_dep_args)
      )

      # if model didn't optimise with Keef 2013 constrains, return NA
      ll <- mod$dependence$loglik
      if (any(is.na(ll)) || any(abs(ll) > 1e9)) {
        system(sprintf(
          'echo "%s"', "Model not fitting properly under Keef constraints"
        ))
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
