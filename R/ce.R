#### Functions for fitting Conditional Extremes model ####

# TODO Write tests for functions
# TODO Add silence option for messages
# TODO Separate/decouple marginal and dependence fitting funs as in `texmex`

#' @title Fit Conditional Extremes model
#' @description Fit Conditional Extremes model using marginal fits form `evgam`.
#' @param data List of dataframes at each location containing data for each
#' variable.
#' @param vars Variable names for each location.
#' @param marg_prob Can be a list of arguments to `quantile_thresh` if a
#' variable quantile is desired, or a numeric value for simple quantile
#' thresholding.
#' TODO Expand marg_val to work with list of lists for each location.
#' @param marg_val Explicit value for marginal thresholds for each variable, 
#' only specified if marg_prob is NULL. 
#' @param cond_prob Conditional quantile for dependence modelling.
#' @param f Formula for `evgam` model.
#' @param start Starting values for dependence model, can be either a vector or
#' list of lists of matrices for each conditioning variable,
#' Default: `c("a" = 0.01, "b" = 0.01)`.
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
  marg_val    = NULL,
  cond_prob   = 0.9,
  f           = list(excess ~ name, ~ 1), # keep shape constant for now
  start       = c("a" = 0.01, "b" = 0.01),
  ncores      = 1,
  fit_no_keef = FALSE,
  output_all  = FALSE
) {
  
  ## Setup ##

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
  
  # check marginal thresholds specified correctly
  if (is.null(marg_val) && is.null(marg_prob)) { # must provide one
    stop("you must provide one of mth or mqu")
  }
  if (!is.null(marg_val) && !is.null(marg_prob)) { # must provide only one
    stop("you must provide precisely one of marg_val or marg_prob")
  }
  # must have marginal values for each variable
  # TODO Expand so marg_val can be a list of locations like `start`
  # TODO Expand so marg_val can be a list of dataframes with varying thresh
  if (!is.null(marg_val) && is.numeric(marg_val)) { 
    stopifnot(length(marg_val) == length(vars) && all(names(marg_val) == vars))
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
    data_df <- as.data.frame(data)
    # convert matrix names if required
    names(data_df)[names(data_df) %in% paste0("V", seq_len(nvars))] <- vars
  }

  # must have names column, and all of vars must be columns in data_df
  stopifnot("Must have a `name` column" = "name" %in% names(data_df))
  stopifnot("All of `vars` must be in data" = all(vars %in% names(data_df)))
  
  # locations
  locs <- unique(data_df$name)
  
  ## Threshold Selection ##

  # If threshold isn't modelled:
  # If marg_val not specified, calculate thresh as quantile across all locs
  # TODO Allow marg_val to change by location
  if (is.null(marg_prob) || (!is.list(marg_prob) && is.numeric(marg_prob))) {
    # TODO Could also calculate for each location/name??
    if (is.null(marg_val)) {
      marg_val <- apply(
        data_df[, c(vars)], 2, stats::quantile, marg_prob, na.rm = TRUE
      )
    }
    # for each variable, calculate excess over threshold
    data_thresh <- lapply(vars, \(x) {
      data_df |>
        # remove other responses, will be joined together after
        dplyr::select(-dplyr::all_of(vars[vars != x])) |>
        dplyr::mutate(
          thresh = marg_val[x],
          excess = !!rlang::sym(x) - marg_val[x]
        ) |>
        dplyr::filter(excess > 0)
    })
    # If thresh is a list, assume it is arguments to quantile_thresh
    # TODO: May be easier to just copy each vars column as response in data_df
    # Would allow for simpler formula specification
  } else if (is.list(marg_prob)) {
    data_thresh <- lapply(vars, \(x) {
      # Change formula to include response in question
      spec_params <- marg_prob
      spec_params$f <- lapply(marg_prob$f, \(f_spec) {
        stats::formula(stringr::str_replace_all(f_spec, "response", x))
      })
      
      # Run `quantile_thresh` for each response with specified args
      ret <- do.call(
        quantile_thresh,
        args = c(
          list(data = data_df, response = x), # data args
          spec_params
        )
      )
      
      # return message where no exceedances are observed for any locations
      thresh_locs <- unique(ret$name)
      if (length(thresh_locs) < length(locs)) {
        loc_missing <- setdiff(locs, thresh_locs)
        message(paste0(
          "No exceedances for variable ", x, " at: ",
          paste(loc_missing, collapse = ", "), 
          ", removing for all variables"
        ))
      }
      return(ret)
    })
  } else {
    stop("marg_prob must be numeric or arguments to `quantile_thresh`")
  }
  # Only keep locs with exceedances for all vars, otherwise can't do CE!
  # TODO evc only works for locations, add error if data not in this form!
  locs_keep <- Reduce(intersect, lapply(data_thresh, \(x) unique(x$name)))
  data_df <- dplyr::filter(data_df, name %in% locs_keep)
  data_thresh <- lapply(data_thresh, \(x)
    dplyr::filter(x, name %in% locs_keep)
  )
  
  ## Marginal Model ##

  # If f NULL, fit ordinary marginal models with `ismev::gpd.fit` for each loc
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
        # texmex::migpd(as.matrix(x[, vars]), mth = mth)
        gpd_fits <- lapply(seq_along(vars), \(i) {
          fit <- ismev::gpd.fit(x[[vars[i]]], threshold = mth[i], show = FALSE)
          return(list(
            "sigma"     = fit$mle[1], 
            "xi"        = fit$mle[2], 
            "thresh"    = fit$threshold[[1]]
          )) 
        })
        names(gpd_fits) <- vars
        return(gpd_fits)
      })
    
  # fit evgam model for each marginal
  } else {
    evgam_fit <- loop_fun(data_thresh, \(x) {
      fit_evgam(
        data      = x,
        pred_data = data_df,
        f         = f 
      )
    })
    
    # pull scale, shape and threshold for each variable 
    marginal <- lapply(seq_along(evgam_fit), \(i) {
      
      # pull for each location (or predictor combo)
      params <- dplyr::distinct(
        evgam_fit[[i]]$predictions, sigma = scale, xi = shape
      ) |>
        # TODO Assumes threshs in same order as evgam preds, may lead to bugs
        dplyr::mutate(thresh = unique(data_thresh[[i]]$thresh)) |>
        dplyr::group_split(dplyr::row_number(), .keep = FALSE) |>
        lapply(as.vector, mode = "list")
    })
    names(marginal) <- vars
    
    # transpose list from variables -> locations to locations -> variables
    marginal <- purrr::transpose(marginal)
  }
  names(marginal) <- locs_keep
  
  ## Dependence ##
  
  # Test that start values for dependence parameters (a and b)  
  test_start <- \(start, marginal) {
    # check same locations
    test_loc <- length(start) == length(marginal)
    # check same variables
    test_var <- all(unlist(lapply(seq_along(start), \(i) {
      length(start[[i]]) == length(marginal[[i]]) && 
        all(names(start[[i]]) == names(marginal[[i]]))
    })))
    # check start values for each variable against each conditioning variable
    test_dim <- all(unlist(lapply(seq_along(start), \(i) {
      lapply(seq_along(start[[i]]), \(j) {
        all(colnames(start[[i]][[j]]) == names(marginal[[i]])[-j]) && 
          nrow(start[[i]][[j]]) == 2
      })
    })))
    stopifnot(
      "Start values not correct" = all(c(test_loc, test_var, test_dim))
    )
  }

  # check if start values for dependence is a list
  is_list_start <- FALSE
  if (is.list(start)) {
    is_list_start <- TRUE
    test_start(start, marginal)
  }
  
  # Calculate dependence from marginals (default output object)
  # first, transform margins to Laplace
  marginal_trans <- lapply(seq_along(marginal), \(i) {
    # semi-parametric CDF
    # TODO: More efficient to also split data_df by name and subset with i
    F_hat <- data_df |>
      dplyr::filter(name == locs_keep[i]) |>
      dplyr::select(dplyr::all_of(vars)) |>
      semi_par_cdf(marginal[[i]])
    # Laplace transform
    Y <- laplace_trans(F_hat)
    colnames(Y) <- vars
    return(Y)
  })
  names(marginal_trans) <- locs_keep
  
  # now fit dependence model
  ret <- loop_fun(seq_along(marginal_trans), \(i) {
    if (is_list_start) {
      start_spec <- start[[i]]
    } else start_spec <- start
    o <- ce_optim(
      Y         = marginal_trans[[i]],
      dqu       = cond_prob,
      control   = list(maxit = 1e6),
      constrain = !fit_no_keef, 
      start     = start_spec
    )
    return(o)
  })
  names(ret) <- locs_keep
  
  # check that all dependence models have run successfully, message if not
  locs_fail <- locs_keep[vapply(ret, \(x) any(is.na(unlist(x))), logical(1))]
  if (length(locs_fail) > 0) {
    message(paste0(
      length(locs_fail), 
      " locations failed to fit CE model for at least one variable: ",
      paste(locs_fail, collapse = ", ")
    ))  
  }

  # output more than just dependence object, if desired
  if (output_all) {
    arg_vals <- as.list(match.call())[-1]  
    arg_vals <- arg_vals[
      !names(arg_vals) %in% c("data", "start", "ncores", "output_all")
    ]
    arg_vals <- lapply(arg_vals, eval, envir = parent.frame())
    
    # original data
    orig_dat <- dplyr::group_split(data_df, name, .keep = FALSE)
    names(orig_dat) <- locs_keep
    
    # TODO Output object with all these as attributes
    ret <- list(
      "arg_vals"    = arg_vals,
      "thresh"      = data_thresh,
      "marginal"    = marginal,
      "original"    = orig_dat,
      "transformed" = marginal_trans,
      "dependence"  = ret
    )
    if (exists("evgam_fit", envir = environment())) {
      ret$evgam_fit <- evgam_fit
    }
  }
  return(ret)
}
