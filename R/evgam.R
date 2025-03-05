# Functions for fitting marginal models and estimating thresholds with evgam

#' @title `qgam` varying threshold
#' @description Fit varying threshold using quantile regression via `qgam`.
#' @param data Dataframe for one location which we wish to
#' threshold.
#' @param response Name of variable to threshold.
#' @param f Formula for `evgam` model.
#' @param tau Quantile to threshold at (see \link[evgam]{evgam} for details)
#' @param jitter Add jitter to data to remove 0s.
#' @param thresh return thresholded (i.e. filtered) data if TRUE.
#' @return Dataframe with `thresh` and `excess` columns, optionally thresholded.
#' @rdname qgam_thresh
#' @export
qgam_thresh <- function(
    data,
    response,
    f,
    tau    = .95,
    jitter = TRUE,
    thresh = TRUE
) {
  
  excess <- NULL
  
  # jitter, if specified, to remove 0s when calculating quantiles
  # TODO Change jitter to match magnitude of data
  if (jitter == TRUE) {
    data <- data |>
      dplyr::mutate(dplyr::across(
        dplyr::all_of(response), ~ . + abs(
          stats::rnorm(dplyr::n(), 0, 1e-6)
        )
      ))
  }
  
  # fit the quantile regression model at tau'th percentile
  qgam_fit <- qgam::qgam(
    f, data, qu = tau
  )
  # print(summary(qgam_fit))
  
  # add threshold to data
  predictors <- predictors <- names(qgam_fit$var.summary)
  predictions <- data |> 
    dplyr::mutate(thresh = qgam_fit$fitted.values) |> 
    dplyr::distinct(across(c(all_of(predictors), thresh)))
  
  data_thresh <- data |>
    dplyr::left_join(predictions, by = predictors) |>
    dplyr::mutate(excess = !!rlang::sym(response) - thresh)
  # threshold if desired
  if (thresh == TRUE) {
    data_thresh <- dplyr::filter(data_thresh, excess > 0)
  }
  
  # return model fit and thresholded data 
  return(list(
    "m"           = qgam_fit,
    "data_thresh" = data_thresh
  ))
}


#' @title `evgam` varying threshold
#' @description Fit varying threshold using quantile regression on 
#' asymmetric Laplace distribution, as in
#' https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R.
#' @param data Dataframe for one location which we wish to
#' threshold.
#' @param response Name of variable to threshold.
#' @param f Formula for `evgam` model.
#' @param tau Quantile to threshold at (see \link[evgam]{evgam} for details)
#' @param jitter Add jitter to data to remove 0s.
#' @param thresh return thresholded (i.e. filtered) data if TRUE.
#' @return Dataframe with `thresh` and `excess` columns, optionally thresholded.
#' @rdname evgam_ald_thresh 
#' @export
evgam_ald_thresh <- function(
    data,
    response,
    # f      = list(
    #   stats::formula(paste(response, " ~ s(lon, lat, k = 50)")), # location
    #   ~ s(lon, lat, k = 40)                                      # logscale
    # ),
    f,
    tau    = .95,
    jitter = TRUE,
    thresh = TRUE
) {
  
  excess <- NULL
  
  # jitter, if specified, to remove 0s when calculating quantiles
  if (jitter == TRUE) {
    data <- data |>
      dplyr::mutate(dplyr::across(
        dplyr::all_of(response), ~ . + abs(stats::rnorm(dplyr::n(), 0, 1e-6))
      ))
  }
  
  # fit the quantile regression model at tau'th percentile
  # f <- list(rain ~ s(lon, lat), ~s(lon, lat))
  ald_fit <- evgam::evgam(
    f, data, family = "ald", ald.args = list(tau = tau)
  )
  # print(summary(ald_fit))
  
  # make predictions on distinct predictors
  predictors <- ald_fit$predictor.names
  data_distinct <- data |>
    dplyr::distinct(dplyr::across(dplyr::all_of(predictors)))
  predictions <- cbind(
    data_distinct, 
    "thresh" = stats::predict(ald_fit, data_distinct)$location
  )
  
  # add threshold to data
  data_thresh <- data |>
    dplyr::left_join(predictions, by = predictors) |>
    dplyr::mutate(excess = !!rlang::sym(response) - thresh)
  # threshold if desired
  if (thresh == TRUE) {
    data_thresh <- dplyr::filter(data_thresh, excess > 0)
  }
  
  # return model fit and thresholded data 
  return(list(
    "m"           = ald_fit,
    "data_thresh" = data_thresh
  ))
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
    excess ~ s(lon, lat), # increase smoothing on scale parameter
    ~ s(lon, lat) # shape parameter
  )
) {
  # ensure f is a formula
  f <- lapply(f, stats::formula)
  # fit evgam model
  m <- evgam::evgam(f, data = data, family = "gpd")
  
  # create predictions for unique rows in pred_data (ensures one pred per loc)
  predictors <- m$predictor.names
  pred_dat_distinct <- pred_data |>
    dplyr::distinct(name, dplyr::across(dplyr::all_of(predictors)))
  predictions <- cbind(
    pred_dat_distinct, 
    stats::predict(m, pred_dat_distinct, type = "response")
  )
  
  # return model fit and predictions
  return(list(
    "m"           = m,
    "predictions" = predictions
  ))
}
