# Functions for fitting marginal models and estimating thresholds with evgam

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

