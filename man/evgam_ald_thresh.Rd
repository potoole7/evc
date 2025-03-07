% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evgam.R
\name{evgam_ald_thresh}
\alias{evgam_ald_thresh}
\title{\code{evgam} varying threshold}
\usage{
evgam_ald_thresh(data, response, f, tau = 0.95, jitter = TRUE, thresh = TRUE)
}
\arguments{
\item{data}{Dataframe for one location which we wish to
threshold.}

\item{response}{Name of variable to threshold.}

\item{f}{Formula for \code{evgam} model.}

\item{tau}{Quantile to threshold at (see \link[evgam]{evgam} for details)}

\item{jitter}{Add jitter to data to remove 0s.}

\item{thresh}{return thresholded (i.e. filtered) data if TRUE.}
}
\value{
Dataframe with \code{thresh} and \code{excess} columns, optionally thresholded.
}
\description{
Fit varying threshold using quantile regression on
asymmetric Laplace distribution, as in
https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R.
}
