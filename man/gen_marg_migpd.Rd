% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ce.R
\name{gen_marg_migpd}
\alias{gen_marg_migpd}
\title{Generate marginal \code{migpd} objects}
\usage{
gen_marg_migpd(data_gpd, data, vars, loop_fun = lapply)
}
\arguments{
\item{data_gpd}{Scale and shape parameters for each location.}

\item{data}{Data for each location.}

\item{vars}{Variable names for each location.}

\item{loop_fun}{Function to loop through each location, Default: \code{lapply}.}
}
\value{
List of \code{migpd} objects for each location.
}
\description{
Generate marginal \code{migpd} objects from \code{evgam} objects for each
site.
}
