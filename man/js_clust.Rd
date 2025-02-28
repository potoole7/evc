% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/js_clust.R
\name{js_clust}
\alias{js_clust}
\title{Cluster based on Jensen-Shannon divergence for conditional extremes
model}
\usage{
js_clust(
  dependence,
  k = NULL,
  dist_mat = NULL,
  cluster_mem = NULL,
  dat_max_mult = 2,
  n_dat = 10,
  ncores = 1,
  par_dist = FALSE,
  scree_k = 1:5,
  return_dist = FALSE
)
}
\arguments{
\item{dependence}{List of \code{mexDependence} objects for each location.}

\item{k}{Number of clusters to fit, set to NULL to produce scree plot and
return distance matrix.}

\item{dist_mat}{Distance matrix, set to NULL to calculate.}

\item{cluster_mem}{Optional known cluster memberships for each location,
which can be used to evaluate the quality of the clustering.}

\item{dat_max_mult}{Multiplier for the maximum multiple of the largest
threshold across all locations.}

\item{n_dat}{Number of data points to use in the Jensen-Shannon divergence.}

\item{ncores}{Number of cores to use for parallel computation, Default: 1.}

\item{par_dist}{Logical, whether to parallelise the distance computation,
which is useful for a high number of locations, but has significant
overheads versus embarrassing parallelism over each variable, Default: FALSE.}

\item{scree_k}{Vector of k values to use in scree plot, if \code{k} is NULL.}

\item{return_dist}{Logical, whether to return distance matrix after
clustering.}
}
\value{
List containing the clustering results and, if \code{cluster_mem} is
provided, the adjusted Rand index.
}
\description{
Function to calculate and cluster on the Jensen-Shannon
Divergence for the conditional extremes model.
When comparing conditional extremes fits for a single variable using the JS
between any two locations, we need to use the same "data" for each location.
Therefore, we look at values from the maximum of the thresholds
at each location, to some multiple \code{dat_max_mult} of this, with
\code{n_dat} equally spaced points.
\code{cluster:pam} is used to cluster the Jensen-Shannon divergence distance
matrix between all locations, summed across all variables.
}
