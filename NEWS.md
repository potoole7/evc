# evc 0.2.1

- `fit_ce` and `js_clust` now can optionally be run in parallel using 
`parallel::mcapply` (on macOS and Linux). 

# evc 0.2.0

- `fit_ce` and `js_clust` can now work for data with more than two variables on 
which to condition and cluster on.

# evc 0.1.7

- Add plots for silhouette width associated with clustering solution:
  - `sil_boxplot` to plot a boxplot of silhouette widths for different choices 
  of cluster number, 
  - `plt_sil_map` to plot a map of silhouette widths for each location.

# evc 0.1.6

- `marg_prob` argument to `fit_ce` can now also be arguments (other than data and response) to `evgam::evgam` , as well as numeric quantile to threshold at.

# evc 0.1.5

- Add utility and plotting functions often used in conjunction with this package.

# evc 0.1.4

- Add functions to cluster based Vignotto 2021 method using empirical KL divergence between locations with bivariate extremal data.

# evc 0.1.3

- Add functions to cluster based on Jensen-Shannon divergence calculated from the estimated conditional extremes model.

# evc 0.1.2 

- Add functions to fit conditional extremes model with a given threshold using `evgam::evgam` to produce marginal distributions.
- Document functions and clean code so that it passes `devtools::check()`.

# evc 0.1.1

- Add a `NEWS.md` file to track changes to the package.
