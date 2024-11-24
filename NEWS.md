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
