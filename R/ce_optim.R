#### Functions for custom MLE of conditional extremes model ####

#### Marginal Transformation ####

# TODO: Document these functions
# Semiparametric CDF
semi_par <- \(dat, gpd, n = nrow(dat)) {
  # As in Heff & Tawn '04, semiparametric mod uses ecdf below thresh, GPD above
  return(vapply(seq_along(gpd), \(i) {
    dat_spec <- dat[, i, drop = TRUE]
    stopifnot(names(gpd[[i]]) == c("sigma", "xi", "thresh"))
    spec_sigma <- gpd[[i]][[1]]
    spec_xi <- gpd[[i]][[2]]
    spec_loc <- gpd[[i]][[3]]
    
    # TODO: Replace with ecdf fun from evc
    # order and sort data
    dat_spec_ord <- order(dat_spec)
    dat_spec_sort <- dat_spec[dat_spec_ord]
    
    # calculate ECDF
    m <- length(dat_spec)
    ecdf_vals <- (seq_len(m)) / (m + 1)
    # convert back to original order
    ecdf_dat_ord <- numeric(m)
    ecdf_dat_ord[dat_spec_ord] <- ecdf_vals
    
    # initialise
    cdf <- numeric(n)
    # ecdf (i.e. non-parametric) below threshold
    cdf[dat[, i] <= spec_loc] <- ecdf_dat_ord[dat_spec <= spec_loc]
    # GPD above threshold
    # parametric part of model
    para <- pmax(
      0,
      1 + spec_xi * (dat_spec[dat_spec > spec_loc] - spec_loc) / spec_sigma
    )^(-1 / spec_xi)
    # TODO: texmex uses different formulation to H&T, investigate
    # cdf[dat[, i] > spec_loc] <- 1 - (1 - ecdf_spec_loc) * para
    cdf[dat[, i] > spec_loc] <- 1 - (mean(dat_spec > spec_loc) * para)
    return(cdf)
  }, FUN.VALUE = numeric(n)))
}

# transform to Laplace margins
laplace_trans <- \(F_hat, tol = .Machine$double.eps) {
  apply(F_hat, 2, \(x) {
    y <- pmin(pmax(x, tol), 1 - tol) 
    return(ifelse(y < 0.5, log(2 * y), -log(2 * (1 - y))))
  })
}


##### MLE ####

# functions for MLE largely taken from `texmex::mexDependence`

# check constraints on parameters under constrained Laplace estimation
ConstraintsAreSatisfied <- \(a, b, z, zpos, zneg, v) {
  C1e <- a <= min(
    1, 1 - b * min(z) * v^(b - 1), 1 - v^(b - 1) * min(z) + min(zpos) / v
  ) &
    a <= min(
      1, 1 - b * max(z) * v^(b - 1), 1 - v^(b - 1) * max(z) + max(zpos) / v
    )

  C1o <- a <= 1 &
    a > 1 - b * min(z) * v^(b - 1) &
    a > 1 - b * max(z) * v^(b - 1) &
    (1 - 1 / b) * (b * min(z))^(1 / (1 - b)) * 
    (1 - a)^(-b / (1 - b)) + min(zpos) > 0 &
    (1 - 1 / b) * (b * max(z))^(1 / (1 - b)) * 
    (1 - a)^(-b / (1 - b)) + max(zpos) > 0

  C2e <- -a <= min(
    1, 1 + b * v^(b - 1) * min(z), 1 + v^(b - 1) * min(z) - min(zneg) / v
    ) &
    -a <= min(
      1, 1 + b * v^(b - 1) * max(z), 1 + v^(b - 1) * max(z) - max(zneg) / v
    )

  C2o <- -a <= 1 &
    -a > 1 + b * v^(b - 1) * min(z) &
    -a > 1 + b * v^(b - 1) * max(z) &
    (1 - 1 / b) * (-b * min(z))^(1 / (1 - b)) * 
    (1 + a)^(-b / (1 - b)) - min(zneg) > 0 &
    (1 - 1 / b) * (-b * max(z))^(1 / (1 - b)) * 
    (1 + a)^(-b / (1 - b)) - max(zneg) > 0

  if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
    message("Strayed into impossible area of parameter space")
    C1e <- C1o <- C2e <- C2o <- FALSE
  }

  (C1e | C1o) && (C2e | C2o)
}

# function to evaluate (negative log) likelihood
laplace_nll <- \(yex, ydep, a, b, m, s, constrain, v, aLow) {
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)

  # give large value if parameters are out of bounds
  if (a < aLow | s < WeeNumber | a > 1 - WeeNumber | b > 1 - WeeNumber) {
    res <- BigNumber
  } else {
    # Assuming normal distribution for excesses, calculate mean & sd
    mu <- a * yex + m * yex^b
    sig <- s * yex^b

    # calculate log likelihood
    res <- sum(0.5 * log(2 * pi) + log(sig) + 0.5 * ((ydep - mu) / sig)^2)

    if (is.infinite(res)) {
      if (res < 0) {
        res <- -BigNumber
      } else {
        res <- BigNumber
      }
      warning("Infinite value of Q in mexDependence")
      # Apply Keef, Papastathopoulos constraints if specified
    } else if (constrain) {
      zpos <- range(ydep - yex) # q0 & q1
      z <- range((ydep - yex * a) / (yex^b)) # q0 & q1
      zneg <- range(ydep + yex) # q0 & q1

      if (!ConstraintsAreSatisfied(a, b, z, zpos, zneg, v)) {
        res <- BigNumber
      }
    }
  }
  res
}

# function to evaluate (negative) profile (log) likelihood and optimise over
laplace_npll <- function(yex, ydep, a, b, constrain, v, aLow) {
  # first, estimate Z by rearanging the conditional extremes equation
  Z <- (ydep - yex * a) / (yex^b)
  stopifnot(
    "NaNs in Z, conditional quantile may be negative" = all(!is.nan(Z))
  )
  # estimate nuisance parameters
  m <- mean(Z)
  s <- stats::sd(Z)

  # now estimate a and b
  res <- laplace_nll(
    yex, ydep, a, b,
    m = m, s = s, constrain, v, aLow = aLow
  )
  res <- list(profLik = res, m = m, s = s)
  res
}

# function to evaluate profile likelihood and optimise over
Qpos <- function(param, yex, ydep, constrain, v, aLow) {
  a <- param[1]
  b <- param[2]

  res <- laplace_npll(yex, ydep, a, b, constrain, v, aLow)
  res$profLik
}

# fit CE for all pairs of variables
ce_optim <- \(
  Y,
  dqu,
  start = c("a" = 0.01, "b" = 0.01),
  control = list(maxit = 1e6),
  constrain = TRUE,
  v = 10,
  aLow = -1
) {
  
  # object to return if we find an error
  err_obj <- as.matrix(c("a" = NA, "b" = NA, "m" = NA, "s" = NA))
  
  # optimise for a single variable vs another
  single_optim <- \(yex, ydep) {
    if (any(is.infinite(yex))) {
      message("Inf values in Laplace transformed data, optimisation failed")
      return(err_obj)
    }
    # threshold data
    thresh <- stats::quantile(yex, dqu)
    wch <- yex > thresh
    o_single <- try(stats::optim(
      par = start,
      fn = Qpos,
      control = control,
      yex = yex[wch],
      ydep = ydep[wch],
      constrain = constrain,
      v = v,
      aLow = aLow
    ), silent = TRUE)
    if (inherits(o_single, "try-error")) {
      message(paste("optimisation failed for constrain =", constrain))
      return(err_obj)
    }
    # set to NA if no change from starting values
    if (all(o_single$par[1:2] == start)) {
      message("No change from starting values, optimisation failed")
      return(err_obj)
    }
    # calculate nuisance parameters from a and b estimates
    if (all(!is.na(o_single$par))) {
      Z <- (ydep[wch] - yex[wch] * o_single$par[1]) /
        (yex[wch]^o_single$par[2])
      o_single$par <- c(o_single$par[1:2], "m" = mean(Z), "s" = stats::sd(Z))
    }
    # TODO: Add checks afterwards on o
    return(c(o_single$par, "ll" = o_single$value, "dth" = thresh[[1]]))
  }
  names_y <- colnames(Y)
  # loop through variables, fit CE model against other variables
  ret <- lapply(seq_len(ncol(Y)), \(i) {
    o_yex <- apply(Y[, -i, drop = FALSE], 2, \(x) single_optim(Y[, i], x))
    if (!is.null(names_y)) {
      colnames(o_yex) <- names_y[-i]
    }
    return(o_yex)
  })
  if (!is.null(names_y)) {
    names(ret) <- names_y
  }
  return(ret)
}
