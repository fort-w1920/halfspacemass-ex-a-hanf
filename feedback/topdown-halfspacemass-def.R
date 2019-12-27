# Implements "Halfspace Mass" (a data depth) as well as
# a Monte Carlo-approximate halfspace (Tukey) depth
#
# Reference:
#   Chen, B., Ting, K.M., Washio, T. et al., Mach Learn (2015): 100(2):677--699
#   <Half-space mass: a maximally robust and efficient data depth method>
#   https://doi.org/10.1007/s10994-015-5524-x
# Author: Fabian Scheipl, 2017-11-20

################################################################################
# implements Algorithm 1 of Chen et al.

# inputs:
#   data: a matrix: **rows** are vectors / points to project
#   n_halfspace: how many halfspaces to draw
#   subsample: what proportion of data to use for each halfspace computation
#   scope: >= 1, controls size of region of convexity for halfspace mass
#       ($\lambda$ in the paper) i.e., how far outside of sampled data range
#       the sampled hyperplanes can lie...
#   seed: optional RNG seed
# output:
#   a list of n_halfspace halfspaces, defined by their normal vector and offset
#   from origin, with estimated data frequencies above/below halfspace boundary
train_depth <-
  function(data, n_halfspace = 1e3, subsample = 1, scope = 1, seed = NULL) {
    if (inherits(data, "data.frame")) data <- as.matrix(data)
    checkmate::assert_matrix(data,
      mode = "numeric", any.missing = FALSE, min.cols = 2,
      min.rows = 1
    )
    checkmate::assert_integerish(n_halfspace, lower = 1)
    checkmate::assert_number(subsample, lower = 1e-3, upper = 1)
    checkmate::assert_number(scope, lower = 1, na.ok = FALSE, finite = TRUE)
    if (!is.null(seed)) {
      checkmate::assert_integerish(seed, lower = 1L)
      set.seed(as.integer(seed))
    }
    dims <- ncol(data)
    normals <- get_directions(n_halfspace, dims = dims)
    # see end of file for alternative definition of get_direction
    halfspaces <- apply(normals, 2, get_halfspace,
      data = data, subsample = subsample, scope = scope
    )
    # save seed, training data for reproducibility
    structure(halfspaces, seed = seed, train_data = data)
  }

# return normal direction & location of halfspace boundary and relative
#   frequencies of data above/below it
# inputs: see subroutines
# output: a list with entries
#   normal: (as given),
#   split: the split point along normal defining the halfspace,
#   mass: a named vector giving the relative frequencies of subsampled data that
#     lie above and below the halfspace boundary
get_halfspace <- function(normal, data, subsample, scope) {
  sampledata <- subsample(data, subsample)
  projections <- project_scalar(sampledata, normal)
  split <- sample_split(projections, scope)
  list(
    normal = normal,
    split = split,
    mass = compute_masses(projections, split)
  )
}

# uniformly sample n d-dimensional direction vectors (i.e., points on the unit
#   sphere in d dimensions, see
#   http://mathworld.wolfram.com/SpherePointPicking.html., eq. [16])
# used as normal vectors of the (hyper-)planes defining the halfspaces
# inputs:
#   n: how many
#   dims: required dimension
# output:
#   a <dims> x <n> matrix of <n> <dims>-dimensional directions
get_directions <- function(n_halfspace, dims = 2) {
  checkmate::assert_integerish(dims, lower = 2)
  # no need to scale to length 1 here since later steps adjust for length anyway
  matrix(rnorm(dims * n_halfspace), nrow = dims)
}

# compute projection(s) of point(s) on a vector
#   Equations: https://en.wikipedia.org/wiki/Scalar_projection
#   used to project data points onto the normal vector defining the orientation
#   of a halfspace
# inputs:
#   data: a matrix: rows(!) are vectors / points to project
#   direction: numeric ncol(data) vector: the direction on which to project
# output:
#   numeric nrow(data) vector: projections of <data> on <direction> in
#   units of "length of <direction>"
project_scalar <- function(data, direction) {
  checkmate::assert_numeric(direction, any.missing = FALSE, min.len = 2)
  checkmate::assert_matrix(data,
    any.missing = FALSE, ncols = length(direction),
    mode = "numeric"
  )
  # would be crossprod(x, d)/sqrt(crossprod(d)) for column vectors but
  # no need to scale with sqrt(crossprod(d)) since split is drawn from
  # observed projection values
  data %*% direction
}

# sample a number in [midrange - scope/2 (max-min), midrange + scope/2 (max-min)]
#   i.e. from [min, max] with scope = 1, from [mid, mid] with  scope = 0
# inputs:
#   projections: numeric vector: projections of data on a halfspace normal vector
#   scope: scalar: parameter controlling sampled range
# output:
#   "offset" of the halfspace boundary along its normal vector
sample_split <- function(projections, scope) {
  checkmate::assert_numeric(projections, any.missing = FALSE, finite = TRUE)
  minmax <- range(projections)
  span <- minmax[2] - minmax[1]
  mid <- mean(minmax)
  runif(1, min = mid - scope / 2 * span, max = mid + scope / 2 * span)
}

# get proportion of observations above/below boundary separating 2 halfspaces
# inputs  :
#   projections: numeric vector: projections of data on normal vector of halfspace
#   split: "offset" of the halfspace boundary along its normal vector
# output:
#   a vector with relative data frequencies "below" and "above" the boundary
compute_masses <- function(projections, split) {
  checkmate::assert_numeric(projections, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(split, finite = TRUE)
  c(
    below = mean(projections < split),
    above = mean(projections >= split)
  )
}

# subsample a proportion <subsample> from <data>
subsample <- function(data, subsample = 1) {
  nrows <- NROW(data)
  subsample_size <- round(subsample * seq_len(nrows))
  subsample_size <- min(max(1, subsample_size), nrows)
  if (subsample_size == nrows) return(data)
  use <- sample(seq_len(nrows), subsample_size)
  data[use, ]
}

################################################################################
# functions to implement Algorithm 2 of Chen et al.

# evaluate halfspace-mass or (approximate) -depth for points in data based
#   halfspaces returned by get_halfspaces()
# inputs:
#   data: numeric data.frame or matrix containing data points in rows
#   halfspaces: list, return object of get_halfspaces()
#   metric: "mass" for Halfspace mass, "depth" for HS depth, defaults to "mass"
# output:
#   numeric vector with depth metric values for data: either HS mass or HS depth
evaluate_depth <- function(data, halfspaces, metric = c("mass", "depth")) {
  if (inherits(data, "data.frame")) data <- as.matrix(data)
  checkmate::assert_matrix(data,
    mode = "numeric", any.missing = FALSE, min.cols = 2,
    min.rows = 1
  )
  check_halfspaces(halfspaces, data)
  metric <- match.arg(metric)

  if (metric == "mass") {
    get_mass(data, halfspaces)
  } else {
    get_depth(data, halfspaces)
  }
}

# check correct structure of halfspaces w.r.t. data
check_halfspaces <- function(halfspaces, data) {
  names <- c("normal", "split", "mass")
  lengths <- c(NCOL(data), 1, 2)
  check_halfspace <- function(halfspace, names, lengths) {
    is.list(halfspace) & all(names(halfspace) %in% names) &
      all(sapply(halfspace, length) == lengths) &
      sum(halfspace$mass) == 1
  }
  wrong <- !sapply(halfspaces, check_halfspace, names = names, lengths = lengths)
  if (all(wrong)) stop("Halfspaces are not in correct format.")
  if (any(wrong)) {
    stop(
      "Halfspaces ", paste(which(wrong), collapse = ", "),
      " not in correct format."
    )
  }
}

# computes approximate halfspace mass
# inputs/outputs: see evaluate_depth()
get_mass <- function(data, halfspaces) {
  result <- numeric(NROW(data))
  for (halfspace in halfspaces) {
    projections <- project_scalar(data, halfspace$normal)
    # for each combination of a halfspace and a data point, use the mass on the
    # side of the split which the data point lies on...
    result <- result + ifelse(projections < halfspace$split,
      yes = halfspace$mass["below"],
      no = halfspace$mass["above"]
    )
  }
  # ... and take the mean:
  result / length(halfspaces)
}

# computes approximate halfspace depth ("Tukey depth")
# inputs/outputs: see evaluate_depth()
get_depth <- function(data, halfspaces) {
  # init in 1 to iteratively find minimal masses
  result <- rep(1, NROW(data))
  for (halfspace in halfspaces) {
    projections <- project_scalar(data, halfspace$normal)
    result <- pmin(result, ifelse(projections < halfspace$split,
      yes = halfspace$mass["below"],
      no = halfspace$mass["above"]
    ))
  }
  # Tukey depth is usually given in [0, nrow(training data)/2] not [0, 0.5]:
  result * nrow(attr(halfspaces, "train_data"))
}

################################################################################
# DEPRECATED ALTERNATIVE FUNCTIONS:
if (FALSE) {
  # Draw regular "grid" of <n> directions in <dims> dimensions
  #   See http://mathworld.wolfram.com/HyperspherePointPicking.html; eq. (5).
  # inputs:
  #   n: how many
  #   dims: required dimension
  # output:
  #   a <dims> x <n> matrix of <n> <dims>-dimensional directions
  get_directions <- function(n, dims = 2) {
    checkmate::assert_integerish(n, lower = 1)
    checkmate::assert_integerish(dims, lower = 2)

    n_each <- ceiling(n^(1 / dims))
    # use grid of gaussian quantiles along each axis instead of random
    # draws to set up a regular grid of points on <dims>-dimensional sphere:
    gauss_quantiles <- qnorm(p = seq(0, 1, l = n_each + 2)[-c(1, n_each + 2)])
    directions <- do.call(
      expand.grid,
      replicate(dims, gauss_quantiles, simplify = FALSE)
    )
    # drop 0-vector
    directions <- subset(directions, rowSums(directions^2) != 0)
    # reduce back to n vectors if n ^ (1 / d) is not integer
    directions <- directions[sample(1:nrow(directions), min(n, nrow(directions))), ]
    t(directions / sqrt(rowSums(directions^2)))
  }
}
