#' Generate a vector from a uniform distribution on an n-dimensional sphere
#'
#' @param dimension: dimension of the sphere
#'
#' @return an n-dimensional numeric vector
#' @export
pick_direction <- function(dimension) {
  direction_vector <- rnorm(dimension)
  direction_vector / sqrt(sum(direction_vector^2))
}

#' Sample a fraction of rows from a dataframe
#'
#' @param data: dataframe to be sampled
#' @param fraction: fraction of rows to be sampled
#'
#' @return dataframe with fraction * rows of the input data (rounded up)
#' @export
generate_subsample <- function(data, fraction) {
  sample_size <- ceiling(nrow(data) * fraction)
  sample_indices <- sample(seq_len(nrow(data)), sample_size)
  data[sample_indices, , drop = FALSE]
}

#' Perform a scalar projection of the columns of a dataframe onto a vector
#'
#' @param sample: dataframe containing the observations
#' @param direction: numeric vector to project the columns on
#'
#' @return dataframe with one column containing the scalar projections
#' @export
project_subsample <- function(sample, direction) {
  result <- matrix(0, nrow = nrow(sample), ncol = 1)
  for (row in seq_len(nrow(sample))) {
    projected_vector <- as.numeric(sample[row, ]) %*% direction
    result[row, 1] <- projected_vector
  }
  as.data.frame(result, row.names = rownames(sample))
}

#' Checks inputs for train_depth() function
#'
#' @param data: expected to be a data frame with non-missing numeric entries
#' @param n_halfspace: positive integer
#' @param subsample: float between 0 and 1
#' @param scope: positive number
#' @param seed: positive integer
#'
#' @export
check_train_inputs <- function(data, n_halfspace, subsample, scope, seed) {
  checkmate::assert_data_frame(data, types = c("integer", "integerish", 
                                               "double", "numeric"), 
                               any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_count(n_halfspace)
  checkmate::assert_count(seed)
  checkmate::assert_number(subsample, lower = 0, upper = 1)
  checkmate::assert_number(scope, lower = 0)
}

#' Training algorithm for halfspace mass, for details see 
#' http://scheipl.userweb.mwn.de//downloads/fortprog/ChenEtAl-HalfspaceMass-MachLearn2015.pdf
#'
#' @param data: dataframe, for which the half space mass should be computed
#' @param n_halfspace: number of halfspaces to fit
#' @param subsample: fraction of data to be used for fitting in each iteration
#' @param scope: scope parameter, see linked paper
#' @param seed: for reproducability of results 
#'
#' @return named list containing for each halfspace:
#' direction: sampled unit vector l_i
#' score: score value s_i
#' mass_left: mass "left" of hyperplan m^l_i
#' mass_right: mass "right" of hyperplan m^r_i
#' @export
train_depth <- function(data, n_halfspace, subsample = 1, scope = 1, seed) {
  check_train_inputs(data, n_halfspace, subsample, scope, seed)
  dimension <- ncol(data)
  directions <- list()
  score <- list()
  mass_left <- list()
  mass_right <- list()
  for (iteration in seq_len(n_halfspace)) {
    direction <- pick_direction(dimension)
    current_data <- generate_subsample(data, subsample)
    projected_data <- project_subsample(current_data, direction)
    current_max <- max(projected_data)
    current_min <- min(projected_data)
    current_mid <- mean(c(current_max, current_min))
    lower <- current_mid - scope / 2 * (current_max - current_min)
    upper <- current_mid + scope / 2 * (current_max - current_min)
    current_s <- runif(1, min = lower, max = upper)
    left <- projected_data < current_s
    right <- projected_data >= current_s
    directions[[iteration]] <- direction
    score[[iteration]] <- current_s
    mass_left[[iteration]] <- sum(left) / nrow(current_data)
    mass_right[[iteration]] <- sum(right) / nrow(current_data)
  }
  list(directions = directions, score = score, 
       mass_left = mass_left, mass_right = mass_right)
}

#' initialize the halfspace measure according to the method used
#'
#' @param used_method: the method used
#'
#' @return 0 for method "mass", Inf otherwise
#' @export
setup_hm <- function(used_method) {
  if (used_method == "mass") {
    return(0)
  }
  return(Inf)
}

#' update the value of the halfspace measure for the "mass" method
#'
#' @param hm_prev: numeric value of halfspace measure in previous iteration
#' @param current_projection: projected values of current data sample
#' @param halfspaces: list of halfspaces fitted with train_depth()
#' @param iteration: integer describing current iteration
#'
#' @return updated value for the halfspace measure as numeric
#' @export
update_hm_mass <- function(hm_prev, current_projection, halfspaces, iteration) {
  if (current_projection < halfspaces[['score']][[iteration]]) {
    hm_new <- hm_prev + halfspaces[['mass_left']][[iteration]]
  }
  else {
    hm_new <- hm_prev + halfspaces[['mass_right']][[iteration]]
  }
  hm_new
}

#' update the value of the halfspace measure for the "depth" method
#'
#' @param hm_prev: numeric value of halfspace measure in previous iteration
#' @param current_projection: projected values of current data sample
#' @param halfspaces: list of halfspaces fitted with train_depth()
#' @param iteration: integer describing current iteration
#'
#' @return updated value for the halfspace measure as numeric
#' @export
update_hm_depth <- function(hm_prev, current_projection, halfspaces, iteration){
  if (current_projection < halfspaces[['score']][[iteration]]) {
    hm_new <- min(hm_prev, halfspaces[['mass_left']][[iteration]])
  }
  else {
    hm_new <- min(hm_prev, halfspaces[['mass_right']][[iteration]])
  }
  hm_new
}

#' Checks inputs for evaluate_depth() function
#'
#' @param data: data frame or matrix with non-missing numeric entries
#' @param halfspaces: list returned by train_depth() function
#'
#' @export
check_test_inputs <- function(data, halfspaces) {
  # either dataframe OR matrix
  checkmate::assert(
    checkmate::check_data_frame(
    data, types = c("integer", "integerish", "double", "numeric"), 
    any.missing = FALSE, null.ok = FALSE
    ),
    checkmate::check_matrix(
    data, any.missing = FALSE, null.ok = FALSE)
    )
  checkmate::assert_list(halfspaces, len = 4, any.missing = FALSE)
}

#' Testing algorithm for halfspace mass, for details see 
#' http://scheipl.userweb.mwn.de//downloads/fortprog/ChenEtAl-HalfspaceMass-MachLearn2015.pdf
#'
#' @param data: dataframe, for which the half space mass should be computed
#' @param halfspaces: list of halfspaces fitted with train_depth
#' @param metric: either "mass" (Halfspace mass) or "depth" (Tukey depth)
#'
#' @return numeric vector containing the halfspace depth of each data point
#' according to the chosen metric
#' @export
evaluate_depth <- function(data, halfspaces, metric = c("mass", "depth")) {
  check_test_inputs(data, halfspaces)
  result <- numeric()
  n_halfspace <- length(halfspaces[["score"]])
  tryCatch({ used_method <- match.arg(metric) },
           error = function() {
             stop("Could not match input for metric. 
                  Available options are 'mass' und 'depth'.")
           })
  for (row in seq_len(nrow(data))) {
    halfspace_mass <- setup_hm(used_method)
    for (iteration in seq_len(n_halfspace)) {
      current_projection <- project_subsample(
        data[row, , drop = FALSE], halfspaces[['directions']][[iteration]])
      
      if (used_method == "mass") {
        halfspace_mass <- update_hm_mass(halfspace_mass, current_projection, 
                                         halfspaces, iteration)
      }
      if (used_method == "depth") {
        halfspace_mass <- update_hm_depth(halfspace_mass, current_projection, 
                                          halfspaces, iteration)
      }
    }
    result[row] <- halfspace_mass
  }
  result
}
