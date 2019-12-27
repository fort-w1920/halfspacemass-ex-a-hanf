# Testing approximate vs exact Halfspace Depth in 2D

# make depth::depth accept matrices containing points u from test as rows:
hsdepth <- function (test, train) {
  apply(test, 1, depth::depth, x = train)
}

# wrap train and test steps into a single function:
my_hsdepth <- function (test, train, ...) {
  evaluate_depth(
    data = test,
    halfspaces= train_depth(train, ...),
    metric = "depth")
}


# generate list of test cases:
set.seed(121133)
angles <- seq(0, 2*pi, length = 101)[-101]

data_list <- list(
  circle = cbind(cos(angles), sin(angles)),
  gaussian = cbind(rnorm(200), rnorm(200)),
  clustered = cbind(rnorm(100, mean = rep(c(-1, 0, 1), l = 100), sd = .2),
                  rnorm(100, mean = rep(c(-1, 0.5, -1), l = 100), sd = .2)),
  grid = expand.grid(seq(-2, 2, l = 15), seq(-2, 2, l = 15))
)

layout(matrix(1:4, 2, 2)); lapply(data_list, plot)

# test approximate equivalence for the 4 datasets above:
# (values are scaled differently, so compare ranks and/or check correlation)
for (train in data_list) {
  approx_depth <- my_hsdepth(test = data_list$grid, train = train,
                             n_halfspace = 1e4, scope = 1, seed = 23)
  exact_depth <- hsdepth(test = data_list$grid, train = train)
  testthat::expect_true(
    cor(approx_depth,
        exact_depth) > 0.99
  )
  testthat::expect_equivalent(
    rank(approx_depth),
    rank(exact_depth),
    tol = nrow(data_list$grid)/20
  )
}





