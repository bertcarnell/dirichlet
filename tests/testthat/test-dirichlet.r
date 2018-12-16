context("test-dirichlet")

describe("ddirichlet", {
  it("A two parameter Dirichlet is equivalent to a beta", {
    expect_equal(ddirichlet(c(0.3, 0.7), c(2, 3)), dbeta(0.3, 2, 3))
  })
  it("Dirichlet which is uniform across the parameter space has the appropriate density", {
    # uniform on 0,1 with height 2
    expect_equal(ddirichlet(c(.1,.8,.1), c(1,1,1)), 2)
    expect_equal(ddirichlet(c(.1,.6,.3), c(1,1,1)), 2)
  })
  it("That the function returns values equal to an independent calculation", {
    expect_equal(ddirichlet(c(.1,.4,.5), c(1,2,3)),
                       gamma(1 + 2 + 3)/gamma(1)/gamma(2)/gamma(3) *
                         .1^(1 - 1)*.4^(2 - 1)*.5^(3 - 1))
    expect_equal(ddirichlet(c(.1,.2,.7), c(4,5,8)),
                       gamma(4 + 5 + 8)/gamma(4)/gamma(5)/gamma(8) *
                         .1^(4 - 1)*.2^(5 - 1)*.7^(8 - 1))
  })
  it("That values of a outside the Dirichlet parameters yield values of zero", {
    expect_equal(ddirichlet(c(1,1,1), c(1,1,1)), 0)
    expect_equal(ddirichlet(c(.1,.1,.1), c(1,1,1)), 0)
  })
})

describe("rdirichlet", {
  set.seed(1976)
  n <- 100000
  p <- c(.4, .3, .2, .1)
  Y <- rdirichlet(n, p)
  expect_true(all(dim(Y) == c(n, length(p))))

  alpha <- c(4,8,2)
  Y2 <- rdirichlet(n, alpha)
  expect_true(all(dim(Y2) == c(n, length(alpha))))
  it("The returned vector is on (0,1)",{
    expect_true(all(Y >= 0 && Y <= 1))
    expect_true(all(Y2 >= 0 && Y2 <= 1))
  })
  it("The sum of the values in each draw is one.", {
    expect_equal(rowSums(Y), rep(1, n))
    expect_equal(rowSums(Y2), rep(1, n))
  })
  it("The mean of all the draws for each parameter is equal to the desired mean probabilities", {
    expect_equal(apply(Y, 2, mean), p, tol = 0.05)
    expect_equal(apply(Y2, 2, mean), alpha/sum(alpha), tol = 0.05)
  })
})

describe("qdirichlet",{
  set.seed(1976)
  X <- matrix((1:12)/13, ncol = 3, nrow = 4)
  X2 <- X
  X2[2,3] <- NA
  A <- matrix(
    c(0.01605770, 0.2673003, 0.7166420,
      0.02892119, 0.2698036, 0.7012752,
      0.03886065, 0.2672579, 0.6938815,
      0.04514385, 0.2554320, 0.6994241), nrow = 4, ncol = 3, byrow = TRUE)

  B <- matrix(
    c(0.02191580, 0, 0.9780842,
      0.03960742, 0, 0.9603926,
      0.05303454, 0, 0.9469655,
      0.06063094, 0, 0.9393691), nrow = 4, ncol = 3, byrow = TRUE)
  assert_that(all(abs(rowSums(A) - rep(1,4)) < 1E-7))
  assert_that(all(abs(rowSums(B) - rep(1,4)) < 1E-7))
  it("The result of the function matches known output within 1E-7", {
    expect_equal(qdirichlet(X, c(1,2,3)), A, tolerance = 1E-7)
    expect_equal(qdirichlet(X, c(1,0,3)), B, tolerance = 1E-7)
  })
  it("Errors are generated for illegal input", {
    expect_error(qdirichlet(c(2,3), c(1,2,3)), silent = TRUE)
    expect_error(qdirichlet(matrix(nrow = 3, ncol = 2), c(1,2,3)), silent = TRUE)
    expect_error(qdirichlet(X2, c(1,2,3)), silent = TRUE)
    expect_error(qdirichlet(X, c(1,NA,3)), silent = TRUE)
  })
  it("Other", {
    Z <- qdirichlet(rdirichlet(100, c(3,4,5)), c(3,4,5))
    expect_true(all( Z <= 1 & Z >= 0))
    expect_equal(rowSums(Z), rep(1, 100))

    # test for numerical instability in the qgamma calculation
    expect_warning(Z <- qdirichlet(matrix(c(.001, .1, .899,
                             .001, .1, .899), nrow = 2, ncol = 3, byrow = TRUE),
                    c(.001, 2, 2)))
    expect_false(any(is.nan(Z)))

    A <- matrix(runif(20*1000), nrow = 1000, ncol = 20)
    A <- A / rep(rowSums(A), 20)
    expect_warning(Z <- qdirichlet(A, c(.001, .01, .1, 1, rep(2, 16))))
    expect_false(any(is.nan(Z)))

    expect_warning(qdirichlet(matrix(c(.001, .1, .899,
                                       .001, .1, .899), nrow = 2, ncol = 3, byrow = TRUE),
                              c(.001, 2, 2)), silent = TRUE)
  })
})

test_that("fit.dirichlet", {
  tempTest <- function(desired_p, desired_k, type="mm")
  {
    set.seed(1976)
    Y <- rdirichlet(1000, desired_k*desired_p)
    fit.dirichlet(Y, type)
  }

  desired_p <- c(.4, .3, .2, .1)
  temp <- tempTest(desired_p, 5)
  expect_equal(temp$most.likely.k, 5, tolerance = 0.5)
  expect_equal(temp$weighted.k, 5, tolerance = 0.5)
  expect_equal(temp$p, desired_p, tolerance = 0.05)

  desired_p <- c(.5, .3, .2)
  temp <- tempTest(desired_p, 100)
  expect_equal(temp$most.likely.k, 100, tolerance = 0.5)
  expect_equal(temp$weighted.k, 100, tolerance = 0.5)
  expect_equal(temp$p, desired_p, tolerance = 0.05)

  desired_p <- c(.4, .3, .2, .1)
  temp <- tempTest(desired_p, 5, "ml")
  expect_equal(temp$k, 5, tolerance = 0.5)
  expect_equal(temp$p, desired_p, tolerance = 0.05)

  desired_p <- c(.5, .3, .2)
  temp <- tempTest(desired_p, 100, "ml")
  expect_equal(temp$k, 100, tolerance = .1)
  expect_equal(temp$p, desired_p, tolerance = 0.05)

  expect_error(fit.dirichlet(4), silent = TRUE)
  expect_error(fit.dirichlet(matrix(1,2,3,4, nrow = 2), "not there yet"), silent = TRUE)
})
