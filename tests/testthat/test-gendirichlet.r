context("test-gendirichlet")

test_that("dGenDirichlet", {
  expect_equal(0, dGenDirichlet(c(0.1, 0.1, 0.1), c(0.7, 0.2, 0.1), c(2, 3, 4)))
  expect_equal(0, dGenDirichlet(c(0.1, 0.1, 0.1), c(1.1, 0.2, 0.1), c(2, 3, 4)))
  expect_equal(0, dGenDirichlet(c(0.1, 0.1, 0.1), c(0.7, 0.2, -0.1), c(2, 3, 4)))
  expect_equal(0, dGenDirichlet(c(0.1, 0.1, 0.1), c(0.7, 0.2, 0.1), c(2, -3, 4)))
})

test_that("dGenDirichletStd", {
  expect_equal(0, dGenDirichletStd(c(0.1, 0.1, 0.1), c(0.7, 0.2, 0.1), c(2, 3, 4)))
})

describe("rGenDirichlet", {
  set.seed(1976)
  n <- 10000
  p <- c(.4, .3, .2, .1)
  k <- c(10, 20, 30, 40)
  Y <- rGenDirichlet(n, p, k)
  it("The returned vector is on (0,1)", {
    expect_true(all(dim(Y) == c(n, length(p))))
    expect_true(all(Y >= 0 && Y <= 1))
  })
  it("The sum of the values in each draw is one.", {
    expect_equal(rowSums(Y), rep(1, n))
  })
  it("The mean of all the draws for each parameter is equal to the desired mean probabilities", {
    expect_equal(apply(Y, 2, mean), p, tol = 0.05)
  })
  it("Errors are generated on illegal inputs", {
    p <- c(.1, .7, .05, .15)
    expect_error(rGenDirichlet(-n, p, k))
    expect_warning(rGenDirichlet(10, c(.5, .5, .5), c(2, 3, 4)))
  })
  it("That rdirichlet and rGenDirichlet are equivalent", {
    Y <- matrix(c(.5, .3, .2), nrow = 1, ncol = 3)
    p <- c(.4, .35, .25)
    A <- rGenDirichlet(n, p, c(2,1.2,0))
    B <- rdirichlet(n, p*2)
    expect_equal(apply(A,2,var), apply(B,2,var), tolerance = 0.01)
  })
})

describe("qGenDirichlet", {
  set.seed(1976)
  n <- 100
  p <- c(.4, .3, .2, .1)
  k <- c(10, 20, 30, 40)
  Y <- rGenDirichlet(n, p, k)
  Z <- qGenDirichlet(Y, p, k)
  Z2 <- qGenDirichlet(Y, c(.4, 0, .35, .25), k)
  Z3 <- qGenDirichlet(Y, c(.4, .35, 0, .25), c(2, 3, 4, 5))
  W <- qGenDirichlet(Y, c(.35, .4, 0, .25), c(3, 2, 4, 5))
  it("The quantiles are on (0,1)", {
    expect_true(all( Z <= 1 & Z >= 0))
    expect_equal(rowSums(Z), rep(1, n))
    expect_true(all(Z2 <= 1 & Z2 >= 0))
    expect_equal(rowSums(Z2), rep(1, n))
    expect_true(all(Z2[,2] == 0))
  })
  it("Errors are generated for illegal input", {
    expect_error(rGenDirichlet(-1, p, k))
    expect_warning(qGenDirichlet(Y, c(0.5, 0.5, 0.5, 0.5), c(2, 3, 4, 5)))
  })
  expect_equal(apply(Z3, 2, mean), apply(W, 2, mean)[c(2,1,3,4)])
})

test_that("fit_genDirichlet", {
  tempTest <- function(desired_p, desired_k, type="mm", N)
  {
    set.seed(1976)
    Y <- rGenDirichlet(N, desired_p, desired_k)
    fit.genDirichlet(Y, type)
  }

  desired_p <- c(.4, .3, .2, .1)
  desired_k <- c(5, 6, 7, 10)
  temp <- tempTest(desired_p, desired_k, N = 10000)
  expect_equal(temp$k[1:3], desired_k[1:3], tolerance = 0.05*max(desired_k))
  expect_equal(temp$p, desired_p, tolerance = 0.01)

  desired_p <- c(.5, .3, .2)
  desired_k <- c(50, 60, 70)
  temp <- tempTest(desired_p, desired_k, N = 10000)
  expect_equal(temp$k[1:2], desired_k[1:2], tolerance = 0.05*max(desired_k))
  expect_equal(temp$p, desired_p, tolerance = 0.01)

  # this test failed before adding the sumx numerical accuracy code
  a <- c(1, 10, 20, 40, 100, 500)
  p1 <- 1/a / sum(1/a)
  a <- c(1, 10, 10, 50, 200, 200)
  p2 <- 1/a / sum(1/a)
  set.seed(1876)
  Y <- rdirichlet(100, p1*3)
  Z <- rdirichlet(100, p2*4)
  X <- as.matrix(rbind(Y, Z))
  temp <- fit.genDirichlet(X)
  expect_true(all(is.finite(temp$k)))
  expect_true(all(is.finite(temp$p)))

  desired_p <- c(.4, .3, .2, .1)
  desired_k <- c(5, 6, 7, 10)
  temp <- tempTest(desired_p, desired_k, "mm", N = 10000)
  expect_equal(temp$k[1:3], desired_k[1:3], tolerance = 0.05*max(desired_k))
  expect_equal(temp$p, desired_p, tolerance = 0.1)

  desired_p <- c(.5, .3, .2)
  desired_k <- c(50, 60, 70)
  temp <- tempTest(desired_p, desired_k, "mm", N = 10000)
  expect_equal(temp$k[1:2], desired_k[1:2], tolerance = 0.05*max(desired_k))
  expect_equal(temp$p, desired_p, tolerance = 0.1)

  expect_error(fit.dirichlet(4))
  expect_error(fit.dirichlet(matrix(c(1,2,3,4), nrow = 2), "not there yet"))

  desired_p <- c(.5, .3, .2)
  desired_k <- c(5, 6, 7)
  temp <- tempTest(desired_p, desired_k, "ml", N = 100)
  expect_equal(temp$k[1:2], desired_k[1:2], tolerance = 0.05*max(desired_k))
  expect_equal(temp$p, desired_p, tolerance = 0.05)
})

test_that("qApproxMarginalGenDirichlet", {
  X <- qApproxMarginalGenDirichlet(c(0.05,0.5,0.95), c(.4,.3,.2,.1), c(3,4,5,6))
  expect_equal(dim(X)[1], 3)
  expect_equal(dim(X)[2], 4)
  expect_true(all(X < 1 & X > 0))
  expect_equal(X[2,1], qbeta(0.5, 3*.4, 3*.6))
  expect_error(qApproxMarginalGenDirichlet(.5, .5, 2))
  expect_error(qApproxMarginalGenDirichlet(.5, c(.6,.4), c(3,4,5)))
  expect_error(qApproxMarginalGenDirichlet(.5, c(.6,.4), c(3,4)))
  expect_error(qApproxMarginalGenDirichlet(.5, c(.2,.3,.5), c(3,4,5)))
  expect_warning(qApproxMarginalGenDirichlet(c(0.05,0.5,0.95), c(.5,.5,.5,.5), c(3,4,5,6)))
})

test_that("calculateGenDirichletCV", {
  expect_equal(1, calculateGenDirichletCV(c(0.5, 0.3, 0.2), c(1, 1, 1))[2])
})

test_that("calculateConstantCVGenDirichletK", {
  expect_warning(k <- calculateConstantCVGenDirichletK(c(.5,.3,.2), 10))
  expect_true(length(k) == 3)
  expect_true(all(k[1:2] > 0))
  expect_equal(k[1], 10)

  a <- c(1, 10, 20, 40, 100, 500)
  p1 <- 1/a / sum(1/a)
  expect_warning(k <- calculateConstantCVGenDirichletK(p1, 3))
  expect_true(length(k) == 6)
  expect_true(all(k[1:5] > 0))
  expect_equal(k[1], 3)

  expect_warning(k <- calculateConstantCVGenDirichletK(c(.4,.3,.2,.1), 1))
  expect_equal(k, c(1, 2.2, .Machine$double.xmax, .Machine$double.xmax))

  k <- calculateConstantCVGenDirichletK(c(.3,.2,.15,.10,.05, .05, .05, .05, .05), 20)
  expect_equal(k[c(1,9)], c(20,  32.02777777777888))

  expect_warning(calculateConstantCVGenDirichletK(c(.5,.5,.5), 2))
})
