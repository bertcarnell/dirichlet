################################################################################
# Program Name:  runit_testGenDirichlet.R
# Purpose:       Test dirichlet functions
# Charge:        G882602
# Author:        Rob Carnell
# Date:          December 2007
#
# R version:     >= 2.4.1
#
################################################################################

#/**
# /testedfunction rGenDirichlet
# /testedpurpose  a random draw from the multi-variate generalized Dirichlet distribution
# /test           1. The returned vector is on (0,1)
# /test           2. The sum of the values in each draw is one.
# /test           3. The mean of all the draws for each parameter is equal to the desired mean probabilities
#*/
test.rGenDirichlet <- function()
{
  set.seed(1976)
  n <- 100000
  p <- c(.4, .3, .2, .1)
  k <- c(10, 20, 30, 40)
  Y <- rGenDirichlet(n, p, k)
  checkTrue(all(dim(Y) == c(n, length(p))))
  checkTrue(all(Y >= 0 && Y <= 1))
  checkEqualsNumeric(rowSums(Y), rep(1, n))
  checkEqualsNumeric(apply(Y, 2, mean), p, tol=.05)
  p <- c(.1, .7, .05, .15)
  checkException(rGenDirichlet(Y, p, k), silent=TRUE)

  # check that it is the same as rdirichlet
  Y <- matrix(c(.5, .3, .2), nrow=1, ncol=3)
  p <- c(.4, .35, .25)
  A <- rGenDirichlet(n, p, c(2,1.2,0))
  B <- rdirichlet(n, p*2)
  checkEqualsNumeric(apply(A,2,var), apply(B,2,var), tolerance=0.01)
}

#/**
# /testedfunction qGenDirichlet
# /testedpurpose  Find the marginal quantiles of a generalized Dirichlet distribution
# /test           1. The quantiles are on (0,1)
# /test           2. Errors are generated for illegal input
#*/
test.qGenDirichlet <- function()
{
  set.seed(1976)
  n <- 100
  p <- c(.4, .3, .2, .1)
  k <- c(10, 20, 30, 40)
  Y <- rGenDirichlet(n, p, k)
  Z <- qGenDirichlet(Y, p, k)
  checkTrue(all( Z <= 1 & Z >= 0))
  checkEqualsNumeric(rowSums(Z), rep(1, n))

  Z <- qGenDirichlet(Y, c(.4, 0, .35, .25), k)
  checkTrue(all( Z <= 1 & Z >= 0))
  checkEqualsNumeric(rowSums(Z), rep(1, n))
  checkTrue(all(Z[,2] == 0))
  
  Z <- qGenDirichlet(Y, c(.4, .35, 0, .25), c(2, 3, 4, 5))
  W <- qGenDirichlet(Y, c(.35, .4, 0, .25), c(3, 2, 4, 5))
  checkEqualsNumeric(apply(Z, 2, mean), apply(W, 2, mean)[c(2,1,3,4)])
}

test.fit_genDirichlet <- function()
{
  tempTest <- function(desired_p, desired_k, type="mm", N)
  {
    set.seed(1976)
    Y <- rGenDirichlet(N, desired_p, desired_k)
    fit.genDirichlet(Y, type)
  }

  desired_p <- c(.4, .3, .2, .1)
  desired_k <- c(5, 6, 7, 10)
  temp <- tempTest(desired_p, desired_k, N=10000)
  checkEqualsNumeric(temp$k[1:3], desired_k[1:3], tolerance=0.05*max(desired_k))
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.01)

  desired_p <- c(.5, .3, .2)
  desired_k <- c(50, 60, 70)
  temp <- tempTest(desired_p, desired_k, N=10000)
  checkEqualsNumeric(temp$k[1:2], desired_k[1:2], tolerance=0.05*max(desired_k))
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.01)

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
  checkTrue(all(is.finite(temp$k)))
  checkTrue(all(is.finite(temp$p)))

  desired_p <- c(.4, .3, .2, .1)
  desired_k <- c(5, 6, 7, 10)
  temp <- tempTest(desired_p, desired_k, "mm", N=10000)
  checkEqualsNumeric(temp$k[1:3], desired_k[1:3], tolerance=0.05*max(desired_k))
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.1)

  desired_p <- c(.5, .3, .2)
  desired_k <- c(50, 60, 70)
  temp <- tempTest(desired_p, desired_k, "mm", N=10000)
  checkEqualsNumeric(temp$k[1:2], desired_k[1:2], tolerance=0.05*max(desired_k))
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.1)

  checkException(fit.dirichlet(4), silent=TRUE)
  checkException(fit.dirichlet(matrix(1,2,3,4, nrow=2), "not there yet"), silent=TRUE)

  desired_p <- c(.5, .3, .2)
  desired_k <- c(5, 6, 7)
  temp <- tempTest(desired_p, desired_k, "ml", N=100)
  checkEqualsNumeric(temp$k[1:2], desired_k[1:2], tolerance=0.05*max(desired_k))
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.05)
}

test.qApproxMarginalGenDirichlet <- function()
{
  X <- qApproxMarginalGenDirichlet(c(0.05,0.5,0.95), c(.4,.3,.2,.1), c(3,4,5,6))
  checkEqualsNumeric(dim(X)[1], 3)
  checkEqualsNumeric(dim(X)[2], 4)
  checkTrue(all(X < 1 & X >0))
  checkEqualsNumeric(X[2,1], qbeta(0.5, 3*.4, 3*.6))
  checkException(qApproxMarginalGenDirichlet(.5, .5, 2), silent=TRUE)
  checkException(qApproxMarginalGenDirichlet(.5, c(.6,.4), c(3,4,5)), silent=TRUE)
  checkException(qApproxMarginalGenDirichlet(.5, c(.6,.4), c(3,4)), silent=TRUE)
  checkException(qApproxMarginalGenDirichlet(.5, c(.2,.3,.5), c(3,4,5)), silent=TRUE)
}

test.calculateConstantCVGenDirichletK <- function()
{
  #expected warning
  suppressWarnings(
    k <- calculateConstantCVGenDirichletK(c(.5,.3,.2), 10)
  )
  checkTrue(length(k) == 3)
  checkTrue(all(k[1:2] > 0))
  checkEqualsNumeric(k[1], 10)

  a <- c(1, 10, 20, 40, 100, 500)
  p1 <- 1/a / sum(1/a)
  #expected warning
  suppressWarnings(
    k <- calculateConstantCVGenDirichletK(p1, 3)
  )
  checkTrue(length(k) == 6)
  checkTrue(all(k[1:5] > 0))
  checkEqualsNumeric(k[1], 3)
  
  #expected warning
  suppressWarnings(
    k <- calculateConstantCVGenDirichletK(c(.4,.3,.2,.1), 1)
  )
  checkEqualsNumeric(k, c(1, 2.2, .Machine$double.xmax, .Machine$double.xmax))

  k <- calculateConstantCVGenDirichletK(c(.3,.2,.15,.10,.05, .05, .05, .05, .05), 20)
  checkEqualsNumeric(k[c(1,9)], c(20,  32.02777777777888))
}
