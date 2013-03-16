################################################################################
# Program Name:  runit_dirichlet.R
# Purpose:       Test dirichlet functions
# Charge:        G882602
# Author:        Rob Carnell
# Date:          December 2006
#
# R version:     >= 2.4.1
#
################################################################################

#/**
# /testedfunction ddirichlet
# /testedpurpose  The density of the Dirichlet distribution
# /test           1. That a two parameter Dirichlet is equivalent to a beta
# /test           2. That a Dirichlet which is uniform across the parameter space has the appropriate density
# /test           3. That the function returns values equal to an independent calculation
# /test           4. That values of a outside the Dirichlet parameters yield values of zero
#*/
test.ddirichlet <- function()
{
  # generalization of the beta
  checkEqualsNumeric(ddirichlet(c(0.3, 0.7), c(2, 3)), dbeta(0.3, 2, 3))
  # uniform on 0,1 with height 2
  checkEqualsNumeric(ddirichlet(c(.1,.8,.1), c(1,1,1)), 2)
  checkEqualsNumeric(ddirichlet(c(.1,.6,.3), c(1,1,1)), 2)
  # non-uniform
  checkEqualsNumeric(ddirichlet(c(.1,.4,.5), c(1,2,3)),
                     gamma(1+2+3)/gamma(1)/gamma(2)/gamma(3)*
                     .1^(1-1)*.4^(2-1)*.5^(3-1))
  checkEqualsNumeric(ddirichlet(c(.1,.2,.7), c(4,5,8)),
                     gamma(4+5+8)/gamma(4)/gamma(5)/gamma(8)*
                     .1^(4-1)*.2^(5-1)*.7^(8-1))
  checkEqualsNumeric(ddirichlet(c(1,1,1), c(1,1,1)), 0)
  checkEqualsNumeric(ddirichlet(c(.1,.1,.1), c(1,1,1)), 0)
}

#/**
# /testedfunction rdirichlet
# /testedpurpose  a random draw from the multi-variate Dirichlet distribution
# /test           1. The returned vector is on (0,1)
# /test           2. The sum of the values in each draw is one.
# /test           3. The mean of all the draws for each parameter is equal to the desired mean probabilities
#*/
test.rdirichlet <- function()
{
  set.seed(1976)
  n <- 100000
  p <- c(.4, .3, .2, .1)
  Y <- rdirichlet(n, p)
  checkTrue(all(dim(Y) == c(n, length(p))))
  checkTrue(all(Y >= 0 && Y <= 1))
  checkEqualsNumeric(rowSums(Y), rep(1, n))
  checkEqualsNumeric(apply(Y, 2, mean), p, tol=.05)

  alpha <- c(4,8,2)
  Y <- rdirichlet(n, alpha)
  checkTrue(all(dim(Y) == c(n, length(alpha))))
  checkTrue(all(Y >= 0 && Y <= 1))
  checkEqualsNumeric(rowSums(Y), rep(1, n))
  checkEqualsNumeric(apply(Y, 2, mean), alpha/sum(alpha), tol=.05)
}

#/**
# /testedfunction qdirichlet
# /testedpurpose  Find the marginal quantiles of a dirichlet distribution
# /test           1. The result of the function matches known output within 1E-7
# /test           2. Errors are generated for illegal input
#*/
test.qdirichlet <- function()
{
  set.seed(1976)
  X <- matrix((1:12)/13, ncol=3, nrow=4)
  X2 <- X; X2[2,3] <- NA
  A <- matrix(
       c(0.01605770, 0.2673003, 0.7166420,
         0.02892119, 0.2698036, 0.7012752,
         0.03886065, 0.2672579, 0.6938815,
         0.04514385, 0.2554320, 0.6994241), nrow=4, ncol=3, byrow=TRUE)

  B <- matrix(
       c(0.02191580, 0, 0.9780842,
         0.03960742, 0, 0.9603926,
         0.05303454, 0, 0.9469655,
         0.06063094, 0, 0.9393691), nrow=4, ncol=3, byrow=TRUE)
         
  checkEqualsNumeric(rowSums(A), rep(1,4), tolerance=1E-7)
  checkEqualsNumeric(rowSums(B), rep(1,4), tolerance=1E-7)

  checkException(qdirichlet(c(2,3), c(1,2,3)), silent=TRUE)
  checkException(qdirichlet(matrix(nrow=3, ncol=2), c(1,2,3)), silent=TRUE)
  checkException(qdirichlet(X2, c(1,2,3)), silent=TRUE)
  checkException(qdirichlet(X, c(1,NA,3)), silent=TRUE)
  checkEqualsNumeric(qdirichlet(X, c(1,2,3)), A, tolerance=1E-7)
  checkEqualsNumeric(qdirichlet(X, c(1,0,3)), B, tolerance=1E-7)

  Z <- qdirichlet(rdirichlet(100, c(3,4,5)), c(3,4,5))
  checkTrue(all( Z <= 1 & Z >= 0))
  checkEqualsNumeric(rowSums(Z), rep(1, 100))

  # test for numerical instability in the qgamma calculation
  Z <- qdirichlet(matrix(c(.001, .1, .899,
                      .001, .1, .899), nrow=2, ncol=3, byrow=TRUE),
                      c(.001, 2, 2))
  checkTrue(any(is.nan(Z)) == FALSE)

  A <- matrix(runif(20*1000), nrow=1000, ncol=20)
  A <- A / rep(rowSums(A), 20)
  Z <- qdirichlet(A, c(.001, .01, .1, 1, rep(2, 16)))
  checkTrue(any(is.nan(Z)) == FALSE)
  
  options(warn=2)
  checkException(qdirichlet(matrix(c(.001, .1, .899,
                      .001, .1, .899), nrow=2, ncol=3, byrow=TRUE),
                      c(.001, 2, 2)), silent=TRUE)
  options(warn=1)

}

test.fit_Dirichlet <- function()
{
  tempTest <- function(desired_p, desired_k, type="mm")
  {
    set.seed(1976)
    Y <- rdirichlet(1000, desired_k*desired_p)
    fit.dirichlet(Y, type)
  }
  
  desired_p <- c(.4, .3, .2, .1)
  temp <- tempTest(desired_p, 5)
  checkEqualsNumeric(temp$most.likely.k, 5, tolerance=0.5)
  checkEqualsNumeric(temp$weighted.k, 5, tolerance=0.5)
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.05)

  desired_p <- c(.5, .3, .2)
  temp <- tempTest(desired_p, 100)
  checkEqualsNumeric(temp$most.likely.k, 100, tolerance=0.5)
  checkEqualsNumeric(temp$weighted.k, 100, tolerance=0.5)
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.05)

  desired_p <- c(.4, .3, .2, .1)
  temp <- tempTest(desired_p, 5, "ml")
  checkEqualsNumeric(temp$k, 5, tolerance=0.5)
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.05)

  desired_p <- c(.5, .3, .2)
  temp <- tempTest(desired_p, 100, "ml")
  checkEqualsNumeric(temp$k, 100, tolerance=.1)
  checkEqualsNumeric(temp$p, desired_p, tolerance=0.05)

  checkException(fit.dirichlet(4), silent=TRUE)
  checkException(fit.dirichlet(matrix(1,2,3,4, nrow=2), "not there yet"), silent=TRUE)
}

