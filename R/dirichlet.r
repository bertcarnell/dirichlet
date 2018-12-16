# Copyright 2018 Robert Carnell


#' Dirichlet Functions
#'
#' @description  pdf and random deviates of the Dirichlet distribution.
#'
#' @param x vector of quantiles.  For \code{x < 0 || x > 1} or if the sum of \code{x} is not equal to 1, zero will be returned for the density
#' @param alpha The Dirichlet parameters.  \code{alpha > 0} unless \code{allowZero=TRUE}
#' @param n number of observations
#' @param allowZero Indicator that zeros are allowed in the alpha vector and should be handled appropriately
#'
#' @return ddirichlet gives the density and rdirichlet generates random deviates.
#' @name dirichlet
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' rdirichlet(10, c(4,3,2))
#' ddirichlet(c(.2, .7, .1), c(4, 3, 2))
#' rdirichlet(5, c(4,3,0,2), allowZero=TRUE)
#' (ddirichlet(c(-0.2, 1.2), c(2,2)) == 0)
ddirichlet <- function(x, alpha) {
  ## probability density for the Dirichlet function
  ## x = vector of values at which the density is desired
  ## alpha = vector of dirichlet parameters ( k * p )

  assert_that(length(alpha) == length(x), msg = "the length of x and alpha must be equal")
  assert_that(all(alpha > 0), msg = "All alphas must be > 0")
  # the density outside the support space for x is zero
  if (any(x <= 0))
    return(0)
  if (any(x > 1))
    return(0)
  if (abs(sum(x) - 1) > 10*.Machine$double.eps)
    return(0)

  s <- sum((alpha - 1 ) * log(x))
  return(exp(s - (sum(lgamma(alpha)) - lgamma(sum(alpha)))))
}

################################################################################

#' @rdname dirichlet
#' @export
rdirichlet <- function(n, alpha, allowZero=FALSE) {
  ## pick n random deviates from the Dirichlet function with shape
  ## parameters alpha.  alpha = 0 is allowed and produces 0 for those options

  lena <- length(alpha)
  assert_that(lena > 1 && n > 0, msg = "The length of alpha must be > 1 and n > 0")
  if (!allowZero)
  {
    assert_that(all(alpha > 0), msg = "All alpha must be >0 when allowZer0 = FALSE")
    x <- matrix( rgamma(lena*n, alpha), ncol = lena, byrow = TRUE)
    sm <- x %*% rep(1, lena)
    return(x / as.vector(sm))
  } else
  {
    assert_that(all(alpha >= 0), msg = "All alpha must be >= 0")
    ind <- which(alpha != 0)
    X <- sapply(alpha[ind], function(x) rgamma(n, x, 1))
    Xsum <- apply(X, 1, sum)
    Y <- apply(X, 2, "/", Xsum)
    Z <- matrix(0, nrow = n, ncol = length(alpha))
    Z[,ind] <- Y
    return(Z)
  }
}

################################################################################


#' Approxminate Dirichlet Quantiles
#'
#' qdirichlet is not an exact quantile function since the quantile of a
#' multivariate distribtion is not unique.
#' qdirichlet is also not the quantiles of the marginal distributions since
#' those quantiles do not sum to one
#' qdirichlet is the quantile of the underlying gamma functions, normalized
#' This has been tested to show that qdirichlet approximates the dirichlet
#' distribution well and creates the correct marginal means and variances
#' when using a latin hypercube sample
#' @param X a matrix of marginal probabilities
#' @param alpha The Dirichlet parameters.  \code{alpha > 0}
#'
#' @return approximate dirichlet quantiles
#' @export
#'
#' @examples
#' qdirichlet(matrix(c(0.1, 0.1, 0.3), nrow=1, byrow=TRUE), c(2,3,4))
qdirichlet <- function(X, alpha) {
  lena <- length(alpha)
  assert_that(is.matrix(X), msg = "X must be a matrix")
  sims <- dim(X)[1]
  assert_that(dim(X)[2] == lena, msg = "The number of columns of X must equal the length of alpha")
  if (any(is.na(alpha)) || any(is.na(X)))
    stop("NA values not allowed in qdirichlet")

  Y <- matrix(0, nrow = sims, ncol = lena)
  ind <- which(alpha != 0)
  for (i in ind)
  {
    # add check to trap numerical instability in qgamma
    # start to worry if alpha is less than 1.0
    if (alpha[i] < 1)
    {
      # look for places where NaN will be returned by qgamma
      nanind <- which(pgamma(.Machine$double.xmin, alpha[i], 1) >= X[,i])
      # if there are such places
      if (length(nanind) > 0)
      {
        # set the output probability to near zero
        Y[nanind,i] <- .Machine$double.xmin
        # calculate the rest
        Y[-nanind,i] <- qgamma(X[-nanind,i], alpha[i], 1)
        warning("at least one probability set to the minimum machine double")
      } else {
        Y[,i] <- qgamma(X[,i], alpha[i], 1)
      }
    } else {
      Y[,i] <- qgamma(X[,i], alpha[i], 1)
    }
  }
  Y <- Y / rowSums(Y)
  return(Y)
}

################################################################################

#' Fit Dirichlet parameters
#'
#' @param X samples from a dirichlet distribution
#' @param type \code{mm} for Method of Moments estimation or \code{ml} for Maximum Likelihood estimation.
#'
#' @return A list giving the \code{k} value or vector of \code{k} values for the Generalized
#' Dirichlet and vector of \code{p} estimates.  The normal Dirichlet parameterization can be obtained by \code{k*p}.
#' The \code{fit.dirichlet} function returns two estimates of \code{k}.  One based on the most likely
#' parameter in the distribution, the other based on a weighted mean of the \code{k} estimates for each parameter.
#'
#' @name fit_dirichlet
#' @export
#'
#' @examples
#' fit.dirichlet(matrix(c(.1, .9, .2, .8, .3, .7), nrow=3, ncol=2, byrow=TRUE), type="mm")
#' Z <- rdirichlet(100, c(4,3,2))
#' fit.dirichlet(Z)
#' fit.dirichlet(Z, "ml")
fit.dirichlet <- function(X, type="mm")
{
  assert_that(is.matrix(X), msg = "X must be a matrix")

  # get method of moments estimates to return to use as starting values for the
  #  maximum likelihood method
  p <- apply(X, 2, mean)
  ind <- which.max(p)
  v <- apply(X, 2, var)
  temp <- p*(1 - p)/v - 1
  #alpha <- p[ind]*temp
  #beta <- (1-p[ind])*temp
  #k <- alpha + beta

  if (type == "mm")
  {
    return(list(most.likely.k = temp[ind], weighted.k = weighted.mean(temp, v), p = p))
  } else if (type == "ml")
  {
    f <- function(a)
    {
      # pass in a on the log scale so that it is guaranteed to be a strictly
      #  positive number when transformed back
      a <- exp(a)
      return(-sum(log(apply(X, 1, ddirichlet, alpha = a))))
    }

    #  use starting values on the log scale
    o <- optim(log(temp[ind]*p), f)

    if (o$convergence != 0)
      stop(o$message)

    # the probabilities are the normalized a values
    p <- exp(o$par)
    p <- p / sum(p)

    k <- exp(o$par)[ind]/p[1]

    return(list(k = k, p = p))
  } else stop("type not recognized")
}

################################################################################

#' Quantiles of the marginal beta distributions in a dirichlet distribution
#'
#' @param p the probability of the desired quantile
#' @param alpha the vector of dirichlet parameters
#'
#' @return the vector of marginal quantiles
#' @export
#'
#' @examples
#' qMarginalDirichlet(c(0.1, 0.2, 0.4), c(2,3,4))
qMarginalDirichlet <- function(p, alpha) {
  # the marginals of the Dirichlet distribution are beta distributions
  #  with paramters alpha = k*p, beta = k(1-p)
  #  where the Dirichlet alpha = k*p
  #  p is the vector of desired percentiles
  # returns a length p x length alpha matrix

  assert_that(all(p > 0 & p < 1), msg = "p must be > 0 and < 1")
  assert_that(all(alpha > 0), msg = "all alpha must be > 0")

  p_dirichlet <- alpha / sum(alpha)
  k <- alpha[1] / p[1]

  sapply(p_dirichlet, function(x) qbeta(p, k*x, k*(1 - x)))
}

################################################################################

#' The coefficient of variation of the marginal betas
#'
#' @param alpha the dirichlet parameters
#'
#' @return the marginal CVs
#' @export
#'
#' @examples
#' calculateDirichletCV(c(5, 4, 2))
calculateDirichletCV <- function(alpha) {
  p <- alpha / sum(alpha)
  k <- sum(alpha)
  margCV <- sqrt(p*(1 - p)/(k + 1))/p

  return(margCV)
}
