################################################################################
# Program Name:         Dirichlet
# Purpose:              package of tools for using dirichlet and generalized
#                       dirichlet distributions
# Charge:               G882602
# Author:               carnellr
# Date:                 July 2007
#
# R version:            2.4.1
#
################################################################################


ddirichlet <- function(x, alpha) {
  ## probability density for the Dirichlet function
  ## x = vector of values at which the density is desired
  ## alpha = vector of dirichlet parameters ( k * p )

  stopifnot( length(alpha) == length(x) )
  stopifnot( all(alpha > 0) )
  # the density outside the support space for x is zero
  if (any(x <= 0)) return(0)
  if (any(x > 1)) return(0)
  if ( abs(sum(x)-1) > 10*.Machine$double.eps ) return(0)

  s <- sum( ( alpha - 1 ) * log(x) )
  exp( s - ( sum( lgamma(alpha) ) - lgamma( sum(alpha) ) ) )
}

################################################################################

rdirichlet <- function(n, alpha, allowZero=FALSE) {
  ## pick n random deviates from the Dirichlet function with shape
  ## parameters alpha.  alpha = 0 is allowed and produces 0 for those options

  lena <- length(alpha)
  stopifnot( lena > 1 && n > 0 )
  if (!allowZero)
  {
    stopifnot(all(alpha > 0))
    x <- matrix( rgamma(lena*n, alpha), ncol=lena, byrow=TRUE)
    sm <- x %*% rep(1, lena)
    return(x / as.vector(sm))
  } else
  {
    stopifnot(all(alpha >= 0))
    ind <- which(alpha != 0)
    X <- sapply(alpha[ind], function(x) rgamma(n, x, 1))
    Xsum <- apply(X, 1, sum)
    Y <- apply(X, 2, "/", Xsum)
    Z <- matrix(0, nrow=n, ncol=length(alpha))
    Z[,ind] <- Y
    return(Z)
  }
}

################################################################################

qdirichlet <- function(X, alpha) {
  # qdirichlet is not an exact quantile function since the quantile of a
  #  multivariate distribtion is not unique
  # qdirichlet is also not the quantiles of the marginal distributions since
  #  those quantiles do not sum to one
  # qdirichlet is the quantile of the underlying gamma functions, normalized
  # This has been tested to show that qdirichlet approximates the dirichlet
  #  distribution well and creates the correct marginal means and variances
  #  when using a latin hypercube sample
  lena <- length(alpha)
  stopifnot(is.matrix(X))
  sims <- dim(X)[1]
  stopifnot(dim(X)[2] == lena)
  if(any(is.na(alpha)) || any(is.na(X)))
    stop("NA values not allowed in qdirichlet")

  Y <- matrix(0, nrow=sims, ncol=lena)
  ind <- which(alpha != 0)
  for(i in ind)
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

fit.dirichlet <- function(X, type="mm")
{
  stopifnot(is.matrix(X))

  # get method of moments estimates to return to use as starting values for the
  #  maximum likelihood method
  p <- apply(X, 2, mean)
  ind <- which.max(p)
  v <- apply(X, 2, var)
  temp <- p*(1-p)/v - 1
  #alpha <- p[ind]*temp
  #beta <- (1-p[ind])*temp
  #k <- alpha + beta

  if( type == "mm" )
  {
    return(list(most.likely.k=temp[ind], weighted.k=weighted.mean(temp, v), p=p))
  } else if ( type == "ml" )
  {
    f <- function(a)
    {
      # pass in a on the log scale so that it is guaranteed to be a strictly
      #  positive number when transformed back
      a <- exp(a)
      -sum(log(apply(X, 1, ddirichlet, alpha=a)))
    }

    #  use starting values on the log scale
    o <- optim(log(temp[ind]*p), f)
    
    if ( o$convergence != 0 ) stop(o$message)

    # the probabilities are the normalized a values
    p <- exp(o$par)
    p <- p / sum(p)
    
    k <- exp(o$par)[ind]/p[1]

    return(list(k=k, p=p))
  } else stop("type not recognized")
}

################################################################################

qMarginalDirichlet <- function(p, alpha) {
  # the marginals of the Dirichlet distribution are beta distributions
  #  with paramters alpha = k*p, beta = k(1-p)
  #  where the Dirichlet alpha = k*p
  #  p is the vector of desired percentiles
  # returns a length p x length alpha matrix

  stopifnot(all(p > 0 & p < 1))
  stopifnot(all(alpha > 0))
  
  p_dirichlet <- alpha / sum(alpha)
  k <- alpha[1] / p[1]
  
  sapply(p_dirichlet, function(x) qbeta(p, k*x, k*(1-x)))
}

################################################################################

calculateDirichletCV <- function(p, k) {
  margCV <- sqrt(p*(1-p)/(k+1))/p

  return(margCV)
}

################################################################################
