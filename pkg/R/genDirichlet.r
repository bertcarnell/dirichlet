################################################################################
# Program Name:         Generalized Dirichlet
# Purpose:              package of tools for using dirichlet and generalized
#                       dirichlet distributions
# Charge:               G882602
# Author:               carnellr
# Date:                 December 2007
#
# R version:            2.6.0
#
################################################################################

dGenDirichlet <- function(x, p, k)
{
  ## probability density for the Generalized Dirichlet function
  ## x = vector of values at which the density is desired
  len <- length(x)
  stopifnot( length(p) == len & length(k) == len )
  stopifnot( len >= 2 )
  if ( abs(sum(x)-1) > .Machine$double.eps ) return(0)
  if ( any(p <= 0) || any(p >= 1) ) return(0)
  if ( any(k < 0) ) return(0)

  # put k and p in standard form for a and b
  theta <- numeric(len)
  a <- numeric(len)
  b <- numeric(len)
  # a_1 = p_1 * k_1
  # b_1 = (1-p_1)*k_1
  a[1] <- p[1]*k[1]
  b[1] <- (1-p[1])*k[1]
  theta[1] <- p[1]
  # a_2 = k_2*p_2/(1-p_1)
  # b_2 = k_2*(1-p_2/(1-p_1))
  # a_3 = k_3*p_3/(1-p_1)/(1-p_2)
  # b_3 = k_3*(1-p_3/(1-p_1)/(1-p_2))
  for ( i in 2:(len-1) )
  {
    temp <- p[i]
    theta[i] <- p[i]
    for ( j in 1:(i-1) )
    {
      temp <- temp / ( 1-theta[j])
      theta[i] <- theta[i] / (1-theta[j])
    }
    a[i] <- temp*k[i]
    b[i] <- (1-temp)*k[i]
  }
  if(any(c(a, b) < 0))
  {
    return(.Machine$double.xmax)
  }
  dGenDirichletStd(x, a, b)
}

################################################################################

dGenDirichletStd <- function(x, a, b)
{
  ## generalized dirichlet density function using standard notation
  len <- length(x)
  stopifnot( length(a) == len & length(b) == len )
  stopifnot( len >= 2 )
  if ( abs(sum(x) - 1) > 10*.Machine$double.eps ) return(0)

  # normalizing constant = Prod(i=1, len-1) 1 / Beta(a_i, b_i)
  norm <- sum(mapply(lbeta, a, b)[1:(len-1)])
  # part of the density
  #  for k=3, ... = x1^(a1-1) x2^(a2-2)
  #  = Prod(i=1, len-1) x_i^(a_i-1)
  temp <- sum( ( a - 1 ) * log(x) )
  # to ensure that b[len-1] - a[len] - b[len] = b[len-1] - 1
  a[len] <- 1
  b[len] <- 0
  # the rest of the density
  # for k=3, ... = (1-x_1)^(b_1-a_2-b_2) (1-x_1-x_2)^(b_2-1)
  for ( i in 1:(len-1) )
  {
    temp2 <- 1
    for ( j in 1:i )
    {
      temp2 <- temp2 - x[j]
    }
    temp <- temp + (b[i]-a[i+1]-b[i+1]) * log(temp2)
  }
  # return from the log scale
  exp(temp - norm)
}

################################################################################

rGenDirichlet <- function(n, p, k)
{
  # parameterization according to Rust's "Generalized Dirichlet Distribution"
  # n is the number of samples
  # p is the vector of m probabilities
  # k is the vector of m k values

  lenp <- length(p)
  stopifnot(lenp == length(k))
  if(abs(sum(p) - 1) > 1E-12)
  {
    warning("renormalizing p")
    p <- p / sum(p)
  }
  ind <- order(p, decreasing=TRUE)
  if(!all(ind==1:lenp))
  {
    stop("p must be in decreasing order")
  }

  # theta is the individual p_i divided by the sum of that p_i and smaller p_i
  theta <- p / (1 + p - cumsum(p))

  B <- mapply(function(TH,K) rbeta(n, K*TH, K*(1-TH)),
              theta[1:(lenp-1)], k[1:(lenp-1)])

  X <- matrix(0, nrow=n, ncol=lenp)
  for(i in 1:(lenp-1))
  {
    if(i == 1)
    {
      X[,i] = B[,i]
    } else
    {
      X[,i] <- (1 - apply(X[,1:i], 1, sum))*B[,i]
    }
  }
  X[,lenp] <- 1 - apply(X, 1, sum)
  X
}

################################################################################

qGenDirichlet <- function(X, p, k)
{
  stopifnot(all(p >= 0 ))
  lenp <- length(p)
  lenk <- length(k)
  stopifnot(lenp == lenk)
  if(abs(sum(p) - 1) > 1E-12)
  {
    warning("renormalizing p")
    p <- p / sum(p)
  }
  colX <- dim(X)[2]
  sims <- dim(X)[1]
  stopifnot( sims > 0 & lenp == colX )
  stopifnot(all(!is.na(X)) & all(!is.na(p)) & all(!is.na(k)))

  ind <- which(p > 0)
  lenind <- length(ind)

  stopifnot(lenind > 1)

  indo <- order(p[ind], decreasing=TRUE)
  p[ind] <- p[ind][indo]
  k[ind] <- k[ind][indo]

  theta <- p / rev( cumsum( rev(p) ) )

  Y <- matrix(0, nrow=sims, ncol=lenind)

  for( i in 1:(lenind-1) )
  {
    if ( i == 1 )
    {
      Y[,i] <- qbeta(X[,ind[i]], k[ind[i]]*theta[ind[i]],
                     k[ind[i]]*(1 - theta[ind[i]]))
    } else if ( i == 2 )
    {
      Y[,2] <- (1 - Y[,1]) *
               qbeta(X[,ind[i]], k[ind[i]]*theta[ind[i]],
                     k[ind[i]]*(1 - theta[ind[i]]))
    } else
    {
      Y[,i] <- (1 - rowSums( Y[,1:(i-1)] )) *
               qbeta(X[,ind[i]], k[ind[i]]*theta[ind[i]],
                     k[ind[i]]*(1 - theta[ind[i]]))
    }
  }
  Y[,lenind] <- 1 - rowSums(Y)

  X[,] <- 0
  X[,ind[indo]] <- Y
  X
}

################################################################################

fit.genDirichlet <- function(X, type="mm")
{
  stopifnot(is.matrix(X))

  # first get method of moments estimates for the mm return or ml starting values
  n <- nrow(X)
  nparms <- ncol(X)

  sumx <- 0
  sump <- 0
  p <- numeric(nparms-1)
  k <- numeric(nparms-1)
  for (i in 1:(nparms-1))
  {
    # check for numerical accuracy in sumx
    ind <- which( abs(sumx-1) < .Machine$double.eps )
    sumx[ind] <- 1 - .Machine$double.eps
    # Convert data to beta form
    b    <- X[,i] / (1-sumx)
    sumx <- sumx + X[,i]

    m <- mean(b)
    v  <- var(b)
    # Calculate method-of-moments estimate of k parameter
    # k = p(1-p) / ((p(1-p)/(k+1)) - 1
    k[i] <- (m*(1-m)/v) - 1
    # Calculate method-of-moments estimate of p parameter
    p[i] <- (1-sump) * m
    sump <- sump + p[i]
  }
  k <- c(k, k[nparms-1])
  ind <- which(k < 0)
  if (length(ind) > 0)
  {
    warning("k estimates may not be accurate")
    k[ind] <- .Machine$double.eps
  }
  p <- c(p, 1-sump)
  ind <- which(p < 0)
  if (length(ind) > 0)
  {
    warning("p estimates may not be accurate")
    p[ind] <- .Machine$double.eps
  }
  if ( type == "mm" )
  {
    return(list(k=k, p=p))
  } else if ( type == "ml" )
  {
    # optimization function
    g <- function(pk)
    {
      # bring in parameters in a vector of p's then k's
      p_g <- pk[1:nparms]
      k_g <- pk[(nparms+1):(2*nparms)]
      # bring in parameters on a -Inf to Inf scale for the optimizer
      #  logit transform back to a probability (guaranteed to be on [0,1]
      p_g <- exp(p_g)/(1+exp(p_g))
      #  log transform back to a strictly positive number
      k_g <- exp(k_g)
      # check to ensure that the probability vector sums to one.
      if ( abs(sum(p_g) - 1 ) > .Machine$double.eps ) return(.Machine$double.xmax)
      # return the sum of the negative log likelihood
      -sum(log(apply(X, 1, dGenDirichlet, p=p_g, k=k_g)))
    }

    # start the optimizer at the mean values of the paramters
    m <- apply(X, 2, mean)
    # use the method of moments k estimate
    # both parameter vectors start out on the transformed scale
    o <- optim(c(log(m/(1-m)), log(k)), g, control=list(maxit=10000))

    if ( o$convergence != 0 ) stop(o$message)

    # extract the paramters and transform them to the proper scales
    p <- o$par[1:nparms]
    p <- exp(p)/(1+exp(p))
    k <- o$par[(nparms+1):(2*nparms)]
    k <- exp(k)
    # ensure the probability is strictly normalized
    p <- p / sum(p)

    return(list(k=k, p=p))
  } else stop("type not recognized")
}

################################################################################

qApproxMarginalGenDirichlet <- function(pquant, p, k)
{
  # pquants - vector of probabilities for desired quantiles
  # p - vector of probabilities for the dirichlet distribution
  # k - vector of k parameters

  lenp <- length(p)

  if (lenp == 1) stop("more than one probability required")
  if (lenp != length(k)) stop("length p and length k must be equal")
  stopifnot(lenp > 2)

  # if the p's are not probabilities, change them into probabilities
  if (sum(p) != 1)
  {
    warning("normalizing p in qApproxMarginalGenDirichlet")
    p <- p / sum(p)
  }

  # if the p's are not in order, error
  ind <- order(p, decreasing=TRUE)
  stopifnot(all(p == p[ind]))

  # theta is the individual p_i divided by the sum of that p_i and smaller p_i
  theta <- p / (1 + p - cumsum(p))

  margVariance <- numeric(lenp)
  tempVector <- (k*(1-theta) + 1) / (k+1)

  margVariance[1] <- p[1] * (1 - p[1]) / (k[1] + 1)
  for (i in 2:(lenp-1))
  {
    temp <- prod(tempVector[1:(i-1)])
    margVariance[i] <- p[i] * ((k[i] * theta[i] + 1) / (k[i] + 1) * temp - p[i])
  }

  margVariance[lenp] <- p[lenp] * (prod(tempVector[1:(lenp-1)]) - p[lenp])

  ind <- which(margVariance <= 0)
  if(length(ind) > 0) margVariance[ind] <- .Machine$double.eps

  # compute beta parameters using method of moments
  a <- p * (p * (1 - p) / margVariance - 1)
  b <- (1 - p) * (p * (1 - p) / margVariance - 1)

  quants <- mapply(function(x,y) qbeta(pquant, x, y), a, b)

  return(quants)
}

################################################################################

calculateGenDirichletCV <- function(p, k)
{
  # calculate the coefficient of variation for the generalized dirichlet
  # not to be used as a standalone function
  lenp <- length(p)

  stopifnot(lenp == length(k))

  margVariance <- numeric(lenp)
  margCV <- numeric(lenp)

  theta <- p / (1 + p - cumsum(p))
  tempVector <- (k * (1 - theta) + 1) / (k + 1)
  margVariance[1] <- p[1]*(1-p[1])/(k[1]+1)

  for (i in 2:(lenp-1))
  {
    temp <- prod(tempVector[1:(i-1)])
    margVariance[i] <- p[i]*((k[i]*theta[i] + 1)/(k[i]+1)*temp - p[i])
  }

  margVariance[lenp] <- p[lenp] * (prod(tempVector[1:(lenp-1)]) - p[lenp])

  ind <- which(margVariance <= 0)
  if (length(ind) > 0) margVariance[ind] <- .Machine$double.eps

  margCV <- sqrt(margVariance) / p

  # so that the optimizer does not return k < 0
  ind <- which(k < 0)
  if (length(ind) > 0) margCV[ind] <- .Machine$double.xmax

  return(margCV)
}

################################################################################

calculateConstantCVGenDirichletK <- function(p, most.likely.k)
{
  lenp <- length(p)
  stopifnot(lenp > 2)

  # if the p's are not probabilities, change them into probabilities
  if (sum(p)!= 1)
  {
    warning("Renormalizing p in calculateConstantCVGenDiricletK")
    p <- p/sum(p)
  }

  # if the p's are not in order, order them
  ind <- order(p, decreasing=TRUE)
  p <- p[ind]

  # calculate CV of most likely branch
  varMostLikely <- p[1] * (1 - p[1]) / (most.likely.k + 1)
  CVMostLikely <- sqrt(varMostLikely) / p[1]

  k <- rep(1, lenp)
  k[1] <- most.likely.k
  theta <- p / (1 + p - cumsum(p))
  tempVector <- (k * (1 - theta) + 1) / (k + 1)

  for (i in 2:(lenp-1))
  {
    # must re-calculate temp based on new k's
    tempVector <- (k * (1 - theta) + 1) / (k + 1)
    temp <- prod(tempVector[1:(i-1)])
    k[i] <- (p[i] * (CVMostLikely^2 + 1) - temp) /
            (temp * theta[i] - p[i] * (CVMostLikely^2 + 1))
    if (k[i] < 0 || !is.finite(k[i]))
    {
      warning("Negative or Infinite k at index ", i, " k = ", k)
      # in this situation, set the k to the largest machine number
      k[i] <- .Machine$double.xmax
    }
  }
  k[lenp] <- k[lenp-1]
  
  ktemp <- numeric(length(k))
  ktemp[ind] <- k
  ktemp
}

################################################################################

calculateSlopeGenDirichletK <- function(p, k0, k_m)
{
  # Bio method of a k and a slope for k to calculate the k vector
  lenp <- length(p)
  stopifnot(lenp > 2)

  # if the p's are not probabilities, change them into probabilities
  if(sum(p)!= 1) p <- p / sum(p)

  # if the p's are not in order, order them
  ind <- order(p, decreasing=TRUE)
  p <- p[ind]

  k <- numeric(lenp)

  k[1] <- k0 + k_m

  for (i in 2:lenp)
  {
    k[i] <- k0 + k_m * 1 / (1 - sum(p[1:(i-1)]))
  }

  return(k)
}


