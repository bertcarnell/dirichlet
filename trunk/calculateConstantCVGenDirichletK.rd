\name{calculateConstantCVGenDirichletK}
\alias{calculateConstantCVGenDirichletK}
\alias{calculateGenDirichletCV}
\title{Calculate the K Parameter for a Generalized Dirichlet Distribution \\
       using a Constant Coefficient of Variation Procedure
}
\description{
  Given the k value of the most likely probability in a generalzed Dirichlet
  distribution and the probabilties of the branches, this function will
  calculate the remaining generalized Dirichlet k values such that each of the
  branches has a constant coefficent of variation.
}
\usage{
calculateConstantCVGenDirichletK(p, most.likely.k)
calculateGenDirichletCV(p, k)
}
\arguments{
  \item{p}{A vector of probabilities}
  \item{k}{A vector of k values for the generalized Dirichlet.  Note that the last k in the vector is not used.}
  \item{most.likely.k}{The k value for the most likely branch}
}
\value{
  \code{calculateConstantCVGenDirichletK} returns a vector of k values.  A warning
  is issued if one of the k values must be set to the maximum machine double
  \code{.Machine$double.xmax}.  This does not affect future sampling from
  the generalized Dirichlet distribution.
  \code{calculateGenDirichletCV} returns a vector of coefficients of varaiation
  for each of the branches.
}
\author{Rob Carnell}
\seealso{\code{\link{rGenDirichlet}} for taking random samples from the generalized Dirichlet distribution.}
\examples{
  p <- c(.3,.2,.15,.10,.05, .05, .05, .05, .05)
  k <- calculateConstantCVGenDirichletK(p, 20)
  k
  # should return a constant cv
  cv <- calculateGenDirichletCV(p, k)
  cv

  # warning expected
  p <- c(.4,.3,.2,.1)
  suppressWarnings(
    k <- calculateConstantCVGenDirichletK(p, 1)
  )
  # will include .Machine$double.xmax as noted by the suppressed warning
  k
  # should return a constant cv
  cv <- calculateGenDirichletCV(p, k)
  cv
  
  X <- rGenDirichlet(1000, p, k)
  # should return values closed to the initial p values
  apply(X, 2, mean)
}
\keyword{distribution}

