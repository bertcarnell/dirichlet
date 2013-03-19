\name{dirichlet}
\alias{rdirichlet}
\alias{ddirichlet}
\title{Dirichlet functions}
\description{
  pdf and random deviates of the Dirichlet distribution.
}
\usage{
rdirichlet(n, alpha, allowZero=FALSE)
ddirichlet(x, alpha)
}
\arguments{
  \item{alpha}{The Dirichlet parameters.  \code{alpha > 0} unless \code{allowZero=TRUE}}
  \item{n}{number of observations}
  \item{x}{vector of quantiles.  For \code{x < 0 || x > 1} or if the sum of \code{x} is not equal to 1, zero will be returned for the density}
  \item{allowZero}{Indicator that zeros are allowed in the alpha vector and should be handled appropriately}
}
\value{
  ddirichlet gives the density and rdirichlet generates random deviates.  
}
\author{Rob Carnell}
\seealso{\code{\link{fit.dirichlet}}}
\examples{
  rdirichlet(10, c(4,3,2))
  ddirichlet(c(.2, .7, .1), c(4, 3, 2))
  rdirichlet(5, c(4,3,0,2), allowZero=TRUE)
  (ddirichlet(c(-0.2, 1.2), c(2,2)) == 0)
}
\keyword{distribution}

