\name{dirichlet}
\alias{rdirichlet}
\alias{ddirichlet}
\title{Dirichlet functions}
\description{
  pdf and random deviates of the Dirichlet distribution.
}
\usage{
rdirichlet(n, alpha)
ddirichlet(x, alpha)
}
\arguments{
  \item{alpha}{The Dirichlet parameters.  The Dirichlet parameters may contain zeros.}
  \item{n}{number of observations. If length(n) > 1, the length is taken to be the number required.}
  \item{x}{vector of quantiles.}
}
\value{
  ddirichlet gives the density and rdirichlet generates random deviates.  
}
\author{Rob Carnell}
\seealso{\code{\link{fit.dirichlet}}}
\examples{
  require(lhs)
  rdirichlet(10, c(4,3,2))
  ddirichlet(c(.2, .7, .1), c(4, 3, 2))
}
\keyword{distribution}

