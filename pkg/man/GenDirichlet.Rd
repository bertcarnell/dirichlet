\name{GenDirichlet}
\alias{rGenDirichlet}
\title{Generalized Dirichlet functions}
\description{
  Random deviates of the Generalized Dirichlet distribution.
}
\usage{
rGenDirichlet(n, p, k)
}
\arguments{
  \item{p}{The Dirichlet normalized parameters.  The Dirichlet parameters may contain zeros if the branch to which a parameter is referring is not allowed.}
  \item{k}{k * p = alpha (the Dirichlet parameters.}
  \item{n}{number of observations. If length(n) > 1, the length is taken to be the number required.}
}
\value{
  rGenDirichlet generates random deviates.
}
\author{Rob Carnell}
\seealso{\code{\link{rdirichlet}}}
\examples{
  rGenDirichlet(10, c(4,3,2)/9, c(2, 2, 2))
}
\keyword{distribution}

