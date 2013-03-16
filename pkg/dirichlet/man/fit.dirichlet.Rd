\name{fit.dirichlet}
\alias{fit.dirichlet}
\alias{fit.genDirichlet}
\title{Estimate Dirichlet Parameters}
\description{
  Estimate Dirichlet and Generalized Dirichlet distribution parameters.
}
\usage{
fit.dirichlet(X, type)
fit.genDirichlet(X, type)
}
\arguments{
  \item{X}{Matrix with columns corresponding to the entries of \code{a} and rows corresponding to the number of simulations in the hypercube.}
  \item{type}{\code{mm} for Method of Moments estimation or \code{ml} for Maximum Likelihood estimation.}
}
\value{
  A list giving the \code{k} value or vector of \code{k} values for the Generalized
  Dirichlet and vector of \code{p} estimates.  The normal Dirichlet parameterization can be obtained by \code{k*p}.
  The \code{fit.dirichlet} function returns two estimates of \code{k}.  One based on the most likely
  parameter in the distribution, the other based on a weighted mean of the \code{k} estimates for each parameter.
}
\author{Rob Carnell}
\seealso{\code{\link{rdirichlet}, \code{\link{rGenDirichlet}}}}
\examples{
  Z <- rdirichlet(100, c(4,3,2))
  fit.dirichlet(Z)
  fit.dirichlet(Z, "ml")
  
  Z2 <- rGenDirichlet(100, c(0.7, 0.2, 0.1), c(10, 20, 20))
  fit.genDirichlet(Z2)
}
\keyword{distribution}


