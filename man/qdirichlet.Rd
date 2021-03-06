% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dirichlet.r
\name{qdirichlet}
\alias{qdirichlet}
\title{Approxminate Dirichlet Quantiles}
\usage{
qdirichlet(X, alpha)
}
\arguments{
\item{X}{a matrix of marginal probabilities}

\item{alpha}{The Dirichlet parameters.  \code{alpha > 0}}
}
\value{
approximate dirichlet quantiles
}
\description{
qdirichlet is not an exact quantile function since the quantile of a
multivariate distribtion is not unique.
qdirichlet is also not the quantiles of the marginal distributions since
those quantiles do not sum to one
qdirichlet is the quantile of the underlying gamma functions, normalized
This has been tested to show that qdirichlet approximates the dirichlet
distribution well and creates the correct marginal means and variances
when using a latin hypercube sample
}
\examples{
qdirichlet(matrix(c(0.1, 0.1, 0.3), nrow=1, byrow=TRUE), c(2,3,4))
}
