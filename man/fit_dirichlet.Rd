% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dirichlet.r, R/genDirichlet.r
\name{fit_dirichlet}
\alias{fit_dirichlet}
\alias{fit.dirichlet}
\alias{fit.genDirichlet}
\title{Fit Dirichlet parameters}
\usage{
fit.dirichlet(X, type = "mm")

fit.genDirichlet(X, type = "mm")
}
\arguments{
\item{X}{samples from a dirichlet distribution}

\item{type}{\code{mm} for Method of Moments estimation or \code{ml} for Maximum Likelihood estimation.}
}
\value{
A list giving the \code{k} value or vector of \code{k} values for the Generalized
Dirichlet and vector of \code{p} estimates.  The normal Dirichlet parameterization can be obtained by \code{k*p}.
The \code{fit.dirichlet} function returns two estimates of \code{k}.  One based on the most likely
parameter in the distribution, the other based on a weighted mean of the \code{k} estimates for each parameter.
}
\description{
Fit Dirichlet parameters
}
\examples{
fit.dirichlet(matrix(c(.1, .9, .2, .8, .3, .7), nrow=3, ncol=2, byrow=TRUE), type="mm")
Z <- rdirichlet(100, c(4,3,2))
fit.dirichlet(Z)
fit.dirichlet(Z, "ml")
  Z2 <- rGenDirichlet(100, c(0.7, 0.2, 0.1), c(10, 20, 20))
  fit.genDirichlet(Z2)
}
