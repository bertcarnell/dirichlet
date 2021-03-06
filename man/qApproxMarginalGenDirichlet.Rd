% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genDirichlet.r
\name{qApproxMarginalGenDirichlet}
\alias{qApproxMarginalGenDirichlet}
\title{Approximate Marginal Quantiles for the Generalized Dirichlet}
\usage{
qApproxMarginalGenDirichlet(pquant, p, k)
}
\arguments{
\item{pquant}{vector of probabilities for desired quantiles}

\item{p}{vector of probabilities for the dirichlet distribution}

\item{k}{vector of k parameters}
}
\value{
blah
}
\description{
Approximate Marginal Quantiles for the Generalized Dirichlet
}
\examples{
qApproxMarginalGenDirichlet(c(.1, .2, .3), c(0.5, 0.3, 0.2), c(5,5,3))
}
