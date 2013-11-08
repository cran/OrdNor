\name{valid.limits}
\alias{valid.limits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Computes the lower and upper bounds of correlation in the form of two matrices
}
\description{
The function computes the lower and upper bounds for the target correlation based on the marginal probabilities.  
}
\usage{
valid.limits(plist, no.ord, no.norm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plist}{
A list of probability vectors corresponding to each ordinal variable. The i-th element of \code{plist} is a vector of the cumulative probabilities defining the marginal distribution of the i-th ordinal component of the multivariate variables. If the i-th ordinal variable has k categories, the i-th vector of the \code{plist} will contain k-1 probability values. The k-th element is implicitly 1.
}
  \item{no.ord}{
Number of ordinal variables.
}
  \item{no.norm}{
Number of normal variables.
}
}
\details{
The function returns a list of two matrices. The \code{lower} contains the lower bounds and the \code{upper} contains the upper bounds of the feasible correlations.
}

\examples{
marginal = list( c(0.2, 0.5), c(0.4, 0.7, 0.9))
valid.limits (marginal, 2,2)
}

