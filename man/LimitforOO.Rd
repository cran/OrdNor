\name{LimitforOO}
\alias{LimitforOO}
\alias{Limit_forOO}
\title{
Finds the feasible correlation range for a pair of ordinal variables
}
\description{
The function computes the lower and upper bounds of correlation between two ordinal variables via the method of Demirtas and Hedeker (2011).
}
\usage{
LimitforOO(pvec1, pvec2)
Limit_forOO(pvec1, pvec2) #deprecated
}
\arguments{
  \item{pvec1}{
A vector of marginal probabilities for the first ordinal variable.
}
  \item{pvec2}{
A vector of marginal probabilities for the second ordinal variable.
}
}
\value{
A vector of two elements. The first element is the lower bound and the second element is the upper bound.
}
\references{
Demirtas, H., Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. The American Statistician, \bold{65(2)}, 104-109.
}
\examples{
	pvec1 = cumsum( c(0.30, 0.40) )
	pvec2=cumsum(c(0.4, 0, 0.3) ) # The second category is skipped in this setting
 	LimitforOO(pvec1, pvec2)
}