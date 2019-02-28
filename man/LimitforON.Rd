\name{LimitforON}
\alias{LimitforON}
\alias{Limit_forON}
\title{
Finds the feasible correlation range for a pair of ordinal and normal variable
}
\description{
The function computes the lower and upper bounds of a pairwise correlation between an ordinal and a normal variable via the method of Demirtas and Hedeker (2011).
}
\usage{
LimitforON(pvec1)
Limit_forON(pvec1) #deprecated
}
\arguments{
  \item{pvec1}{
A vector of marginal probabilities for an ordinal variable of the pair.
}
}
\value{
A vector of two elements. The first element is the lower bound and the second element is the upper bound.
}
\references{
Demirtas, H., Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. The American Statistician, \bold{65(2):104-109}.
}
\examples{
	pvec = cumsum( c(0.30, 0.40) )
 	LimitforON(pvec)
}