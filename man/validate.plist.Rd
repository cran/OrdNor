\name{validate.plist}
\alias{validate.plist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Checks the validity of ordinal probabilities
}
\description{
The function checks the validity of the probability vectors of the ordinal variables. It verifies that the elements in the vectors are cumulative probabilities and the values are between 0 and 1. It also checks a number of vectors within the list matches the specified number of ordinal variables.
}
\usage{
validate.plist(plist, no.ord)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	  \item{plist}{
A list of probability vectors corresponding to each ordinal variable. The i-th element of \code{plist} is a vector of the cumulative probabilities defining the marginal distribution of the i-th ordinal component of the multivariate variables. If the i-th ordinal variable has k categories, the i-th vector of the \code{plist} will contain k-1 probability values. The k-th element is implicitly 1.
	}
	  \item{no.ord}{
		Number of ordinal variables.
	}
}
\details{
The function returns error message if there are any violations. No message is displayed for the correct specifications.
}

\examples{
	marginal = list( c(0.2, 0.5), c(0.4, 0.7, 0.9))
	validate.plist(marginal, 2)
}

