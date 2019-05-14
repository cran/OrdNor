\name{ordinalize}
\alias{ordinalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordinalizes the standard normal variable 
}
\description{
The function transforms the standard normal variable to an ordinal variable with a specified probability for each category.
}
\usage{
ordinalize(pvec, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvec}{
	A vector of probabilities for an ordinal variable. The i-th element of the \code{pvec} is the cumulative probability defining the marginal distribution of the ordinal variable. If the variable has k categories, the i-th element of pvec will contain k-1 probabilities. The k-th element is implicitly 1.
	}
	  \item{z}{
	A vector of samples from the standard normal distribution.
	}
}
\value{
A vector of ordinalized variates.
}

