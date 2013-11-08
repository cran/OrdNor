\name{cmat.star}
\alias{cmat.star}

\title{
Computes the correlation of intermediate multivariate normal data before subsequent ordinalization
}
\description{
The function computes an intermediate correlation matrix which leads to the target correlation matrix after ordinalization of
the samples generated from a multivariate normal distribution with the intermediate correlation matrix. 
}
\usage{
cmat.star(plist, CorrMat, no.ord, no.norm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plist}{
A list of probability vectors corresponding to each ordinal variable. The i-th element of \code{plist} is a vector of the cumulative probabilities defining the marginal distribution of the i-th ordinal component of the multivariate variables. If the i-th ordinal variable has k categories, the i-th vector of the \code{plist} will contain k-1 probability values. The k-th element is implicitly 1.
}
  \item{CorrMat}{
The target correlation matrix which must be positive definite and within the valid limits.
}
  \item{no.ord}{
Number of ordinal variables in the data.
}
  \item{no.norm}{
Number of normal variables in the data.
}
}

\value{
An intermediate correlation matrix of size \eqn{(no.ord + no.norm) \times (no.ord + no.norm) }

}



\seealso{
 \code{\link{validate.target.cormat}}
}
\examples{
	Sigma = diag(4)
	Sigma[lower.tri(Sigma)] = c(0.42, 0.78, 0.29, 0.37, 0.14, 0.26)
	Sigma = Sigma + t(Sigma)
	diag(Sigma)=1

	marginal = list( c(0.2, 0.5), c(0.4, 0.7, 0.9))
	cmat.star(marginal, Sigma, 2, 2) 
}

