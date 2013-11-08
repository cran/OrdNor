\name{IntermediateOO}
\alias{IntermediateOO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes intermediate correlations for ordinal-ordinal pairs before ordinalization
}
\description{
This function computes the correlation of normal-normal pairs before ordinalizing both components.
}
\usage{
IntermediateOO(plist, OOCorrMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plist}{
A list of probability vectors corresponding to each ordinal variable. The i-th element of \code{plist} is a vector of the cumulative probabilities defining the marginal distribution of the i-th ordinal component of the multivariate variables. If the i-th ordinal variable has k categories, the i-th vector of the \code{plist} will contain k-1 probability values. The k-th element is implicitly 1.
}
  \item{OOCorrMat}{
A matrix of pairwise target correlations between ordinal variables. It is a symmetric square matrix whose diagonal elements are 1. 
}
}

\value{
A pairwise correlation matrix of intermediate correlations for ordinal variables.
}


\seealso{
\code{\link{IntermediateON}}, \code{\link{cmat.star}}
}

\examples{
	no.ord=3
	no.norm =4
	n = 200
	q=no.ord + no.norm
	set.seed(12345)
	
	Sigma = diag(q)
	Sigma[lower.tri(Sigma)] = runif( (q*(q-1)/2),-0.4,0.4 )
	Sigma = Sigma + t(Sigma)
	diag(Sigma)=1
	Sigma=as.matrix( nearPD(Sigma,corr = TRUE, keepDiag = TRUE)$mat )

	marginal = list( 0.3, cumsum( c(0.30, 0.40) ), cumsum(c(0.4, 0.2, 0.3) ) )
 
	OOCorrMat = Sigma[1:3, 1:3]
 	IntermediateOO(marginal, OOCorrMat)
}