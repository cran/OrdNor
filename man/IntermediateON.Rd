\name{IntermediateON}
\alias{IntermediateON}

\title{
 Computes intermediate correlations for ordinal-normal pairs before ordinalization
}
\description{
The function computes the intermediate correlation values of pairwise correlations between ordinal and normal variables.
}
\usage{
IntermediateON(plist, ONCorrMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plist}{
A list of probability vectors corresponding to each ordinal variable. The i-th element of \code{plist} is a vector of the cumulative probabilities defining the marginal distribution of the i-th ordinal component of the multivariate variables. If the i-th ordinal variable has k categories, the i-th vector of the \code{plist} will contain k-1 probability values. The k-th element is implicitly 1.
}
  \item{ONCorrMat}{
A matrix of pairwise target correlations between ordinal and normal variables. This is a submatrix of the overall correlation matrix, and it is pertinent to the ordinal-normal part. Hence, the matrix may or may not be square. Even when it is square, it may not be symmetric.
}
}

\value{
A pairwise correlation matrix of intermediate correlations.
}
\seealso{
\code{\link{IntermediateOO}}, \code{\link{cmat.star}}
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
 
	ONCorrMat = Sigma[4:7, 1:3]
 	IntermediateON(marginal, ONCorrMat)
}