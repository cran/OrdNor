\name{genOrdNor}
\alias{genOrdNor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a data set with ordinal and normal variables
}
\description{
The function simulates a sample of size \code{n} from a multivariate ordinal and normal variables with correlation matrix \code{CorrMat} and pre-specified marginals.
}
\usage{
genOrdNor(n, plist, cmat.star, mean.vec, sd.vec, no.ord, no.norm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of rows
}
  \item{plist}{
A list of probability vectors corresponding to each ordinal variable. The i-th element of \code{plist} is a vector of the cumulative probabilities defining the marginal distribution of the i-th ordinal component of the multivariate variables. If the i-th ordinal variable has k categories, the i-th vector of the \code{plist} will contain k-1 probability values. The k-th element is implicitly 1.
}
  \item{cmat.star}{
The intermediate correlation matrix obtained from \code{\link{cmat.star}} function.
}
  \item{mean.vec}{
A vector of means for the normal variables.
}
  \item{sd.vec}{
A vector of standard deviations for the normal variables.
}
  \item{no.ord}{
Number of ordinal variables.
}
  \item{no.norm}{
Number of normal variables.
}
}


\value{
A matrix of size \eqn{n \times (no.ord + no.norm)}, of which first \code{no.ord} are ordinal variables.
}

\author{
Anup Amatya, Hakan Demirtas
}


\seealso{
\code{\link{cmat.star}}, \code{\link{validate.target.cormat}}, \code{\link{validate.plist}}
}
\examples{
	Sigma = diag(4)
	Sigma[lower.tri(Sigma)] = c(0.42, 0.78, 0.29, 0.37, 0.14, 0.26)
	Sigma = Sigma + t(Sigma)
	diag(Sigma)=1

	marginal = list( c(0.2, 0.5), c(0.4, 0.7, 0.9))
	cmat= cmat.star(marginal, Sigma, 2, 2)  
	mean.vec = c(2,4)
	sd.vec = c(0.5, 1.5)
	Y=genOrdNor(10000,marginal, cmat, mean.vec, sd.vec, 2, 2)
	cor(Y)
}

