\name{validate.target.cormat}
\alias{validate.target.cormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks the target correlation matrix
}
\description{
The function checks the validity of the values of pairwise correlations. Additionally, it checks positive definitiveness, symmetry and correct dimension.
}
\usage{
validate.target.cormat(plist, CorrMat, no.ord, no.norm)
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
Number of ordinal variables.
}
  \item{no.norm}{
Number of normal variables.
}
}
\details{
In addition to being positive definite and symmetric, the values of pairwise correlations in the target correlation matrix must also fall within the limits imposed by the marginal distributions of the ordinal variables.The function ensures that the supplied correlation matrix is valid for simulation. If a violation occurs, an error message is displayed that identifies the violation. The function returns a logical value \code{TRUE} when no such violation occurs.
}



\examples{

Sigma = diag(4)
Sigma[lower.tri(Sigma)] = c(0.42, 0.78, 0.29, 0.37, 0.14, 0.26)
Sigma = Sigma + t(Sigma)
diag(Sigma)=1

marginal = list( c(0.2, 0.5), c(0.4, 0.7, 0.9))
validate.target.cormat(marginal, Sigma, 2,2)
}

