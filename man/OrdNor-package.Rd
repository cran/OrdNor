\name{OrdNor-package}
\alias{OrdNor-package}
\alias{OrdNor}
\docType{package}
\title{
Concurrently generates ordinal and normal variables with a specified correlation matrix and marginal distributions
}
\description{
The package implements a procedure for concurrently generating samples from ordinal and normal random variables with a pre-specified correlation matrix and marginal distributions. It is accomplished by first calculating an intermediate correlations (\code{\link{cmat.star}}) which is used to generate a sample from multivariate normal distribution. Then, the first few components (corresponding to number of ordinal variables) are ordinalized. The resulting data are composed of a mixture of ordinal and normal variables that conform with a pre-specified marginal distributions and correlation structure. The function \code{\link{valid.limits}} returns the lower and upper bounds of the correlation coefficients of ordinal-ordinal (OO) and ordinal-normal (ON) pairs given their marginal distributions, i.e. returns the range of feasible pairwise correlations. The function \code{\link{validate.target.cormat}} checks the validity of the values of pairwise correlations. Additionally, it checks positive definitiveness, symmetry and correct dimension. The engine function \code{\link{genOrdNor}} generates mixed data in accordance with the specified marginal and correlational quantities.
}

\details{
\tabular{ll}{
Package: \tab OrdNor\cr
Type: \tab Package\cr
Version: \tab 2.2.3\cr
Date: \tab 2021-03-05\cr
License: \tab GPL\cr
}

}
\author{
Anup Amatya, Hakan Demirtas, Ran Gao

Maintainer: Ran Gao <rgao8@uic.edu>

}

