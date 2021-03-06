% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment.R
\name{imputeMat}
\alias{imputeMat}
\title{Impute missing values on an matrix set from an ptrSet}
\usage{
imputeMat(X, ptrSet, quantiUnit)
}
\arguments{
\item{X}{the peak table matrix with missing values}

\item{ptrSet}{processed by detectPeak function}

\item{quantiUnit}{the unit of the quantities in the matrix \code{X} (ppb, 
cps or ncps)}
}
\value{
the same matrix as in input, with missing values imputing
}
\description{
Imputing missing values by returning back to the raw data and fitting the 
peak shape function on the noise / residuals
}
\examples{
data(exhaledPtrset)
eSet <- alignSamples(exhaledPtrset,pValGreaterThres=0.05,fracGroup=0)
X <-Biobase::exprs(eSet)
X <- imputeMat(X,exhaledPtrset,quantiUnit='ppb')
}
