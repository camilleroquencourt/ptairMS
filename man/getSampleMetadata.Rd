% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetMethods.R
\name{getSampleMetadata}
\alias{getSampleMetadata}
\title{get sampleMetadata of a ptrSet}
\usage{
getSampleMetadata(set)
}
\arguments{
\item{set}{a ptrSet object}
}
\value{
a data.frame
}
\description{
get sampleMetadata of a ptrSet
}
\examples{
data(exhaledPtrset )
SMD<-getSampleMetadata(exhaledPtrset )
}
