% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment.R
\name{makeSubGroup}
\alias{makeSubGroup}
\title{Use in align function. return a peak group thanks to kernel gaussian density 
in a peak matrix.}
\usage{
makeSubGroup(subpeakl, den, plim)
}
\arguments{
\item{subpeakl}{a matrix with mz, ppb, background and group column.}

\item{den}{the kernel gaussian density estimated on subpeakl}

\item{plim}{the limit of a peak in the density of the group who will pe fromed}
}
\value{
the sub peakgroup
}
\description{
Use in align function. return a peak group thanks to kernel gaussian density 
in a peak matrix.
}
\keyword{internal}
