% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/alignment.R
\name{alignSamples}
\alias{alignSamples}
\alias{alignSamples,ptrSet-method}
\title{Alignment between samples}
\usage{
alignSamples(
  X,
  ppmGroup = 70,
  fracGroup = 0.8,
  group = NULL,
  fracExp = 0.3,
  pValGreaterThres = 0.001,
  pValLessThres = 0,
  quantiUnit = c("ppb", "ncps", "cps")[1],
  bgCorrected = FALSE,
  dmzGroup = 0.001
)

\S4method{alignSamples}{ptrSet}(
  X,
  ppmGroup = 70,
  fracGroup = 0.8,
  group = NULL,
  fracExp = 0.3,
  pValGreaterThres = 0.001,
  pValLessThres = 0,
  quantiUnit = c("ppb", "ncps", "cps")[1],
  bgCorrected = FALSE,
  dmzGroup = 0.001
)
}
\arguments{
\item{X}{ptrSet already processed by the \code{\link[ptairMS]{detectPeak}} 
function}

\item{ppmGroup}{ppm maximal width for an mz group}

\item{fracGroup}{only variables present in \code{fracGroup} 
percent of at least one \code{group} will be kept (if 0 the filter is 
not applied)}

\item{group}{character: sampleMetadata data column name. If \code{NULL},
variables not present in \code{fracGroup} percent of samples will be deleted. 
Else, variables not present in \code{fracGroup} percent in in at least one 
group group will 
be removed.}

\item{fracExp}{fraction of samples which must have their p-value less than 
\code{pValGreaterThres} and \code{pValLessThres}}

\item{pValGreaterThres}{threshold of the p-value for the unilateral
testing that quantification (in cps) of expiration points are greater than 
the intensities in the background.}

\item{pValLessThres}{threshold of the p-value for the unilateral
testing that quantification (in cps) of expiration points are less than the 
intensities of the background.}

\item{quantiUnit}{ppb, ncps or cps}

\item{bgCorrected}{logical: should the peak table contain the background 
corrected values?}

\item{dmzGroup}{minimum mz width to be used for grouping the features 
(required for low masses)}
}
\value{
an \code{\link[Biobase]{ExpressionSet}} (Biobase object)
}
\description{
\code{AlignSamples} performs alignment between samples (i.e. the matching of 
variables between the peak lists within the \code{ptrSet} object) by using a 
kernel gaussian density (Delabriere et al, 2017).
This function returns an \code{\link[Biobase]{ExpressionSet}}, which contains 
the matrix of peak intensities, the sample metadata (borrowed from the
input ptrSet) and the variable metadata which contains the peak intensities in 
the background. 
Two filters may be applied to:
\itemize{
\item  keep only variables with a significant higher intensity in the 
expirations compared to the background (i.e., a p-value less than 
\code{pValGreaterThres}) for at least \code{fracExp} % of the samples
\item keep only variables which are detected in more 
than \code{fracGroup} percent of the samples (or \code{group})
}
If you do not want to apply those filters, set \code{fracGroup} to 0 and 
\code{pValGreaterThres} to 1.
}
\examples{
library(ptairData)
dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
exhaledPtrset <- createPtrSet(dir=dirRaw, 
setName="exhaledPtrset", mzCalibRef = c(21.022, 60.0525),
fracMaxTIC = 0.7, saveDir = NULL )
exhaledPtrset<-detectPeak(exhaledPtrset,mzNominal=c(21,60,79))
eset <- alignSamples(exhaledPtrset,pValGreaterThres=0.05,bgCorrected=TRUE)
Biobase::exprs(eset)
Biobase::fData(eset)
Biobase::pData(eset)
}
\references{
Delabriere et al., 2017
}
