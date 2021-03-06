% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetMethods.R
\name{exportSampleMetada}
\alias{exportSampleMetada}
\title{export sampleMetadata}
\usage{
exportSampleMetada(set, saveFile)
}
\arguments{
\item{set}{a ptrSet object}

\item{saveFile}{a file path in tsv extension where the data.frame will be 
exported}
}
\value{
nothing
}
\description{
export sampleMetadata
}
\examples{
data(exhaledPtrset )
saveFile<-file.path(getwd(),'sampleMetadata.tsv')
#exportSampleMetada(exhaledPtrset ,saveFile)
}
