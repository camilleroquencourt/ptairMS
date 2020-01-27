##detectPeak----
#' @rdname detectPeak
#' @export
setGeneric("detectPeak",
           function(x, ...){
             standardGeneric("detectPeak")
           })

## Calibration ----
#' @rdname calibration
#' @export
setGeneric(name = "calibration", 
           function(x,mzCalibRef = c(21.022, 29.013424,41.03858,59.049141,75.04406, 
                                     203.943, 330.8495), tol=70) {
             standardGeneric("calibration")
             })

## plotRaw -----
#' @rdname plotRaw
#' @export
setGeneric("plotRaw",
           function(object,
                    mzRange ,
                    timeRange = c(NA, NA),
                    type = c("classical", "plotly")[1],
                    ppm = 2000,
                    palette = c("heat",
                                "revHeat",
                                "grey",
                                "revGrey",
                                "ramp")[1],
                    showVocDB = TRUE,
                    figure.pdf = "interactive", ...) 
           {
             standardGeneric("plotRaw")
           })

## timeLimit ----
#' @rdname timeLimits
#' @export
setGeneric("timeLimits",
           function(object,fracMaxTIC=0.5, traceMasses=NULL, minPoints=2 , plotDel=FALSE ) {
             standardGeneric("timeLimits")
           })

## PeakList ----
#' @rdname PeakList
#' @export
## SetMethod
setGeneric("PeakList",
           function(raw,
                    mzNominal = unique(round(raw@mz)), 
                    ppm = 130, 
                    minIntensity=5, fctFit=c("Sech2","average")[1], maxIter=2, autocorNoiseMax = 0.3,
                    plotFinal=FALSE, plotAll=FALSE, thNoiseRate=1.1, thIntensityRate = 0.01,
                    countFacFWHM=10, daSeparation=0.005, d=3, windowSize=0.4){
             standardGeneric("PeakList")})


## plot TIC----
#' @rdname plotTIC
#' @export
setGeneric("plotTIC",
           function(object, type = c("plotly","ggplot")[1], baselineRm=FALSE, 
                    showLimits=FALSE ,...){
             standardGeneric("plotTIC")} )

## plotFeatures----
#' plot signal for all files around a mzRange
#' 
#' Plot the raw data spectrum for several files in a ptrSet object around the \code{mz} masses.
#' @param set a ptrSet object
#' @param mz the mz values to plot
#' @param typePlot set "plotly" to get an interactive plot, and "ggplot" for classical plot.    
#' @param addFeatureLine boolean. If TRUE a vertical line at the mz masses is plotted
#' @param ppm windows wize in ppm
#' @param pdfFile a file path to save a pdf with a individual plot per file
#' @param fileNames vector of character. The file names you want to plot. If \code{NULL}, it plot all files
#' @param colorBy character. A column name of sample metadata by which the line are colored. 
#' @return a plotly or ggplot2 object 
#' @examples 
#' library(ptairData)
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' ptrSet <- createPtrSet(directory,setName="ptrSet",mzCalibRef=c(21.022,59.049))
#' plotF <- plotFeatures(ptrSet,mz=59.049)
#' print(plotF)
#' @rdname plotFeatures
#' @export
setGeneric("plotFeatures",
           function(set, mz, typePlot="plotly", addFeatureLine =TRUE,ppm=2000,
                    pdfFile=NULL, fileNames=NULL, colorBy="rownames"){standardGeneric("plotFeatures")})

## ploCalib----
#'@rdname plotCalib
#'@export 
setGeneric("plotCalib",
           function(object,ppm=2000, ...){standardGeneric("plotCalib")})

## getFileNames----
#'@rdname getFileNames
#'@export
setGeneric("getFileNames", function(object,fullNames=FALSE){
  standardGeneric("getFileNames")
})

## annotateVOC ----

#' @rdname annotation
#' @export 
setGeneric("annotateVOC",
           function(x, ...){ standardGeneric("annotateVOC") })

##aligneSamples----
#' @rdname alignSamples
#' @export
setGeneric("alignSamples",function(X, ppmGroup = 70, fracGroup = 0.8, group=NULL,bgThreshold=2,
                                   dmzGroup = 0.001,...){standardGeneric("alignSamples")})

##writte----
#' @rdname writeEset
#' @export
setGeneric("writeEset",function(x,
                                dirName,
                                overwrite = FALSE,
                                verbose = TRUE){standardGeneric("writeEset")})


