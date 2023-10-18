##detectPeak----
#' @rdname detectPeak
#' @export
setGeneric("detectPeak",
           function(x, 
                    ppm=130, 
                    minIntensity=10, 
                    minIntensityRate = 0.01,
                    mzNominal=NULL, 
                    resolutionRange=NULL,
                    fctFit=NULL,
                    smoothPenalty=NULL,
                    parallelize=FALSE,
                    nbCores=2,
                    saving=TRUE,
                    saveDir=getParameters(x)$saveDir,...){
             standardGeneric("detectPeak")
           })

## Calibration ----
#' @rdname calibration
#' @export
setGeneric(name = "calibration", 
           function(x,mzCalibRef = c(21.022, 29.013424,41.03858, 
                                     60.0525,203.943, 330.8495), 
                    calibrationPeriod=60,tol=70,...) {
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
           function(object,fracMaxTIC=0.5,fracMaxTICBg=0.2, 
                    derivThresholdExp=1,derivThresholdBg=0.05,
                    mzBreathTracer=NULL, minPoints=2 ,degreeBaseline=1, 
                    baseline=TRUE,
                    redefineKnots=TRUE,plotDel=FALSE ) {
             standardGeneric("timeLimits")
           })


## defineKnots -----
#' @param ... not used
#' @rdname defineKnots
#' @export
setGeneric("defineKnots",
           function(object, knotsPeriod=3, 
                    method=c("expiration","uniform","manual")[1],
                    knotsList=NULL,...) {
             standardGeneric("defineKnots")
           })

## PeakList ----
#' @rdname PeakList
#' @export
## SetMethod
setGeneric("PeakList",
           function(raw,
                    mzNominal = unique(round(getRawInfo(raw)$mz)), 
                    ppm = 130, resolutionRange=c(3000,5000,8000),
                    minIntensity=5, fctFit=c("sech2","averagePeak")[1], 
                    peakShape=NULL,maxIter=1, R2min=0.995,autocorNoiseMax = 0.3,
                    plotFinal=FALSE, plotAll=FALSE, thNoiseRate=1.1, 
                    minIntensityRate = 0.01,
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
#' Plot raw average spectrum around a mzRange
#' 
#' Plot the raw data spectrum for several files in a ptrSet object around 
#' the \code{mz} masses. The expiration average spectrum is in full lines, 
#' and background in dashed lines.
#' @param set a \code{\link[ptairMS]{ptrSet-class}} object
#' @param mz the mz values to plot
#' @param typePlot set "plotly" to get an interactive plot, or "ggplot"
#' @param addFeatureLine boolean. If TRUE a vertical line at the mz masses 
#' is plotted
#' @param ppm windows size of the plot round \code{mz} in ppm
#' @param pdfFile a file path to save a pdf with a individual plot per file
#' @param fileNames vector of character. The file names you want to plot. 
#' If \code{NULL}, it plot all files
#' @param colorBy character. A column name of sample metadata by which 
#' the line are colored. 
#' @return a plotly or ggplot2 object 
#' @examples 
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, 
#' setName="exhaledPtrset", mzCalibRef = c(21.022, 60.0525),
#' fracMaxTIC = 0.7, saveDir = NULL )
#' plotF<-plotFeatures(exhaledPtrset ,mz=59.049,type="ggplot",colorBy="subfolder")
#' print(plotF)
#' @rdname plotFeatures
#' @export
setGeneric("plotFeatures",
           function(set, mz, typePlot="plotly", 
                    addFeatureLine =TRUE,ppm=2000,
                    pdfFile=NULL, fileNames=NULL, 
                    colorBy="rownames"){standardGeneric("plotFeatures")})

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

setGeneric("getParameters", function(object){
    standardGeneric("getParameters")
})

setGeneric("setParameters", function(object,newparam){
    standardGeneric("setParameters")
})

setGeneric("getDate", function(object){
    standardGeneric("getDate")
})



setGeneric("getCalibrationInfo", function(object){
    standardGeneric("getCalibrationInfo")
})

setGeneric("getPeaksInfo", function(object){
    standardGeneric("getPeaksInfo")
})


setGeneric("setPeakShape", function(object,peakShape){
    standardGeneric("setPeakShape")
})


setGeneric("setCalibration", function(object,newCalibrationInfo,...){
    standardGeneric("setCalibration")
})


setGeneric("getTimeInfo", function(object){
    standardGeneric("getTimeInfo")
})


setGeneric("getPTRInfo", function(object){
    standardGeneric("getPTRInfo")
})

setGeneric("getRawInfo", function(object){
    standardGeneric("getRawInfo")
})

setGeneric("setTimeLimits", function(object,newTimeLimites,index=NULL){
    standardGeneric("setTimeLimits")
})

setGeneric("setKnots", function(object,newKnots,index=NULL){
    standardGeneric("setKnots")
})
setGeneric("deleteFilePtrSet", function(object,deletedFiles){
    standardGeneric("deleteFilePtrSet")
})

setGeneric("mergedPtrSet", function(object,ptrSetNewfile,orderedFile){
    standardGeneric("mergedPtrSet")
})



## annotateVOC ----

#' Putative annotation of VOC mz by using the reference compilation from the 
#' literature
#'
#' Putatively annotate VOC mz by using the reference compilation from the 
#' literature, and porpose an isotope detetcion.
#'
#' @param x Expression set object (resp. data.frame) (resp. numeric vector) 
#' containing
#' the PTR-MS processed data (resp. containing a column with the ion mass values) 
#' (resp. containing the ion mass values)
#' @param ionMassColname Character: column name from the fData 
#' (resp. from the data.frame) containing
#' the ion mass values; [default: 'ion_mass']; this argument is not used when x 
#' is a numeric vector
#' @param ppm Numeric: tolerance
#' @param prefix Character: prefix for the new 'annotation' 
#' columns [default: 'vocDB_']
#' @param fields Characer vector: fields of the 'vocDB' database to be 
#' queried among:
#' 'ion_mass' [default], 'ion_formula' [default], 'formula', 'mass_monoiso',
#' 'name_iupac' [default], 'pubchem_cid', 'inchi', 'inchikey', 'ref_year',
#' 'ref_pmid', 'disease_name', 'disease_meshid'
#' @return Returns the data.frame with additional columns containing the vocDB 
#' informations
#' for the matched ion_mass values as well as the detected isotopes
#' @examples
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, 
#' setName="exhaledPtrset", mzCalibRef = c(21.022, 60.0525),
#' fracMaxTIC = 0.7, saveDir = NULL )
#' exhaledPtrset<-detectPeak(exhaledPtrset,mz=c(21,59,77))
#' exhaled.eset <-alignSamples(exhaledPtrset ,pValGreaterThres=0.05)
#' # Expression Set
#' exhaled.eset <- annotateVOC(exhaled.eset)
#' head(Biobase::fData(exhaled.eset))
#' # Data frame
#' exhaled_fdata.df <- Biobase::fData(exhaled.eset)
#' exhaled_fdata.df <- annotateVOC(exhaled_fdata.df)
#' head(exhaled_fdata.df)
#' # Numeric
#' ionMass.vn <- as.numeric(Biobase::featureNames(exhaled.eset))
#' annotated_ions.df <- annotateVOC(ionMass.vn)
#' head(annotated_ions.df)
#' @rdname annotation
#' @export 
setGeneric("annotateVOC",
           function(x,
                    ionMassColname = "ion_mass",
                    ppm = 20,
                    prefix = "vocDB_",
                    fields = c("ion_mass",
                               "ion_formula",
                               "formula",
                               "mass_monoiso",
                               "name_iupac",
                               "pubchem_cid",
                               "inchi",
                               "inchikey",
                               "ref_year",
                               "ref_pmid",
                               "disease_name",
                               "disease_meshid")[c(1,2,5)]){ 
             standardGeneric("annotateVOC") })

##alignSamples----
#' @rdname alignSamples
#' @export
setGeneric("alignSamples",function(X, 
                                   ppmGroup = 70, 
                                   fracGroup = 0.8, 
                                   group=NULL,
                                   fracExp=0.3,
                                   pValGreaterThres= 0.001,pValLessThres=0,
                                   quantiUnit=c("ppb","ncps","cps")[1],
                                   bgCorrected=TRUE,
                                   dmzGroup = 0.001){
  standardGeneric("alignSamples")})

##writte----
#' @rdname writeEset
#' @export
setGeneric("writeEset",function(x,
                                dirName,
                                overwrite = FALSE,
                                verbose = TRUE){standardGeneric("writeEset")})


