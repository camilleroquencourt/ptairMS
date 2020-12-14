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
           function(x,mzCalibRef = c(21.022, 29.013424,41.03858, 60.0525,203.943, 330.8495), tol=70) {
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
           function(object,fracMaxTIC=0.5,fracMaxTICBg=0.2, derivThresholdExp=1,derivThresholdBg=0.05,
                    mzBreathTracer=NULL, minPoints=2 ,degreeBaseline=1, baseline=TRUE,plotDel=FALSE ) {
             standardGeneric("timeLimits")
           })

## PeakList ----
#' @rdname PeakList
#' @export
## SetMethod
setGeneric("PeakList",
           function(raw,
                    mzNominal = unique(round(raw@mz)), 
                    ppm = 130, resMinMeanMax=c(3000,5000,8000),
                    minIntensity=5, fctFit=c("sech2","average")[1], maxIter=2, R2min,autocorNoiseMax = 0.3,
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
#' plotF <- plotFeatures(ptrSet,mz=59.049,type="ggplot)
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

#' Putatively annotate VOC mz by using the reference compilation from the literature
#'
#' Putatively annotate VOC mz by using the reference compilation from the literature, 
#' and detect isotope thanks to the \code{findIsotope} function. 
#'
#' @param x Expression set object (resp. data.frame) (resp. numeric vector) containing
#' the PTR-MS processed data (resp. containing a column with the ion mass values) (resp. containing the ion mass values)
#' @param ionMassColname Character: column name from the fData (resp. from the data.frame) containing
#' the ion mass values; [default: 'ion_mass']; this argument is not used when x is a numeric vector
#' @param ppm Numeric: tolerance
#' @param prefix Character: prefix for the new 'annotation' columns [default: 'vocDB_']
#' @param fields Characer vector: fields of the 'vocDB' database to be queried among:
#' 'ion_mass' [default], 'ion_formula' [default], 'formula', 'mass_monoiso',
#' 'name_iupac' [default], 'pubchem_cid', 'inchi', 'inchikey', 'ref_year',
#' 'ref_pmid', 'disease_name', 'disease_meshid'
#' @return Returns the data.frame with additional columns containing the vocDB informations
#' for the matched ion_mass values as well as the detected isotopes
#' @examples
#' library(ptairData)
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' bacteria.ptrset <- createPtrSet(directory, setName = "bacteria",
#' mzCalibRef = c(21.022,59.049))
#' bacteria.ptrset <- detectPeak(bacteria.ptrset,mzNominal =c(59,79))
#' bacteria.eset <-alignSamples(bacteria.ptrset,pValGreaterThres=0.05)
#' # Expression Set
#' bacteria.eset <- annotateVOC(bacteria.eset)
#' head(Biobase::fData(bacteria.eset)[, c("vocDB_ion_mass", "vocDB_ion_formula")])
#' # Data frame
#' bacteria_fdata.df <- Biobase::fData(bacteria.eset)
#' bacteria_fdata.df <- annotateVOC(bacteria_fdata.df)
#' head(bacteria_fdata.df[, c("vocDB_ion_mass", "vocDB_ion_formula")])
#' # Numeric
#' ionMass.vn <- as.numeric(Biobase::featureNames(bacteria.eset))
#' annotated_ions.df <- annotateVOC(ionMass.vn)
#' head(annotated_ions.df[, c("vocDB_ion_mass", "vocDB_ion_formula")])
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
                               "disease_meshid")[c(1,2,5)]){ standardGeneric("annotateVOC") })

##aligneSamples----
#' @rdname alignSamples
#' @export
setGeneric("alignSamples",function(X, ppmGroup = 70, fracGroup = 0.8, group=NULL,
                                   bgThreshold=2,pValGreaterThres= 2e-26,pValLessThres=0,fracExp=0.3,
                                   dmzGroup = 0.001,...){standardGeneric("alignSamples")})

##writte----
#' @rdname writeEset
#' @export
setGeneric("writeEset",function(x,
                                dirName,
                                overwrite = FALSE,
                                verbose = TRUE){standardGeneric("writeEset")})


