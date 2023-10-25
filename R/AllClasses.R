#' PTR-TOF-MS raw data from a rhdf5 file
#' 
#' A ptrRaw object contains PTR-TOF-MS raw data from one rhdf5 file. It is 
#' created with the \code{\link[ptairMS]{readRaw}} function.
#' 
#' @slot name the file name 
#' @slot rawM the raw intensities matrix
#' @slot mz array of the m/z axis 
#' @slot time numeric vector of acquisition time (in seconds)
#' @slot calibMzToTof function to convert m/z to Tof
#' @slot calibToftoMz function to convert tof to m/z
#' @slot calibCoef calibration coefficients (a,b) such that: mz= ((tof-b)/a)^2 
#' for each calibration period
#' @slot indexTimeCalib index time of each calibration period
#' @slot calibMassRef the reference masses used for the calibration
#' @slot calibError the shift error in ppm at the reference masses
#' @slot calibSpectr the spectrum of calibration reference masses
#' @slot peakShape average normalized peak shape of the calibration peak 
#' @slot ptrTransmisison matrix with transmission values
#' @slot prtReaction a list containing PTR reaction information: drift 
#' temperature, pressure and voltage
#' @slot date acquisition date and hour 
#' @slot peakList individual peak list in \code{\link[Biobase]{eSet}} 
#' @slot fctFit the peak function used for peak deconvolution for each file
#' @slot resolution estimation of the resolution for each file based on the 
#' calibration reference masses
#' @slot primaryIon the quantity in count per acquisition time of the isotope 
#' of primary ion for each file
#' @name ptrRaw-class
#' @rdname ptrRaw
#' @docType class
#' @references https://www.hdfgroup.org
#' @exportClass ptrRaw
setClass(
  Class = "ptrRaw",
  representation = representation(
    name="character",
    rawM= "matrix",
    mz = "numeric",
    time = "numeric",
    calibMzToTof = "function",
    calibToftoMz = "function",
    calibCoef="list",
    indexTimeCalib="list",
    calibMassRef = "numeric",
    calibError="numeric",
    calibSpectr = "list",
    peakShape = "list",
    ptrTransmisison = "matrix",
    prtReaction = "list",
    date="character",
    fctFit="character",
    peakList ="ExpressionSet",
    resolution ="numeric",
    primaryIon="numeric"
  )
)

#' A set of PTR-TOF-MS raw data informations
#' 
#' A ptrSet object is related to a directory that contains several PTR-TOF-MS 
#' raw data in rhdf5 format. It is created 
#' with the \code{\link[ptairMS]{createPtrSet}} function. This object could be
#' updated when new files are added with the \code{\link[ptairMS]{updatePtrSet}}
#' function.
#' 
#' @slot parameter the input parameters value of the function 
#' \code{\link[ptairMS]{createPtrSet}} and \code{\link[ptairMS]{detectPeak}}
#' @slot sampleMetadata dataframe of sample metadata, with file names 
#' in row names, suborders names and acquisition date in columns
#' @slot date acquisition date for each file
#' @slot mzCalibRef the masses used for calibration for each file
#' @slot signalCalibRef the spectrum of mass calibration for each file
#' @slot errorCalibPpm the calibration error for each file
#' @slot coefCalib the coefficients of mass axis calibration of each calibration 
#' periods for each file
#' @slot indexTimeCalib index time of each calibration period for each file
#' @slot primaryIon the quantity in count per acquisition time of the isotope 
#' of primary ion for each file
#' @slot resolution estimation of the resolution for each file based on the 
#' calibration reference masses
#' @slot prtReaction drift information (temperature, pressure and voltage)
#' @slot ptrTransmisison  transmission curve for each file
#' @slot TIC the total ion current (TIC) for each file
#' @slot breathTracer the tracer for expiration/head spaces detection
#' @slot timeLimit the index of time limit for each file 
#' @slot knots numeric vector correspond to the knot that will be used for the 
#' two dimensional regression for each file
#' @slot fctFit the peak function used for peak deconvolution for each file
#' @slot peakShape average normalized peak shape of the calibration peak for 
#' each file
#' @slot peakList individual peak list in \code{\link[Biobase]{eSet}}
#' @name ptrSet-class
#' @docType class
#' @rdname ptrSet
#' @exportClass ptrSet
setClass(
  Class = "ptrSet",
  representation = representation(
    parameter = "list",
    sampleMetadata = "data.frame",
    date="list",
    mzCalibRef= "list",
    signalCalibRef = "list",
    errorCalibPpm ="list",
    coefCalib="list",
    indexTimeCalib="list",
    primaryIon="list",
    ptrTransmisison = "list",
    prtReaction = "list",
    resolution ="list",
    TIC = "list",
    breathTracer="list",
    timeLimit="list",
    knots="list",
    fctFit="list",
    peakShape="list",
    peakList ="list"
  )
)

