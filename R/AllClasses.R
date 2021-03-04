#' Object contains PTR-TOF-MS raw data frome a file. 
#' 
#' A ptrRaw object contains PTR-TOF-MS raw data from one h5 file. It is created with the \code{readRaw} function.
#' 
#' @slot name the file name 
#' @slot rawM the intensity raw matrix
#' @slot mz array of mz axis 
#' @slot time numeric vector of acquisition time
#' @slot calibMzToTof function to convert mz to Tof
#' @slot calibToftoMz function to convert tof to mz
#' @slot calibCoef calibration coefficients (a,b) such that: mz= ((tof-b)/a)^2 for each calibration periods
#' @slot indexTimeCalib index time of each calibration periods
#' @slot calibMassRef the reference masses used for the calibration
#' @slot calibError the shift error in ppm at the reference masses
#' @slot calibSpectr the spectrum of calibration reference masses
#' @slot peakShape average normalized peak shape of the calibration peak 
#' @slot ptrTransmisison matrix with transmission values
#' @slot prtReaction a list containing PTR reaction information: drift temperature, pressure and voltage
#' @slot date acquisition date and hour 
#' @name ptrRaw-class
#' @rdname ptrRaw-class
#' @docType class
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
    date="character"
  )
)

#' Object contains a set of PTR-TOF-MS raw data
#' 
#' A ptrSet object is related to a directory that contains several PTR-TOF-MS raw data in h5 format. It is created 
#' with the \code{createPtrSet} function. 
#' 
#' @slot parameter the input parameters value of the function createPtrSet
#' @slot sampleMetadata data frame of sample metadata, with file names in row names 
#' @slot date acquisition date for each file
#' @slot mzCalibRef the masses uses for calibration for each file
#' @slot signalCalibRef the spectrum of mass calibration for each file
#' @slot errorCalibPpm the calibration error for each file
#' @slot coefCalib the coefficients of mass axis calibration of each calibration periods for each file
#' @slot indexTimeCalib index time of each calibration periods for each file
#' @slot primaryIon the quantity in count per acquisition time of the isotope of primary ion for each file
#' @slot resolution estimation of the resolution for each file based on the calibration reference masses
#' @slot prtReaction drift information (temperature, pressure and voltage)
#' @slot ptrTransmisison  transmision curve for each file
#' @slot TIC the TIC for each file
#' @slot breathTracer EIC for expirations/head spaces detection
#' @slot timeLimit the index of time limit for each file 
#' @slot knots numeric vector correspond to the knot that will be use for the two dimensional 
#' regression for each file
#' @slot fctFit the peak function use for peak deconvolution for each file
#' @slot peakShape average normalized peak shape of the calibration peak 
#' @slot peakList individual peak list in expression set (Biobase)
#' @name ptrSet-class
#' @docType class
#' @rdname ptrSet-class
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

