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
#' @slot calibCoef calibration coefficients (a,b) such that: mz= ((tof-b)/a)^2 
#' @slot calibMassRef the reference masses used for the calibration
#' @slot calibError the shift error in ppm at the reference masses
#' @slot calibSpectr the spectrum of calibration reference masses
#' @slot ptrTransmisison matrix with transmission values
#' @slot prtReaction a list containing PTR reaction information: drift temperature, pressure and voltage 
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
    calibCoef="matrix",
    calibMassRef = "numeric",
    calibError="numeric",
    calibSpectr = "list",
    ptrTransmisison = "matrix",
    prtReaction = "list"
  )
)

#' Object contains a set of PTR-TOF-MS raw data
#' 
#' A ptrSet object is related to a directory that contains several PTR-TOF-MS raw data in h5 format. It is created 
#' with the \code{createPtrSet} function. 
#' 
#' @slot parameter the input parameters value of the function createPtrSet
#' @slot sampleMetadata data frame of sample metadata, with file names in row names 
#' @slot mzCalibRef the masses uses for calibration for each file
#' @slot signalCalibRef the spectrum of mass calibration for each file
#' @slot errorCalibPpm the calibration error for each file
#' @slot coefCalib the coefficients of mass axis calibration for each file
#' @slot primaryIon the quantity in count per acquisition time of the isotope of primary ion (H30_18, mz 21.022) for each file
#' @slot resolution estimation of the resolution for each file based on the calibration reference masses
#' @slot TIC the TIC for each file
#' @slot timeLimit the index of time limit for each file 
#' @slot peakListRaw the peak list before a potential aggregation  of expirations 
#' @slot peakListAligned the peakList with expirations aggregation  for each file 
#' @name ptrSet-class
#' @docType class
#' @rdname ptrSet-class
#' @exportClass ptrSet
setClass(
  Class = "ptrSet",
  representation = representation(
    parameter = "list",
    sampleMetadata = "data.frame",
    mzCalibRef= "list",
    signalCalibRef = "list",
    errorCalibPpm ="list",
    coefCalib="list",
    primaryIon="list",
    resolution ="list",
    TIC = "list",
    timeLimit="list",
    peakListRaw ="list",
    peakListAligned = "list"
  )
)

