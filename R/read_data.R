#' Read a h5 file of PTR-TOF-MS data
#'
#' \code{readRaw} reads a h5 file with rhdf5 library, and calibrates the 
#' mass axis with \code{mzCalibRef} masses each \code{calibrationPeriod} 
#' seconds. It returns a \code{\link[ptairMS]{ptrRaw-class}} S4 object, 
#' that contains raw data.
#' 
#' @param filePath h5 absolute file path full name.
#' @param calib boolean. If true, an external calibration is performed on 
#' the \code{calibrationPeriod} sum spectrum with mzCalibRef reference masses. 
#' @param mzCalibRef calibration parameter. Vector of exact theoretical masses 
#' values of an intensive peak without overlapping.
#' @param calibrationPeriod in second, coefficient calibration are estimated 
#' for each sum spectrum of 
#' \code{calibrationPeriod} seconds
#' @param tolCalibPpm calibration parameter. The maximum error tolerated in ppm. 
#' A warning appears for error greater than \code{tolCalibPpm}.
#' @param maxTimePoint number maximal of time point to read
#' @return a ptrRaw object, including slot \itemize{
#' \item rawM the data raw matrix, in count of ions  
#' \item mz the mz axis
#' \item time time acquisition in second 
#' } 
#' @examples
#' library(ptairData)
#' filePathRaw <- system.file('extdata/exhaledAir/ind1', 'ind1-1.h5', 
#' package = 'ptairData')
#' raw <- readRaw(filePath=filePathRaw, mzCalibRef=c(21.022, 60.0525), calib=FALSE)
#' @import bit64
#' @export
readRaw <- function(filePath, calib = TRUE, mzCalibRef = c(21.022, 29.013424, 41.03858, 
    60.0525, 203.943, 330.8495), calibrationPeriod = 60, tolCalibPpm = 70, maxTimePoint = 900) {
    
    if (is.null(filePath) | filePath == "") 
        stop("filePath is empty")
    if (is.na(filePath)) 
        stop("filePath is empty")
    if (!is.character(filePath)) 
        stop("filePath must be a character")
    if (length(filePath) > 1) 
        stop("Only one filePath is required")
    if (!grepl("*h5$", filePath)) 
        stop("The file is not in h5 format")
    
    # find the longer of acquisitionTime
    
    name<-rhdf5::h5ls(filePath)
    
    if(name$group[2] == "/AcquisitionLog"){
        file <- rhdf5::H5Fopen(filePath)
        time <- rhdf5::H5Dopen(file, "/TimingData/BufTimes")
        acTime <- rhdf5::H5Sget_simple_extent_dims(rhdf5::H5Dget_space(time))$size
        nbrBuf <- acTime[1]
        nbrWrite <- acTime[2]
        
        if (nbrWrite == 0) {
            stop("The file is empty (0 seconde of acquisition)")
        }
        
        # acquiqision time limit to 900 spectra, indeed the rawAn is more than 3.4 GB
        NbrWriteMax <- ceiling(maxTimePoint/nbrBuf)
        rhdf5::h5closeAll()
        
        # read information needed
        date_heure <- rhdf5::h5read(filePath, "/AcquisitionLog", bit64conversion = "bit64")$Log$timestring[1]
        
        timVn <- rhdf5::h5read(filePath, "/TimingData/BufTimes", bit64conversion = "bit64", 
                               index = list(NULL, seq_len(min(nbrWrite, NbrWriteMax))))
        rawAn <- rhdf5::h5read(filePath, "/FullSpectra/TofData", index = list(NULL, NULL, 
                                                                              NULL, seq_len(min(nbrWrite, NbrWriteMax))), bit64conversion = "bit64")
        mzVn <- rhdf5::h5read(filePath, "FullSpectra/MassAxis", bit64conversion = "bit64")
        reaction <- try(rhdf5::h5read(filePath, "AddTraces/PTR-Reaction"))
        transmission <- try(rhdf5::h5read(filePath, "PTR-Transmission"))
        calibCoef <- try(rhdf5::h5read(filePath, "FullSpectra/MassCalibration", index = list(NULL, 
                                                                                             1)))
        attributCalib <- try(rhdf5::h5readAttributes(filePath, "/FullSpectra"))
        if (!is.null(attr(calibCoef, "condition")) & is.null(attr(attributCalib, "condition"))) {
            calibCoef <- matrix(c(attributCalib$`MassCalibration a`, attributCalib$`MassCalibration b`), 
                                ncol = 1)
        }
        
        #singleIon<-attributCalib$`Single Ion Signal` #2.9
        #sampleInterval<-attributCalib$SampleInterval #2e-10
        
        if (!is.null(attr(transmission, "condition"))) 
            
            transmission <- matrix(0) else transmission <- transmission$Data
        
        if (!is.null(attr(reaction, "condition"))){
            reaction <- list()
        } else {
            index<-grep(pattern = "Data",names(reaction))
            reaction[[index]] <-  matrix(reaction[[index]], nrow = dim(reaction[[index]])[1], ncol = prod(utils::tail(dim(reaction[[index]]), 
                                                                                           2)))
            
        }
        
        # re format
        timVn <- c(timVn)
        mzVn <- c(mzVn)
        rawMn <- matrix(rawAn, nrow = dim(rawAn)[1], ncol = prod(utils::tail(dim(rawAn), 
                                                                             2)), dimnames = list(mzVn, timVn))
        # remove index where timVn =0 except the first time
        index_zero <- which(timVn == 0)[-1]
        
        rm(rawAn)
        if (length(index_zero) != 0) {
            timVn <- timVn[-index_zero,drop=FALSE]
            rawMn <- rawMn[, -index_zero,drop=FALSE]
        }
        #timVn<-sort(timVn)
        #colnames(rawMn)<-timVn
        
        # count ion conversion 
        #factor<- sampleInterval*1e9 /singleIon
        #*(SampleInterval,ns) / (single ion signal mV.ns) for convert to number of ion
        #rawMn <- rawMn*(factor)
        
        # calibration infomration
        if (is.null(attr(calibCoef, "condition"))) {
            calibCoef<-calibCoef[c(1,2),1,drop=FALSE]
            rownames(calibCoef) <- c("a", "b")
            calib_formula <- function(tof, calibCoef) ((tof - calibCoef["b", ])/calibCoef["a", 
            ])^2
            calib_invFormula <- function(m, calibCoef) sqrt(m) * calibCoef["a", ] + calibCoef["b", 
            ]
            calibMassRef = c(attributCalib$`MassCalibration m1`, attributCalib$`MassCalibration m2`)
            if (is.null(calibMassRef)) 
                calibMassRef <- 0
        } else {
            calibCoef <- matrix(0)
            calib_formula <- function(tof, calibCoef) NULL
            calib_invFormula <- function(m, calibCoef) NULL
            calibMassRef <- 0
        }
        calibCoef<-list(calibCoef)
        
        
    } else if (name$group[2] == "/AddTraces"){
        
        file <- rhdf5::H5Fopen(filePath)
        time <- rhdf5::H5Dopen(file, "/SPECdata/Times")
        acTime <- rhdf5::H5Sget_simple_extent_dims(rhdf5::H5Dget_space(time))$size
        nbrWrite <- acTime[2]
        
        if (nbrWrite == 0) {
            stop("The file is empty (0 seconde of acquisition)")
        }
        
        # acquiqision time limit to 900 spectra, indeed the rawAn is more than 3.4 GB
        NbrWriteMax <- maxTimePoint
        rhdf5::h5closeAll()
        
       
        date_heure <- rhdf5::h5readAttributes(filePath,"/")$FileCreatedTimeSTR_LOCAL[1]
        timVn <- rhdf5::h5read(filePath, "/SPECdata/Times", bit64conversion = "bit64", 
                               index = list(NULL, seq_len(min(nbrWrite, NbrWriteMax))))[4,]
        rawMn <- rhdf5::h5read(filePath, "/SPECdata/Intensities", bit64conversion = "bit64", 
                               index = list(NULL, seq_len(min(nbrWrite, NbrWriteMax))))
      
        reaction <- try(rhdf5::h5read(filePath, "/AddTraces/PTR-Instrument"))
        
        index <- try(sapply(c("Udrift_Act","p-Drift_Act","T-Drift_Act","E/N_Act","Udrift[Act]","p-Drift","Temp: T-Drift[Act]","E/N[Act]"), function(x) grep(x = reaction$Info,pattern = x,fixed = TRUE))) 
        index<-Reduce(c,index)
        
        if(length(index)==0){
            reaction <- try(rhdf5::h5read(filePath, "/AddTraces/PTR-Reaction"))
            
            index <- try(sapply(c("Udrift","p-Drift","T-Drift","E/N"), function(x) grep(x = reaction$Info,pattern = x,fixed = TRUE))) 
            index<-Reduce(c,index)
            }
        
        if(length(index)==0){
            warning("missing reaction info")    
            reaction<- list()
        }else {
            reaction$Data<-reaction$Data[index, ]
            
            reaction$Info<-  c("Udrift[V]",         
                               "p-Drift[mbar]",    
                               "T-Drift[C]" ,      
                               "E/N[Td]")
            
            names(reaction)<-c("TwData", "TwInfo")
        }
        
        
        transmission <- try(t(rhdf5::h5read(filePath, "PTR-Transmission")$Masses_Factors[,,1]))
      
        CalibInfo <-rhdf5::h5read(filePath, "/CALdata", index = list(NULL,1))
        if(dim(CalibInfo$Spectrum)[1]!=0) calibCoefFirt<-as.matrix(CalibInfo$Spectrum[,1]) else{
            mzRef <- CalibInfo$Mapping[1,]
            tofRef <- CalibInfo$Mapping[2,]
            a <- (tofRef[2] - tofRef[1]) / (sqrt(mzRef[2])-sqrt(mzRef[1]))
            b<- tofRef[1] - sqrt(mzRef[1])*a
            calibCoefFirt<-as.matrix(c(a,b))
            
        }
        
        rownames(calibCoefFirt)<-c("a","b")
        mzVn <- ptairMS:::tofToMz(seq(0,(dim(rawMn)[1]-1)),calibCoef = calibCoefFirt)
        colnames(rawMn)<-timVn
        rownames(rawMn)<-mzVn
        
        if(dim(CalibInfo$Spectrum)[1]!=0)  calibCoef<-lapply(apply(CalibInfo$Spectrum,2,function(x){
          y<-as.matrix(x,ncol=1,nrow=2)
          rownames(y)<-c("a","b")
           list(y)
       } ),function(x) x[[1]]) else calibCoef<-list(calibCoefFirt)
  
        # singleIon<-attributCalib$`Single Ion Signal`
        # sampleInterval<-attributCalib$SampleInterval
        
        if (!is.null(attr(transmission, "condition"))) 
            transmission <- matrix(0) 
       
        # calibration infomration
     
            calib_formula <- function(tof, calibCoef) ((tof - calibCoef["b", ])/calibCoef["a", 
            ])^2
            calib_invFormula <- function(m, calibCoef) sqrt(m) * calibCoef["a", ] + calibCoef["b", 
            ]
            calibMassRef = CalibInfo$Mapping[1,]

    }
    
    # write ptrRaw objet
    raw <- methods::new(Class = "ptrRaw", name = filePath, rawM = rawMn, 
        mz = mzVn, time = timVn, calibCoef = calibCoef, calibMzToTof = calib_invFormula, 
        calibToftoMz = calib_formula, calibError = 0, calibMassRef = calibMassRef, 
        calibSpectr = list(NULL), peakShape = list(NULL), ptrTransmisison = transmission, 
        prtReaction = reaction, date = date_heure,fctFit="", peakList= Biobase::ExpressionSet(), 
        resolution= c(0,0,0),primaryIon=0)
    rm(rawMn)
    if (calib) {
        raw <- calibration(x = raw, mzCalibRef, calibrationPeriod = calibrationPeriod, 
            tol = tolCalibPpm)
    }
    
    return(raw)
}

#' Creates a ptrSet object form a directory
#'
#' This function creates a \code{\link[ptairMS]{ptrSet-class}} S4 object. It 
#' opens each file and: 
#' \itemize{
#' \item performs an external calibration by using the \code{mzCalibRef} 
#' reference masses on the sum spectra every \code{calibrationPeriod} second
#' \item quantifies the primary ion (H3O+ isotope by default)  on the average 
#' total ion spectrum. 
#' \item calculates expiration on the \code{mzBreathTracer} trace. The part of 
#' the trace where  the intensity is higher than \code{fracMaxTIC * max(trace)} 
#' is considered as expiration. 
#' If \code{fracMaxTIC} is different to zero, this step is skipped
#' \item defines the set of knots for the peak analysis 
#' (see \code{\link[ptairMS]{detectPeak}})
#' \item provides a default sampleMetadata based on the tree structure of the 
#' directory and the acquisition date (a \code{data.frame} with file names as 
#' row names) 
#' \item If \code{saveDir} is not \code{NULL}, the returned object will be saved 
#' as a  \code{.Rdata} in \code{saveDir} with the \code{setName} as name}
#'
#' @param dir character. a directory path which contains several h5 files, 
#' possibly organized in subfolder 
#' @param setName character. name of the ptrSet object. If `saveDir` is 
#' not null, the object will be saved with this name.
#' @param mzCalibRef vector of the reference mass values; those masses should be 
#' accurate, and the corresponding peaks should be of high intensity and 
#' 'unique' in a nominal mass interval (without overlapping peaks) to performs 
#' calibration. See \code{ ?calibration}.
#' @param calibrationPeriod in second, coefficient calibration are estimated 
#' for each sum spectrum of 
#' \code{calibrationPeriod} seconds
#' @param fracMaxTIC Fraction (between 0 and 1) of the maximum of the 
#' Total Ion Current (TIC) amplitude after baseline removal. 
#' Only the part of the spectrum where the TIC intensity is higher than 
#' `fracMaxTIC * max(TIC) ` will be analyzed. If you want to analyze the entire 
#' spectrum,  set this parameter to 0. 
#' @param mzBreathTracer integer: nominal mass of the 
#' Extracted Ion Current (EIC) used to compute the expiration time limits. 
#' If \code{NULL}, the limits will be computed on the Total Ion Current (TIC).
#' @param knotsPeriod period in second (time scale) between two knots for 
#' the two dimensional modeling
#' @param mzPrimaryIon Exact mass of the primary ion isotope
#' @param saveDir Directory where the ptrSet object will be saved as 
#' .RData. If \code{NULL}, nothing will be saved.
#' @return a ptrSet object with slots :
#' \itemize{
#' \item Parameter: list containing the parameters used for 
#' \code{createPrtSet}, \code{detectPeak} and  \code{alignTimePeriods} functions.  
#' \item sampleMetadata: data frame containing information about the data, 
#' with file names in row names
#' \item mzCalibRef: list containing for each file the masses used for the 
#' calibration  (see \code{?ptairMS::calibration} for more details)
#' \item signalCalibRef: mz and intensity +- 0.4Da around the calibration 
#' masses
#' \item errorCalibPpm: list containing for each file the accuracy error 
#' in ppm at each calibration masses
#' \item coefCalib: list containing the calibration coefficients 'a' and 
#' 'b' which enable to convert tof to mz for each file 
#' (see \code{\link[ptairMS]{calibration}} function for more details.  
#' \item resolution: estimated resolution \eqn{m / \Delta m} for each 
#' calibration masses within each file
#' \item TIC: The Total Ion Current for each file
#' \item timeLImit: list containing, for each file, a list of two element: 
#' the matrix of time limit for each file  (if \code{fracMaxTIC} is different 
#' to zero), and the background index. See \code{\link[ptairMS]{timeLimits}}
#' for more details
#' \item peakList: list containing for each file an expression set 
#' \code{\link[Biobase]{eSet}}, with m/z peak center, quantification for background 
#' and exhaled air in cps, ppb and ncps, and quantity for each time points. See 
#' \code{\link[ptairMS]{getPeakList}} for more details. 
#' }
#' @export 
#' @examples
#' library(ptairData)
#' directory <- system.file('extdata/mycobacteria',  package = 'ptairData')
#' ptrSet<-createPtrSet(dir=directory,setName='ptrSet'
#' ,mzCalibRef=c(21.022,59.049),
#' fracMaxTIC=0.9,saveDir= NULL)
createPtrSet <- function(dir, setName, mzCalibRef = c(21.022, 29.013424, 41.03858, 
    60.0525, 203.943, 330.8495), calibrationPeriod = 60, fracMaxTIC = 0.8, mzBreathTracer = NULL, 
    knotsPeriod = 3, mzPrimaryIon = 21.022, saveDir = NULL) {
    
    # test on parameter dir exist
    if (dir == "") 
        stop("dir is empty")
    if (!dir.exists(dir)) 
        stop(dir, " does not exist")
    if (!is.null(saveDir)) {
        if (!dir.exists(saveDir)) 
            stop(saveDir, " does not exist")
    }
    # setName a character
    if (!is.character(setName)) 
        stop("setName must be a character")
    # fracMaxTIC must be between 0 and 1
    if (fracMaxTIC < 0 | fracMaxTIC > 1) 
        stop("fracMaxTIC must be between zero 
                                        and 1, see ?createPtrSet")
    # mzCalibRef numeric
    if (!is.numeric(mzCalibRef)) 
        stop("mzCalibRef must be a numeric vector")
    
    # get h5 files name
    filesFullName <- list.files(dir, recursive = TRUE, pattern = "\\.h5$", full.names = TRUE)
    fileDir <- dirname(list.files(dir, recursive = TRUE, pattern = "\\.h5$"))
    fileName <- basename(filesFullName)
    
    # save parameters
    parameter <- list(dir = dir, 
                      listFile = filesFullName, 
                      name = setName, 
                      mzCalibRef = mzCalibRef, 
        timeLimit = list(fracMaxTIC=fracMaxTIC,
                         fracMaxTICBg = 0.5, 
                         derivThresholdExp = 0.5, 
                         derivThresholdBg = 0.01, 
                         minPoints = 3, 
                         degreeBaseline = 1), 
        mzBreathTracer = mzBreathTracer, 
        knotsPeriod = knotsPeriod, 
        mzPrimaryIon=mzPrimaryIon,
        calibrationPeriod=calibrationPeriod,
        saveDir = saveDir)
    
    # create sampleMetadata test if there is subfolder if there is no subfolder
    if (all(fileDir == ".")) {
        sampleMetadata <- data.frame(row.names = fileName)
    } else {
        subfolder <- strsplit(fileDir, "/")
        # max depth of subfolder
        nSubfolder <- max(vapply(subfolder, length, 0))
        # set every element at the same length
        subfolder <- lapply(subfolder, function(x) {
            length(x) <- nSubfolder
            return(x)
        })
        group <- do.call(rbind, subfolder)
        sampleMetadata <- data.frame(subfolder = group, row.names = fileName)
    }
    
    # checkSet
    check <- checkSet(files = filesFullName,mzCalibRef =  mzCalibRef, 
                      fracMaxTIC, mzBreathTracer, calibrationPeriod, 
        knotsPeriod, mzPrimaryIon)
    
  

    if (length(check$failed) > 0) {
        parameter$listFile <- parameter$listFile[-which(basename(parameter$listFile) %in% 
            check$failed)]
        sampleMetadata <- sampleMetadata[-which(row.names(sampleMetadata) %in% check$failed), 
            , drop = FALSE]
    }
    
    sampleMetadata$date <- Reduce(c, check$date)
    # create ptrSet object
    ptrSet <- methods::new(Class = "ptrSet", parameter = parameter, sampleMetadata = sampleMetadata, 
        mzCalibRef = check$mzCalibRefList, timeLimit = check$timeLimit, knots = check$knots, 
        signalCalibRef = check$signalCalibRef, errorCalibPpm = check$errorCalibPpm, 
        coefCalib = check$coefCalibList, indexTimeCalib = check$indexTimeCalib, primaryIon = check$primaryIon, 
        resolution = check$resolution, prtReaction = check$prtReaction, ptrTransmisison = check$ptrTransmisison, 
        TIC = check$TIC, breathTracer = check$breathTracer, fctFit = check$fctFit, 
        peakShape = check$peakShape, peakList = check$peakList, date = check$date)
    
    # save in Rdata with the name chosen
    if (!is.null(saveDir)) {
        changeName <- parse(text = paste0(setName, "<- ptrSet "))
        eval(changeName)
        eval(parse(text = paste0("save(", setName, ",file= paste0( saveDir,'/', '", 
            setName, ".RData '))")))
    }
    
    return(ptrSet)
}

#' update a ptrSet object 
#'
#' When new files are added to a directory which has already a ptrSet object 
#' associated, run \code{updatePtrSet} to add the new files in the object. 
#' The information on the new files are added to object with the same parameter 
#' used for the function \code{createPtrSet} who has created the object.
#' \code{updatePtrSet} also delete from the ptrSet deleted files in the directory.
#' @param ptrset a ptrset object 
#' @return teh same ptrset object than ininput, but completed with new files 
#' and without deleted files in the directory
#' @export
#' @examples
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' ##add or delete files in the directory 
#' # exhaledPtrset<- updatePtrSet(exhaledPtrset)
updatePtrSet <- function(ptrset) {
    
    if (!methods::is(ptrset, "ptrSet")) 
        stop("ptrset must be a ptrSet object.
                                          Use createPtrSet function")
    
    # get information
    parameter <- getParameters(ptrset)
    sampleMetadata <- getSampleMetadata(ptrset)
    
    # files in the directory
    
    if (methods::is(parameter$dir, "expression")) {
        parameter$dir <- eval(parameter$dir)
        filesDirFullName <- eval(parse(text = "list.files(parameter$dir, 
                                   recursive = TRUE, pattern=\".h5$\",
                                   full.names = TRUE)"))
        filesDirFullNameParam <- parse(text = "list.files(parameter$dir, 
                                 recursive = TRUE, pattern=\".h5$\",
                                 full.names = TRUE)")
    } else {
        filesDirFullName <- list.files(parameter$dir, recursive = TRUE, pattern = "\\.h5$",
            full.names = TRUE)
        filesDirFullNameParam <- filesDirFullName
       
        
    }
    
    # files already processed
    filesProcessed <- rownames(sampleMetadata)
   
    # new files
    newFilesIndex <- which(!basename(filesDirFullName) %in% filesProcessed)
    newFilesFullNames <- filesDirFullName[newFilesIndex]
    
    # files deleted
    deletedFilesIndex <- which(!filesProcessed %in% basename(filesDirFullName))
    deletedFiles <- filesProcessed[deletedFilesIndex]
    
    if (length(newFilesIndex) == 0 & length(deletedFilesIndex) == 0) {
        cat("nothing to update")
        return(ptrset)
    }
    
    # if there is deleted files
    if (length(deletedFilesIndex) > 0) {
        # deleted in sample meta data
        sampleMetadata <- data.frame(sampleMetadata[-deletedFilesIndex, , drop = FALSE])
        
        # deleted in ptrSet
        parameter<-getParameters(ptrset)
        parameter$listFile <- filesDirFullNameParam[]
        ptrset<- setParameters(ptrset,parameter)
        ptrset<-deleteFilePtrSet(ptrset,deletedFiles)
        ptrset<-setSampleMetadata(ptrset,sampleMetadata)
        message(deletedFiles, " deleted \n")
    }
    
    # if there is new files
    if (length(newFilesIndex) > 0) {
        
        # complete sample metaData
        newFilesName <- basename(newFilesFullNames)
        sampleMetadata[newFilesName, ] <- NA
        
        # same order as list.file
        sampleMetadata <- sampleMetadata[basename(filesDirFullName), ,drop=FALSE]
        
        
        # checkset
        check <- checkSet(files = newFilesFullNames, 
                          mzCalibRef = parameter$mzCalibRef, 
                          fracMaxTIC = parameter$timeLimit$fracMaxTIC, 
                          calibrationPeriod = parameter$calibrationPeriod,
                          knotsPeriod =  parameter$knotsPeriod,
                          mzPrimaryIon = parameter$mzPrimaryIon,
                          mzBreathTracer = parameter$mzBreathTracer)
        if (length(check$failed) > 0) {
            filesDirFullName <- filesDirFullName[-which(basename(filesDirFullName) %in% 
                check$failed)]
            filesDirFullNameParam<-filesDirFullName
            sampleMetadata <- sampleMetadata[-which(row.names(sampleMetadata) %in% 
                check$failed), ,drop=FALSE]
        }
        
        # add new file to ptrset and in same order as list.file
        parameter<-getParameters(ptrset)
        parameter$listFile <- filesDirFullNameParam
        ptrset<- setParameters(ptrset,parameter)
        ptrset<- setSampleMetadata(ptrset, as.data.frame(sampleMetadata))
        check$mzCalibRef<-check$mzCalibRefList
        ptrset <-mergedPtrSet(ptrset, check ,orderedFile = filesDirFullName)
       
    }
    
    saveDir <- parameter$saveDir
    objName <- parameter$name
    if (!is.null(saveDir)) {
        if (!is.null(objName)) {
            changeName <- parse(text = paste0(objName, "<- ptrset "))
            eval(changeName)
            eval(parse(text = paste0("save(", objName, ",file= paste0(
                                 saveDir,'/', '", 
                objName, ".RData '))")))
        } else save(ptrset, file = paste0(saveDir, "/ptrSet.RData"))
    }
    
    return(ptrset)
}

checkSet <- function(files, 
                     mzCalibRef, 
                     fracMaxTIC, 
                     mzBreathTracer, 
                     calibrationPeriod, 
                     knotsPeriod, mzPrimaryIon) {
    
    # init output list
    mzCalibRefList <- list()
    signalCalibRef <- list()
    TIC <- list()
    breathTracer <- list()
    timeLimit <- list()
    knots <- list()
    resolution <- list()
    errorCalibPpm <- list()
    peakList <- list()
    coefCalibList <- list()
    indexTimeCalib <- list()
    primaryIon <- list()
    fctFit <- list()
    peakShape <- list()
    transmisison <- list()
    reaction <- list()
    date <- list()
    fileName <- basename(files)
    failed <- c()
    
    # loop over files
    for (j in seq_along(files)) {
        
        # check reading and calibration of file
        raw <- try(readRaw(filePath = files[j], mzCalibRef = mzCalibRef, calibrationPeriod = calibrationPeriod))
        if (attr(raw, "class") == "try-error") {
            message(fileName[j], " opening or calibration failed \n")
            failed <- c(failed, fileName[j])
            next
        }
        
        # calibration infomration
        mzCalibRefList[[fileName[j]]] <- getCalibrationInfo(raw)$calibMassRef
        # signalCalibRef[[ fileName[j] ]] <- getCalibrationInfo(raw)$calibSpectr
        errorCalibPpm[[fileName[j]]] <- getCalibrationInfo(raw)$calibError
        coefCalibList[[fileName[j]]] <- getCalibrationInfo(raw)$calibCoef
        indexTimeCalib[[fileName[j]]] <- getCalibrationInfo(raw)$indexTimeCalib
        
        # estimate the resolution
        calibSpectr <- alignCalibrationPeak(getCalibrationInfo(raw)$calibSpectr, getCalibrationInfo(raw)$calibMassRef, 
                                            length(getRawInfo(raw)$time))
        signalCalibRef[[fileName[j]]] <- calibSpectr
        resolutionEstimated <- try(estimateResol(calibMassRef = getCalibrationInfo(raw)$calibMassRef, calibSpectr = calibSpectr))
        
        if (!is.null(attr(resolutionEstimated, "class"))) {
            message(fileName[j], "estimated resolution failed \n")
            failed <- c(failed, fileName[j])
            mzCalibRefList[[fileName[j]]] <- NULL
            errorCalibPpm[[fileName[j]]] <- NULL
            coefCalibList[[fileName[j]]] <- NULL
            indexTimeCalib[[fileName[j]]] <- NULL
            signalCalibRef[[fileName[j]]] <- NULL
            next
        }
        
        resolution[[fileName[j]]] <- resolutionEstimated
        reaction[[fileName[j]]] <- getPTRInfo(raw)$prtReaction
        transmisison[[fileName[j]]] <- getPTRInfo(raw)$ptrTransmisison
        date[[fileName[j]]] <- getDate(raw)
        resolutionRange <- c(min(resolutionEstimated) * 0.8, mean(resolutionEstimated), 
            max(resolutionEstimated) * 1.2)
        
        # time unit BlockPeriodNS <- try(rhdf5::h5readAttributes(files[j],
        # /TimingData))$BlockPeriod #nano seconde
        
        # TIC
        TIC[[fileName[j]]] <- colSums(getRawInfo(raw)$rawM)
        if (is.null(mzBreathTracer)) {
            breathTracer[[fileName[j]]] <- colSums(getRawInfo(raw)$rawM)
        } else {
            index <- lapply(mzBreathTracer, function(x) {
                th <- 350 * x/10^6
                which(x - th < getRawInfo(raw)$mz & getRawInfo(raw)$mz < x + th)
            })
            breathTracer[[fileName[j]]] <- colSums(getRawInfo(raw)$rawM[unlist(index), ])
        }
        # timeLimit
        
        indLim <- timeLimits(raw, fracMaxTIC = fracMaxTIC, mzBreathTracer = mzBreathTracer)
        timeLimit[[fileName[j]]] <- indLim
        t <- getRawInfo(raw)$time
        
        # default knots location every 3 second on expiration
        if (knotsPeriod == 0) {
            knots <- list(NULL)
        } else {
            background <- indLim$backGround
            knots[[fileName[j]]] <- try(defineKnotsFunc(t, background, knotsPeriod, 
                method = "expiration", fileName[j]))
        }
        
        l.shape <- determinePeakShape(raw)$peakShapemz
        peakShape[[fileName[j]]] <- l.shape
        
        # check best fit
        sech2 <- mean(PeakList(raw, mzNominal = getCalibrationInfo(raw)$calibMassRef, 
                               fctFit = "sech2", 
            maxIter = 1, ppm = 500, minIntensityRate = 0.2, 
            windowSize = 0.2, resolutionRange = resolutionRange, 
            peakShape = l.shape)$peak$R2)
        
        averagePeak <- mean(PeakList(raw, mzNominal = getCalibrationInfo(raw)$calibMassRef, fctFit = "averagePeak", 
            resolutionRange = resolutionRange, maxIter = 1, peakShape = l.shape, 
            ppm = 500, minIntensityRate = 0.2, windowSize = 0.2)$peak$R2)
        asymGauss <- mean(PeakList(raw, mzNominal = getCalibrationInfo(raw)$calibMassRef, fctFit = "asymGauss", 
            resolutionRange = resolutionRange, maxIter = 1, peakShape = l.shape, 
            ppm = 500, minIntensityRate = 0.2, windowSize = 0.2)$peak$R2)
        
        fctFit[[fileName[j]]] <- c("sech2", "averagePeak", "asymGauss")[which.max(c(sech2, 
            averagePeak, asymGauss))]
        
        # check if mass 21 contains in mass axis
        if (!21 %in% unique(round(getRawInfo(raw)$mz))) {
            warning("mass 21 not in mass Axis,the ppb quantification can not be done.")
            primaryIonV <- NA
        } else {
            p <- try(PeakList(raw, mzNominal = round(mzPrimaryIon), ppm = 700, minIntensityRate = 0.5, 
                minIntensity = 0, maxIter = 1, thNoiseRate = 0, fctFit = fctFit[[fileName[j]]], 
                peakShape = l.shape, windowSize = 0.2))
            
            if(is.null(attr(p,"condition"))){
                primaryIndex <- which.min(abs(p$peak$Mz - mzPrimaryIon))
                if (length(primaryIndex) )
                    primaryIonV <- p$peak$quanti_cps[primaryIndex] else primaryIonV <- NA
            } else primaryIonV <- NA
            #primaryIndex <- which(abs(p$peak - mzPrimaryIon) * 10^6/21 < 200)
        }
        if (!38 %in% unique(round(getRawInfo(raw)$mz))) {
            waterCluster <- NA
        } else {
            p <- try(PeakList(raw, mzNominal = c(38), ppm = 700, minIntensityRate = 0.5, 
                minIntensity = 0, maxIter = 1, thNoiseRate = 0, fctFit = fctFit[[fileName[j]]], 
                peakShape = l.shape, windowSize = 0.2))
            
            if(is.null(attr(p,"condition"))){
                clusterIndex <- which(abs(p$peak$Mz - 38.03) * 10^6/21 < 200)
                if (length(clusterIndex) != 0) 
                    waterCluster <- p$peak$quanti_cps[clusterIndex] else waterCluster <- NA
            } else  waterCluster <- NA
            
        }
        primaryIon[[fileName[j]]] <- list(primaryIon = primaryIonV, waterCluster = waterCluster)
        
        message(fileName[j], " check")
        rm(raw)
    }
    
    return(list(mzCalibRefList = mzCalibRefList, signalCalibRef = signalCalibRef, 
        TIC = TIC, breathTracer = breathTracer, timeLimit = timeLimit, knots = knots, 
        indexTimeCalib = indexTimeCalib, resolution = resolution, errorCalibPpm = errorCalibPpm, 
        fctFit = fctFit, peakShape = peakShape, peakList = peakList, coefCalibList = coefCalibList, 
        primaryIon = primaryIon, prtReaction = reaction, ptrTransmisison = transmisison, 
        date = date, failed = failed))
}

#' Convert a h5 file to mzML
#'
#' convert_to_mzML create a mzML file from a h5 file in the same directory 
#' with the \code{writeMLData} of the \code{MSnbase} package
#' @param file A .h5 file path
#' @return create a mzML file in the same directory of the h5 input file
#' @examples 
#' library(ptairData)
#' filePathRaw <- system.file('extdata/exhaledAir/ind1', 'ind1-1.h5', 
#' package = 'ptairData')
#' # write a mzml file in the same directory
#' convert_to_mzML(filePathRaw)
#' file_mzML <- gsub(".h5", ".mzML", filePathRaw)
#' file.remove(file_mzML)
#' @export
convert_to_mzML <- function(file) {
    data <- readRaw(file)
    mzVn <- getRawInfo(data)$mz
    timVn <- getRawInfo(data)$time
    rawMN <- getRawInfo(data)$rawM
    mzML <- apply(rawMN, 2, function(x) matrix(c(mzVn[x != 0], x[x != 0]), ncol = 2))
    pk_count <- vapply(mzML, function(x) dim(x)[1], FUN.VALUE = 1)
    
    hdr <- data.frame(matrix(0, ncol = 26, nrow = dim(rawMN)[2]))
    names(hdr) <- c("seqNum", "acquisitionNum", "msLevel", "polarity", "peaksCount", 
        "totIonCurrent", "retentionTime", "basePeakMZ", "basePeakIntensity", "collisionEnergy", 
        "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum", "precursorMZ", 
        "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
        "mergedResultStartScanNum", "mergedResultEndScanNum", "injectionTime", "filterString", 
        "spectrumId", "centroided", "ionMobilityDriftTime")
    hdr$seqNum <- seq(1, dim(rawMN)[2])
    hdr$acquisitionNum <- seq(1, dim(rawMN)[2])
    hdr$msLevel <- rep(1, dim(rawMN)[2])
    hdr$polarity <- rep(1, dim(rawMN)[2])
    hdr$peaksCount <- pk_count  #nombre de masse detecter par rt
    hdr$retentionTime <- timVn  #retention time
    hdr$totIonCurrent <- colSums(rawMN)  #TIC
    # les mz des plus grande intensitees par rt
    hdr$basePeakMZ <- apply(rawMN, 2, function(x) mzVn[which.max(x)])
    hdr$basePeakIntensity <- apply(rawMN, 2, max)  # max des intensite par rt
    # plus petite detecte mz par rt
    hdr$lowMZ <- vapply(mzML, function(x) min(x[, 1]), FUN.VALUE = 1)
    # plus grande mz detecte par rt
    hdr$hightMZ <- vapply(mzML, function(x) max(x[, 1]), FUN.VALUE = 1)
    hdr$filterString <- rep("NA", dim(rawMN)[2])
    hdr$centroided <- rep(FALSE, dim(rawMN)[2])
    hdr$ionMobilityDriftTime <- rep(1, dim(rawMN)[2])
    hdr$ionMobilityDriftTime <- rep(as.numeric(NA), dim(rawMN)[2])
    for (j in seq_len(ncol(rawMN))) {
        hdr$spectrumId[j] <- paste("spectrum=", j, sep = "")
    }
    file_mzML <- gsub(".h5", ".mzML", file)
    MSnbase::writeMSData(object = mzML, file = file_mzML, header = hdr)
}




# data("exhaledPtrset")
# exhaledPtrset@prtReaction<-lapply(exhaledPtrset@prtReaction,function(x) {
#     x$TwInfo[3]<-"T-Drift[degree C]"
#     x})
# save(exhaledPtrset,version = 2,
#      file = "C:/Users/CR258086/Documents/Source/ptairMS/data/exhaledPtrset.RData")

