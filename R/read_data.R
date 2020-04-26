#' Read a h5 file of PTR-TOF-MS data
#'
#' \code{readRaw} reads a h5 file with rhdf5 library, and calibrates the mass axis with \code{\link[ptairMS]{calibration}} function.
#' It returns a \code{\link[ptairMS]{ptrRaw-class}} S4 object that contains raw data.
#' 
#' @param filePath h5 file path full name.
#' @param calibTIS boolean. If true, the function \code{calibration} is apply on the Total Ion Spectrum
#' with the \code{mzCalibRef} reference masses. 
#' @param mzCalibRef calibration parameter. Vector of exact theoritical mass values of a intensive peak of moleculs
#'  and 'unique' in a nominal mass interval
#' @param tolCalibPpm calibration parameter. The maximum error tolerated in ppm. A warning appears for 
#' error graeter than \code{tolCalibPpm}.
#' @return a ptrRaw object, including slot \itemize{
#' \item rawM the data raw matrix, in count of ions  
#' \item mz the mz axis
#' \item time time acquisition in second 
#' } 
#' @examples
#' library(ptairData)
#' filePathRaw <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
#' raw <- readRaw(filePathRaw,mzCalibRef=c(21.022, 41.03858, 60.0525))
#' @import bit64
#' @export
readRaw <- function(filePath, calibTIS=TRUE, 
                    mzCalibRef = c(21.022, 29.013424,41.03858, 60.0525,203.943, 330.8495),
                    tolCalibPpm=70){

  if(is.null(filePath) | filePath=="") stop("filePath is empty")
  if(is.na(filePath)) stop("filePath is empty")
  if(!is.character(filePath)) stop("filePath must be a character")
  if(length(filePath) >1) stop("Only one filePath is required")
  if(!grepl("*h5$",filePath)) stop("The file is not in h5 format")
  
  #find the longer of acquisitionTime
  file <- rhdf5::H5Fopen(filePath)
  time <- rhdf5::H5Dopen(file,"/TimingData/BufTimes")
  acTime <- rhdf5::H5Sget_simple_extent_dims(rhdf5::H5Dget_space(time))$size
  nbrBuf<-acTime[1]
  nbrWrite <- acTime[2]
  
  if(nbrWrite == 0) { stop("The file is empty (0 seconde of acquisition)")}
  
  # acquiqision time limit to 900 spectra, indeed the rawAn is more than 3.4 GB
  NbrWriteMax <- ceiling(3000/nbrBuf)
  rhdf5::h5closeAll()
  
  #read information needed
  date_heure<-rhdf5::h5read(filePath, "/AcquisitionLog", bit64conversion='bit64')$Log$timestring[1]
  
  timVn <- rhdf5::h5read(filePath, "/TimingData/BufTimes", bit64conversion='bit64',
                         index = list(NULL, seq_len(min(nbrWrite,NbrWriteMax))))
  rawAn <- rhdf5::h5read(filePath, "/FullSpectra/TofData",  
                         index = list(NULL, NULL,NULL, seq_len(min(nbrWrite,NbrWriteMax))), bit64conversion='bit64')
  mzVn <- rhdf5::h5read(filePath, "FullSpectra/MassAxis", bit64conversion='bit64')
  reaction<- try(rhdf5::h5read(filePath,"AddTraces/PTR-Reaction"))
  transmission <- try(rhdf5::h5read(filePath,"PTR-Transmission"))
  calibCoef <- try(rhdf5::h5read(filePath,"FullSpectra/MassCalibration",index=list(NULL,1)))
  
  attributCalib <- try(rhdf5::h5readAttributes(filePath,"/FullSpectra"))
  #singleIon<-attributCalib$`Single Ion Signal`
  #sampleInterval<-attributCalib$SampleInterval
  
  if(!is.null(attr(transmission,'condition'))) transmission <- matrix(0) else transmission<-transmission$Data
  if(!is.null(attr(reaction,'condition'))) reaction <- list()
  
  # re format
  timVn <-c(timVn)
  mzVn <- c(mzVn)
  rawMn <- matrix(rawAn,
                   nrow = dim(rawAn)[1],
                   ncol = prod(utils::tail(dim(rawAn),2)),
                   dimnames = list(mzVn, timVn))  
  index_zero<-which(timVn==0)[-1] # remove index where timVn =0 except the first time
  rm(rawAn)
  if(length(index_zero)!=0){
    timVn<- timVn[-index_zero]
    rawMn<- rawMn[,-index_zero]
  }
  
  #count ion conversion
  #factor<-sampleInterval*1e9 /singleIon #* (SampleInterval,ns) / (single ion signal mV.ns) for convert to number of ion
  #rawMn<-rawMn*as.numeric(factor) 
  
  # calibration infomration
  if(is.null(attr(calibCoef,'condition'))){
    rownames(calibCoef)<-c("a","b")
    calib_formula <- function(tof,calibCoef) ((tof - calibCoef['b',]) / calibCoef['a',]) ^ 2
    calib_invFormula <- function(m,calibCoef) sqrt(m)*calibCoef['a',] + calibCoef['b',]
    calibMassRef = c(attributCalib$`MassCalibration m1`, attributCalib$`MassCalibration m2`)
    if(is.null(calibMassRef)) calibMassRef <- 0 } else { 
      calibCoef <- matrix(0)
      calib_formula <- function(x) NULL
      calib_invFormula <- function(x) NULL
      calibMassRef <- 0
  }
  
  # write ptrRaw objet
  raw <- methods::new(Class = "ptrRaw", name= basename(filePath),rawM= rawMn, mz=mzVn, time=timVn,calibCoef=calibCoef,
            calibMzToTof = calib_invFormula, calibToftoMz = calib_formula, calibError=0,
            calibMassRef= calibMassRef, calibSpectr= list(NULL),
            ptrTransmisison = transmission, prtReaction= reaction,date = date_heure)
  
  if(calibTIS){
    raw <- calibration(raw, mzCalibRef , tol= tolCalibPpm)
  }
  
  return(raw)
}

#' Create a ptrSet object form a directoy
#'
#' This function creates a \code{\link[ptairMS]{ptrSet-class}} S4 object. It opens each file and: 
#' \itemize{
#' \item performs a calibration with \code{\link[ptairMS]{calibration}} function
#' \item quantify the primary ion (H30+) with the isotope H3180+ mz 21.022 on the average total ion spectrum. 
#' \item if \code{fracMaxTIC} is different to zero, calculates boundaries on the 
#' Total Ion Spectrum (TIC) part of the spectrum where 
#' the TIC intensity is higher than \code{fracMaxTIC * max(TIC)} to identify expirations or headspace analyze
#' \item provides a default sampleMetadata based on the tree structure of the directory 
#' (a \code{data.frame} with file names in rowname, and subfolder in column) 
#' \item If \code{saveDir} is not null, the returned object will be saved in \code{.Rdata} in \code{saveDir} with 
#' \code{setName} as name}
#'
#' @param dir character. a directory path which contains several h5 files, possibly organized in subfolder 
#' @param setName character. The name of the ptrSet object. If `saveDir` is not null, the object will be saved with this name.
#' @param mzCalibRef Vector of accurate mass values of intensive peaks and 'unique' in a 
#' nominal mass interval (without overlapping peaks) to performs calibration. See \code{ ?calibration}.
#' @param fracMaxTIC Percentage (between 0 and 1) of the maximum of the Total Ion Chromatogram (TIC) 
#' amplitude with baseline removal. We will analyze only the part of the spectrum where 
#' the TIC intensity is higher than `fracMaxTIC * max(TIC) `. If you want to analyze the entire spectrum, 
#' set this parameter to 0. 
#' @param saveDir The directory where the ptrSet object will be saved in .RData. If NULL, nothing will be saved.
#' @return a ptrSet objet with slot :
#' \itemize{
#' \item Parameter: a list containing the parameters used for \code{createPrtSet}, \code{detectPeak} and 
#' \code{alignTimePeriods} functions.  
#' \item sampleMetadata: a data.frame containing infomration on the data, with file names in row names
#' \item mzCalibRef: a list containing for each file the masses used for the calibration 
#' (see \code{?ptairMS::calibration} for more details)
#' \item signalCalibRef: the mz and intensity +- 0.4Da around the calibration masses
#' \item errorCalibPpm: a list containing for each file the accuracy error in ppm at each calibration masses
#' \item coefCalib: a list containing the calibration coeficients 'a' and 'b' that permit to convert tof to mz
#' for each file (see \code{?ptairMS::calibration} for more details.  
#' \item resolution: the estimated resolution \eqn{m / \Delta m} for each calibration masses within each file
#' \item TIC: The Total Ion Chromatogram for each file
#' \item timeLImit: a list containing list of two element per file: the matrix of time limit for each file 
#' (if \code{fracMaxTIC} is different to zero), and the background index. See ?ptairMS::timeLimits for more details
#' \item peakListRaw: a list containing for each file a peak list per time interval and background. 
#' (the output of \code{detectPeak} function for each file, see \code{?ptairMS::detectPeak})
#' \item peakListAligned: one peak list for each file, obtain after alignment betewwen all time period and backgroud.
#' (the output of \code{alignTimePeriods} function)
#' }
#'
#' @export 
#' @examples
#' library(ptairData)
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' ptrSet<-createPtrSet(directory,setName="ptrSet",mzCalibRef=c(21.022,59.049),
#' fracMaxTIC=0.9,saveDir= NULL)
createPtrSet<-function(dir, setName,
                       mzCalibRef= c(21.022, 29.013424, 41.03858, 59.049141,75.04406, 203.943, 330.8495),
                       fracMaxTIC=0.5,mzBreathTracer=NULL,
                       saveDir=NULL){
  
  # test on parameter
  # dir exist
  if(dir=="") stop("dir is empty")
  if(!dir.exists(dir)) stop(paste( dir ," does not exist"))
  if(!is.null(saveDir)) {if(!dir.exists(saveDir)) stop(paste(saveDir ," does not exist"))}
  #setName a character
  if(!is.character(setName)) stop("setName must be a character")
  #fracMaxTIC must be between 0 and 1
  if(fracMaxTIC<0 | fracMaxTIC >1) stop("fracMaxTIC must be between zero and 1, see ?createPtrSet")
  #mzCalibRef numeric
  if(!is.numeric(mzCalibRef)) stop("mzCalibRef must be a numeric vector")
  
  # get h5 files name
  filesFullName <- list.files(dir, recursive = TRUE, pattern="\\.h5$",full.names = TRUE)
  fileDir <- dirname(list.files(dir, recursive = TRUE, pattern="\\.h5$"))
  fileName <- basename(filesFullName)
  
  # save parameters
  parameter <- list(dir = dir, listFile= filesFullName, name = setName, mzCalibRef = mzCalibRef,  
                    timeLimit= list(fracMaxTIC = fracMaxTIC), mzBreathTracer=mzBreathTracer,saveDir = saveDir)
  
  # create sampleMetadata 
  # test if there is subfolder
  if( all(fileDir == ".")){ # if there is no subfolder
    sampleMetadata <- data.frame( row.names = fileName)
  } else {
    subfolder <- strsplit(fileDir, "/")
    # max depth of subfolder
    nSubfolder <- max ( vapply(subfolder,length,0)) 
    # set every element at the same length
    subfolder<- lapply(subfolder,
                       function(x) {
                         length(x) <-nSubfolder
                         return(x)} 
    )  
    group <- do.call(rbind, subfolder)
    sampleMetadata <- data.frame(subfolder=group, row.names = fileName)
  }
  
  # checkSet
  check <- checkSet(filesFullName, mzCalibRef, fracMaxTIC,mzBreathTracer) 
  
  if(length(check$failed) > 0) {
    parameter$listFile <- parameter$listFile[ - which(
    basename(parameter$listFile) %in% check$failed )]
    sampleMetadata<-sampleMetadata[ - which(
     row.names(sampleMetadata) %in% check$failed ),,drop=FALSE]
  }
  
  # create ptrSet object
  ptrSet <- methods::new(Class = "ptrSet", parameter = parameter, sampleMetadata = sampleMetadata,
              mzCalibRef= check$mzCalibRefList, timeLimit = check$timeLimit, signalCalibRef = check$signalCalibRef, 
              errorCalibPpm= check$errorCalibPpm, coefCalib= check$coefCalibList, primaryIon=check$primaryIon,
              resolution = check$resolution, prtReaction= check$prtReaction,ptrTransmisison=check$ptrTransmisison,
              TIC = check$TIC, breathTracer=check$breathTracer,peakListRaw = check$peakListRaw, 
              peakListAligned = check$peakListAligned,date=check$date)
  
  #save in Rdata with the name choosen 
  if(!is.null(saveDir)){
    changeName <- parse(text=paste0(setName,"<- ptrSet "))
    eval(changeName)
    eval(parse(text =  paste0( "save(" ,setName ,",file= paste0( saveDir,'/', '",setName,".RData '))")))
  }
 
  return(ptrSet)
}

#' update a ptrSet object 
#'
#' When new files are added to a directory which has already a ptrSet object associated, run \code{updatePtrSet}
#' for add the new files to the object. The information on the new files are added to
#' object with the same parameter used for the function \code{createPtrSet} who has created the object.
#' \code{updatePtrSet} also delete from the ptrSet deleted files in the directory.
#' @param ptrset a ptrset object 
#' @return teh same ptrset object than ininput, but completed with new files and without deleted files in the directory
#' @export
#' @examples
#' library(ptairData)
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' mycobacteria <- createPtrSet(dir= directory, setName="mycobacteria",mzCalibRef= c(21.022, 59.049141))
#' ##add or delete files in the directory 
#' mycobacteria<- updatePtrSet(mycobacteria)
updatePtrSet<-function(ptrset){

  if(! methods::is(ptrset,"ptrSet")) stop("ptrset must be a ptrSet object. Use createPtrSet function")
  
  # get information 
  parameter <- ptrset@parameter
  sampleMetadata <- ptrset@sampleMetadata
  
  # files in the diretcory 
  filesDirFullName <- list.files(parameter$dir, recursive = TRUE, pattern="\\.h5$",full.names = TRUE)
  
  # files already processed
  filesProcessed <- rownames(sampleMetadata)
  
  # new files
  newFilesIndex <- which(! basename(filesDirFullName) %in% filesProcessed) 
  newFilesFullNames <-filesDirFullName[newFilesIndex]
 
  # files deleted 
  deletedFilesIndex <-  which(!filesProcessed %in%  basename(filesDirFullName))
  deletedFiles <- filesProcessed[deletedFilesIndex]
  
  if(length(newFilesIndex) ==0 & length(deletedFilesIndex)==0) { 
    cat("nothing to update")
    return(ptrset)
    }
  
  #if there is deleted files
  if(length(deletedFilesIndex) >0) {
    #deleted in sample meta data
    sampleMetadata<-data.frame(sampleMetadata[-deletedFilesIndex,,drop=FALSE] )
    
    #deleted in ptrSet
    ptrset@sampleMetadata <- sampleMetadata
    ptrset@mzCalibRef[deletedFiles] <-  NULL
    ptrset@signalCalibRef[deletedFiles] <- NULL
    ptrset@errorCalibPpm[deletedFiles] <-  NULL
    ptrset@coefCalib[deletedFiles] <-  NULL
    ptrset@primaryIon[deletedFiles]<-NULL
    ptrset@resolution[deletedFiles] <-  NULL
    ptrset@TIC[deletedFiles] <-  NULL
    ptrset@timeLimit[deletedFiles] <-  NULL
    ptrset@peakListRaw[deletedFiles] <-  NULL
    ptrset@peakListAligned[deletedFiles] <-  NULL
    ptrset@parameter$listFile <-filesDirFullName
    ptrset@breathTracer[deletedFiles]<-NULL
    ptrset@ptrTransmisison[deletedFiles]<-NULL
    ptrset@prtReaction[deletedFiles]<-NULL
    ptrset@date[deletedFiles]<-NULL
    
    message(paste(deletedFiles," deleted \n"))
  }
  
  #if there is new files
  if(length(newFilesIndex)>0){
    
    # complete sample metaData
    newFilesName<-basename(newFilesFullNames)
    sampleMetadata[newFilesName,]<-NA
    
    # same order as list.file
    sampleMetadata<-sampleMetadata[basename(filesDirFullName),]
    
    
    #checkset 
    check <- checkSet(newFilesFullNames,mzCalibRef = parameter$mzCalibRef,fracMaxTIC = parameter$timeLimit$fracMaxTIC,
                      mzBreathTracer = parameter$mzBreathTracer)
    if(length(check$failed) > 0) {
      filesDirFullName <- filesDirFullName[ - which(
      basename(filesDirFullName) %in% check$failed )]
      sampleMetadata<-sampleMetadata[ - which(
        row.names(sampleMetadata) %in% check$failed ),]}
    
    # add new file to ptrset and in same order as list.file
    ptrset@parameter$listFile <- filesDirFullName
    ptrset@sampleMetadata <- as.data.frame(sampleMetadata)
    ptrset@mzCalibRef <-  c(ptrset@mzCalibRef,check$mzCalibRefList)[basename(filesDirFullName)]
    ptrset@signalCalibRef <- c(ptrset@signalCalibRef,check$signalCalibRef)[basename(filesDirFullName)]
    ptrset@errorCalibPpm <- c(ptrset@errorCalibPpm ,check$errorCalibPpm)[basename(filesDirFullName)]
    ptrset@coefCalib <- c(ptrset@coefCalib, check$coefCalibList)[basename(filesDirFullName)]
    ptrset@resolution <- c(ptrset@resolution, check$resolution)[basename(filesDirFullName)]
    ptrset@TIC <- c(ptrset@TIC, check$TIC)[basename(filesDirFullName)]
    ptrset@timeLimit <-c(ptrset@timeLimit,check$timeLimit)[basename(filesDirFullName)]
    ptrset@primaryIon<-c( ptrset@primaryIon,check$primaryIon)[basename(filesDirFullName)]
    ptrset@breathTracer<-c(ptrset@breathTracer,check$breathTracer)[basename(filesDirFullName)]
    ptrset@ptrTransmisison<-c(ptrset@ptrTransmisison,check$ptrTransmisison)[basename(filesDirFullName)]
    ptrset@prtReaction <-c(ptrset@prtReaction,check$prtReaction)[basename(filesDirFullName)]
    ptrset@date <-c(ptrset@date,check$date)[basename(filesDirFullName)]
  }
  
  saveDir<-parameter$saveDir
  objName<-parameter$name
  if(!is.null(saveDir)){
    if(!is.null(objName)){
      changeName <- parse(text=paste0(objName,"<- ptrset "))
      eval(changeName)
      eval(parse(text =  paste0( "save(" ,objName ,",file= paste0( saveDir,'/', '",objName,".RData '))")))
    } else save(ptrset, file=paste0(saveDir,"/ptrSet.RData"))
  }
  
  return(ptrset)
}

#' check all files in a directory
#' 
#' main function of createPtrSet: opens each files, performs calibration, calculates TIC boundaries and
#' quantify primary ion . 
#' @param files list of full names files to chexk
#' @param mzCalibRef Vector of accurate mass values of intensive peaks and 'unique' in a 
#' nominal mass interval (without overlapping)
#' @param fracMaxTIC timeLimits argument
#' @return list containing in ptrSet object slot
checkSet <- function(files, mzCalibRef , fracMaxTIC, mzBreathTracer){
  
  # init output list
  mzCalibRefList <- list()
  signalCalibRef <- list()
  TIC <- list()
  breathTracer<-list()
  timeLimit <- list()
  resolution <- list()
  errorCalibPpm<-list()
  peakListRaw <- list()
  peakListAligned<-list()
  coefCalibList<-list()
  primaryIon<-list()
  transmisison<-list()
  reaction<-list()
  date<-list()
  fileName<-basename(files)
  failed<-c()
  
  # loop over files
  for (j in seq_along(files)){
   
    # check reading and calibration of file
    raw <- try(readRaw(files[j], mzCalibRef = mzCalibRef) ) ## files en full name
    if(attr(raw,"class")== "try-error" ){
      message( paste(fileName[j]," opening or calibration failed"))
      failed<-c(failed,fileName[j])
      next 
    }
    
    # calibration infomration
    mzCalibRefList[[ fileName[j] ]] <- raw@calibMassRef
    signalCalibRef[[ fileName[j] ]] <- raw@calibSpectr
    errorCalibPpm[[ fileName[j] ]] <- raw@calibError
    coefCalibList[[ fileName[j] ]] <-raw@calibCoef
    resolution[[ fileName[j] ]] <- estimateResol(raw@calibMassRef,raw@calibSpectr)
    reaction[[ fileName[j] ]] <- raw@prtReaction
    transmisison[[ fileName[j]  ]]<- raw@ptrTransmisison
    date[[ fileName[j]  ]] <- raw@date
    
    #primary ion quantification
    
    #time unit 
    #BlockPeriodNS <- try(rhdf5::h5readAttributes(files[j],"/TimingData"))$BlockPeriod #nano seconde
    
    # TIC
    TIC[[fileName[j]]] <- colSums(raw@rawM)
    if(is.null(mzBreathTracer)){ 
      breathTracer[[fileName[j]]]<-colSums(raw@rawM)} 
    else {
      index<- lapply(mzBreathTracer, function(x) {
        th<-350*x/10^6
        which( x - th < raw@mz & raw@mz < x + th)
      })
      breathTracer[[fileName[j]]] <- colSums(raw@rawM[unlist(index),]) 
    }
    # timeLimit
    indLim <- timeLimits(raw, fracMaxTIC = fracMaxTIC,traceMasses = mzBreathTracer)
    timeLimit[[ fileName[j] ]] <- indLim
    
    # check if mass 21 contains in mass axis
    if( ! 21 %in% unique(round(raw@mz)) ) {
      warning("mass 21 not in mass Axis,the ppb quantification can not be done.")
      primaryIon[[ fileName[j] ]] <- NA
    } else {
      ## sur le background
      if(! is.null(indLim$backGound)) indBg<- indLim$backGound else indBg<-seq_along(raw@time)
      raw.bg<-raw
      raw.bg@rawM<-raw.bg@rawM[,indBg]
      raw.bg@time<-raw.bg@time[indBg]
      p<-PeakList(raw.bg,mz=21,ppm = 500,thIntensityRate = 0.5)
      primaryIon[[ fileName[j] ]]<- sum(p$peak$quanti_cps)
    }
    
    
    message( paste(fileName[j]," check"))
  }  

  return(list(mzCalibRefList=mzCalibRefList,
              signalCalibRef=signalCalibRef,
              TIC=TIC ,
              breathTracer=breathTracer,
              timeLimit=timeLimit,
              resolution=resolution,
              errorCalibPpm=errorCalibPpm,
              peakListRaw=peakListRaw ,
              peakListAligned=peakListAligned,
              coefCalibList=coefCalibList,
              primaryIon=primaryIon, 
              prtReaction= reaction,
              ptrTransmisison= transmisison,
              date=date,
              failed=failed))
}

#' Convert a h5 file to mzML
#'
#' convert_to_mzML create a mzML file from a h5 file in the same directory with the \code{writeMLData} of
#' the \code{MSnbase} package
#' @param file A .h5 file path
#' @return cretae a mzML file in the same directory of the h5 iput file
#' @examples 
#' library(ptairData)
#' filePathRaw <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
#' \dontrun{convert_to_mzML(filePathRaw)}
#' @export
convert_to_mzML<-function(file){
  data<-readRaw(file)
  mzVn<-data$mz
  timVn<-data$time
  rawMN<-data$rawdata
  mzML<-apply(rawMN,2,function(x) matrix(c(mzVn[x!=0],x[x!=0]),ncol=2))
  pk_count<- vapply(mzML, function(x) dim(x)[1],FUN.VALUE = 1)

  hdr<-data.frame(matrix(0,ncol=26,nrow=dim(rawMN)[2]))
  names(hdr)<-c("seqNum","acquisitionNum","msLevel","polarity","peaksCount",
                "totIonCurrent","retentionTime","basePeakMZ","basePeakIntensity",
                "collisionEnergy","ionisationEnergy","lowMZ", "highMZ","precursorScanNum"
                ,"precursorMZ","precursorCharge","precursorIntensity","mergedScan"
                ,"mergedResultScanNum","mergedResultStartScanNum","mergedResultEndScanNum","injectionTime"
                ,"filterString","spectrumId","centroided","ionMobilityDriftTime")
  hdr$seqNum<-seq(1,dim(rawMN)[2])
  hdr$acquisitionNum<-seq(1,dim(rawMN)[2])
  hdr$msLevel<-rep(1,dim(rawMN)[2])
  hdr$polarity<-rep(1,dim(rawMN)[2])
  hdr$peaksCount<-pk_count #nombre de masse detecter par rt
  hdr$retentionTime<-timVn #retention time
  hdr$totIonCurrent<-colSums(rawMN) #TIC
  hdr$basePeakMZ<-apply(rawMN,2,function(x) mzVn[which.max(x)]) #les mz des plus grande intensitées par rt
  hdr$basePeakIntensity<-apply(rawMN,2,max) # max des intensité par rt
  hdr$lowMZ<-vapply(mzML,function(x) min(x[,1]),FUN.VALUE = 1) # plus petite détécté mz par rt
  hdr$hightMZ<-vapply(mzML,function(x) max(x[,1]),FUN.VALUE = 1) # plus grande mz détécté par rt
  hdr$filterString<-rep("NA",dim(rawMN)[2])
  hdr$centroided<-rep(FALSE,dim(rawMN)[2])
  hdr$ionMobilityDriftTime<-rep(1,dim(rawMN)[2])
  hdr$ionMobilityDriftTime<-rep(as.numeric(NA),dim(rawMN)[2])
  for (j in seq_len(ncol(rawMN))){
    hdr$spectrumId[j]<-paste("spectrum=",j,sep="")
  }
  file_mzML<-gsub(".h5",".mzML",file)
  MSnbase::writeMSData(object = mzML, file = file_mzML, header = hdr)
}
