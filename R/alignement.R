#' Alignment with kernel gaussian density
#'
#' @param peakTab table with comlumn : mass, quantification, and groups number to aligned
#' @param ppmGroup width of sub group created beafore density estimation in ppm
#' @param dmzGroup width of sub group created beafore density estimation in Da
#' @return A list containing groups formed by alignment.
align <- function(peakTab, ppmGroup=70, dmzGroup=0.001){
  ### faire un test si c'est null renvoyer nul 
  porder <- order( peakTab[ , "mz"] )
  peaksL <- peakTab[ porder , ,drop=FALSE]
  
  #delete the unit in col names
  colnames(peaksL)[grep("quanti",colnames(peaksL))]<-"quanti"
  legList <- colnames(peaksL)
  
  # Get the number of different values in group
  n.group <- nlevels( as.factor( peaksL[, "group"] ) ) 
  if(n.group==0) return(list(peakTab)) #there is no peak 
  
  # Cuts the mass axis into 1 Da size intervals around the nominal mass
  
  nPoints<-1024
  mzrange <- range(peaksL[, 1])
  #inter <- mzrange[2] * ppmGroup * nPoints / (20 * 10 ^ 6) ###### voir profia
  #massInter <- seq( mzrange[1], mzrange[2], inter/2 )
  inter=1
  massInter <- seq( floor(mzrange[1]) - 0.5, ceiling(mzrange[2]) + 0.5, inter )
  
  # Find the position of peaksL in massInter
  masspos <- findEqualGreaterM(peaksL[, 1], massInter)

  # Initialize variables
  groupPeakList <- list() # The list who will be return, contains each group of same masses 
  numGroupPeak <- 0 # the number of group formed

  # Loop on each group
  index_group <- which( diff(masspos,lag=1) >= 1 )
  for (i in index_group) {

    start <- masspos[i] # Get start and end of the sub group
    end <- masspos[i + 1] - 1
    #end <- masspos[i + 2] - 1
    
    subpeakl <- peaksL[start:end, , drop=FALSE] # Take submatrix

    # calculated the density 
    bw <- max(massInter[i] * ppmGroup / (10 ^ 6), dmzGroup) # Determining the base of the signal to skip for a resolution.
    den <- density(subpeakl[, "mz"], bw,
                   from = massInter[i] -  bw , to = massInter[i + 1] + bw, n = 1024)
    maxy <- max(den$y)
    den$y[which(den$y <= bw * maxy)] <- 0
    plim <- list(linflex=-1, rinflex=0, state=0)
    
    
    # While we detect peak in the density 
    repeat {
      plim <- findLimDensity(den$y, plim$rinflex + 1, plim$state)
      if (plim$linflex == plim$rinflex)  break ###In this case the last peak is found.
      
      subpeakl.group <- makeSubGroup (subpeakl, den, plim)
      if(nrow(subpeakl.group)!=0){
        numGroupPeak <- numGroupPeak + 1
        groupPeakList[[numGroupPeak]] <-subpeakl.group
      }
        
      
    } # end of repeat
  } # end of error
  
  return(groupPeakList)
  
  }

#' Use in align function. return a peak group thanks to kernel gaussian density in a peak matrix.
#' @param subpeakl a matrix with mz, ppb, delta_mz, resolution and group column. 
#' @param den the kernel gaussian density estimated on subpeakl
#' @param plim the limit of a peak in the density of the group who will pe fromed
#' @return the sub peakgroup  
makeSubGroup <- function(subpeakl,den,plim) {
  
  selectedPGroup <- which(subpeakl[, "mz"] >= den$x[plim$linflex] &
                            subpeakl[, "mz"] <= den$x[plim$rinflex])
  
  subpeakl.group<-subpeakl[selectedPGroup,,drop=FALSE]
  
  ## if two varaibles of same level group are confused we keep the more intense one for mz center and sum theme
  if( ( length(unique(subpeakl.group[,"group"])) < dim(subpeakl.group)[1] ) ){
    subpeakl.group <- data.table::as.data.table(subpeakl.group)
    if("background" %in% colnames(subpeakl.group)){
      subpeakl.group<-as.matrix(subpeakl.group[,.(mz=mz[which.max(quanti)], 
                                                  quanti=sum(quanti),
                                                  background =  sum(background),
                                                  delta_mz=max(delta_mz)), 
                                               by=group] )
    } else subpeakl.group<-as.matrix(subpeakl.group[,.(mz=mz[which.max(quanti)], 
                                                       quanti=sum(quanti),
                                                       delta_mz=max(delta_mz)), 
                                                    by=group] )
    
    
  }
  return(subpeakl.group)
  }


#' Align mass tables of different expirations 
#' 
#' @param PeakMat A peakMatrix, with at least column named 
#' "mz" correspond to masses, and an other named "group" correspond to the expiration number.
#' @param ppmGroup the ppm larger of a group
#' @param fracGroup We will keep variables only present in \code{fracGroup} of expirations
#' @param dmzGroup difference of mz of a group for little mz
#' @return A list of data.table with same length as PeakList input, but the alignment between 
#' expirations is done: each matrix contains columns "mz", "delta_mz" and "ppb", "fracExp" and "background"
alignExpirations <- function(PeakMat, ppmGroup = 70, fracGroup=0.8, 
                            dmzGroup = 0.001 ){
  
  if( is.null(PeakMat) ) stop("PeakMat is null")

  #test if there is several group
  if(length(unique(PeakMat$group))==1) {
  PeakMat<-data.table(PeakMat[,c("mz","quanti","delta_mz")])
  PeakMat$fracExp<-rep(1,nrow(PeakMat))
  PeakMat$background<-rep(NA,nrow(PeakMat))
  cat(nrow(PeakMat),"peaks detected \n")
  return(PeakMat)
  }
  
  # we apply at each matrix the function align 
  alignExpList <- align(PeakMat, ppmGroup = ppmGroup)
    
  #Get the number of expiration 
  dim.group.exp <- vapply(alignExpList, function(x) length(which(x[,"group"]!=0)), 
                            FUN.VALUE=0L)
  n.exp <- max(dim.group.exp)

  # aggregate groupe and concatenate the resulted list
  aligExpMat <- do.call(rbind, lapply(alignExpList, function(x) aggregate(x, n.exp)))
 
  #filter on the reproductibility 
  aligExpMat <- aligExpMat[ aligExpMat$fracExp >= fracGroup, ]
 
  cat(paste(nrow(aligExpMat),"peaks detected \n"))
  return(aligExpMat)
}

#' aggregate peakgroup for alignExpirations function  
#' @param subGroupPeak teh group tp aggregate
#' @param n.exp number of expirations done in the file
#' @return a mtrix with the mediane of mz, mean of ppb, ppb in background, 
#'  and percentage of expiration where th epeak is detected
#' @import data.table
aggregate <- function(subGroupPeak, n.exp) {
  
  if(nrow(subGroupPeak) == 0) return(subGroupPeak) #empty data frame
  subGroupPeak <- data.table::data.table(subGroupPeak)
  subGroupPeak<-subGroupPeak[,.(mz=mean(mz),quanti=mean(quanti[group!=0]),
                  delta_mz=mean(delta_mz[group!=0]),fracExp= round( sum(group!=0) /n.exp, 1 ),
                  background=mean(quanti[group==0]))]
  if(is.na(subGroupPeak[,quanti])) return(NULL)
  
  return(subGroupPeak)
}

### Alignsamples method -----

#' Align mass tables of different samples 
#'
#' @param X ptrSet or data.frame contains peakList with at least column mz, quanti and group 
#' @param ppmGroup the ppm larger of a group
#' @param fracGroup We will keep variables only present in \code{fracGroup} 
#' percent of at least one \code{group}, if 0 the filter is not applied
#' @param group character. A sampleMetadata data colnames. If is \code{NULL},variables not presents
#' in \code{fracGroup} percent of samples will be deleted. Else, variables not presents
#' in \code{fracGroup} percent in each group will be deleted. 
#' @param bgThreshold if the background quantity is analyzed, only variables where
#' quantity > bgThreshold x background or quantity < bgThreshold x background for all
#' samples are keeped. 
#' @param dmzGroup difference of mz of a group for little mz
#' @param normalize normalization method 
#' @param ... not used
#' @return an expressionSet (Biobase object), with annotaTion in features Data
#' @examples 
#'library(ptairData)
#'directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' dirSet <- createPtrSet(directory,setName="test",mzCalibRef =c(21.022,59.049))
#' dirSet <- detectPeak(dirSet,mzNominal=c(21,59))
#' getSampleMetadata(dirSet)
#' eset <- alignSamples(dirSet,bgThreshold=0)
#' Biobase::exprs(eset)
#' Biobase::fData(eset)
#' Biobase::pData(eset)
#' @rdname alignSamples
#'@export
setMethod(f="alignSamples",signature = "ptrSet",
            function(X, ppmGroup = 70, fracGroup = 0.8, group=NULL, bgThreshold=2,
                         dmzGroup = 0.001,normalize=FALSE,...
                         ){

              if(normalize){
                if(length(raw@prtReaction)!=0 & nrow(transmission) > 1){
                  U <- c(reaction$TwData[1,,])
                  Td <-c(reaction$TwData[3,,])
                  pd <- c(reaction$TwData[2,,])
                  list_peak$peak[,"quanti"] <- ppbConvert(peakList = list_peak$peak,
                                                          primaryIon = primaryIon,
                                                          transmission = transmission,
                                                          U = U[indBg] , 
                                                          Td = Td[ indBg], 
                                                          pd = pd[ indBg])
                  namesQuanti<-"quanti (ppb)"
                  
                } else {
                  #normalize by primary ions
                  list_peak$peak[,"quanti"] <- list_peak$peak[,"quanti"]/(primaryIon*4.9*10^6)
                  namesQuanti<-"quanti (ncps)"
                }
              }
              
            sampleMetadata <- X@sampleMetadata
            peakList <- X@peakListAligned
            eSet<- alignSamplesFunc(peakList,sampleMetadata,ppmGroup, fracGroup , group,
                          dmzGroup,bgThreshold)
  
  return(eSet)
})

#' @rdname alignSamples
#' @param sampleMetadata if X is a ptrSet, not used, else a data.frame containing sample Metadata
#'@export
setMethod(f="alignSamples",signature = "data.frame",
          function(X, ppmGroup = 70, fracGroup = 0.8, group=NULL,bgThreshold=2,
                   dmzGroup = 0.001,sampleMetadata
          ){
            
            eSet<- alignSamplesFunc(X,sampleMetadata,ppmGroup, fracGroup , group,
                                    dmzGroup,bgThreshold)
            
            return(eSet)
          })


alignSamplesFunc <- function(peakList,sampleMetadata, ppmGroup=100, fracGroup=0.8, group=NULL, 
                 dmzGroup=0.001,bgThreshold=2){
  
  ## add column group with Samples group number
  mat<-NULL
  for(i in 1:length(peakList)){
    mat_temp <- cbind(peakList[[i]], group = i)
    mat <- rbind(mat,mat_temp)
  }
  
  
  
  # make subgroup of peak
  groupList <- align(as.matrix(mat), ppmGroup, dmzGroup)
  
  # Get number of samples
  nSample<-length(peakList)
  
  # aggregate
  groupMat <- do.call(rbind, lapply( groupList, 
                                     function(y){
                                       y<-data.table::as.data.table(y)
                                       y <- y[,.(mz= median(mz), 
                                                 quanti= paste(quanti,collapse = "_"),
                                                 background = if("background" %in% colnames(y)) paste(background,collapse = "_"),
                                                 Samples=paste(group,collapse="_"),
                                                 nSamples= round(length(group)/nSample,1))]
                                       y
                                     }
  ) 
  )
  
  # fomratting the final matrix 
  mat.final.Exp<-apply(groupMat[,c("quanti","Samples")], MARGIN=1, function(x){
    output<-rep(NA,nSample)
    ch.Area <- x[1]
    ch.Area.split <- strsplit(ch.Area, split = "_")[[1]]
    vec.Area <- as.numeric(ch.Area.split)
    
    ch.samples <- x[2]
    vec.samples <- as.numeric(strsplit(ch.samples, split = "_")[[1]])
    
    output[vec.samples]<-vec.Area
    output
  })
  
  X <- t(matrix(mat.final.Exp, nrow = nSample)) # matrix function for the case of one sample and mat.final is an array
  rownames(X) <- round(groupMat[,mz],4)
  colnames(X) <- names(peakList)
  
  #filter on samples frenquency
  if(!is.null(group)){
    groupFac <- as.factor(sampleMetadata[,group])
    test <- apply(X,1,
                  function(x) sapply(levels(groupFac),
                                     function(y) 
                                       sum(!is.na(x[groupFac==y]))/length(x[groupFac==y]) >= fracGroup ))
    keepVar<-apply(test,2,any)
    
  }else {
    keepVar<- apply(X, 1, function(y) sum(!is.na(y))/length(y) >= fracGroup)
    
  }
  
  X<- X[keepVar,,drop=FALSE]
  
  if( nrow(X) ==0 ) { 
    warning(paste("peakList is empty, there is no peak presents in more thant ", fracGroup, "% of samples"))
    return(NULL)
  }
  
  if("background" %in% colnames(groupMat) & bgThreshold!=0 ){
    mat.final.bg<-apply(groupMat[,c("background","Samples")], MARGIN=1, function(x){
      output<-rep(NA,nSample)
      ch.Area <- x[1]
      ch.Area.split <- strsplit(ch.Area, split = "_")[[1]]
      vec.Area <- as.numeric(ch.Area.split)
      
      ch.samples <- x[2]
      vec.samples <- as.numeric(strsplit(ch.samples, split = "_")[[1]])
      
      output[vec.samples]<-vec.Area
      output
    })
    Xbg <- t(matrix(mat.final.bg, nrow = nSample))
    Xbg <- Xbg[keepVar,,drop=FALSE]
    rownames(Xbg) <- rownames(X)
    colnames(Xbg) <- paste("Background -",names(peakList))

    # backgroud discarding
    
    ## wilcoxon test
    #test<-rep(0,nrow(X))
    #for(i in 1:nrow(X)){
    #  test[i]<-wilcox.test(x=X[i,],y=Xbg[i,],
      #       alternative="two.sided",paired=T)$p.value

    ## fixed threshold
    if(bgThreshold!=0){
      rap<-Xbg/X
      delete <- which(apply(rap,1,function(x) all( 1/bgThreshold <=x & x<= bgThreshold , na.rm=T) ))
      X <- X[-delete,,drop=F]
      Xbg<-Xbg[-delete,,drop=F]
      if( nrow(X)==0 ) {
        warning(paste("peakList is empty, there is no peak presents in more thant ", bgThreshold, "* background"))
        return(NULL)
      }
    }
      
      
     
      
  } else Xbg<-data.frame(row.names = rownames(X) )
 

  
  
  cat(paste(nrow(X),"peaks aligned \n"))
 
  featureData <- Biobase::AnnotatedDataFrame(as.data.frame(Xbg))
  sampleMetadata<-as.data.frame(sampleMetadata[colnames(X),,drop=F])
  return(Biobase::ExpressionSet(assayData=X,
                                phenoData = Biobase::AnnotatedDataFrame(sampleMetadata),
                                featureData = featureData)
  )
}





#'Impute missing values on an expression set from an ptrSet
#'
#'Imputing missing values by returning back to the raw data and fitting the 
#'peak shape function on the noise / residuals
#'@param eSet an expression set retun by alignSamples function 
#'@param ptrSet processed by detectPeak function
#'@return the same expression set as in input, with missing values imputing
#'@export 
#'@examples
#'library(ptairData)
#'directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' dirSet <- createPtrSet(directory,setName="test",mzCalibRef =c(21.022,59.049))
#' dirSet <- detectPeak(dirSet,mzNominal=c(21,63))
#' getSampleMetadata(dirSet)
#' eSet <- alignSamples(dirSet,fracGroup=0)
#' Biobase::exprs(eSet)
#' eSet <- impute(eSet,dirSet)
#' Biobase::exprs(eSet)
impute <- function(eSet,ptrSet){
  
  #get peak list 
  peakList <- getPeakList(ptrSet)$aligned
  
  #get index of missing values
  missingValues <-which(is.na(Biobase::exprs(eSet)),arr.ind=TRUE)
  indexFilesMissingValues <- unique(missingValues[,"col"])
  namesFilesMissingValues <- colnames(Biobase::exprs(eSet))[indexFilesMissingValues]
  
  #get primry ion quantity
  primaryIon<-ptrSet@primaryIon
  
  #get files full names in ptrSet object
  filesFullName<-ptrSet@parameter$listFile
  for (file in namesFilesMissingValues){
    j<-which(file==colnames(Biobase::exprs(eSet)))
    filesFullName.j<-filesFullName[which(basename(filesFullName)==file)]
    
    mzMissing <- as.numeric(rownames(missingValues[missingValues[,"col"]==j,,drop=F]))
    
    #opense mz Axis 
    mzAxis <- rhdf5::h5read(filesFullName.j,name="FullSpectra/MassAxis")
    indexMzList <- lapply(unique(round(mzMissing)),function(m) which( m - 0.6 < mzAxis & 
                                                                        mzAxis < m+ 0.6))
    names(indexMzList)<-unique(round(mzMissing))
    indexMz<-Reduce(union,indexMzList)
    
    #get index time
    indexLim <- ptrSet@timeLimit[[file]]
    indexTime<-Reduce(c,apply(indexLim,2,function(x) seq(x["start"],x["end"])))
    nbExp<-ncol(indexLim)
    
    #open raw data
    raw <- rhdf5::h5read(filesFullName.j, name = "/FullSpectra/TofData",
                         index=list(indexMz,NULL,NULL,NULL))
    
    rawMn <- matrix(raw,
                    nrow = dim(raw)[1],
                    ncol = prod(tail(dim(raw),2))) #* 0.2 ns / 2.9 (single ion signal) if convert to cps
    
    # information for ppb convertion
    reaction <-  try(reaction<- rhdf5::h5read(filesFullName.j,"AddTraces/PTR-Reaction"))
    transmission <-try(rhdf5::h5read(filesFullName.j,"PTR-Transmission"))
    
    #calibrate mass axis
    FirstcalibCoef <- rhdf5::h5read(filesFullName.j,"FullSpectra/MassCalibration",index=list(NULL,1))
    tof <- sqrt(mzAxis)*FirstcalibCoef[1,1] + FirstcalibCoef[2,1]
    coefCalib<-ptrSet@coefCalib[[file]]
    mzAxis <- ((tof-coefCalib['b',])/coefCalib['a',])^2
    
    #peak list raw
    peakListRaw.j<-ptrSet@peakListRaw[[file]]
    
    for (m in unique(round(mzMissing))){
      #exact missing mz
      mz <- mzMissing[round(mzMissing)==m]
      
      #mzAxis around m
      mzAxis.m <- mzAxis[indexMzList[[as.character(m)]]]
      
      indexExp<-Reduce(c,apply(indexLim,2,function(x) seq(x[1],x[2])))
      length.exp<-length(indexExp)
      quantiMat<-matrix(0,ncol=length(mz),nrow=nbExp)
        spectrum <- rowSums(rawMn[which(indexMz %in% indexMzList[[as.character(m)]]),
                                  indexExp]) /length.exp
        
        #substract fitted peak also find
        peakAlsoDetected <- peakListRaw.j[round(peakListRaw.j$mz)==m & peakListRaw.j$group==1,]
        if(nrow(peakAlsoDetected)!=0){ 
          # cumFitPeak
          
          fitPeaks <- apply(peakAlsoDetected,1,
                            function(x) sech2(x["mz"],x["parameter.1"],
                                              x["parameter.2"],x["parameter.3"],x = mzAxis.m)
          )
          if(nrow(peakAlsoDetected)>1) cumFitPeak <- rowSums(fitPeaks) else cumFitPeak <- c(fitPeaks)
          spectrum<- spectrum-cumFitPeak
        }
        
        #fit on the missing values
        resolution_upper<-8000
        resolution_mean<- 5000
        resolution_lower<-3500
        
        n.peak<-length(mz)
        delta<-rep(m/resolution_mean,2*n.peak)
        h<- vapply(mz, function(m) max(max(spectrum[which(abs(mzAxis.m-m)<(m*50/10^6))]),1),0)
        
        initMz <- matrix(c(mz,log(sqrt(2)+1)*2/delta,
                           h),nrow=n.peak)
        colnames(initMz)<-c("m","delta1","delta2","h")
        
        lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak),
                                            rep(0, n.peak*2),
                                            rep(0.1, n.peak)),ncol = 4) 
                          -
                            matrix(c(initMz[,"m"]/(resolution_mean*100),
                                     -resolution_lower*log(sqrt(2)+1)*2/initMz[,"m"],
                                     -resolution_lower*log(sqrt(2)+1)*2/initMz[,"m"],
                                     rep(0, n.peak)),ncol = 4)))
        
        upper.cons <- c(t( initMz * matrix(c(rep(1, n.peak),
                                             rep(0, n.peak*2),
                                             rep(Inf, n.peak)),ncol = 4) 
                           +
                             matrix(c(initMz[,"m"]/(resolution_mean*100),
                                      resolution_upper*log(sqrt(2)+1)*2/initMz[,"m"],
                                      resolution_upper*log(sqrt(2)+1)*2/initMz[,"m"],
                                      rep(0, n.peak)),ncol = 4)))
        
        
        fit <- fit_sech2(initMz, spectrum, mzAxis.m, lower.cons, upper.cons)
        
        fit.peak <- fit$fit.peak
        par_estimated<-fit$par_estimated
        
        quanti.m <- apply(par_estimated,2, function(x){
          th<-10*0.5*(log(sqrt(2)+1)/x[2]+log(sqrt(2)+1)/x[3])
          mz.x <- mzAxis.m[ x[1] - th < mz & mz < x[1]+th ]
          sum(sech2(x[1],x[2],x[3],x[4],mz.x),na.rm =T)}) 
        
        list_peak<-cbind(mz,quanti=quanti.m)
        
        # convert to ppb or ncps
        #if there is reaction ans transmission information
        if(ptrSet@parameter$detectPeakParam$normalize){
          if(is.null(attr(transmission,'condition')) & is.null(attr(reaction,'condition'))){
            U <- c(reaction$TwData[1,,])
            Td <-c(reaction$TwData[3,,])
            pd <- c(reaction$TwData[2,,])
            quanti.m <- ppbConvert(peakList = list_peak,
                                   primaryIon = primaryIon[[file]],
                                   transmission = transmission$Data,
                                   U = U[ indexExp] , 
                                   Td = Td[ indexExp], 
                                   pd = pd[indexExp])
            
          } else {
            #normalize by primary ions
            quanti.m <- quanti.m/(primaryIon[[basename(fullNamefile)]]*4.9*10^6)
          }
        }


      # add to peak table
      Biobase::exprs(eSet)[as.character(mz),file] <- quanti.m
    }
    message(basename(file)," done")
  }
  return(eSet)
}


# imputeOld <- function(eSet,ptrSet){
#   
#   #get peak list 
#   peakList <- getPeakList(ptrSet)$aligned
#   
#   #get index of missing values
#   missingValues <-which(is.na(Biobase::exprs(eSet)),arr.ind=TRUE)
#   indexFilesMissingValues <- unique(missingValues[,"col"])
#   namesFilesMissingValues <- colnames(Biobase::exprs(eSet))[indexFilesMissingValues]
#   
#   #get primry ion quantity
#   primaryIon<-ptrSet@primaryIon
#   
#   #get files full names in ptrSet object
#   filesFullName<-ptrSet@parameter$listFile
#   for (file in namesFilesMissingValues){
#     j<-which(file==colnames(Biobase::exprs(eSet)))
#     filesFullName.j<-filesFullName[which(basename(filesFullName)==file)]
#     
#     mzMissing <- as.numeric(rownames(missingValues[missingValues[,"col"]==j,,drop=F]))
#    
#     #opense mz Axis 
#     mzAxis <- rhdf5::h5read(filesFullName.j,name="FullSpectra/MassAxis")
#     indexMzList <- lapply(unique(round(mzMissing)),function(m) which( m - 0.6 < mzAxis & 
#                                                                 mzAxis < m+ 0.6))
#     names(indexMzList)<-unique(round(mzMissing))
#     indexMz<-Reduce(union,indexMzList)
#     
#     #get index time
#     indexLim <- ptrSet@timeLimit[[file]]
#     indexTime<-Reduce(c,apply(indexLim,2,function(x) seq(x["start"],x["end"])))
#     nbExp<-ncol(indexLim)
#     
#     #open raw data
#     raw <- rhdf5::h5read(filesFullName.j, name = "/FullSpectra/TofData",
#                          index=list(indexMz,NULL,NULL,NULL))
#     
#     rawMn <- matrix(raw,
#                     nrow = dim(raw)[1],
#                     ncol = prod(tail(dim(raw),2))) #* 0.2 ns / 2.9 (single ion signal) if convert to cps
#     
#     # information for ppb convertion
#     reaction <-  try(reaction<- rhdf5::h5read(filesFullName.j,"AddTraces/PTR-Reaction"))
#     transmission <-try(rhdf5::h5read(filesFullName.j,"PTR-Transmission"))
#     
#     #calibrate mass axis
#     FirstcalibCoef <- rhdf5::h5read(filesFullName.j,"FullSpectra/MassCalibration",index=list(NULL,1))
#     tof <- sqrt(mzAxis)*FirstcalibCoef[1,1] + FirstcalibCoef[2,1]
#     coefCalib<-ptrSet@coefCalib[[file]]
#     mzAxis <- ((tof-coefCalib['b',])/coefCalib['a',])^2
#     
#     #peak list raw
#     peakListRaw.j<-ptrSet@peakListRaw[[file]]
#     
#     for (m in unique(round(mzMissing))){
#       #exact missing mz
#       mz <- mzMissing[round(mzMissing)==m]
#       
#       #mzAxis around m
#       mzAxis.m <- mzAxis[indexMzList[[as.character(m)]]]
#       
#       #loop over expirations
#       quantiMat<-matrix(0,ncol=length(mz),nrow=nbExp)
#       for (n.exp in 1:nbExp ){
#         length.exp<-indexLim["end",n.exp] - indexLim["start",n.exp] + 1
#         indexExp<-seq( indexLim["start",n.exp],indexLim["end",n.exp] )
#         spectrum <- rowSums(rawMn[which(indexMz %in% indexMzList[[as.character(m)]]),
#                                   indexExp]) /length.exp
#         
#         #substract fitted peak also find
#         peakAlsoDetected <- peakListRaw.j[round(peakListRaw.j$mz)==m & peakListRaw.j$group==n.exp,]
#         if(nrow(peakAlsoDetected)!=0){ 
#           # cumFitPeak
# 
#           fitPeaks <- apply(peakAlsoDetected,1,
#                           function(x) sech2(x["mz"],x["parameter.1"],
#                                                       x["parameter.2"],x["parameter.3"],x = mzAxis.m)
#                             )
#           if(nrow(peakAlsoDetected)>1) cumFitPeak <- rowSums(fitPeaks) else cumFitPeak <- c(fitPeaks)
#           spectrum<- spectrum-cumFitPeak
#         }
#        
#         #fit on the missing values
#         resolution_upper<-8000
#         resolution_mean<- 5000
#         resolution_lower<-3500
#         
#         n.peak<-length(mz)
#         delta<-rep(m/resolution_mean,2*n.peak)
#         h<- vapply(mz, function(m) max(max(spectrum[which(abs(mzAxis.m-m)<(m*50/10^6))]),1),0)
#         
#         initMz <- matrix(c(mz,log(sqrt(2)+1)*2/delta,
#                            h),nrow=n.peak)
#         colnames(initMz)<-c("m","delta1","delta2","h")
#         
#         lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak),
#                                             rep(0, n.peak*2),
#                                             rep(0.1, n.peak)),ncol = 4) 
#                           -
#                             matrix(c(initMz[,"m"]/(resolution_mean*100),
#                                      -resolution_lower*log(sqrt(2)+1)*2/initMz[,"m"],
#                                      -resolution_lower*log(sqrt(2)+1)*2/initMz[,"m"],
#                                      rep(0, n.peak)),ncol = 4)))
#         
#         upper.cons <- c(t( initMz * matrix(c(rep(1, n.peak),
#                                              rep(0, n.peak*2),
#                                              rep(Inf, n.peak)),ncol = 4) 
#                            +
#                              matrix(c(initMz[,"m"]/(resolution_mean*100),
#                                       resolution_upper*log(sqrt(2)+1)*2/initMz[,"m"],
#                                       resolution_upper*log(sqrt(2)+1)*2/initMz[,"m"],
#                                       rep(0, n.peak)),ncol = 4)))
#         
#         
#         fit <- fit_sech2(initMz, spectrum, mzAxis.m, lower.cons, upper.cons)
#         
#         fit.peak <- fit$fit.peak
#         par_estimated<-fit$par_estimated
#         
#         quanti.m <- apply(par_estimated,2, function(x){
#           th<-10*0.5*(log(sqrt(2)+1)/x[2]+log(sqrt(2)+1)/x[3])
#           mz.x <- mzAxis.m[ x[1] - th < mz & mz < x[1]+th ]
#           sum(sech2(x[1],x[2],x[3],x[4],mz.x),na.rm =T)}) 
#         
#         list_peak<-cbind(mz,quanti.m)
# 
#         # convert to ppb or ncps
#         #if there is reaction ans transmission information
#         if(ptrSet@parameter$detectPeakParam$normalize){
#         if(is.null(attr(transmission,'condition')) & is.null(attr(reaction,'condition'))){
#           U <- c(reaction$TwData[1,,])
#           Td <-c(reaction$TwData[3,,])
#           pd <- c(reaction$TwData[2,,])
#           quanti.m <- ppbConvert(peakList = list_peak,
#                                                             primaryIon = primaryIon[[file]],
#                                                             transmission = transmission$Data,
#                                                             U = U[ indexExp] , 
#                                                             Td = Td[ indexExp], 
#                                                             pd = pd[indexExp])
#           
#         } else {
#           #normalize by primary ions
#           quanti.m <- quanti.m/(primaryIon[[basename(fullNamefile)]]*4.9*10^6)
#         }
#         }
#         quantiMat[n.exp,]<-quanti.m
#       }
#       
#           # aggregation 
#       quantiMissing<-apply(quantiMat,2,mean)
#       
#          
#           # add to peak table
#       Biobase::exprs(eSet)[as.character(mz),file] <- quantiMissing
#     }
#     message(basename(file)," done")
#   }
#   return(eSet)
# }

