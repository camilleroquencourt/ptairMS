utils::globalVariables("::<-")

## detectPeak----

#' Detection and quantification of peaks for ptrSet object. 
#' 
#' \code{detectPeak} function detect peak in all file present in ptrSet who have not alredy been processed, in the signal 
#' and the background. To see the peakList use getPeakList() method. It give a list of data.frame, whith column:
#' peak center in mz, peak quanitfication in ppb, peak width and background intenisty. 
#' which each element correspond to the Peak List of one file.
#' For each file the following steps are done: \cr
#' for each expiration and ambiant air :
#' \itemize{
#' \item calibration 
#' \item peak detection and quantification 
#' \item alignment between expirations and ambiant air
#' }
#' \code{calibration} \cr
#' \code{timeLimits} \cr
#' for each time periods and background : \code{PeakList} \cr
#' \code{alignExpirations}  
#' @param x ptrSet object 
#' @param mzNominal nominal mass who peak will be detected, by default all 
#' nominal mass present in the mass axis
#' @param ppm the minimum distance in ppm betewwen two peaks
#' @param ppmGroupBkg the ppm width of a group made by the alignment between signal and background peaks
#' @param fracGroup numeric, between 0 and 1. defining the minimum fraction of expiration (if there are) 
#' which the peak have to be present to be considered as a peak group (feature)
#' @param minIntensity the minimum intenisty for peak detection. The threshold for peak detection
#' will be : max ( \code{minPeakDetect} , threshold noise ). The theshold noise correspond to
#'  max(\code{thNoiseRate} * max( noise around the nominal mass), \code{thIntensityRate} * 
#'  max( intenisty in the nominal mass)
#' @param fctFit the function for the quantification of Peak, should be sech2 or Average
#' @param parallelize boolean. If \code{TRUE} loop aver files will be parallelized
#' @param nbCores number of cluster to use for parrallel computation.
#' @param normalize boolean. if TRUE result in ppb or ncps, normalized by primary ions
#' @param fracMaxTIC if x is a ptrRaw, the same paramter as \code{createPtrSet} function:Percentage 
#' (between 0 and 1) of the maximum of the Total Ion Chromatogram (TIC) amplitude with baseline removal. 
#' We will analyze only the part of the spectrum where the TIC intensity is higher than 'fracMaxTIC * max(TIC) '. 
#' If you want to analyze the entire spectrum, set this parameter to 0. 
#' @param saving boolean. If TRUE, the object will be saved in saveDir with the
#' \code{setName} paramter of \code{createPtrSet} function
#' @param saveDir The directory where the ptrSet object will be saved in .RData. If NULL, nothing will be saved
#' @param processFun the function used to process a file. Can not be changed for the moent 
#' @param ... not used
#' @return a S4 object ptrSet, that contains ptrset and a list of data.frame (data.table), which contains : \cr
#' mass of peak center, peak quantifictaion in ppb if it is possible, peak width, percentage of expiration where the peak 
#' have been detected, and background quantification. If the peak has not been deteted in 
#' the background, last column is NA.
#' @examples 
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' dirSet <- createPtrSet(directory,setName="test")
#' dirSet <- detectPeak(dirSet , mzNominal=59)
#' getPeakList(dirSet)$aligned
#' @rdname detectPeak
#' @import doParallel foreach parallel
#' @export
setMethod(f="detectPeak",
            signature = "ptrSet",
          function(x, 
                   mzNominal=NULL , ppm=130, ppmGroupBkg=50, fracGroup=0.8, minIntensity=10, 
                   fctFit=c("Sech2","average")[1],parallelize=FALSE,nbCores=2,normalize=TRUE,
                   saving=TRUE,saveDir=x@parameter$saveDir,processFun=processFileSepExp,...)
          {
            ptrset<-x
            
            #get infomration
            massCalib<-ptrset@mzCalibRef
            primaryIon<-ptrset@primaryIon
            indTimeLim<-ptrset@timeLimit
            parameter <- ptrset@parameter
            dir <- parameter$dir
            peakListRaw <-ptrset@peakListRaw
            paramOld <- parameter$detectPeakParam
            if(is.null(mzNominal)) mzNominalParam <- "NULL" else mzNominalParam<- mzNominal
            
            #files already check by chexkSet 
            files <- parameter$listFile 
            fileName <- basename(files)
           
            paramNew <- list(mzNominal=mzNominalParam, ppmGroupBkg=ppmGroupBkg,fracGroup=fracGroup,
                          ppm=ppm, minIntensity=minIntensity, fctFit=fctFit,normalize=normalize)

           # keep files that not alreday processed
            fileDone <- names(peakListRaw)
            fileToProcess <- which(!(fileName %in% fileDone))
            fileName <- fileName[ fileToProcess ]
            allFilesName<-basename(files) #to save the oder
            
            files <- files[ fileToProcess ]
            
            if(length(files)==0) {
              message("All files have already been processed")
              return(ptrset)
            }
            
            if( is.null(paramOld)){
              ptrset@parameter$detectPeakParam <- paramNew
            }  else {
              mzNominal<-paramOld$mzNominal
              if(mzNominal[1]=="NULL") {
                mzNominal<-NULL} else paramOld$mzNominal<-paste(range(paramOld$mzNominal),collapse="-")
              ppmGroupBkg<-paramOld$ppmGroupBkg
              ppm<-paramOld$ppm
              minIntensity<-paramOld$minIntensity
              fctFit<-paramOld$fctFit
              paramOldMat<-Reduce(cbind,paramOld)
              colnames(paramOldMat)<-names(paramOld)
              message("the peak list will be calculated with the same parameters as the other files :\n" )
              #print(paramOldMat)
              message("\n If you want to change theme, remove the peak list before with rmPeakList() function and restart detectPeak() \n")
            }
            
            
            FUN <- function(x) { 
              test<-try(processFun(x,massCalib[[basename(x)]],
                                                  primaryIon[[basename(x)]],
                                                  indTimeLim[[basename(x)]], mzNominal,
                                                      ppm, ppmGroupBkg, fracGroup,minIntensity, 
                                                   fctFit,normalize))
              if(!is.null(attr(test,'condition'))){
                return(list(raw=0,aligned=0))
              } else return(test)
              
              }

            #parallel 1
            #peakLists<-BiocParallel::bplapply(files,FUN = FUN)
            
            #parrallel 2 
            if(parallelize){
              cl <- parallel::makeCluster(nbCores)
              doParallel::registerDoParallel(cl)
              `%dopar%`<-foreach::`%dopar%`
              peakLists <- foreach::foreach(file=files, .packages = c("data.table")) %dopar% {
                FUN(file)}
              parallel::stopCluster(cl)
              
            } else  peakLists<-lapply(files,FUN)
            
            #separate raw and aligned peak list
            peakListRaw<-lapply(peakLists,function(x) x$raw)
            names(peakListRaw)<-basename(files)
            peakListAligned<-lapply(peakLists,function(x) x$aligned)
            names(peakListAligned)<-basename(files)
            
            # add list and order as in list.file
            ptrset@peakListRaw<-c(ptrset@peakListRaw,peakListRaw)[allFilesName]
            ptrset@peakListAligned<-c(ptrset@peakListAligned,peakListAligned)[allFilesName ]
            
            #delete processFile failed
            failed<-which(Reduce(c,lapply(ptrset@peakListRaw, function(x) is.null(dim(x)))))
            if( length(failed) ) {
              ptrset@peakListAligned<- ptrset@peakListAligned[-failed]
              ptrset@peakListRaw<- ptrset@peakListRaw[-failed]
              warning(allFilesName[failed],"failed")
            }
            
            # save
            if(!is.null(saveDir) & saving ){
              if(!dir.exists(saveDir)){
                warning("saveDir does not exist, object not saved")
                return(ptrset)
              } 
              changeName <- parse(text=paste0(ptrset@parameter$name,"<- ptrset "))
              eval(changeName)
              eval(parse(text =  paste0( "save(" ,ptrset@parameter$name ,",file = paste0(saveDir,'/', '",ptrset@parameter$name,".RData '))")))
            }
            
            return(ptrset)
          } )

#'extract quantity over time of VOC in a pre-defined peak list
#'@param raw ptrRaw object 
#'@param peak a data.frame with a column named "Mz". The Mz of the VOC detected
#'@param timeCalib duration in second of every calibration
#'@param deconvMethod the deconvolution 2d method. Should be "deconv2d2linearIndependant",
#'"deconv2d2LinearCoupled" or "deconv2d2NonlinearIndependant"
#'@return list containing the matrix with all temporal profile (matPeak) 
#'of VOC and a list of all raw EIC (EIClist)
computeTemporalFile<-function(raw,peak,timeCalib=20,deconvMethod=deconv2d2linearIndependant){
  
  listCalib<-multiCalib(raw,timeCalib)
  
  EICex<-extractEIC(raw,peak)
  
  borne<-EICex$borne
  EIC<-EICex$EIC
  
  nbPeak <- sum(borne[,"overlap"])
  
  # compute EIC
  
  # without overlap
  borne<-cbind(borne,matrix(0,nrow = nrow(borne),ncol=length(raw@time))) #ajoter colonne de time
  colnames(borne)[5:ncol(borne)]<-raw@time
  if(nrow(borne)>1) borneUnique<-borne[!duplicated(data.frame(borne[,2:3])),,drop=FALSE] else
    borneUnique <- borne
  borne[borne[,"overlap"]==0,5:ncol(borne)] <- Reduce(rbind,lapply(EIC[borneUnique[,"overlap"]==0], colSums))
  
  #with overlapp
  XICdeconv2<-matrix(0,ncol=length(raw@time),nrow=nbPeak)
  c<-1
  for(j in which(borneUnique[,"overlap"]==1)){
    mz<-borne[which(borne[,"lower"]==borneUnique[j,"lower"]),"Mz"]
    peak.detect<-peak[peak[,"Mz"]%in% mz,,drop=F]
    
    rawM<-EIC[[j]]
    deconv2<-deconvMethod(rawM,raw@time,peak.detect,raw,listCalib)
    
    for (i in seq_len(nrow(peak.detect))){
      XICdeconv2[c,]<-deconv2$predPeak[,i]
      c<-c+1
    }
  }
  XICdeconv<-XICdeconv2
  borne[borne[,"overlap"]==1,5:ncol(borne)]<-XICdeconv
 
  
  # test significativité
  # comparaison de moyenne normal 
  indLim<-timeLimits(raw,fracMaxTIC = 0.7)
  indExp<-Reduce(c,apply(indLim,2,function(x) seq(from=x[1],to=x[2])))
  indLimBg<- bakgroundDetect(colSums(raw@rawM))
  indBg<-Reduce(c,apply(indLimBg,2,function(x) seq(from=x[1],to=x[2])))
  
  XIC<-borne[,5:ncol(borne),drop=F]
  p.greater<-apply(XIC,1,function(x) t.test(x[indExp],x[indBg],alternative="greater")$p.value)
  p.less<-apply(XIC,1,function(x) t.test(x[indExp],x[indBg],alternative="less")$p.value)
  p.greater<-p.adjust(p.greater,method = "fdr")
  p.less<-p.adjust(p.less,method = "fdr")
  
  borne<-cbind(borne, data.frame(p.less=p.less,p.upper=p.greater))
  
  return(list(matPeak=borne,EIClist=EIC))
}

#'extract all raw EIC from a pre-definied peak List
#'@param raw ptrRaw object
#'@param peak a data.frame with a column named "Mz". The Mz of the VOC detected
#'@param peakQuantil the quantile of the peak shape to determine the borne of the EIC
#'@return list containing all EIC and the mz borne for all peak
extractEIC<-function(raw,peak,peakQuantil=0.05){
  #borne integration
  borne<-apply(peak,1,function(x) sechInv(x["Mz"],
                                          x["parameter.1"],
                                          x["parameter.2"],
                                          x["parameter.3"],
                                          x["parameter.3"]*0.05) )
  borne<-cbind(peak$Mz,t(borne))
  colnames(borne)<- c("Mz","lower","upper")
  
  #overlap detection and fusion
  overlap<-list(NULL)
  c<-0
  
  o<-1
  for (i in seq_len(nrow(borne)-1)){
    if(borne[i+1,"lower"]< borne[i,"upper"]){
      #save the overlap
      o<-c(o,i+1)
    } else {
      if(length(o)>1){
        #change
        borne[o,"lower"]<-borne[o[1],"lower"]
        borne[o,"upper"]<-borne[tail(o,1),"upper"]
        c<-c+1
        overlap[[c]]<-o
      }
      o <-i+1 
    }
  }
  if(length(o)>1){
    #change
    borne[o,"lower"]<-borne[o[1],"lower"]
    borne[o,"upper"]<-borne[tail(o,1),"upper"]
    c<-c+1
    overlap[[c]]<-o
  }
  
  borne<-cbind(borne,overlap=0)
  borne[Reduce(c,overlap),"overlap"]<-1
  if(nrow(borne)>1) borneUnique<-borne[!duplicated(data.frame(borne[,2:3])),,drop=FALSE] else
    borneUnique <- borne
  
 
  #extract EIC
  EIC<-list(NULL)
  for (j in seq_len(nrow(borneUnique))){
    EIC[[j]]<-raw@rawM[raw@mz>borneUnique[j,"lower"] & raw@mz<borneUnique[j,"upper"],]
  }
  
  return(list(EIC=EIC,borne=borne))
}

#'Performs a calibration every x time step
#'@param raw a ptrRaw object already calibrate
#'@param step number of time step  
#'@return a list that each element contains: index time, coef calibration, and the function to convert mz to tof
multiCalib<-function(raw,step){
  t<-raw@time
  mzCalibRef<-raw@calibMassRef
  calib_invformulaLIst<-list(NULL)
  mzRef<-raw@mz
  
  indexTime<-round(diff(t)[1],3)
  nbIndex<- round(step/indexTime)
  for (i in seq(0,(floor(length(t)/nbIndex)-1))) {
    index <- seq(from=i*nbIndex+1,
                 to=min((i*nbIndex+nbIndex),length(raw@time)))
    if(i == (floor(length(t)/nbIndex)-1) & utils::tail(index,1) < length(t)) index<-c(index,seq(utils::tail(index,1)+1,length(raw@time)))
    sp.i <- rowSums(raw@rawM[,index])
    
    width.window<-0.4
    
    # calculate tof axis
    if(is.null(mzToTof(1,raw@calibCoef))) tof<-seq(0,length(mzRef)-1) else tof <- mzToTof(mzRef,raw@calibCoef)
    
    #spectrum of mass calib 
    spTronc <- lapply(mzCalibRef,function(m) {
      index<- which((mzRef < m + width.window) & (mzRef > m - width.window))
      tof<-tof[index]
      mzRef <- mzRef[index]
      signal <- sp.i[index] 
      return(list(signal=signal,mz=mzRef,tof=tof))})
    
    # calculate the tof of reference masses
    tofMax <- vapply(spTronc, function(x) {
      tofrange <- x$tof
      sp<- x$signal
      t<-tofrange[which.max(sp)]
      delta<-10000 * log(sqrt(2)+1)*2/ t
      init<-list(m=t,d1=delta,d2=delta,h=max(sp))
      fit<-suppressWarnings(minpack.lm::nls.lm(par=init, 
                                               fn =function(par,x,y) y- sech2(
                                                 par$m,par$d1,par$d2,par$h,x),
                                               x= tofrange , y = sp))
      tofMax<-fit$par$m
      
      return(tofMax)
    },FUN.VALUE = 0.1)
    
    # re estimated calibration coefficient with reference masses
    regression <- stats::nls( rep(1,length(mzCalibRef))  ~  I( ( (tofMax - b) / a) ^ 2 /mzCalibRef ),
                              start = list(a = raw@calibCoef["a",], b= raw@calibCoef["b",] ), algorithm = "port")
    
    coefs <- stats::coefficients(regression)
    calib_invformula <- function(m,coefs) sqrt(m)*coefs['a'] + coefs['b']
    calib_invformulaLIst[[i+1]]<-list(index=index,coef=coefs,funcMzToTof=calib_invformula)
    
  }
  
  return(calib_invformulaLIst)
  
}

#'Estimate the calibration for precise mz
#'@param raw a ptrRaw object
#'@param peak a data.frame with a column named "Mz". The Mz of the VOC detected
#'@param listCalib a list containing calibration coeficient for different time periods:
#'the output of multicalib function
#'@return a vector with mz shift values in Da for eaf mz of peak list
calibShitEstimate<-function(raw,peak,listCalib){
  
  shift<-matrix(0,ncol=length(listCalib),nrow=nrow(peak))
  for(k in seq_along(listCalib)){
    index<-listCalib[[k]]$index
    mz<- tofToMz(listCalib[[k]]$funcMzToTof(peak$Mz,listCalib[[k]]$coef),raw@calibCoef)
    shift[,k]<-mz-peak$Mz
  }
  
  return(shift)
}

#'corretcion of the tof shift
#'@param rawM an EIC
#'@param raw ptrRaw object (replaec by the coef calib) 
#'@param peak mz peak in rawM
#'@param listCalib mutlicalib list
#'@return the rawM corrected 
calibShiftCorr<-function(rawM,raw,peak,listCalib){
  calibShit<-calibShitEstimate(raw,peak,listCalib)
  nbIndex<-length(listCalib[[1]]$index)
  size<-utils::tail(listCalib[[length(listCalib)]]$index,1)
  time<-seq(1, size, by=nbIndex)[1:ncol(calibShit)]+nbIndex/2
  
  # image(x=raw@time,y=as.numeric(rownames(rawM)),z=t(rawM))
  # lines(time,calibShit[1,]+peak$Mz[1],pch=19)
  # lines(time,calibShit[2,]+peak$Mz[2],pch=19)
  
  rawMCorr<-matrix(0,ncol=ncol(rawM),nrow=nrow(rawM))
  rawMCorr[,1]<-rawM[,1]
  nr<-nrow(rawM)
  for (i in seq_along(listCalib)){
    index<-listCalib[[i]]$index
    rawMCorr[,index]<-vapply(index,function(j){
      approx(x=as.numeric(rownames(rawM)) ,
             y=rawM[,j],
             xout=as.numeric(rownames(rawM))+ calibShit[1,i])$y
    },FUN.VALUE = rep(0,nrow(rawM)))
    
  }
  
  rownames(rawMCorr)<-rownames(rawM)
  rawMCorr<-rawMCorr[! apply(rawMCorr,1,function(x) any(is.na(x))),]
  # image(x=raw@time,y=as.numeric(rownames(rawMCorr)),t(rawMCorr))
  # abline(h=peak$Mz)
  return(rawMCorr)
}

#t<-raw@time

#process<-deconv2dLinearCoupled(rawM,t,peak,raw,listCalib,smoothParam=1)
#apply(process$predPeak,2,plot)
#plot(colSums(rawM))
#lines(colSums(process$predRaw))
#image(t(rawM))
#image(t(process$predRaw))
deconv2dLinearCoupled<-function(rawM,t,peak.detect,raw,listCalib,K=round(length(t)/2),smoothParam=100){
  
  #mass shifed correction
  mzNom<-round(peak.detect$Mz)[1]
  rawM<-calibShiftCorr(rawM,raw,peak.detect,listCalib)
  m<-as.numeric(row.names(rawM))
  
  #mass function
  n.peak<-nrow(peak.detect)
  sp <-function(m,par) {
    1/((cosh(par["lf"]*(m-par["Mz"]))^2)*(m<=par["Mz"]) + (cosh(par["lr"]*(m-par["Mz"]))^2)*(m >par["Mz"]))}
  
  #time function 
  tic<- function(t,K,tRef,smp){
    smooth<-mgcv::smooth.construct(mgcv::s(t,k=K,bs="cr",sp=smp),data=list(t=tRef),
                                   list(t=quantile(tRef,probs = seq(0,1,length=K))))
    return(list(TIC=mgcv::Predict.matrix(smooth,data.frame(t=t)),S=smooth$S))
  }
  
  #fit 
  #g<-ggplot2::ggplot()+ggplot2::geom_point(mapping = ggplot2::aes(time,EIC,colour=param),
  #                                          data=data.frame(time=raw@time,EIC=colSums(rawM),param="raw"))
  
  
  par<-peak.detect[,c("Mz","parameter.1","parameter.2")]
  colnames(par)<-c("Mz","lf","lr")
  SP <- apply(par,1,function(z) sp(m,z) )
  TICconstruc<- tic(t,K,t,smoothParam)
  TIC<-TICconstruc$TIC
  S<-  diag(rep(1,n.peak)) %x% (smoothParam *TICconstruc$S[[1]])
  X<-SP %x% TIC
  Y<-matrix(t(rawM),ncol=1)
  param<- solve(t(X)%*%X+S) %*% t(X) %*% Y
  predRaw <- t(matrix(X %*% param, ncol=nrow(rawM),nrow=length(t)))
  
  # g<-g+ggplot2::geom_line(mapping = ggplot2::aes(time,EIC,colour=param),
  #                         data=data.frame(time=raw@time,EIC=colSums(predRaw),param=as.character(smoothParam)))
  # plot(colSums(rawM),main=paste(K,smoothParam))
  # lines(colSums(predRaw),col="red")
  
  #deconvolution in ppb 
  predPeak<-matrix(0,ncol=n.peak,nrow=length(t))
  for (i in seq_len(nrow(peak.detect))){
    pred <- t(matrix(SP[,i] %x% TIC %*% param[((i-1)*K+1):(i*K)],nrow=length(t)) )
    # quanti<- ppbConvert(peakList = cbind(Mz=rep(peak.detect$Mz[i]),quanti=colSums(pred)),
    #                               primaryIon = primaryIon,
    #                               transmission = transmission,
    #                               U = U, 
    #                               Td = Td, 
    #                               pd = pd)
    #plot(colSums(pred),main=peak.detect$Mz[i],pch=19)
    predPeak[,i]<-colSums(pred)
  }
  
  # image(t(rawM))
  # image(t(predRaw))
  return(list(predRaw=predRaw,predPeak=predPeak))
}


#process<-deconv2d2linearIndependant(rawM,t,peak,raw,listCalib)
#apply(process$predPeak,2,plot)
#plot(colSums(rawM))
#lines(colSums(process$predRaw))
#image(t(rawM))
#image(t(process$predRaw))
deconv2d2linearIndependant<-function(rawM,time,peak.detect,raw,listCalib){
  
  shift<-calibShitEstimate(raw,peak.detect,listCalib)
  
  Indextime<-lapply(listCalib,function(x) {x$index})
  Indextime<-Reduce(rbind,lapply(seq_len(length(listCalib)),function(i) cbind(Indextime[[i]],i)))
  
  mzAxis<-as.numeric(row.names(rawM))
  mzNom<-round(peak.detect$Mz)[1]
  
  n.peak <- nrow(peak.detect)
  matRaw<-matrix(0,nrow=nrow(rawM),ncol=ncol(rawM))
  matPeak<-matrix(0,ncol=n.peak,nrow=ncol(rawM))
  #mzlargeIndex<-which(raw@mz > mzNom-0.4 & raw@mz < mzNom+0.4)
  for (i in seq_along(time)){
    
    #base line correction 
    #spectrum.large<-raw@rawM[mzlargeIndex,i]
    #bl<- MALDIquant::estimateBaseline(spectrum.large, method="SNIP", iterations=75)
    #spectrum.large <- spectrum.large - snipBase(spectrum.large)
    
    # first<-which( raw@mz[mzlargeIndex]>= mzAxis[1])[1]
    # spectrum.m <-spectrum.large[first:(nrow(rawM)+first-1)]
    spectrum.m <- rawM[,i]
    
    #initialisation
    mz<-peak.detect$Mz+ shift[,Indextime[i,2]]
    lf<-peak.detect$parameter.1
    lr<-peak.detect$parameter.2
    param<-cbind(mz,lf,lr)
    
    model<- apply(param,1,
                  function(x) 1/(cosh(x["lf"]*(mzAxis-x["mz"]))^2*(mzAxis<=x["mz"])+
                                   cosh(x["lr"]*(mzAxis-x["mz"]))^2*(mzAxis>x["mz"])))
    
    fit<-lm(spectrum.m ~ model-1)
    
    par_estimated<-fit$coefficients
    fit.peak<- fit$fitted.values
    
    matRaw[,i]<-fit.peak
    
    #quantification 
    quanti <- apply(rbind(mz,par_estimated,lf,lr),2, function(x){
      sum(sech2(x[1],x[3],x[4],x[2],mzAxis))})
    
    ## normalisation
    #if there is reaction ans transmission information
    # if(normalize){
    #   #normalize by primary ions
    #   quanti <- quanti/(primaryIon*4.9*10^6)
    #   namesQuanti<-"quanti (ncps)"
    # }
    
    matPeak[i,]<-quanti
    
    
  }
  return(list(predRaw=matRaw,predPeak=matPeak))
}


#process<-deconv2d2NonlinearIndependant(rawM,t,peak,raw,listCalib)
#apply(process$predPeak,2,plot)
#plot(colSums(rawM))
#lines(colSums(process$predRaw))
#image(t(rawM))
#image(t(process$predRaw))
deconv2d2NonlinearIndependant<-function(rawM,time,peak.detect,raw,listCalib){
  
  shift<-calibShitEstimate(raw,peak.detect,listCalib)
  
  Indextime<-lapply(listCalib,function(x) {x$index})
  Indextime<-Reduce(rbind,lapply(seq_len(length(listCalib)),function(i) cbind(Indextime[[i]],i)))
  
  mzNom<-round(peak.detect$Mz)[1]
  n.peak <- nrow(peak.detect)
  matRaw<-matrix(0,nrow=nrow(rawM),ncol=ncol(rawM))
  matPeak<-matrix(0,ncol=n.peak,nrow=ncol(rawM))
  for (i in seq_along(time)){
    
    #base line correction 
    mzlargeIndex<-which(raw@mz > mzNom-0.4 & raw@mz < mzNom+0.4)
    spectrum.large<-raw@rawM[mzlargeIndex,i]
    #spectrum.large <- spectrum.large - snipBase(spectrum.large)
    
    mzAxis<-as.numeric(row.names(rawM))
    #first<-which( raw@mz[mzlargeIndex]>= mzAxis[1])[1]
    #spectrum.m <-spectrum.large[first:(nrow(rawM)+first-1)]
    spectrum.m<-rawM[,i]
    plot(mzAxis,spectrum.m,type="l")
    
    #initialisation
    mz<-peak.detect$Mz+ shift[,Indextime[i,2]]
    deltal<-peak.detect$parameter.1
    deltar<-peak.detect$parameter.2
    h<- vapply(mz, function(m) max(max(spectrum.m[which(abs(mzAxis-m)<(m*50/10^6))]),0),0)
    initMz <- matrix(c(mz, h),nrow=n.peak)
    lf<-peak.detect$parameter.1
    lr<-peak.detect$parameter.2
    colnames(initMz)<-c("m","h")
    
    ppmCons<-30
    #constrain
    
    lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak),
                                        rep(0, n.peak)),ncol = 2) 
                      -
                        matrix(c(initMz[,"m"]*ppmCons/10^6,
                                 rep(0, n.peak)),ncol = 2)))
    
    upper.cons <- c(t( initMz * matrix(c(rep(1, n.peak),
                                         rep(Inf, n.peak)),ncol = 2) 
                       +
                         matrix(c(initMz[,"m"]*ppmCons/10^6,
                                  rep(0, n.peak)),ncol = 2)))
    
    #fit
    fit_function_str<-"sech2"
    par_var_str<-c("m","h")
    par_fix_str<-c("x")
    formul.character <- ""
    init.names<- ""
    for (j in seq(1, n.peak)){
      par_str.j<- paste(rep("par$",2),par_var_str,rep(j,2),sep = "")
      
      for(n in seq_along(par_var_str)){
        init.names<-c(init.names,paste(par_var_str[n],j,sep=""))
      }
      formul.character <- paste(formul.character,fit_function_str,"(",
                                par_str.j[1],",lf[",j,"],lr[",j,"],", par_str.j[2],",",par_fix_str,")+",sep="")
      
      if (j == n.peak){
        formul.character <- substr(formul.character, start=1, stop=nchar(formul.character) - 1)
      }
    }
    initMz.names<-init.names[-1]
    func.eval <- parse(text = formul.character)
    
    function.char <- function(par, x, y){
      eval(func.eval) - y
    }
    
    initMz<-as.list(t(initMz))
    names(initMz)<-initMz.names
    
    #Regression
    fit<-suppressWarnings(minpack.lm::nls.lm(par=initMz, lower = lower.cons, upper = upper.cons,
                                             fn =function.char,
                                             x= mzAxis , y = spectrum.m))
    
    par_estimated<-matrix(unlist(fit$par),nrow=2)
    fit.peak<- function.char(fit$par,mzAxis,rep(0,length(spectrum.m)))
    
    matRaw[,i]<-fit.peak
    
    #quantification 
    quanti <- apply(rbind(par_estimated,lf,lr),2, function(x){
      lines(mzAxis,sech2(x[1],x[3],x[4],x[2],mzAxis))
      return(sum(sech2(x[1],x[3],x[4],x[2],mzAxis)))})
    
    ## normalisation
    #if there is reaction ans transmission information
    # if(normalize){
    #   #normalize by primary ions
    #   quanti <- quanti/(primaryIon*4.9*10^6)
    #   namesQuanti<-"quanti (ncps)"
    # }
    
    matPeak[i,]<-quanti
    
    
  }
  return(list(predRaw=matRaw,predPeak=matPeak))
}


processFileAvgExp <-function(fullNamefile, massCalib,primaryIon,indTimeLim, mzNominal, ppm, 
                             ppmGroupBkg, fracGroup,
                             minIntensity, 
                             fctFit,normalize){
  
  if(is.character(fullNamefile)){
    cat(basename(fullNamefile),": ")
    # read file
    raw <- readRaw(fullNamefile, calibTIS = FALSE)
  } else raw <- fullNamefile

  # information for ppb convertion
  reaction = raw@prtReaction
  transmission = raw@ptrTransmisison
  
  # time limit and background
  bg<-FALSE
  matPeak.bg<-NULL
  indLim <- indTimeLim
  if(ncol(indLim) == 0 || all(is.na(indLim)) ) {
    indLim <- matrix(c(1,length(raw@time)),
                     dimnames =list(c("start","end"),NULL) )
  } else {
    bg<-TRUE
    indLimExpLarge <- timeLimits(raw,fracMaxTIC = 0.1)
    indExp<-unique(Reduce(c,apply(indLimExpLarge,2,function(x) seq(x["start"]-2,x["end"]+2))))
    indBg<-seq(1,length(raw@time))[-indExp]
    raw.bg <- raw
    raw.bg@rawM <- raw.bg@rawM[, indBg]
    raw.bg@time <-  raw.bg@time[ indBg]
    
    ## calibration
    raw.bg <- calibration(raw.bg, massCalib)
    
    ## peak detection 
    if(is.null(mzNominal)) mzNominal= unique(round(raw@mz))
    list_peak <- PeakList(raw=raw.bg, mzNominal = mzNominal, 
                          ppm=ppm, minIntensity=minIntensity,
                          fctFit=fctFit)
    
    ## normalisation
    #if there is reaction ans transmission information
    namesQuanti<-"quanti (cps)"
    if(normalize){
       
        #normalize by primary ions
        list_peak$peak[,"quanti"] <- list_peak$peak[,"quanti"]/(primaryIon*4.9*10^6)
        namesQuanti<-"quanti (ncps)"
      
    }
    matPeak.bg<- cbind(list_peak$peak , group=rep(0,dim(list_peak$peak)[1]))
  } 
  

    indLimAll<-Reduce(c,apply(indLim,2,function(x) seq(x[1],x[2])))
    ## average on expriration
    raw.i <- raw
    raw.i@rawM <- raw.i@rawM[, indLimAll]
    raw.i@time <-  raw.i@time[ indLimAll]
    
    ## calibration
    raw.i <- calibration(raw.i, massCalib)
    
    ## peak detection 
    if(is.null(mzNominal)) mzNominal= unique(round(raw@mz))
    list_peak <- try(PeakList(raw=raw.i, mzNominal = mzNominal, 
                          ppm=ppm, minIntensity=minIntensity,
                          fctFit=fctFit))
    
    
    ## normalisation
    #if there is reaction ans transmission information
    namesQuanti<-"quanti (cps)"
    if(normalize){
        #normalize by primary ions
        list_peak$peak[,"quanti"] <- list_peak$peak[,"quanti"]/(primaryIon*4.9*10^6)
        namesQuanti<-"quanti (ncps)"
    }
    
    list_peak.exp<- cbind(list_peak$peak , group=rep(1,dim(list_peak$peak)[1]))
    
  
  ## Concatenate all expirations peak lists 
  
  matPeak.j<-rbind(matPeak.bg,list_peak.exp)
  
  ## align expirations and bakcground
  matAligned <- alignExpirations(matPeak.j, ppmGroup=ppmGroupBkg, fracGroup=fracGroup)
  names(matAligned)[grep("quanti",names(matAligned))]<-namesQuanti
  
  #return(matPeak.j)
  return(list(raw=matPeak.j,aligned=matAligned))
}

processFileSepExp <-function(fullNamefile, massCalib,primaryIon,indTimeLim, mzNominal, ppm, 
                      ppmGroupBkg, fracGroup,
                      minIntensity, 
                      fctFit,normalize){
  
  if(is.character(fullNamefile)){
    cat(basename(fullNamefile),": ")
    # read file
    raw <- readRaw(fullNamefile, calibTIS = FALSE)
  } else raw <- fullNamefile
  
  
  # information for ppb convertion
  reaction = raw@prtReaction
  transmission = raw@ptrTransmisison
  
  # time limit and background
  bg <- FALSE
  matPeak.bg <- NULL
  indLim <- indTimeLim
  if(ncol(indLim) == 0 || all(is.na(indLim)) ) {
    indLim <- matrix(c(1,length(raw@time)),
                     dimnames =list(c("start","end"),NULL) ) } 
  if( indLim[1,1] > 1 || utils::tail(c(indLim),1) < length(raw@time) ) {
    bg <- TRUE
    LimBg<-bakgroundDetect(colSums(raw@rawM))
    LimBg<-LimBg[,1] 
    indBg<-seq(LimBg[1],LimBg[2])
    raw.bg <- raw
    raw.bg@rawM <- raw.bg@rawM[, indBg]
    raw.bg@time <-  raw.bg@time[ indBg]
    
    ## calibration
    raw.bg <- calibration(raw.bg, massCalib)
    
    ## peak detection 
    if(is.null(mzNominal)) mzNominal= unique(round(raw@mz))
    list_peak <- PeakList(raw=raw.bg, mzNominal = mzNominal, 
                          ppm=ppm, minIntensity=minIntensity,
                          fctFit=fctFit)
    
    ## normalisation
    #if there is reaction ans transmission information
    namesQuanti<-"quanti (cps)"
    if(normalize){
      
        #normalize by primary ions
        list_peak$peak[,"quanti"] <- list_peak$peak[,"quanti"]/(primaryIon*4.9*10^6)
        namesQuanti<-"quanti (ncps)"
    }
    matPeak.bg<- cbind(list_peak$peak , group=rep(0,dim(list_peak$peak)[1]))
  } 
  
  list_peak.j <- NULL
  n.limit <- dim(indLim)[2]
  for (i in seq_len(n.limit)){
    ## average on selected expiration or background
    raw.i <- raw
    raw.i@rawM <- raw.i@rawM[, seq(indLim["start", i], indLim["end", i])]
    raw.i@time <-  raw.i@time[ seq(indLim["start", i], indLim["end", i])]
    
    ## calibration
    raw.i <- calibration(raw.i, massCalib)
    
    ## peak detection 
    if(is.null(mzNominal)) mzNominal= unique(round(raw@mz))
    list_peak <- PeakList(raw=raw.i, mzNominal = mzNominal, 
                          ppm=ppm, minIntensity=minIntensity,
                          fctFit=fctFit)
    
    ## normalisation
    #if there is reaction ans transmission information
    namesQuanti<-"quanti (cps)"
    if(normalize){
      
        #normalize by primary ions
        list_peak$peak[,"quanti"] <- list_peak$peak[,"quanti"]/(primaryIon*4.9*10^6)
        namesQuanti<-"quanti (ncps)"
      
    }
    
   list_peak.j[[i]]<- cbind(list_peak$peak , group=rep(i,dim(list_peak$peak)[1]))
    
  } # end loop over expirations
  
  ## Concatenate all expirations peak lists 
  matPeak.j <- do.call(rbind, list_peak.j) 
  
  if(bg) matPeak.j <- rbind(matPeak.bg,matPeak.j)
  
  ## align expirations and bakcground
  matAligned <- alignExpirations(matPeak.j, ppmGroup=ppmGroupBkg, fracGroup=fracGroup)
  names(matAligned)[grep("quanti",names(matAligned))]<-namesQuanti
  
  #return(matPeak.j)
  return(list(raw=matPeak.j,aligned=matAligned))
}


processFileTemporal<-function(fullNamefile, massCalib,primaryIon,indTimeLim, mzNominal, ppm, 
                      ppmGroupBkg, fracGroup,
                      minIntensity, 
                      fctFit,normalize,...){
  
  if(is.character(fullNamefile)){
    cat(basename(fullNamefile),": ")
    # read file
    raw <- readRaw(fullNamefile)
  } else raw <- fullNamefile
  
  # information for ppb convertion
  reaction = raw@prtReaction
  transmission = raw@ptrTransmisison
  
  # peak detection on TIS 
  if(is.null(mzNominal)) mzNominal= unique(round(raw@mz))
  list_peak <- PeakList(raw, mzNominal = mzNominal, 
                          ppm=ppm, minIntensity=minIntensity,
                          fctFit=fctFit)$peak

  
  # compute temporal profile
  fileProccess<-ptairMS:::computeTemporalFile(raw,list_peak)
  
  matPeak<-data.table::as.data.table(fileProccess$matPeak)
  #apply(matPeak,1,function(x) plot(x[5:(length(x)-2)],main=round(x[1],2)))
  
  # agregation of expirations and bg

  # adjustement of time limit
  # indLim <- indTimeLim
  # if(ncol(indLim) == 0 || is.na(indLim) ) {
  #   indLim <- matrix(c(1,length(raw@time)),
  #                    dimnames =list(c("start","end"),NULL) )
  # } else if(indLim[1,1] > 4) indLim <- cbind( matrix( c(1, indLim[1, 1] - 3)), indLim)
  # 
  # bg<-raw@time[seq(indLim["start",1],indLim["end",1])]
  # if(ncol(indLim)>1){
  #   exp<-raw@time[Reduce(c,lapply(seq(2,ncol(indLim)),function(i) seq(indLim["start",i],
  #                                                                     indLim["end",i])))]
  #   matPeakAg<-matPeak[ , list(mz=stats:: median(mz), quanti=mean(quanti[t %in% exp]), 
  #                           background=mean(quanti[t %in% bg])),by=group]
  #   #matPeakAg[group `:=` NULL] #delete column group 
  # }else matPeakAg<-matPeak[,list(mz=stats:: median(mz), quanti=mean(quanti)),by=group]
  # 
  # names(matPeakAg)[grep("quanti",names(matPeakAg))]<-namesQuanti
  
  return(list(raw=matPeak,aligned=NULL))
}


#' Detect peak in a single nominal mass, same parameter as peakList
#' @param i the nominal mass
#' @param mz peak mass axis
#' @param sp peak spectrum 
#' @param ppmPeakMinSep the minimum distance between two peeks in ppm 
#' @param calibCoef calibration coeficient a and b
#' @param minPeakDetect the minimum intenisty for peaks detection. The final threshold for peak detection
#' will be : max ( \code{minPeakDetect} , thresholdNoise ). The thresholdNoise correspond to
#'  max(\code{thNoiseRate} * max( noise around the nominal mass), \code{thIntensityRate} * 
#'  max( intenisty in the nominal mass). The noise around the nominal mass correspond : 
#'  \code{[m-windowSize-0.2,m-windowSize]U[m+windowSize,m+WindowSize+0.2]}.
#' @param fitFunc the function for the quantification of Peak, should be average or Sech2
#' @param maxIter maximum ittertion of residual analysis
#' @param autocorNoiseMax the autocorelation threshold for Optimal windows Savitzky Golay 
#' filter in \code{OptimalWindowSG} ptairMS function. See \code{?OptimalWindowSG}
#' @param plotFinal boolean. If TRUE, plot the spectrum for all nominal masses, with the final fitted peaks
#' @param plotAll boolean. Tf TRUE, plot all step to get the final fitted peaks
#' @param thNoiseRate The rate who is multiplie by the max noise intensity
#' @param thIntensityRate The rate who is mutluplie by the max signal intensity
#' @param countFacFWHM integer. We will sum the fitted peaks on a windows's size of countFacFWHM * FWHM, 
#' centered in the mass peak center.
#' @param daSeparation the minimum distance between two peeks in Da for nominal mass < 17.
#' @param d the degree for the \code{Savitzky Golay} filtrer
#' @param windowSize peaks will be detected only around  m - windowSize ; m + windowSize, for all 
#' m in \code{mzNominal}
#' @return a list that contains the peak list of nomminal mass i, and information to plot the peaks detected
peakListNominalMass <- function(i,mz,sp,ppmPeakMinSep=130 ,calibCoef,
                                minPeakDetect=10, fitFunc="Sech2", maxIter=2, autocorNoiseMax=0.3 ,
                                plotFinal=FALSE, plotAll=FALSE, thNoiseRate=1.1, thIntensityRate=0.01 ,
                                countFacFWHM=10, daSeparation=0.001, d=3, windowSize=0.4 ){
  
  emptyData<- data.frame(Mz=double(),quanti=double(),delta_mz=double(),resolution=double())
  warning_mat <- NULL
  infoPlot <- list()
  no_peak_return<-list(emptyData,warning_mat,infoPlot)
  
  # select spectrum around the nominal mass i
  index.large <- which(mz < i + windowSize + 0.2 & mz > i - windowSize-0.2)
  mz.i.large <- mz[index.large]
  sp.i.large  <- sp[index.large]
  
  if( range(mz.i.large)[1] > i - windowSize | range(mz.i.large)[2] < i + windowSize  ) 
    return(no_peak_return)
  # baseline correction
  sp.i.large.corrected <- sp.i.large - snipBase(sp.i.large)
  
  # noise auto correlation estimation
  index.noise <- c (which(i - windowSize > mz & mz > i - windowSize - 0.2),
                 which(i + windowSize < mz & mz < i + windowSize + 0.2) )
  noise <- sp.i.large.corrected[ match( index.noise , index.large )]
  
  real_acf<- stats::acf(noise ,lag.max=1, plot=FALSE)[1]$acf
  if(is.na(real_acf)) return(no_peak_return)
  noiseacf <- min( real_acf , autocorNoiseMax) # max auto correlation at 0.3 indeed OptimnalWinfowSg will overfitted
  thr <- max(noise) # dynamic threshold for peak detetcion
  
  # signal who wil be process
  if(i > 250){
    mz.i <- mz.i.large
    sp.i <- sp.i.large.corrected
  }else{
    index <-which( i - windowSize  <= mz & mz <= i + windowSize )
    mz.i <- mz[ index ]
    sp.i <- sp.i.large.corrected[ match( index , index.large ) ]
  }
  
  infoPlot[[as.character(i)]]<-list(mz=mz.i,sp=sp.i,main=i,pointsPeak=NA)
  no_peak_return<-list(emptyData,warning_mat,infoPlot)
  
  ## fit itterative
  sp.i.fit <- sp.i  # initialization of sp.i.fit
  c=1 # initialize number of itteration
  
  repeat{
    
    # Initialization for regression :   
    
    if(c==1) { minpeakheight <- max( max(thr*thNoiseRate , thIntensityRate*max(sp.i),
                                         minPeakDetect)) 
    } else minpeakheight <- max( minpeakheight*0.8 , 1 ) # minimum intenisty
    
    init <- initializeFit(i,sp.i.fit, sp.i, mz.i, calibCoef,
                          minpeakheight, noiseacf,  ppmPeakMinSep,
                          daSeparation,  d, plotAll,c) # Find local maximum with Savitzky golay filter or wavelet
    
    # Regression :
    if(!is.null(init)){
      
      mz_init <- init$mz[,"m"]

      if(c > 1){
        mz_par <- par_estimated[1,]
        
        # delete peak also find
        to_delete <- NULL
        for (p in seq_along(par_estimated[1,])){
          Da_proxi<- abs(mz_par[p]-mz_init)
          if(i > 17) test_proxi<- Da_proxi*10^6/i > ppmPeakMinSep
          if(i <=17) test_proxi<- Da_proxi > daSeparation
          if(!all(test_proxi)) to_delete<-c(to_delete,p)
        }
        
        n.peak <- nrow(init$mz)
        # if same number of peak detected and all old peak deleted, no new peak are detected
        if( n.peak == dim(par_estimated)[2] & 
            length(to_delete) == dim(par_estimated)[2] ) break
        
        # add new peak
        if( !is.null(to_delete) ) {
          par_estimated<- par_estimated[,-to_delete,drop=FALSE]
        } 
      } # END if c> 1
      
      n.peak <- nrow(init$mz)
      
      if(fitFunc == "average"){
          
          if(c == 1) initTof <- init$tof else initTof <- rbind(init$tof,t(par_estimated))
          n.peak<- nrow(initTof)
          
          resolution_upper<-15000
          resolution_lower<-8000
          
          lower.cons <- c(t(initTof * matrix(c(rep(1, n.peak),
                                               rep(0, n.peak), 
                                               rep(0, n.peak)),ncol = 3)  
                            -
                              matrix(c(initTof[,"delta_mz"]/4,
                                       - initTof[ ,"t"]/resolution_upper,
                                       rep(0, n.peak)),ncol = 3)  ))
          
          upper.cons <- c(t(initTof * matrix(c(rep(1, n.peak),
                                               rep(0, n.peak), 
                                               rep(Inf, n.peak)),ncol = 3)  
                            +
                              matrix(c(initTof[,"delta_mz"]/4,
                                       initTof[ ,"t"]/resolution_lower,
                                       rep(0, n.peak)),ncol = 3)  ))
          
          l.shape<-determinePeakShape(sp,mz)
          fit <- fit_averagePeak(initTof, l.shape, sp.i, mzToTof(mz.i,calibCoef)) 
        }
      if(fitFunc =="Sech2"){
        if(c == 1) initMz <-init$mz else initMz <- rbind(init$mz,t(par_estimated))
        
        n.peak<- nrow(initMz)
        resolution_upper<-8000
        resolution_mean<- 5000
        resolution_lower<-3500
        
        lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak),
                                            rep(0, n.peak*2),
                                            rep(0, n.peak)),ncol = 4) 
                          -
                            matrix(c(initMz[,"m"]/(resolution_mean*4),
                                     -resolution_lower*log(sqrt(2)+1)*2/initMz[,"m"],
                                     -resolution_lower*log(sqrt(2)+1)*2/initMz[,"m"],
                                     rep(0, n.peak)),ncol = 4)))
        
        upper.cons <- c(t( initMz * matrix(c(rep(1, n.peak),
                                             rep(0, n.peak*2),
                                             rep(Inf, n.peak)),ncol = 4) 
                           +
                             matrix(c(initMz[,"m"]/(resolution_mean*4),
                                      resolution_upper*log(sqrt(2)+1)*2/initMz[,"m"],
                                      resolution_upper*log(sqrt(2)+1)*2/initMz[,"m"],
                                      rep(0, n.peak)),ncol = 4)))
        
        fit <- fit_sech2(initMz, sp.i, mz.i, lower.cons, upper.cons)
      }
      
      if(is.na(fit$fit$deviance)){ if(c==1) return(no_peak_return) else break }
      if (fit$fit$niter== 50 ) warning_mat<-rbind(warning_mat,c(i,NA,"max itter lm algo"))
      
      fit.peak <- fit$fit.peak
      par_estimated <- fit$par_estimated
      cum_function.fit.peak <- fit$function.fit.peak
      
      # residual 
      sp.i.fit <- sp.i - fit.peak
      
      #indicators
      R2 <- 1 - sum(sp.i.fit^2)/sum((sp.i-mean(sp.i))^2)
      auto_cor_res <- abs(stats::acf(sp.i.fit,plot=FALSE)[1]$acf[1])
      
      if(plotAll){
        plot(mz.i , sp.i, main=paste(i,"fit itteration :",c),xlab="mz",ylab="intensity",
             type="b",pch=19,cex=0.7)
        graphics::lines(mz.i,fit.peak,col="blue",lwd=2)
        graphics::lines(mz.i,sp.i.fit,col="green3",lwd=2)
        
        if(fitFunc=="average"){
          for(k in seq_len(ncol(par_estimated))){
            graphics::lines(mz.i,
                  fit_averagePeak_function(par_estimated[1,k],par_estimated[2,k],par_estimated[3,k],
                                      l.shape$tofRef,l.shape$peakRef,mzToTof(mz.i,calibCoef)),
                  lwd=2,col="red",lty=2)
          }
          graphics::points(tofToMz(par_estimated[1,],calibCoef),
                 cum_function.fit.peak(fit$fit$par,l.shape$tofRef,l.shape$peakRef,par_estimated[1,],rep(0,length(par_estimated[2,]))),
                 cex=2,col="red",lwd=2)
        }
        if(fitFunc=="Sech2"){
          for(k in seq_len(ncol(par_estimated))){
            graphics::lines(mz.i,
                  sech2(par_estimated[1,k],par_estimated[2,k],
                        par_estimated[3,k],par_estimated[4,k],mz.i),
                  lwd=2,col="red",lty=2)
          }
          graphics::points(par_estimated[1,],cum_function.fit.peak(fit$fit$par,par_estimated[1,],rep(0,length(par_estimated[2,]))),
                 cex=2,col="red",lwd=2)
        }
        
        graphics::legend("topleft",legend=c("Raw","fit sum","fit peak","residual",paste("R2=",round(R2,3)),paste("autocor res=",round(auto_cor_res,3))),
               col=c("black","blue","red","green3"),
               lty=c(NA,1,2,1,NA,NA),pch=c(19,NA,NA,NA,NA,NA))
      }
    } else {
      ## if init is null : no peak find
      if (c==1){
        # if first iteration 
        if(plotFinal){
          plot(mz.i , sp.i, main=paste(i ,c,"fit",sep=" - "), xlab="mz", ylab="intensity",
               type="p",pch=19,cex=0.7)
        }
        return(no_peak_return) } else break}
    
    # condition for break
    if( auto_cor_res< 0.3) break
    if(R2 > 0.995) break
    if(c > 1 ) {
      #degradation of R2
      if(R2 < old_R2) {
        fit<-fit_old
        
        fit.peak<-fit$fit.peak
        par_estimated<-fit$par_estimated
        cum_function.fit.peak<-fit$function.fit.peak
        
        sp.i.fit<-sp.i-fit.peak
        n.peak <- dim(par_estimated)[2]
        
        break
      }
    }
    if(n.peak > 5 ) break
    
    # Iteration limit 
    if(c== maxIter) {
      warning_mat<-rbind(warning_mat,c(i,NA,"residual iteration max"))
      break
    }
    c=c+1
    
    fit_old <- fit
    old_R2 <- R2
    
  } # END OF REPEAT
  
  
  # last fit less constrained
  c=1
  repeat{ # until no peak to close or inside an other peak
    init <- t(par_estimated)
    n.peak <- nrow(init)
    
    if(fitFunc=="average"){
      lower.cons <- c(t(matrix(c(rep(0, n.peak),rep(0, n.peak),init[2,]/15000),ncol = 3)))
      
      upper.cons <- c(t(matrix(c(rep(Inf, n.peak),rep(Inf, n.peak),init[2,]/6000),ncol = 3)))
      
      bin.i<- mzToTof(mz.i,calibCoef)
      fit <- fit_averagePeak(init,l.shape,sp.i,bin.i,lower.cons,upper.cons)
      
      fit.peak<-fit$fit.peak
      par_estimated<-fit$par_estimated
      
      delta_mz <- apply(par_estimated,2, function(x) diff(tofToMz(c(x[2]-x[3]/2,x[2]+x[3]/2),calibCoef)) )
      quanti<- apply(par_estimated,2, function(x){
        th<-countFacFWHM*0.5*diff(tofToMz(c(x[2]-x[3]/2,x[2]+x[3]/2),calibCoef))
        bin.x<- which( tofToMz(x[2]-1,calibCoef) - th < mz & mz < tofToMz(x[2]-1,calibCoef)+th )
        sum(fit_averagePeak_function(x[1],x[2],x[3],l.shape$tofRef,l.shape$peakRef,bin.x))})
      center_peak<- unname(tofToMz(par_estimated[2,]-1,calibCoef))
      
      peaks <- apply(par_estimated,2,function(x)
        fit_averagePeak_function(x[1],x[2],x[3],l.shape$tofRef,l.shape$peakRef,bin.i))
      
    }
    
    if(fitFunc=="Sech2"){
      lower.cons <- c( t(matrix(c(rep(0, n.peak),
                                  rep(3000*log(sqrt(2)+1)*2/init[,1],2),
                                  rep(0, n.peak)),ncol = 4)))
      
      upper.cons <-  c( t(matrix(c(rep(Inf, n.peak)
                                   ,rep(9000*log(sqrt(2)+1)*2/init[,1],2),
                                   rep(Inf, n.peak)),ncol = 4)))
      
      fit <- fit_sech2(init, sp.i, mz.i, lower.cons, upper.cons)
      fit.peak<-fit$fit.peak
      par_estimated<-fit$par_estimated
      
      delta_mz <- apply(par_estimated,2, function(x) log(sqrt(2)+1)/x[2]+log(sqrt(2)+1)/x[3] )
      quanti<- apply(par_estimated,2, function(x){
        th<-countFacFWHM*0.5*(log(sqrt(2)+1)/x[2]+log(sqrt(2)+1)/x[3])
        mz.x <- mz[ x[1] - th < mz & mz < x[1]+th ]
        sum(sech2(x[1],x[2],x[3],x[4],mz.x))})
      center_peak<- par_estimated[1,]
      
      peaks <-apply(par_estimated,2,function(x) sech2(x[1],x[2],x[3],x[4],mz.i))
    }
    
    R2 <-1-sum(sp.i.fit^2)/sum((sp.i-mean(sp.i))^2)
    
    X<-data.frame(Mz=center_peak,
             quanti=quanti,delta_mz=delta_mz, resolution = center_peak/delta_mz, 
             R2=R2 , parameter = t(par_estimated[-1,])  )
    X<-X[quanti>1,,drop=FALSE]
    par_estimated<-par_estimated[,quanti>1,drop=FALSE]
    
    X <- X[order(X[,1]),,drop=FALSE]
    par_estimated<-par_estimated[,order(par_estimated[1,]),drop=FALSE]
    if(dim(X)[1]==1) break
    
    to_delete<-NULL
    #tets if peak is under an other peak
    for(j in seq_len(nrow(X)-1)){
      sign_diff<-levels(factor(sign(peaks[,j+1]-peaks[,j])))
      sign_diff<-sign_diff[sign_diff!=0]
      if(length(sign_diff) <2) {
        if("-1" %in% sign_diff) to_delete <- c(to_delete,j+1) else to_delete <- c(to_delete,j)
      }
    }
    #if no
    if(is.null(to_delete)) {
      #tets proximity
      control_proxi <- c(FALSE,diff(X[,1])*10^6/i < ppmPeakMinSep)
      #if no break
      if(all(!control_proxi)) break
      #else delete
      par_estimated<-par_estimated[,-which(control_proxi),drop=FALSE]
    } else {par_estimated<-par_estimated[,-to_delete,drop=FALSE]
    #if more than one peak
    if(dim(X)[1]>1) {
      #test proximity
      control_proxi <- c(FALSE,diff(X[,1])*10^6/i < ppmPeakMinSep)
      if(any(control_proxi)) par_estimated<-par_estimated[,-which(control_proxi),drop=FALSE]
    } }
    c=c+1
  } # end second repeat
  if(fitFunc=="average") pointsPeak<- list(x=tofToMz(par_estimated[1,],calibCoef),
                                          y=fit$function.fit.peak(fit$fit$par,
                                                                l.shape$tofRef,l.shape$peakRef,
                                                                par_estimated[1,],
                                                                rep(0,length(par_estimated[2,]))))
  if(fitFunc=="Sech2") pointsPeak<- list(x=par_estimated[1,],
                                         y=fit$function.fit.peak(fit$fit$par,par_estimated[1,],
                                                               rep(0,length(par_estimated[2,]))))
  infoPlot[[as.character(i)]]<-list(mz=mz.i,sp=sp.i,main=i,fitPeak=fit.peak,
                 peak=peaks, pointsPeak=pointsPeak)
  if(plotFinal){
    
    plot(mz.i, sp.i, main=i, xlab="mz",ylab="intensity",
         ylim=c(min(sp.i,fit.peak),max(sp.i,fit.peak)),type="b",pch=19,cex=0.7)
    graphics::lines(mz.i,fit.peak,lwd=2,col="blue")
    
    for(k in seq_len(ncol(par_estimated))) graphics::lines(mz.i,peaks[,k] ,lwd=2,col="red",lty=2) 
    
    if(fitFunc=="average") graphics::points(pointsPeak,
                                 cex=2,col="red",lwd=2,pch=19)
    if(fitFunc=="Sech2") graphics::points(pointsPeak,
                                cex=2,col="red",lwd=2,pch=19)
    
    graphics::legend("topleft",legend=c("Raw","fit sum","fit peak","peak center",paste("R2=",round(R2,3)),paste("autocor res=",round(auto_cor_res,3))),
           col=c("black","blue","red","red"),
           lty=c(NA,1,2,NA,NA,NA),pch=c(19,NA,NA,19,NA,NA))
  }
  
  if(!is.null(warning_mat)){
    warning_mat<-data.frame(warning_mat)
    names(warning_mat)<-c("m/z","peak","type")
  }

  return(list(peak=X,warning=warning_mat,plot=infoPlot))
}


## initialization ----
#' initialization for apply fit function in the spectrum
#' @param i the nominal mass
#' @param sp.i.fit the vector who will be fetted (spectrum pf residual)
#' @param sp.i the spectrum around a nominal mass
#' @param mz.i the mass vector around a nominal mass
#' @param calibCoef calibration coeficient
#' @param minpeakheight the minimum peak intensity
#' @param noiseacf aytocorelation of the noise
#' @param ppmPeakMinSep the minimum distance between two peeks in ppm 
#' @param daSeparation the minimum distance between two peeks in da
#' @param d the degree of savitzky golay filter
#' @param plotAll bollean if TRUE, it plot all the initialiation step
#' @param c the number of current itteration
#' @return a list with fit input
initializeFit<-function(i,sp.i.fit, sp.i, mz.i, calibCoef, minpeakheight, noiseacf, ppmPeakMinSep,
                        daSeparation, d, plotAll,c ){
  
  init <- NULL
  
  # find local maxima in the spectrum: return the index
  prePeak <- LocalMaximaSG(sp = sp.i.fit, minPeakHeight = minpeakheight,noiseacf = noiseacf,d = d)
  
  if(!is.null(prePeak) ) {
    # get the mass and the intenisty 
    prePeak<-cbind(mz = mz.i[prePeak],
                   intensity = vapply(prePeak,
                                      function(x) max( sp.i[ (x-1):(x+1) ]),
                                      FUN.VALUE = 1.1)
                   ) #the maximum auround the index find
    
    # delete peak to close
    if(nrow(prePeak) > 1){
      #chexk proximity
      Da_proxi <- diff( prePeak[,"mz"] )
      if(i > 17) test_proxi <- c(TRUE, Da_proxi*10^6/i > ppmPeakMinSep)
      if(i <= 17) test_proxi <- c(TRUE,Da_proxi > daSeparation) ## max 
      prePeak <- prePeak[test_proxi,,drop=FALSE]
    }
    
    # plot 
    if(plotAll){
      plot(mz.i, sp.i, 
           main= paste("initialization iteration :",c), 
           xlab="mz", ylab="intensity", type="b", pch=19, cex=0.7)
      graphics::points(prePeak[,"mz"],prePeak[,"intensity"], col="red", cex=2, lwd=2.5)
      graphics::abline(h = minpeakheight)
      graphics::legend("topleft",legend=c("Raw","Local maximum"), 
             col=c("black","red"), pch=c(19,1), lty=c(1,NA))
    }
    
    # Calculate the initialization
    
      # in tof fot average fit function 
      resolution_mean <- 10000 #in tof : t/delta(t)
      t0 <- unname( mzToTof(prePeak[,"mz"],calibCoef) ) # tof peak center
      delta0 <- t0/resolution_mean # FWHM in tof
      h0 <- unname( prePeak[,"intensity"] )
      initTof<- cbind(t=t0, delta=delta0,h=h0)
    
      #in mz for sech 2 function
      resolution_mean <- 5000 #in mass m / dela(m)
      m0 <- unname( prePeak[,"mz"] ) # peak center
      delta0 <- resolution_mean * log(sqrt(2)+1)*2/m0  #lambda estimation for Sec2 function
      h0<- unname( prePeak[,"intensity"]) #peak height
      initMz <- cbind(m=m0,delta=delta0,delta2=delta0,h=h0)
    
      init<- list(mz=initMz, tof=initTof)

    } #END if prePrek not null
  return(init)
  
  
}
## fit function -----
#' Fit function Sech2
#' @param h peak higth
#' @param p peak center
#' @param lf lambda 1 , peak width left
#' @param lr lambda 2 , peak width right
#' @param x vector values
#' @return numeric value
sech2<-function(p,lf,lr,h,x) { h/(cosh(lf*(x-p))^2*(x<=p)+cosh(lr*(x-p))^2*(x>p))}

sechInv<-function(p,lf,lr,h,x){c(p-1/lf*acosh(sqrt(h/x)),p+1/lr*acosh(sqrt(h/x)))}

#'fit function average
#'@param t tof center of peak
#'@param delta FWHM of peak
#'@param h peak height
#'@param intervRef reference intreval for peak shape
#'@param peakShape peak shape estimated on \code{intervalRef}
#'@param bin bin interval of peak will be fitted
#'@return peak function made on an average of reference peaks normalized
fit_averagePeak_function <- function(t, delta,h, intervRef, peakShape, bin){
  res <- rep(0,length(bin))
  intervFit <- intervRef * delta + t
  interpol_ok <- which( intervFit[1] < bin & bin < utils::tail(intervFit,1) )
  if(length(interpol_ok) !=0 ) res[interpol_ok] <- stats::spline(intervFit,h * peakShape, xout = bin[interpol_ok])$'y'
  res
}


#'Create cumulative function fit
#'
#'@param fit_function_str fit function who will be use in character
#'@param par_var_str parameters of fit function who change with the peak in a vector of character
#'@param par_fix_str parameters of fit function indepedent of the peak in a vector of character
#'@param n.peak number of peak
#'@return a list: \cr
#'   \code{init.names}: names of paramters for the initialization \cr
#'   \code{func.eval}: function who will be fitted
cumulative_fit_function <- function(fit_function_str,par_var_str,par_fix_str,n.peak){
  formul.character <- ""
  init.names<- ""
  for (j in seq(1, n.peak)){
    par_str.j<- ""
    for (n in seq(from=length(par_var_str),to=1)){
      par_str.j<-paste(paste("par$",par_var_str[n],j,sep=""),par_str.j,sep=",")
    }

    for(n in seq_along(par_var_str)){
      init.names<-c(init.names,paste(par_var_str[n],j,sep=""))
    }
    formul.character <- paste(formul.character,fit_function_str,"(",
                              par_str.j,par_fix_str,")+",sep="")

  if (j == n.peak){
    formul.character <- substr(formul.character, start=1, stop=nchar(formul.character) - 1)
    }
  }
  init.names<-init.names[-1]
  func.eval <- parse(text = formul.character)
  return(list(init.names=init.names,func.eval=func.eval))
}

#' Fit peak with Sech2 function
#' @param initMz list of initial parmeter (mz,delta,h) 
#' @param sp spectrum
#' @param mz.i mass spectrum 
#' @param lower.cons lower constrain of fit 
#' @param upper.cons upper constrain of fit 
#' @return list with fit result
fit_sech2<-function(initMz, sp, mz.i, lower.cons, upper.cons){

  n.peak<-nrow(initMz)
  
  fit_function_str<-"sech2"
  par_var_str<-c("m","lf","lr","h")
  par_fix_str<-c("x")
  cum_fit<-cumulative_fit_function(fit_function_str,par_var_str,par_fix_str,n.peak)

  initMz.names<-cum_fit$init.names
  func.eval<-cum_fit$func.eval

  function.char <- function(par, x, y){
    eval(func.eval) - y
  }

  initMz<-as.list(t(initMz))
  names(initMz)<-initMz.names

  #Regression
  fit<-suppressWarnings(minpack.lm::nls.lm(par=initMz, lower = lower.cons, upper = upper.cons,
                               fn =function.char,
                               x= mz.i , y = sp))

  par_estimated<-matrix(unlist(fit$par),nrow=4)
  fit_peak<- function.char(fit$par,mz.i,rep(0,length(sp)))

  return(list(fit.peak=fit_peak,par_estimated=par_estimated,function.fit.peak=function.char,fit=fit))

}

#' Fit peak with average function
#' @param initTof list of initialisation in tof
#' @param l.shape peak shape average
#' @param sp spectrum 
#' @param bin tof axis
#' @param lower.cons lower constrain for fit
#' @param upper.cons upper constrain for fit
#' @return list with fit information 
fit_averagePeak<-function(initTof, l.shape, sp, bin, lower.cons, upper.cons){
  n.peak<-nrow(initTof)
  
  fit_function_str<-"fit_averagePeak_function"
  par_var_str<-c("t","delta","h")
  par_fix_str<-c("intervRef,peakShape,bin")

  cum_fit <- cumulative_fit_function(fit_function_str, par_var_str, par_fix_str, n.peak)

  initTof.names<-cum_fit$init.names
  func.eval<-cum_fit$func.eval

  function.char <- function(par,  intervRef, peakShape, bin, vec.peak){
    eval(func.eval) - vec.peak
  }

  initTof <- as.list(t(initTof))
  names(initTof)<-initTof.names

  fit<-suppressWarnings(minpack.lm::nls.lm(par= initTof, lower= lower.cons, upper=upper.cons,
                          fn = function.char,
                          intervRef=l.shape$tofRef, peakShape = l.shape$peakRef,
                          bin = bin, vec.peak= sp))



  fit.peak<-function.char( fit$par, l.shape$tofRef, l.shape$peakRef, bin, rep(0,length(sp)) )
  par_estimated<-matrix(unlist(fit$par),nrow=3)

  return(list(fit.peak=fit.peak,par_estimated=par_estimated,function.fit.peak=function.char,fit=fit))

  }

#' Determine peak shape from raw data
#'
#'This function use the method descibe by average and al 2013, for determine a peak 
#'shape from the raw data : \cr 
#'$peak_i(Delta_i,A_i,t_i) = interpolation(x= tof.ref * Delta_i + t_i,y = A_i * 
#'peak.ref, xout= TOF_i )$ where peak.ref and tof.ref are peaks reference use for mass calibration.
#' @param sp vector with spectrum values where the peak shape will be extimated
#' @param mz vectr with mass axis
#' @param massRef the accurate mass of the reference peak (for example the same that used for calibration)
#' @param plotShape if true plot each reference peak and the average peak (the peak shape)
#' @return A list of two vectors which are the reference peak normalized tof and intensity.
determinePeakShape <- function(sp, mz, massRef = c(21.022, 29.013424,41.03858,75.04406, 203.943, 330.8495),
                               plotShape=FALSE ){

  massRef.o <- massRef[order(massRef)]

  # test if a peak is present at each massRef
  spTronc <- lapply(massRef.o,function(x) {
    threshold<- 2000*0.5*x/10^6
    interval <- which( (mz < x + threshold) & (mz > x - threshold) )
    sp[interval] })

  ## test if there is only one peak
  badMass<-which(unlist(lapply(spTronc, function(x) length(LocalMaximaSG(sp = x,minPeakHeight = 0.1*max(x))))!=1))

  if(length(badMass)!=0){
    warning(paste("mass ", massRef.o[badMass] ,"is excluded",collapse=""))
    massRef.o<-massRef.o[-badMass]
  }
  if(length(massRef.o)<2) stop("To few references masses for determine peak shape")

  n.mass <- length(massRef.o)
  interval<- lapply(massRef.o,function(x) {
    th<-0.4
      #min(x*2000/10^6,0.4)
  list(which(mz<= x+ th & mz>=x - th),
       sp[mz<= x+ th & mz>=x - th])})


  #calcul t and delta for  each mass
  delta<- vapply(interval, 
                 function(x) width(tof = x[[1]],peak = x[[2]])$delta,
                 FUN.VALUE = 1.1)
  t_centre <- vapply(interval,function(x) x[[1]][which.max(x[[2]])],
                     FUN.VALUE = interval[[1]][[1]][1] )

  #normalization of intreval
  interval.n<-list()
  for (j in seq_len(n.mass)){
    interval.n[[j]]<-(interval[[j]][[1]]-t_centre[[j]])/delta[j]
  }

  interval.n.length <-vapply(interval.n,length,1)
  interval.n.borne <- vapply(interval.n, range,c(1,2))
  index.inter.longest <- which.max(interval.n.length)
  interval.n.longest <- interval.n[[index.inter.longest]]

  interval.ref.borne <- interval.n.borne[,which.min(apply(abs(interval.n.borne),2,min))]
  interval.ref <- interval.n.longest[ interval.ref.borne[1] < interval.n.longest &
                                       interval.n.longest < interval.ref.borne[2]]

  peak<-matrix(0,ncol=n.mass,nrow=length(interval.ref))

  peak[,index.inter.longest]<-interval[[index.inter.longest]][[2]][ interval.ref.borne[1] < interval.n.longest &
                                                                      interval.n.longest < interval.ref.borne[2]]
  peak[,index.inter.longest]<-peak[,index.inter.longest]/max( peak[,index.inter.longest])

  #interpolation from the ref interval
  for (i in seq_len(n.mass[-which.max(interval.n.length)])){
    interpolation<-stats::spline(interval.n[[i]],interval[[i]][[2]],xout = interval.ref)
    peak[,i]<-interpolation$y/max(interpolation$y)
  }

  #average of cumulated signal
  peak_ref<-rowSums(peak)/n.mass
  
  # if(plotShape){
  #   col=grDevicesrainbow(n.mass)
  #   plot(interval.ref,peak[,1],type="p",pch=20,cex=0.9,col=col[n.mass],xlab="interval normalized",main="Peak shape")
  #   for (i in 2:(n.mass)){
  #     points(interval.ref,peak[,i],pch=20,cex=0.9,col=col[i])
  #   }
  #   graphics::lines(interval.ref,peak_ref,type="l",lwd=2)
  #   graphics::lines(interval.ref,sqrt(2*pi)/2*dnorm(interval.ref,0,1/2),lty=2,lwd=2)
  #   legend('topright',legend=c("Average Peak","Gaussian",massRef.o),col=c("black","black",col),lty=c(1,2,rep(NA,n.mass)),pch=c(NA,NA,rep(19,n.mass)))
  # }

  return(list(tofRef = interval.ref, peakRef = peak_ref))
}

# Peak detection -------------------------------------------------

#'Find optimal window's size for Savitzky Golay filter
#'@param sp the array of specrtum values
#'@param noiseacf autocorrelation of the noise
#'@param d the degree of Savitzky Golay filter
#'@return the optimal size of Savitzky Golay filter's windows
OptimalWindowsSG<-function(sp, noiseacf, d=3){
  n=5
  spf <- signal::sgolayfilt(sp, p=d, n=n)
  res<- sp - spf
  acf_res0 <- stats::acf(res,plot=FALSE)[1]$acf[1]
  n=n+2
  repeat {
    spf <-signal::sgolayfilt(sp, p=d, n=n)
    res<- sp - spf
    acf_res1 <- stats::acf(res,plot=FALSE)[1]$acf[1]
    if( abs(acf_res0 - noiseacf) < abs(acf_res1 - noiseacf ) ) break #si res 0 est plus proche que res1
    n=n+2
    if( n >= length(sp)) break
    acf_res0 <- acf_res1
  }
  n
}

#' Find local maxima with Savitzky Golay filter
#'
#' Apply Savitzky Golay filter to the spectrum and find local maxima such that : 
#' second derivate Savitzky Golay filter < 0 and first derivate = 0 and intensity >  minPeakHeight 
#' 
#'@param sp the array of spectrum values
#'@param minPeakHeight minimum intensity of a peak
#'@param noiseacf autocorrelation of the noise
#'@param d the degree of Savitzky Golay filter, defalut 3
#'@return array with peak's index in the spectrum
#' @examples 
#' spectrum<-dnorm(x=seq(-5,5,length.out = 100))
#' index.max<-LocalMaximaSG(spectrum)
#'@export
LocalMaximaSG<-function(sp, minPeakHeight= -Inf, noiseacf=0.1, d=3){

  if(sum(sp==0)/length(sp) > 0.8) return(NULL) ### write a test

  n<-OptimalWindowsSG(sp,noiseacf,d)
  spf<-signal::sgolayfilt(sp, p=d, n=n)
  spf_derivate1<-signal::sgolayfilt(sp, p=d, n=n,m=1)
  spf_derivate2<-signal::sgolayfilt(sp, p=d, n=n,m=2)

  peak<-NULL

  #concavity (second derivate <0) and threshold
  concave <- which(spf_derivate2 < 0 & sp > minPeakHeight)
  if(length(concave)!=0){
    concave.end<- c(0 , which(diff(concave) > 1), length(concave))
    n.peak<- length(concave.end)-1

    #local maxima (first derivate =0 )
    peak<- rep(0,n.peak) #matrix(0,ncol=5,nrow = n.peak)
    for (j in seq_len(n.peak)){
      group_concave<-concave[(concave.end[j]+1):concave.end[j+1]]
      if(length(group_concave) < 3 ) next #|| max(spf[group])-min(spf[group]) < minPeakHeight ) peak[j] <-0
      else {
        peak.j <- group_concave[which.min(abs(spf_derivate1[group_concave]))]
        peak[j]<- peak.j
      }
    }
    peak<-peak[peak!=0]
  }
  if(length(peak)==0) peak<-NULL
  peak
}

# Baseline Correction --------------------------------------------

baselineEstimation<-function(sp,d=2){
  x<-seq(1,length(sp))
  reg0<-stats::lm(sp ~ x)
  a0<-reg0$coefficients
  yhat<- stats::predict(reg0)
  y <- (sp<yhat)*sp + yhat*(sp>=yhat)

  reg1<-stats::lm(y~stats::poly(x, d) )
  a1<-reg1$coefficients
  yhat<-stats::predict(reg1)
  y<-(sp<yhat)*sp + yhat*(sp>=yhat)

  reg2<-stats::lm(y~ stats::poly(x,d))
  a2<-reg2$coefficients
  c=1
  while(sqrt(sum((a1-a2)^2))>1 & c <500){
    a1<-a2
    yhat<-stats::predict(reg2)
    y<-(sp<yhat)*sp + yhat*(sp>=yhat)
    reg2<-stats::lm(y~stats::poly(x,  d))
    a2<-reg2$coefficients
    c=c+1
  }
  return(yhat)
}

#'Baseline estimation
#'
#'@param sp an array with spectrum values
#'@param widthI width of interval
#'@param iteI number of iteration
#'
#'@return baseline estimation of the spectrum
snipBase <- function(sp, widthI = 15, iteI = 3) {

  runave <- function(x, widI) {
    z <- x
    smoVi <- seq(widI + 1, length(x) - widI)
    z[smoVi] <- (z[smoVi - widI] + z[smoVi + widI]) / 2
    z
  }

  GspeVn <- log(log(sp + 1) + 1)

  for (i in seq_len(iteI)) {

    GsmoVn <- runave(GspeVn, widthI)

    GspeVn <- apply(cbind(GspeVn, GsmoVn), 1, min)

  }

  exp(exp(GspeVn) - 1) - 1
}

## Calcul concentration -----

#' Estimation of the transmisison curve
#' @param x masses 
#' @param y transmission data
#' @return a numeric vector
TransmissionCurve <- function(x,y){
  model<-stats::lm(y~I(x^2)+x+I(x^0.5))
  curve<-stats::predict(model,new_data=data.frame(x=seq(1,utils::tail(y,1))))
  extra<-Hmisc::approxExtrap(x, curve, xout=seq(1,1000), method = "linear", n = 50)
  return(extra)
}

#' PPb concentration
#' 
#' @param peakList the peak list return by peakList, must contain the isotope 
#' of primary ion 21.021
#' @param primaryIon primary ion quantification in count per extraction  
#' @param transmission transmission information
#' @param U drift voltage
#' @param Td drift temperature
#' @param pd drift pressure
#' @param k reaction rate constant 
#' @return a vector with ppb concentration
ppbConvert <- function(peakList, primaryIon,transmission, U, Td, pd, k=2*10^-9 ){
  
    Tr_primary<-transmission[2,1]
    x<-transmission[1,][transmission[1,]!=0]
    y<-transmission[2,][transmission[1,]!=0]
    TrCurve<-TransmissionCurve(x, y)
    TrCurve$y<-vapply(TrCurve$y,function(x) max(0.001,x),0.1)
    
    ppb <- peakList[,"quanti"]*
      10^9*
      mean(U)*
      2.8*
      22400*1013^2*
      (273.15+mean(Td))^2*Tr_primary/
      (k*9.2^2*
         primaryIon*488* 
         mean(pd)^2*
         6022*10^23*
         273.15^2 *
         TrCurve$y[peakList[,"Mz"]]
      )
  
    ppb * 1000
}

## Other ----
#' Calculate the FWHM (Full Width at Half Maximum) in raw data
#'
#' @param tof A vector of tof interval
#' @param peak A vector of peak Intensity
#' @param fracMaxTIC the fraction of the maximum intenisty to compute the width
#' @return the delta FWHM in tof 
width <- function(tof, peak, fracMaxTIC=0.5){
  hm<-max(peak)*fracMaxTIC
  sup1<-findEqualGreaterM(peak[seq_len(which.max(peak))],hm)
  sup2<-unname(which.max(peak))+FindEqualLess(peak[(which.max(peak)+1):length(peak)],hm)
  delta_tof<-unname(vapply(c(sup1,sup2),
                           function(x) { 
                             (hm*(tof[x]-tof[x-1])-
                                (tof[x]*peak[x-1]-tof[x-1]*peak[x]))/
                             (peak[x]-peak[x-1])},
                           FUN.VALUE = 0.1 )) # équation intrepolation linéaire entre sup et sup-1 = hm
  delta <- unname(abs(delta_tof[2] - delta_tof[1]))
  list(delta=delta,delta_tof=delta_tof-1)
}


FindEqualLess<-function(x,values){
  idx<-1
  while(x[idx] > values) idx=idx+1
  idx
}
