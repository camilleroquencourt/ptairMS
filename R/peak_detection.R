utils::globalVariables("::<-")
## detectPeak----
#' Detection and quantification of peaks for a ptrSet object. 
#' 
#' The \code{detectPeak} function detects peaks on the average total spectrum 
#' around nominal masses, for all files present in ptrSet which have not already 
#' been processed. The temporal evolution of each peak is then evaluated by using
#' a two-dimensional penalized spline regression. Finally, 
#' the expiration points (if defined in the ptrSet) are averaged, 
#' and a t-test is performed between expiration and ambient 
#' air. The peakList can be accessed with the \code{\link[ptairMS]{getPeakList}} 
#' function, which returns the information about the detected peaks in each file 
#' as a list of ExpressionSet objects.The peak detection steps within each file 
#' are as follows: \cr
#' for each nominal mass:
#' \itemize{
#' \item correction of the calibration drift
#' \item peak detection on the average spectrum
#' \item estimation of temporal evolution 
#' \item t-test between expiration and ambient air
#' }
#' @param x a \code{\link[ptairMS]{ptrSet}} object 
#' @param mzNominal nominal masses at which peaks will be detected; if \code{NULL}, 
#' all nominal masses of the mass axis
#' @param ppm minimum distance in ppm between two peaks
#' @param resolutionRange vector with the minimum, average, and 
#' maximum resolution of the PTR instrument. If \code{NULL}, the values are 
#' estimated by using the calibration peaks.
#' @param minIntensity minimum intensity for peak detection. The final threshold 
#' for peak detection will be: max(\code{minIntensity}, threshold noise ). 
#' The threshold noise corresponds to
#'  max(max(noise around the nominal mass), \code{minIntensityRate} * 
#'  max(intensity within the nominal mass)
#' @param minIntensityRate Fraction of the maximum 
#' intensity to be used for noise thresholding
#' @param smoothPenalty second order penalty coefficient of the p-spline used 
#' for two-dimensional regression. If \code{NULL}, the coefficient is estimated by 
#' generalized cross validation (GCV) criteria 
#' @param fctFit function for the peak quantification: should be sech2 
#' or averagePeak. If \code{NULL}, the best function is selected by using the 
#' calibration peaks
#' @param parallelize Boolean. If \code{TRUE}, loops over files are parallelized
#' @param nbCores number of cluster to use for parallel computation.
#' @param saving boolean. If TRUE, the object will be saved in saveDir with the
#' \code{setName} parameter of the \code{createPtrSet} function
#' @param saveDir The directory where the ptrSet object will be saved as .RData. 
#' If NULL, nothing will be saved
#' @param ... may be used to pass parameters to the processFileTemporal function
#' @return an S4 ptrSet object, that contains the input ptrSet completed with the 
#' peakLists.
#' @references Muller et al 2014, Holzinger et al 2015, Marx and Eilers 1992
#' @examples 
#' 
#' ## For a ptrSet object
#' library(ptairData)
#' directory <- system.file("extdata/exhaledAir",  package = "ptairData")
#' exhaledPtrset<-createPtrSet(dir=directory,setName="exhaledPtrset",
#' mzCalibRef=c(21.022,59.049),
#' fracMaxTIC=0.9,saveDir= NULL)
#' exhaledPtrset  <- detectPeak(exhaledPtrset)
#' peakListEset<-getPeakList(exhaledPtrset)
#' Biobase::fData(peakListEset[[1]])
#' Biobase::exprs(peakListEset[[1]])
#' @rdname detectPeak
#' @import doParallel foreach parallel
#' @export
setMethod(f = "detectPeak", signature = "ptrSet", 
            function(x, ppm = 130, minIntensity = 10, 
                     minIntensityRate = 0.01, mzNominal = NULL, 
                     resolutionRange = NULL,fctFit = NULL, smoothPenalty = 0, 
                     parallelize = FALSE, nbCores = 2, saving = TRUE, 
                     saveDir = getParameters(x)$saveDir, ...){
                
    ptrset <- x
    
    # get infomration
    knots <- ptairMS:::getPeaksInfo(ptrset)$knots
    massCalib <- getCalibrationInfo(ptrset)$mzCalibRef
    primaryIon <- getPTRInfo(ptrset)$primaryIon
    indTimeLim <- getTimeInfo(ptrset)$timeLimit
    parameter <- getParameters(ptrset)
    dir <- parameter$dir
    peakList <- getPeakList(ptrset)
    peakShape <- getPeaksInfo(ptrset)$peakShape
    paramOld <- parameter$detectPeakParam
    if (methods::is(dir, "expression")) 
        dir <- eval(dir)
    if (is.null(mzNominal)) 
        mzNominalParam <- "NULL" else mzNominalParam <- mzNominal
    if (is.null(resolutionRange)) {
        resolutionEstimated <- Reduce(c, getPTRInfo(ptrset)$resolution)
        resolutionRange <- c(floor(min(resolutionEstimated)/1000) * 1000, 
                                round(mean(resolutionEstimated)/1000) * 
            1000, ceiling(max(resolutionEstimated)/1000) * 1000)
    }
    if (is.null(fctFit)) 
        fctFitParam <- "NULL" else fctFitParam <- fctFit
    if (is.null(smoothPenalty)) 
        smoothPenaltyParam <- "NULL" else smoothPenaltyParam <- smoothPenalty
    
    # files already check by checkSet
    files <- parameter$listFile
    if (methods::is(files, "expression")) 
        files <- eval(files)
    fileName <- basename(files)
    paramNew <- list(mzNominal = mzNominalParam, ppm = ppm, 
                        minIntensityRate = minIntensityRate, 
        minIntensity = minIntensity, fctFit = fctFitParam, 
        smoothPenalty = smoothPenalty, 
        resolutionRange = resolutionRange)
    
    # keep files that not alreday processed
    fileDone <- names(peakList)
    fileToProcess <- which(!(fileName %in% fileDone))
    fileName <- fileName[fileToProcess]
    allFilesName <- basename(files) 
    files <- files[fileToProcess]
    
    #to save the oder
    if (length(files) == 0) {
        message("All files have already been processed")
        return(ptrset)
    }
    
    if (is.null(paramOld)) {
        #ptrSet<-setParameters(ptrset,paramNew)
    } else {
        # we take parameters that already processed the other file
        mzNominal <- paramOld$mzNominal
        if (mzNominal[1] == "NULL") {
            mzNominal <- NULL
        } else paramOld$mzNominal <- paste(range(paramOld$mzNominal), 
                                            collapse = "-")
        ppm <- paramOld$ppm
        minIntensity <- paramOld$minIntensity
        fctFit <- paramOld$fctFit
        minIntensityRate <- paramOld$minIntensityRate
        smoothPenalty <- paramOld$smoothPenalty
        if (smoothPenalty == "NULL") 
            smoothPenalty <- NULL
        resolutionRange <- paramOld$resolutionRange
        paramOldMat <- Reduce(cbind, paramOld)
        colnames(paramOldMat) <- names(paramOld)
        message("the peak list will be calculated with the same parameters 
                      as the other files :\n")
        print(paramOldMat)
        message("\n If you want to change theme, remove the peak list 
                      before with rmPeakList() function and restart 
                      detectPeak() \n")
    }
    
    if (!is.null(fctFit)) {
        fctFit <- as.list(rep(fctFit, length(files)))
        names(fctFit) <- basename(files)
    } else fctFit <- getPeaksInfo(ptrset)$fctFit
    
    FUN <- function(x) {
        test <- try(processFileTemporal(fullNamefile = x, 
                                        massCalib = massCalib[[basename(x)]], 
            primaryIon = primaryIon[[basename(x)]], 
            indTimeLim = indTimeLim[[basename(x)]], 
            mzNominal = mzNominal, ppm = ppm, resolutionRange = resolutionRange, 
            minIntensity = minIntensity, fctFit = fctFit[[basename(x)]], 
            minIntensityRate = minIntensityRate, 
            knots = knots[[basename(x)]], smoothPenalty = smoothPenalty, 
            peakShape = peakShape[[basename(x)]]))
        if (!is.null(attr(test, "condition"))) {
            return(list(raw = NULL, aligned = NULL))
        } else return(test)
    }
    
    # parallel peakLists<-BiocParallel::bplapply(files,FUN = FUN)
    if (parallelize) {
        cl <- parallel::makeCluster(nbCores)
        doParallel::registerDoParallel(cl)
        `%dopar%` <- foreach::`%dopar%`
        peakLists <- foreach::foreach(file = files, .packages = c("data.table")) %dopar% 
            {
                FUN(file)
            }
        parallel::stopCluster(cl)
    } else peakLists <- lapply(files, FUN)
    
    
    # delete processFile failed
    failed <- which(Reduce(c, lapply(peakLists, 
                                     function(x) is.null(x$raw))))
    if (length(failed)) {
        peakLists <- peakLists[-failed]
        warning(basename(files)[failed], "failed")
        files<-files[-failed]
    }
    
    
    # create an expression set for each file
    peakListsEset <- lapply(peakLists, function(x) {
        infoPeak <- grep("parameter", colnames(x$raw))
        colnames(x$raw)[infoPeak] <- c("parameterPeak.delta1", 
                                       "parameterPeak.delta2", 
                                       "parameterPeak.height")
        assayMatrix <- as.matrix(x$raw)[, -c(1, 2, 3, 4, infoPeak), drop = FALSE]
        rownames(assayMatrix) <- round(x$raw$Mz, 4)
        featuresMatrix <- data.frame(cbind((as.matrix(x$raw)[, c(1, 2, 3, 4, infoPeak), 
            drop = FALSE]), (x$aligned[, -1])))
        rownames(featuresMatrix) <- rownames(assayMatrix)
        Biobase::ExpressionSet(assayData = assayMatrix, 
                               featureData = Biobase::AnnotatedDataFrame(featuresMatrix))
    })
    names(peakListsEset) <- basename(files)
    ptrset@peakList <-c(peakListsEset,ptrset@peakList)

    # save
    if (!is.null(saveDir) & saving) {
        if (!dir.exists(saveDir)) {
            warning("saveDir does not exist, object not saved")
            return(ptrset)
        }
        changeName <- parse(text = paste0(getParameters(ptrset)$name, "<- ptrset "))
        eval(changeName)
        eval(parse(text = paste0("save(", getParameters(ptrset)$name,
                                    ",file = paste0(saveDir,'/', '", 
                                 getParameters(ptrset)$name, ".RData '))")))
    }
    return(ptrset)
})

processFileTemporal <- function(fullNamefile, massCalib, 
                                primaryIon, indTimeLim, 
                                mzNominal, ppm, resolutionRange, 
                                minIntensity, fctFit, minIntensityRate, 
                                knots, peakShape, smoothPenalty = 0, 
                                funAggreg = mean, ...) {
    if (is.character(fullNamefile)) {
        cat(basename(fullNamefile), ": ")
        # read file
        raw <- readRaw(fullNamefile, mzCalibRef = massCalib)
    } else raw <- fullNamefile
    # processing for each masses
    if (is.null(mzNominal)) 
        mzNominal = unique(round(getRawInfo(raw)$mz))
    if (fctFit == "averagePeak") {
        l.shape <- peakShape
        raw<-setPeakShape(raw,peakShape)
    } else l.shape = list(NULL)
    

    # p<-list()
    # for(m in mzNominal){
    # 
    #     p[[m+1]]<-ptairMS:::processFileTemporalNominalMass(m = m,
    #                                                          raw = raw, mzNominal = mzNominal, ppm = ppm,
    #                                                          resolutionRange = resolutionRange,
    #                                                          minIntensity = minIntensity, fctFit = fctFit,
    #                                                          minIntensityRate = minIntensityRate,
    #                                                          knots = knots, smoothPenalty = smoothPenalty, l.shape = l.shape,
    #                                                          timeLimit = indTimeLim)
    #     print(m)
    # }

    process <- lapply(mzNominal, function(m) ptairMS:::processFileTemporalNominalMass(m = m, 
        raw = raw, mzNominal = mzNominal, ppm = ppm, 
        resolutionRange = resolutionRange, 
        minIntensity = minIntensity, fctFit = fctFit, 
        minIntensityRate = minIntensityRate, 
        knots = knots, smoothPenalty = smoothPenalty, l.shape = l.shape, 
        timeLimit = indTimeLim))
    
    
    matPeak <- Reduce(rbind, process)
    ## agregate
    matPeakAg <- aggregateTemporalFile(time = getRawInfo(raw)$time, indTimeLim = indTimeLim, 
        matPeak = matPeak, funAggreg = funAggreg)
    indLim <- indTimeLim$exp
    indBg <- indTimeLim$backGround
    bg <- FALSE
    if (!is.null(indBg)) 
        bg <- TRUE
    # normalize by primary ions and ppb conversion
    if (!is.na(primaryIon$primaryIon)) {
        matPeakAg[, "quanti_ncps"] <- matPeakAg[, "quanti_cps"]/(
            (primaryIon$primaryIon * 
            488))
        if (bg) {
            matPeakAg[, "background_ncps"] <- matPeakAg[, "background_cps"]/
                ((primaryIon$primaryIon * 
                488))
            matPeakAg[, "diffAbs_ncps"] <- matPeakAg[, "diffAbs_cps"]/
                ((primaryIon$primaryIon *488))
        }
        
        indExp <- Reduce(c, apply(indLim, 2, function(x) seq(x[1], x[2])))
        
        if (length(getPTRInfo(raw)$prtReaction) != 0 & nrow(getPTRInfo(raw)$ptrTransmisison) > 1) {
            matPeakAg[, "quanti_ppb"] <- ppbConvert(peakList = data.frame(Mz = matPeakAg$Mz,
                                                                          quanti = matPeakAg$quanti_ncps), 
                                                    transmission = getPTRInfo(raw)$ptrTransmisison, 
                                                    U = c(getPTRInfo(raw)$prtReaction$TwData[1, ])[indExp], 
                                                    Td = c(getPTRInfo(raw)$prtReaction$TwData[3,])[indExp], 
                                                    pd = c(getPTRInfo(raw)$prtReaction$TwData[2,  ])[indExp])
            if (bg) {
                matPeakAg[, "background_ppb"] <- ppbConvert(peakList = data.frame(Mz = matPeakAg$Mz,
                                                                                  quanti = matPeakAg$background_ncps), 
                                                            transmission = getPTRInfo(raw)$ptrTransmisison, 
                                                            U = c(getPTRInfo(raw)$prtReaction$TwData[1, ])[indBg], 
                                                            Td = c(getPTRInfo(raw)$prtReaction$TwData[3,])[indBg], 
                                                            pd = c(getPTRInfo(raw)$prtReaction$TwData[2, ])[indBg])
                matPeakAg[, "diffAbs_ppb"] <- ppbConvert(peakList = data.frame(Mz = matPeakAg$Mz,
                                                                               quanti = matPeakAg$diffAbs_ncps), 
                                                         transmission = getPTRInfo(raw)$ptrTransmisison, 
                                                         U = c(getPTRInfo(raw)$prtReaction$TwData[1, ])[indBg], 
                                                         Td = c(getPTRInfo(raw)$prtReaction$TwData[3, ])[indBg], 
                                                         pd = c(getPTRInfo(raw)$prtReaction$TwData[2, ])[indBg])
            }
        }
    }
    cat(paste(nrow(matPeakAg), "peaks detected \n"))
    # ordered column
    # matPeakAg<-matPeakAg[,c('Mz','quanti_cps','background_cps','diffAbs_cps',
    # 'quanti_ncps','background_ncps','diffAbs_ncps',
    # 'quanti_ppb','background_ppb','diffAbs_ppb', 'pValGreater','pValLess')]
    return(list(raw = matPeak, aligned = matPeakAg))
}

processFileTemporalNominalMass <- function(m, raw, mzNominal, 
    ppm, resolutionRange, 
    minIntensity, fctFit, minIntensityRate, knots, smoothPenalty, l.shape, 
    timeLimit, 
    ...) {
    # select raw data around the nominal mass
    mz <- getRawInfo(raw)$mz
    time <- getRawInfo(raw)$time
    rawM <- getRawInfo(raw)$rawM
    rawSub <- rawM[mz > m - 0.6 & mz < m + 0.6, ]
    # estimate the calibration shift
    calib_List <- getCalibrationInfo(raw)$calibCoef
    indexTimeCalib <- getCalibrationInfo(raw)$indexTimeCalib
    if (length(calib_List) > 1) {
        shiftm <- rep(0, length(calib_List))
        for (k in seq_along(calib_List)) {
            mz.shift <- tofToMz(mzToTof(m, calib_List[[k]]), calib_List[[1]])
            shiftm[k] <- mz.shift - m
        }
        # correct the shift
        rawMCorr <- matrix(0, ncol = ncol(rawSub), nrow = nrow(rawSub))
        for (i in seq_along(calib_List)) {
            index <- indexTimeCalib[[i]]
            rawMCorr[, index] <- vapply(index, function(j) {
                stats::approx(x = as.numeric(rownames(rawSub)), y = rawSub[, j], 
                                xout = as.numeric(rownames(rawSub),ties=min) + 
                  shiftm[i])$y
            }, FUN.VALUE = rep(0, nrow(rawSub)))
        }
    } else rawMCorr <- rawSub
    
    rownames(rawMCorr) <- rownames(rawSub)
    colnames(rawMCorr) <- colnames(rawSub)
    rawMCorr <- rawMCorr[!apply(rawMCorr, 1, function(x) any(is.na(x))), ]
    # peak detection on the average spectrum
    sp <- rowSums(rawMCorr)/ncol(rawMCorr)
    mz <- as.numeric(rownames(rawMCorr))
    PeakListm <- peakListNominalMass(i = m, mz = mz, sp = sp, 
                                        calibCoef = getCalibrationInfo(raw)$calibCoef[[1]], 
        ppmPeakMinSep = ppm, resolutionRange = resolutionRange, 
        minPeakDetect = minIntensity, 
        fitFunc = fctFit, minIntensityRate = minIntensityRate, l.shape = l.shape)
    peak <- PeakListm$peak
    baseline<- PeakListm$baseline[[1]]
    if (is.null(peak)) 
        return(NULL) else {
        # 2d deconvolution
        if (is.null(knots)) 
            deconvMethod <- deconv2d2linearIndependant else deconvMethod <- deconv2dLinearCoupled
        fileProccess <- computeTemporalFile(raw = raw, peak = peak, 
                                            baseline = baseline, 
                                            deconvMethod = deconvMethod, 
                                            fctFit = fctFit, knots = knots, 
                                            smoothPenalty = smoothPenalty,
                                            timeLimit = timeLimit)
        matPeak <- data.table::as.data.table(fileProccess$matPeak)
        # convert in cps
        matPeak[, (ncol(matPeak) - length(getRawInfo(raw)$time) + 1):ncol(matPeak)] <- matPeak[, 
            (ncol(matPeak) - length(getRawInfo(raw)$time) + 1):ncol(matPeak)]/diff(getRawInfo(raw)$time)[2]
        # change names of quanti colnames(matPeak)[5:(ncol(matPeak)-2)] <-
        # paste('quanti_cps', colnames(matPeak)[5:(ncol(matPeak)-2)] , sep = ' - ')
        return(matPeak)
    }
}
### peak detection for 1D sptectrum ----
peakListNominalMass <- function(i, mz, sp, ppmPeakMinSep = 130, calibCoef, resolutionRange = c(300, 
    5000, 8000), minPeakDetect = 10, fitFunc = "sech2", maxIter = 4, R2min = 0.998, 
    autocorNoiseMax = 0.3, plotFinal = FALSE, plotAll = FALSE, thNoiseRate = 1.1, 
    minIntensityRate = 0.01, countFacFWHM = 10, daSeparation = 0.001, d = 3, windowSize = 0.4, 
    l.shape = NULL, blCor = TRUE) {
    emptyData <- data.frame(Mz = double(), quanti = double(), delta_mz = double(), 
        resolution = double())
    warning_mat <- NULL
    infoPlot <- list()
    baseline <- list()
    no_peak_return <- list(emptyData, warning_mat, infoPlot, baseline)
    # select spectrum around the nominal mass i
    index.large <- which(mz < i + windowSize + 0.2 & mz > i - windowSize - 0.2)
    
    mz.i.large <- mz[index.large]
    sp.i.large <- sp[index.large]
    if (range(mz.i.large)[1] > i - windowSize | range(mz.i.large)[2] < i + windowSize) 
        return(no_peak_return)
    # baseline correction
    if (blCor) {
        bl <- snipBase(sp.i.large)
        sp.i.large.corrected <- sp.i.large - bl
    } else {
        sp.i.large.corrected <- sp.i.large
        bl <- NA
    }
    baseline[[as.character(i)]] <- bl
    # noise auto correlation estimation
    index.noise <- c(which(i - windowSize > mz & mz > i - windowSize - 0.2), which(i + 
        windowSize < mz & mz < i + windowSize + 0.2))
    noise <- sp.i.large.corrected[match(index.noise, index.large)]
    real_acf <- stats::acf(noise, lag.max = 1, plot = FALSE)[1]$acf
    if (is.na(real_acf)) 
        return(no_peak_return)
    noiseacf <- min(real_acf, autocorNoiseMax)  # max auto correlation at 0.3 
    # indeed OptimnalWinfowSg will overfitted
    thr <- max(noise)  # dynamic threshold for peak detetcion
    # signal who wil be process
    if (i > 250) {
        mz.i <- mz.i.large
        sp.i <- sp.i.large.corrected
    } else {
        index <- which(i - windowSize <= mz & mz <= i + windowSize)
        mz.i <- mz[index]
        sp.i <- sp.i.large.corrected[match(index, index.large)]
    }
    infoPlot[[as.character(i)]] <- list(mz = mz.i, sp = sp.i, main = i, pointsPeak = NA)
    no_peak_return <- list(emptyData, warning_mat, infoPlot,baseline)
    ## fit itterative
    sp.i.fit <- sp.i  # initialization of sp.i.fit
    c = 1  # initialize number of itteration
    repeat {
        # Initialization for regression :
        if (c == 1) {
            minpeakheight <- max(max(thr * thNoiseRate, minIntensityRate * max(sp.i), 
                minPeakDetect))
        } else minpeakheight <- max(minpeakheight * 0.8, 0.1)  # minimum intenisty
        # Find local maximum with Savitzky golay filter or wavelet
        init <- initializeFit(i, sp.i.fit, sp.i, mz.i, calibCoef, 
                              resmean = resolutionRange[2], 
                              minpeakheight, noiseacf, ppmPeakMinSep, 
                              daSeparation, d, plotAll, c)   
         
        # Regression :
        if (!is.null(init)) {
            mz_init <- init$mz[, "m"]
            if (c > 1) 
                {
                  mz_par <- par_estimated[1, ]
                  # delete peak also find
                  to_delete <- NULL
                  for (p in seq_along(mz_init)) {
                    Da_proxi <- abs(mz_par - mz_init[p])
                    if (i > 17) 
                      test_proxi <- Da_proxi * 10^6/i > ppmPeakMinSep
                    if (i <= 17) 
                      test_proxi <- Da_proxi > daSeparation
                    if (!all(test_proxi)) 
                      to_delete <- c(to_delete, p)
                  }
                  n.peak <- nrow(init$mz)
                  # if same number of peak detected and all old peak deleted,
                  #  no new peak are detected
                  if (n.peak == dim(par_estimated)[2] & length(to_delete) == dim(par_estimated)[2]) 
                    break
                  # add new peak
                  if (!is.null(to_delete)) {
                      init$mz <- init$mz[ -to_delete, ,drop = FALSE]
                  }
                }  # END if c> 1
            n.peak <- nrow(init$mz)
            if (c == 1) 
                initMz <- init$mz else initMz <- rbind(init$mz, t(par_estimated))
            n.peak <- nrow(initMz)
            resolution_upper <- resolutionRange[3]
            resolution_mean <- resolutionRange[2]
            resolution_lower <- resolutionRange[1]
            lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak), 
                                                rep(0, n.peak * 2),
                                                rep(0, n.peak)), ncol = 4) - 
                                  matrix(c(initMz[, "m"]/(resolution_mean * 4), 
                                           -initMz[, "m"]/(resolution_upper * 2), 
                                           -initMz[, "m"]/(resolution_upper *2), 
                                           rep(0, n.peak)), ncol = 4)))
            upper.cons <- c(t(initMz * matrix(c(rep(1, n.peak), 
                                                rep(0, n.peak * 2),
                                                rep(Inf, n.peak)), ncol = 4) + 
                                  matrix(c(initMz[, "m"]/(resolution_mean * 4), 
                                           initMz[, "m"]/(resolution_lower * 2), 
                                           initMz[, "m"]/(resolution_lower * 2), 
                                           rep(0, n.peak)), ncol = 4)))
            fit <- fitPeak(initMz = initMz, sp = sp.i, mz.i = mz.i, lower.cons, 
                           upper.cons, funcName = fitFunc, l.shape)
            
            if (is.na(fit$fit$deviance)) {
                if (c == 1) 
                  return(no_peak_return) else break
            }
            
            if (fit$fit$niter == 50) 
                warning_mat <- rbind(warning_mat, c(i, NA, "max itter lm algo"))
            
            fit.peak <- fit$fit.peak
            par_estimated <- fit$par_estimated
            cum_function.fit.peak <- fit$function.fit.peak
            
            # residual
            sp.i.fit <- sp.i - fit.peak
            
            # indicators
            R2 <- 1 - sum(sp.i.fit^2)/sum((sp.i - mean(sp.i))^2)
            auto_cor_res <- abs(stats::acf(sp.i.fit, plot = FALSE)[1]$acf[1])
            
            if (plotAll) {
                plot(mz.i, sp.i, main = paste(i, "fit itteration :", c), 
                     xlab = "mz", ylab = "intensity", type = "b", pch = 19, 
                     cex = 0.7)
                graphics::lines(mz.i, fit.peak, col = "blue", lwd = 2)
                graphics::lines(mz.i, sp.i.fit, col = "green3", lwd = 2)
                
                for (k in seq_len(ncol(par_estimated))) {
                  graphics::lines(mz.i, eval(parse(text = fitFunc))(par_estimated[1, 
                    k], par_estimated[2, k], par_estimated[3, k], par_estimated[4, 
                    k], mz.i, l.shape), lwd = 2, col = "red", lty = 2)
                }
                
                graphics::points(par_estimated[1, ], cum_function.fit.peak(fit$fit$par, 
                  par_estimated[1, ], rep(0, length(par_estimated[2, ]))), cex = 2, 
                  col = "red", lwd = 2)
                
                graphics::legend("topleft", legend = c("Raw", "fit sum", "fit peak", 
                  "residual", paste("R2=", round(R2, 3)), paste("autocor res=", 
                                                                round(auto_cor_res, 
                    3))), col = c("black", "blue", "red", "green3"), lty = c(NA, 
                  1, 2, 1, NA, NA), pch = c(19, NA, NA, NA, NA, NA))
            }
        } else {
            ## if init is null : no peak find
            if (c == 1) {
                # if first iteration
                if (plotFinal) {
                  plot(mz.i, sp.i, main = paste(i, c, "fit", sep = " - "), 
                       xlab = "mz", ylab = "intensity", type = "p", pch = 19, 
                       cex = 0.7)
                }
                return(no_peak_return)
            } else break
        }
        # condition for break
        if (auto_cor_res < autocorNoiseMax) 
            break
        if (R2 > R2min) 
            break
        if (c > 1) {
            # degradation of R2
            if (R2 < old_R2) {
                fit <- fit_old
                fit.peak <- fit$fit.peak
                par_estimated <- fit$par_estimated
                cum_function.fit.peak <- fit$function.fit.peak
                sp.i.fit <- sp.i - fit.peak
                n.peak <- dim(par_estimated)[2]
                break
            }
        }
        if (n.peak > 10) 
            break
        # Iteration limit
        if (c == maxIter) {
            warning_mat <- rbind(warning_mat, c(i, NA, "residual iteration max"))
            break
        }
        c = c + 1
        fit_old <- fit
        old_R2 <- R2
    }  # END OF REPEAT
    # last fit less constrained
    c = 1
    repeat {
        # until no peak to close or inside an other peak
        init <- t(par_estimated)
        n.peak <- nrow(init)
        lower.cons <- c(t(matrix(c(rep(range(mz.i)[1], n.peak), rep(init[, 1]/
                                                           (resolution_upper * 
            2), 2), rep(0, n.peak)), ncol = 4)))
        upper.cons <- c(t(matrix(c(rep(range(mz.i)[2], n.peak), rep(init[, 1]/(resolution_lower * 
            2), 2), rep(Inf, n.peak)), ncol = 4)))
        fit <- fitPeak(init, sp.i, mz.i, lower.cons, upper.cons, fitFunc, l.shape)
        fit.peak <- fit$fit.peak
        par_estimated <- fit$par_estimated
        delta_mz <- par_estimated[2, ] + par_estimated[3, ]
        quanti <- apply(par_estimated, 2, function(x) {
            th <- countFacFWHM * 0.5 * (x[2] + x[3])
            mz.x <- mz[x[1] - th < mz & mz < x[1] + th]
            sum(eval(parse(text = fitFunc))(x[1], x[2], x[3], x[4], mz.x, l.shape))
        })
        center_peak <- par_estimated[1, ]
        peaks <- apply(par_estimated, 2, function(x) eval(parse(text = fitFunc))(x[1], 
            x[2], x[3], x[4], mz.i, l.shape))
        X <- data.frame(Mz = center_peak, quanti_cps = quanti, delta_mz = delta_mz, 
            resolution = center_peak/delta_mz, parameter = t(par_estimated[-1, ]))
        
        X <- X[quanti > 0, , drop = FALSE]
        par_estimated <- par_estimated[, quanti > 0, drop = FALSE]
        peaks <- peaks[, quanti > 0, drop = FALSE]
        peaks<-peaks[,order(X[, 1]),drop=FALSE]
        X <- X[order(X[, 1]), , drop = FALSE]
        par_estimated <- par_estimated[, order(par_estimated[1, ]), drop = FALSE]
        # R2
        exprs <- parse(text = paste0(fitFunc, "Inv"))
        borne <- apply(par_estimated, 2, function(x) eval(exprs)(x[1], x[2], x[3], 
            x[4], x[4] * 0.02, l.shape))
        borne <- cbind(par_estimated[1, ], t(borne))
        colnames(borne) <- c("Mz", "lowerMz", "upperMz")
        
        borne <- borne[order(borne[, "Mz"]), , drop = FALSE]
        borne <- overlapDedect(borne)
        R2glob <- 1 - sum(sp.i.fit^2)/sum((sp.i - mean(sp.i))^2)
        R2 <- apply(borne, 1, function(x) {
            1 - sum(sp.i.fit[mz.i > x["lowerMz"] & mz.i < x["upperMz"]]^2)/sum((sp.i[mz.i > 
                x["lowerMz"] & mz.i < x["upperMz"]] - mean(sp.i[mz.i > x["lowerMz"] & 
                mz.i < x["upperMz"]]))^2)
        })
        X$R2 <- R2
        X$overlap <- borne[, "overlap"]
        if (dim(X)[1] == 1) 
            break
        to_delete <- NULL
        # tets if peak is under an other peak
        # for (j in seq_len(nrow(X) - 1)) {
        #     sign_diff <- levels(factor(sign(peaks[, j + 1] - peaks[, j])))
        #     sign_diff <- sign_diff[sign_diff != 0]
        #     if (length(sign_diff) < 2) {
        #         if ("-1" %in% sign_diff) 
        #           to_delete <- c(to_delete, j + 1) else to_delete <- c(to_delete, j)
        #     }else if( abs( mean((peaks[, j + 1] - peaks[, j])[peaks[, j + 1] - peaks[, j] <0])) < 0.1 ){
        #         to_delete <- c(to_delete, j) # if the abs of mean values where peak j+1 < peaks j is very small
        #     } else if( abs( mean((peaks[, j + 1] - peaks[, j])[peaks[, j + 1] - peaks[, j] >0])) < 0.1){
        #         to_delete <- c(to_delete, j + 1) # if the abs of mean values where peak j < peaks j+1 is very small
        #     }
        # }
        
        # if no
        if (is.null(to_delete)) {
            # test proximity
            control_proxi <- c(FALSE, diff(X[, 1]) * 10^6/i < ppmPeakMinSep)
            # if no break
            if (all(!control_proxi)) 
                break
            # else delete
            par_estimated <- par_estimated[, -which(control_proxi), drop = FALSE]
        } else {
            par_estimated <- par_estimated[, -to_delete, drop = FALSE]
            # if more than one peak
            if (dim(X)[1] > 1) {
                # test proximity
                control_proxi <- c(FALSE, diff(X[, 1]) * 10^6/i < ppmPeakMinSep)
                if (any(control_proxi)) 
                  par_estimated <- par_estimated[, -which(control_proxi), 
                                                 drop = FALSE]
            }
        }
        c = c + 1
    }  # end second repeat
    pointsPeak <- list(x = par_estimated[1, ], 
                       y = fit$function.fit.peak(fit$fit$par, 
                                                 par_estimated[1, ], 
                                                 rep(0, length(par_estimated[2, ])
                                                     )))
    infoPlot[[as.character(i)]] <- list(mz = mz.i, sp = sp.i, 
                                        main = i, fitPeak = fit.peak,
                                        peak = peaks, pointsPeak = pointsPeak)
    if (plotFinal) {
        plot(mz.i, sp.i, main = i, xlab = "mz", ylab = "intensity", ylim = c(min(sp.i, 
            fit.peak), max(sp.i, fit.peak)), type = "b", pch = 19, cex = 0.7)
        graphics::lines(mz.i, fit.peak, lwd = 2, col = "blue")
        for (k in seq_len(ncol(par_estimated))) graphics::lines(mz.i, peaks[, k], 
            lwd = 2, col = "red", lty = 2)
        graphics::points(pointsPeak, cex = 2, col = "red", lwd = 2, pch = 19)
        graphics::legend("topleft", legend = c("Raw", "fit sum", "fit peak", "peak center", 
            paste("R2=", round(R2glob, 3)), paste("autocor res=", round(auto_cor_res, 
                3))), col = c("black", "blue", "red", "red"), lty = c(NA, 1, 2, NA, 
            NA, NA), pch = c(19, NA, NA, 19, NA, NA))
    }
    if (!is.null(warning_mat)) {
        warning_mat <- data.frame(warning_mat)
        names(warning_mat) <- c("m/z", "peak", "type")
    }
    return(list(peak = X, warning = warning_mat, plot = infoPlot, 
                baseline = baseline))
}
#' initialization for apply fit function in the spectrum
#' @param i the nominal mass
#' @param sp.i.fit the vector who will be fetted (spectrum pf residual)
#' @param sp.i the spectrum around a nominal mass
#' @param mz.i the mass vector around a nominal mass
#' @param calibCoef calibration coeficient
#' @param resmean resolution m/delta(m) mean
#' @param minpeakheight the minimum peak intensity
#' @param noiseacf aytocorelation of the noise
#' @param ppmPeakMinSep the minimum distance between two peeks in ppm 
#' @param daSeparation the minimum distance between two peeks in da
#' @param d the degree of savitzky golay filter
#' @param plotAll bollean if TRUE, it plot all the initialiation step
#' @param c the number of current itteration
#' @return a list with fit input
#' @keywords internal
initializeFit <- function(i, sp.i.fit, sp.i, mz.i, calibCoef, resmean, 
                          minpeakheight,noiseacf, ppmPeakMinSep, daSeparation, 
                          d, plotAll, c) {
    init <- NULL
    # find local maxima in the spectrum: return the index
    prePeak <- LocalMaximaSG(sp = sp.i.fit, minPeakHeight = minpeakheight, 
                             noiseacf = noiseacf, d = d)
    if (!is.null(prePeak)) 
        {
            # get the mass and the intenisty
            prePeak <- cbind(mz = mz.i[prePeak], 
                             intensity = vapply(prePeak, 
                                                function(x) max(sp.i[(x - 1):(x + 1)]), FUN.VALUE = 1.1))  #the maximum auround the index find
            # delete peak to close
            if (nrow(prePeak) > 1) {
                # chexk proximity
                Da_proxi <- diff(prePeak[, "mz"])
                if (i > 17) 
                  test_proxi <- c(TRUE, Da_proxi * 10^6/i > ppmPeakMinSep)
                if (i <= 17) 
                  test_proxi <- c(TRUE, Da_proxi > daSeparation)  ## max 
                prePeak <- prePeak[test_proxi, , drop = FALSE]
            }
            # plot
            if (plotAll) {
                plot(mz.i, sp.i, main = paste("initialization iteration :", c), 
                     xlab = "mz",ylab = "intensity", type = "b", pch = 19, 
                     cex = 0.7)
                graphics::points(prePeak[, "mz"], prePeak[, "intensity"],
                                 col = "red", cex = 2, lwd = 2.5)
                graphics::abline(h = minpeakheight)
                graphics::legend("topleft", legend = c("Raw", "Local maximum"),
                                 col = c("black", 
                  "red"), pch = c(19, 1), lty = c(1, NA))
            }
            # Calculate the initialization
            # in tof fot average fit function
            resolution_mean <- 10000  #in tof : t/delta(t)
            t0 <- unname(mzToTof(prePeak[, "mz"], calibCoef))  # tof peak center
            delta0 <- t0/resolution_mean  # FWHM in tof
            h0 <- unname(prePeak[, "intensity"])
            initTof <- cbind(t = t0, delta = delta0, h = h0)
            # in mz for sech 2 function
            resolution_mean <- resmean  #in mass m / dela(m)
            m0 <- unname(prePeak[, "mz"])  # peak center
            delta0 <- m0/resolution_mean  #lambda estimation for Sec2 function
            h0 <- unname(prePeak[, "intensity"])  #peak height
            initMz <- cbind(m = m0, delta = delta0/2, delta2 = delta0/2, h = h0)
            init <- list(mz = initMz, tof = initTof)
        }  #END if prePrek not null
    return(init)
}
#'Find optimal window's size for Savitzky Golay filter
#'@param sp the array of specrtum values
#'@param noiseacf autocorrelation of the noise
#'@param d the degree of Savitzky Golay filter
#'@return the optimal size of Savitzky Golay filter's windows
#'@keywords internal
OptimalWindowsSG <- function(sp, noiseacf, d = 3) {
    n = 5
    spf <- signal::sgolayfilt(sp, p = d, n = n)
    res <- sp - spf
    acf_res0 <- stats::acf(res, plot = FALSE)[1]$acf[1]
    n = n + 2
    repeat {
        spf <- signal::sgolayfilt(sp, p = d, n = n)
        res <- sp - spf
        acf_res1 <- stats::acf(res, plot = FALSE)[1]$acf[1]
        if (abs(acf_res0 - noiseacf) < abs(acf_res1 - noiseacf)) 
            break
        # si res 0 est plus proche que res1
        n = n + 2
        if (n >= length(sp)) 
            break
        acf_res0 <- acf_res1
    }
    n
}
#' Find local maxima with Savitzky Golay filter
#'
#' Apply Savitzky Golay filter to the spectrum and find local maxima such that : 
#' second derivate Savitzky Golay filter < 0 and first derivate = 0 and 
#' intensity >  minPeakHeight 
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
LocalMaximaSG <- function(sp, minPeakHeight = -Inf, noiseacf = 0.1, d = 3) {
    if (sum(sp == 0)/length(sp) > 0.8) 
        return(NULL)  ### write a test
    n <- OptimalWindowsSG(sp, noiseacf, d)
    spf <- signal::sgolayfilt(sp, p = d, n = n)
    spf_derivate1 <- signal::sgolayfilt(sp, p = d, n = n, m = 1)
    spf_derivate2 <- signal::sgolayfilt(sp, p = d, n = n, m = 2)
    peak <- NULL
    # concavity (second derivate <0) and threshold
    concave <- which(spf_derivate2 < 0 & sp > minPeakHeight)
    if (length(concave) != 0) {
        concave.end <- c(0, which(diff(concave) > 1), length(concave))
        n.peak <- length(concave.end) - 1
        # local maxima (first derivate =0 )
        peak <- rep(0, n.peak)  #matrix(0,ncol=5,nrow = n.peak)
        for (j in seq_len(n.peak)) {
            group_concave <- concave[(concave.end[j] + 1):concave.end[j + 1]]
            if (length(group_concave) < 3) {
                next  #|| max(spf[group])-min(spf[group]) < #minPeakHeight ) peak[j] <-0
            } else {
                peak.j <- group_concave[which.min(abs(spf_derivate1[group_concave]))]
                peak[j] <- peak.j
            }
        }
        peak <- peak[peak != 0]
    }
    if (length(peak) == 0) 
        peak <- NULL
    peak
}


### two dimensional modelization -----
# quantity over time of VOC in a pre-defined peak list
computeTemporalFile <- function(raw, peak, baseline, 
                                deconvMethod = deconv2d2linearIndependant, 
    knots = NULL, smoothPenalty = 0, fctFit = "sech2", timeLimit = NULL) {
    # compute EIC
    EICex <- extractEIC(raw = raw, peak = peak, peakQuantil = 0.01, 
                        fctFit = fctFit)
    borne <- EICex$borne
    EIC <- EICex$EIC
    individualBorne <- EICex$individualBorne
    nbPeakOvelap <- sum(borne[, "overlap"])
    infoPeak <- peak[, grep("parameter", colnames(peak))]
    borne <- cbind(borne, infoPeak, matrix(0, nrow = nrow(borne), 
                                           ncol = length(getRawInfo(raw)$time)))
    colnames(borne)[(5 + ncol(infoPeak)):ncol(borne)] <- getRawInfo(raw)$time
    if (nrow(borne) > 1) 
        borneUnique <- borne[!duplicated(data.frame(borne[, 2:3])), , drop = FALSE] else borneUnique <- borne
    XICdeconv <- matrix(0, ncol = length(getRawInfo(raw)$time), nrow = nrow(peak))
    c <- 1
    for (j in seq_along(EIC)) {
        mz <- borne[which(borne[, "lowerMz"] == borneUnique[j, "lowerMz"]), "Mz"]
        peak.detect <- peak[peak[, "Mz"] %in% mz, , drop = FALSE]
        # baseline correction (plan)
        rawM <- EIC[[j]]
        borneMz <- range(as.numeric(rownames(rawM)))
        bl1D <- baseline
        bl1D <- bl1D[as.numeric(names(bl1D)) >= borneMz[1] & as.numeric(names(bl1D)) <=
             borneMz[2]]
        BL <- matrix(rep(bl1D,ncol(rawM)), ncol = ncol(rawM), nrow = nrow(rawM))
        rawM <- rawM - BL
        rawM[rawM < 0] <- 0
        # find best smooth param
        if (!is.null(knots) & is.null(smoothPenalty)) {
            smoothPenalty <- GCV(rawM = rawM, knots = knots, t = getRawInfo(raw)$time, 
                                 timeLimit = timeLimit)
        }
        deconv2 <- deconvMethod(rawM = rawM, t = getRawInfo(raw)$time, peak.detect = peak.detect, 
            raw = raw, fctFit = fctFit, knots = knots, smoothPenalty = smoothPenalty)
        for (i in seq_len(nrow(peak.detect))) {
            XICdeconv[c, ] <- deconv2$predPeak[, i]
            c <- c + 1
        }
    }
    borne[, (5 + ncol(infoPeak)):ncol(borne)] <- XICdeconv
    borne[, c("lowerMz", "upperMz")] <- individualBorne[, c("lowerMz", "upperMz")]
    return(list(matPeak = borne, EIClist = EIC))
}

#'extract all raw EIC from a pre-definied peak List
#'@param raw ptrRaw object
#'@param peak a data.frame with a column named 'Mz'. The Mz of the VOC detected
#'@param peakQuantil the quantile of the peak shape to determine the borne of 
#'the EIC
#'@param fctFit function used to fit peak
#'@return list containing all EIC and the mz borne for all peak
#'@keywords internal
extractEIC <- function(raw, peak, peakQuantil = 0.01, fctFit = "sech2") {
    # borne integration
    individualBorne <- apply(peak, 1, function(x) eval(parse(text = paste0(fctFit, 
        "Inv")))(x["Mz"], x["parameter.1"], x["parameter.2"], x["parameter.3"], 
                 x["parameter.3"] * 
        peakQuantil, getPeaksInfo(raw)$peakShape))
    individualBorne <- cbind(peak$Mz, t(individualBorne))
    colnames(individualBorne) <- c("Mz", "lowerMz", "upperMz")
    individualBorne <- individualBorne[order(individualBorne[, "Mz"]), , 
                                       drop = FALSE]
    # overlap detection and fusion
    borne <- overlapDedect(individualBorne)
    if (nrow(borne) > 1) {
        borneUnique <- borne[!duplicated(data.frame(borne[, 2:3])), , 
                             drop = FALSE]
    } else borneUnique <- borne
    # extract EIC
    EIC <- list(NULL)
    for (j in seq_len(nrow(borneUnique))) {
        rawInfo<-getRawInfo(raw)
        EIC[[j]] <- rawInfo$rawM[rawInfo$mz > borneUnique[j, "lowerMz"] & 
                                    rawInfo$mz < borneUnique[j, "upperMz"], ]
    }
    return(list(EIC = EIC, borne = borne, individualBorne = individualBorne))
}

overlapDedect <- function(borne) {
    overlap <- list(NULL)
    c <- 0
    o <- 1
    for (i in seq_len(nrow(borne) - 1)) {
        if (borne[i + 1, "lowerMz"] < borne[i, "upperMz"]) {
            # save the overlap
            o <- c(o, i + 1)
        } else {
            if (length(o) > 1) {
                # change
                borne[o, "lowerMz"] <- borne[o[1], "lowerMz"]
                borne[o, "upperMz"] <- borne[utils::tail(o, 1), "upperMz"]
                c <- c + 1
                overlap[[c]] <- o
            }
            o <- i + 1
        }
    }
    if (length(o) > 1) {
        # change
        borne[o, "lowerMz"] <- borne[o[1], "lowerMz"]
        borne[o, "upperMz"] <- borne[utils::tail(o, 1), "upperMz"]
        c <- c + 1
        overlap[[c]] <- o
    }
    borne <- cbind(borne, overlap = 0)
    borne[Reduce(c, overlap), "overlap"] <- 1
    return(borne)
}
deconv2dLinearCoupled <- function(rawM, t, peak.detect, raw, fctFit, 
                                  smoothPenalty = 0, 
    knots = unique(c(t[1], stats::quantile(t, probs = seq(0, 1, 
                                                          length = (round(length(t)/3)))), 
        utils::tail(t, 1))), d = 3, l.shape = NULL) {
    K = length(knots) + d - 1
    # mass shifed correction
    mzNom <- round(peak.detect$Mz)[1]
    m <- as.numeric(row.names(rawM))
    # mass function
    n.peak <- nrow(peak.detect)
    spasymGauss <- function(x, par, raw) {
        exp(-(x - par["Mz"])^2/(2 * par["parameter.1"]^2)) * (x <= par["Mz"]) + exp(-(x - 
            par["Mz"])^2/(2 * par["parameter.2"]^2)) * (x > par["Mz"])
    }
    spsech2 <- function(m, par, raw) {
        1/(cosh((log(sqrt(2) + 1)/par[["parameter.1"]]) * (m - par[["Mz"]]))^2 * 
            (m <= par[["Mz"]]) + cosh((log(sqrt(2) + 1)/par[["parameter.2"]]) * (m - 
            par[["Mz"]]))^2 * (m > par[["Mz"]]))
    }
    spaveragePeak <- function(m, par, raw) {
        peakShape <- getPeaksInfo(raw)$peakShape$peakRef
        intervRef <- getPeaksInfo(raw)$peakShape$tofRef
        res <- rep(0, length(m))
        intervFit <- intervRef * (par[["parameter.1"]] + par[["parameter.2"]]) + 
            par[["Mz"]]
        interpol_ok <- which(intervFit[1] < m & m < utils::tail(intervFit, 1))
        if (length(interpol_ok) != 0) 
            res[interpol_ok] <- stats::spline(intervFit, peakShape, 
                                              xout = m[interpol_ok])$y
        res
    }
    SP <- apply(peak.detect, 1, function(z) eval(parse(text = paste0("sp", fctFit)))(m, 
        z, raw))
    t[1]<-ceiling(t[1])
    t[length(t)]<- floor(utils::tail(t,1))
    TIC <- splines::spline.des(knots = c(seq(-d, -1), knots, seq(utils::tail(knots, 
        1) + 1, utils::tail(knots, 1) + d)), x = t, ord = d + 1)$design  # add exterior knot
    X <- SP %x% TIC  #tensor product
    D <- diff(diag(K), differences = 2)  #root square o fthe penality matrix of order 2
    Y <- matrix(t(rawM), ncol = 1)
    n <- length(Y)
    Xa <- rbind(X, (diag(rep(1, n.peak)) %x% D) * sqrt(smoothPenalty))
    Y[(n + 1):(n + ncol(SP) * nrow(D))] <- 0
    lm <- lm(Y ~ Xa - 1)  # fit and return the penalized regression spline
    param <- lm$coefficients
    predRaw <- t(matrix(X %*% param, ncol = nrow(rawM), nrow = length(t)))
    mzMEAN <- mean(as.numeric(row.names(rawM)))
    mLarge <- seq(mzMEAN - 500 * round(mzMEAN)/10^6, mzMEAN + 500 * round(mzMEAN)/10^6, 
        by = diff(m)[1])
    SPlarge <- apply(peak.detect, 1, function(z) eval(parse(text = paste0("sp", 
                                                                          fctFit)))(mLarge, 
        z, raw))
    predRawLarge <- t(matrix(SPlarge %x% TIC %*% param, ncol = length(mLarge), nrow = length(t)))
    # g<-ggplot()+ggplot2::geom_line(mapping = ggplot2::aes(time,EIC,colour=param),
    # data=data.frame(time=getTimeInfo(raw)$time,EIC=colSums(predRaw),
    # param=as.character(smoothPenalty)))
    # plot(colSums(rawM),main=paste(K,smoothPenalty))
    # lines(colSums(predRaw),col='red')
    # deconvolution in ppb
    predPeak <- matrix(0, ncol = n.peak, nrow = length(t))
    for (i in seq_len(nrow(peak.detect))) {
        pred <- t(matrix(SP[, i] %x% TIC %*% param[((i - 1) * K + 1):(i * K)], 
                         nrow = length(t)))
        # quanti<- ppbConvert(peakList = cbind(Mz=rep(peak.detect$Mz[i]),
        # quanti=colSums(pred)), primaryIon = primaryIon, transmission = transmission, U
        # = U, Td = Td, pd = pd) plot(colSums(pred),main=peak.detect$Mz[i],pch=19)
        predPeak[, i] <- colSums(pred)
    }
    # image(t(rawM)) image(t(predRaw))
    
    return(list(predRaw = predRaw, predPeak = predPeak, 
                predRawLarge = predRawLarge,param=param))
}
GCV <- function(rawM, knots, t, timeLimit, stepSp = 0.01, d = 3) {
    trace <- colSums(rawM)
    bg <- timeLimit$backGround
    periods <- which(diff(bg) > 1)
    if (length(periods) >= 2) 
        borne <- c(t[1], t[bg[periods[2] - 1]]) else borne <- range(t)
    trace <- trace[t <= borne[2] & t >= borne[1]]
    knotsSub <- sort(unique(c(borne, knots[knots <= borne[2] & knots >= borne[1]])))
    tSub <- t[t <= borne[2] & t >= borne[1]]
    # plot(tSub,trace)
    X <- splines::spline.des(knots = c(seq(-d, -1), knotsSub, seq(utils::tail(knotsSub, 
        1) + 1, utils::tail(knotsSub, 1) + d)), x = tSub, ord = d + 1)$design  # add exterior knot
    K <- length(knotsSub) + d - 1
    D <- diff(diag(K), differences = 2)
    Y <- matrix(trace, ncol = 1)
    n <- length(Y)
    Y[(n + 1):(n + nrow(D))] <- 0
    gcv <- c()
    sp <- 0
    Xa <- rbind(X, D * sqrt(sp))
    fit <- stats::lm(Y ~ Xa - 1)
    diagA <- stats::influence(fit)$hat
    fitted <- fit$fitted.values
    gcv[as.character(sp)] <- n * sum((Y[seq_len(n)] - fitted[seq_len(n)])^2)/(n - 
        sum(diagA[seq_len(n)]))^2
    repeat {
        sp <- sp + stepSp
        Xa <- rbind(X, D * sqrt(sp))
        fit <- stats::lm(Y ~ Xa - 1)
        diagA <- stats::influence(fit)$hat
        fitted <- fit$fitted.values
        gcv[as.character(sp)] <- n * sum((Y[seq_len(n)] - fitted[seq_len(n)])^2)/(n - 
            sum(diagA[seq_len(n)]))^2
        if (utils::tail(diff(gcv), 1) > 0) 
            break
    }
    return(sp - stepSp)
}


deconv2d2linearIndependant <- function(rawM, time, peak.detect, raw, fctFit, 
                                       knots = NULL, 
    smoothPenalty = 0, l.shape = NULL) {
    mzAxis <- as.numeric(row.names(rawM))
    mzNom <- round(peak.detect$Mz)[1]
    n.peak <- nrow(peak.detect)
    matRaw <- matrix(0, nrow = nrow(rawM), ncol = ncol(rawM))
    matPeak <- matrix(0, ncol = n.peak, nrow = ncol(rawM))
    t1 <- Sys.time()
    for (i in seq_along(time)) {
        spectrum.m <- rawM[, i]
        # initialisation
        mz <- peak.detect$Mz
        lf <- peak.detect$parameter.1
        lr <- peak.detect$parameter.2
        param <- cbind(mz, lf, lr)
        if (fctFit == "sech2") {
            model <- apply(param, 1, function(x) 1/(cosh((log(sqrt(2) + 1)/x["lf"]) * 
                (mzAxis - x["mz"]))^2 * (mzAxis <= x["mz"]) + cosh((log(sqrt(2) + 
                1)/x["lr"]) * (mzAxis - x["mz"]))^2 * (mzAxis > x["mz"])))
        } else if (fctFit == "assymGauss") {
            model <- apply(param, 1, function(x) 1 * (exp(-(mzAxis - x["mz"])^2/(2 * 
                x["lf"]^2)) * (mzAxis <= x["mz"]) + exp(-(mzAxis - x["mz"])^2/(2 * 
                x["lr"]^2)) * (mzAxis > x["mz"])))
        } else if (fctFit == "averagePeak") {
            model <- apply(param, 1, function(x) {
                peakShape <- getPeaksInfo(raw)$peakShape$peakRef
                intervRef <- getPeaksInfo(raw)$peakShape$tofRef
                res <- rep(0, length(mzAxis))
                intervFit <- intervRef * (x[["lf"]] + x[["lr"]]) + x[["mz"]]
                interpol_ok <- which(intervFit[1] < mzAxis & mzAxis < utils::tail(intervFit, 
                  1))
                if (length(interpol_ok) != 0) 
                  res[interpol_ok] <- stats::spline(intervFit, peakShape, 
                                                    xout = mzAxis[interpol_ok])$y
                res
            })
        }
        fit <- stats::lm(spectrum.m ~ model - 1)
        par_estimated <- fit$coefficients
        fit.peak <- fit$fitted.values
        matRaw[, i] <- fit.peak
        # quantification
        quanti <- apply(rbind(mz, par_estimated, lf, lr), 2, function(x) {
            sum(eval(parse(text = fctFit))(x[1], x[3], x[4], x[2], mzAxis, 
                                           getPeaksInfo(raw)$peakShape))
        })
        matPeak[i, ] <- quanti
    }
    t2 <- Sys.time()
    return(list(predRaw = matRaw, predPeak = matPeak))
}
# Detect peak in a single nominal mass, same parameter as peakList
#### fit function ------
sech2 <- function(p, lf, lr, h, x, l.shape = NULL) {
    h/(cosh((log(sqrt(2) + 1)/lf) * (x - p))^2 * (x <= p) + cosh((log(sqrt(2) + 1)/lr) * 
        (x - p))^2 * (x > p))
}
sech2Inv <- function(p, lf, lr, h, x, l.shape = NULL) {
    c(p - 1/(log(sqrt(2) + 1)/lf) * acosh(sqrt(h/x)), p + 1/(log(sqrt(2) + 1)/lr) * 
        acosh(sqrt(h/x)))
}
asymGauss <- function(p, lf, lr, h, x, l.shape = NULL) {
    h * (exp(-(x - p)^2/(2 * lf^2)) * (x <= p) + exp(-(x - p)^2/(2 * lr^2)) * (x > 
        p))
}
asymGaussInv <- function(p, lf, lr, h, x, l.shape = NULL) {
    c(p - sqrt(-log(x/h) * 2 * lf^2), p + sqrt(-log(x/h) * 2 * lr^2))
}
averageInv <- function(p, delta, h, x, peakShape, tofRef, bin) {
    peak <- fit_averagePeak_function(t = p, delta = delta, h = h, 
                                     intervRef = tofRef, 
        peakShape = peakShape, bin)
    borneleft <- findEqualGreaterM(peak, x) - 1
    bornerigth <- findEqualGreaterM(peak[order(seq(1, length(peak)), decreasing = TRUE)], 
        x) - 1
    return(matrix(bin[c(max(borneleft, 1), length(peak) - max(bornerigth, 1) + 1)], 
        nrow = 1))
}
averagePeakInv <- function(p, lf, lr, h, x, l.shape) {
    mzAxis <- seq(p - 1000 * p/10^6, p + 1000 * p/10^6, length.out = 1000)
    averageInv(p = p, delta = lf + lr, h = h, x = x, peakShape = l.shape$peakRef, 
        tofRef = l.shape$tofRef, bin = mzAxis)
}
averagePeak <- function(m, lf, lr, h, x, l.shape) {
    fit_averagePeak_function(t = m, delta = lf + lr, h = h, 
                             intervRef = l.shape$tofRef, 
        peakShape = l.shape$peakRef, bin = x)
}
#'fit function average
#'@param t tof center of peak
#'@param delta FWHM of peak
#'@param h peak height
#'@param intervRef reference interval for peak shape
#'@param peakShape peak shape estimated on \code{intervalRef}
#'@param bin bin interval of peak will be fitted
#'@return peak function made on an average of reference peaks normalized
#'@keywords internal
fit_averagePeak_function <- function(t, delta, h, intervRef, peakShape, bin) {
    res <- rep(0, length(bin))
    intervFit <- intervRef * delta + t
    interpol_ok <- which(intervFit[1] < bin & bin < utils::tail(intervFit, 1))
    if (length(interpol_ok) != 0) 
        res[interpol_ok] <- stats::spline(intervFit, h * peakShape, 
                                          xout = bin[interpol_ok])$y
    res
}
#'Create cumulative function fit
#'
#'@param fit_function_str fit function who will be use in character
#'@param par_var_str parameters of fit function who change with the peak in a 
#'vector of character
#'@param par_fix_str parameters of fit function independent of the peak in a 
#'vector of character
#'@param n.peak number of peak
#'@return a list: \cr
#'   \code{init.names}: names of paramters for the initialization \cr
#'   \code{func.eval}: function who will be fitted
#'@keywords internal
cumulative_fit_function <- function(fit_function_str, par_var_str, 
                                    par_fix_str, n.peak) {
    formul.character <- ""
    init.names <- ""
    for (j in seq(1, n.peak)) {
        par_str.j <- ""
        for (n in seq(from = length(par_var_str), to = 1)) {
            par_str.j <- paste(paste("par$", par_var_str[n], j, sep = ""), par_str.j, 
                sep = ",")
        }
        for (n in seq_along(par_var_str)) {
            init.names <- c(init.names, paste(par_var_str[n], j, sep = ""))
        }
        formul.character <- paste(formul.character, fit_function_str, "(", par_str.j, 
            par_fix_str, ")+", sep = "")
        if (j == n.peak) {
            formul.character <- substr(formul.character, start = 1, 
                                       stop = nchar(formul.character) - 
                1)
        }
    }
    init.names <- init.names[-1]
    func.eval <- parse(text = formul.character)
    return(list(init.names = init.names, func.eval = func.eval))
}
fitPeak <- function(initMz, sp, mz.i, lower.cons, upper.cons, funcName, l.shape) {
    n.peak <- nrow(initMz)
    fit_function_str <- funcName
    par_var_str <- c("m", "lf", "lr", "h")
    par_fix_str <- c("x,l.shape")
    cum_fit <- cumulative_fit_function(fit_function_str, par_var_str, par_fix_str, 
        n.peak)
    initMz.names <- cum_fit$init.names
    func.eval <- cum_fit$func.eval
    function.char <- function(par, x, y) {
        eval(func.eval) - y
    }
    initMz <- as.list(t(initMz))
    names(initMz) <- initMz.names
    # Regression
    fit <- suppressWarnings(minpack.lm::nls.lm(par = initMz, lower = lower.cons, 
        upper = upper.cons, fn = function.char, x = mz.i, y = sp))
    par_estimated <- matrix(unlist(fit$par), nrow = 4)
    fit_peak <- function.char(fit$par, mz.i, rep(0, length(sp)))
    return(list(fit.peak = fit_peak, par_estimated = par_estimated, 
                function.fit.peak = function.char, 
        fit = fit))
}
#' Fit peak with average function
#' @param initTof list of initialisation in tof
#' @param l.shape peak shape average
#' @param sp spectrum 
#' @param bin tof axis
#' @param lower.cons lower constrain for fit
#' @param upper.cons upper constrain for fit
#' @return list with fit information 
#' @keywords internal
fit_averagePeak <- function(initTof, l.shape, sp, bin, lower.cons, upper.cons) {
    n.peak <- nrow(initTof)
    fit_function_str <- "fit_averagePeak_function"
    par_var_str <- c("t", "delta", "h")
    par_fix_str <- c("intervRef,peakShape,bin")
    cum_fit <- cumulative_fit_function(fit_function_str, par_var_str, par_fix_str, 
        n.peak)
    initTof.names <- cum_fit$init.names
    func.eval <- cum_fit$func.eval
    function.char <- function(par, intervRef, peakShape, bin, vec.peak) {
        eval(func.eval) - vec.peak
    }
    initTof <- as.list(t(initTof))
    names(initTof) <- initTof.names
    fit <- suppressWarnings(minpack.lm::nls.lm(par = initTof, lower = lower.cons, 
        upper = upper.cons, fn = function.char, intervRef = l.shape$tofRef, 
        peakShape = l.shape$peakRef, 
        bin = bin, vec.peak = sp))
    fit.peak <- function.char(fit$par, l.shape$tofRef, l.shape$peakRef, bin, rep(0, 
        length(sp)))
    par_estimated <- matrix(unlist(fit$par), nrow = 3)
    return(list(fit.peak = fit.peak, par_estimated = par_estimated, 
                function.fit.peak = function.char, 
        fit = fit))
}
#' Determine peak shape from raw data in tof
#'
#'This function use the method descibe by average and al 2013, for determine a 
#'peak shape from the raw data : \cr 
#'$peak_i(Delta_i,A_i,t_i) = interpolation(x= tof.ref * Delta_i + t_i,y = A_i * 
#'peak.ref, xout= TOF_i )$ where peak.ref and tof.ref are peaks reference use 
#'for mass calibration.
#' @param raw a \code{\link[ptairMS]{ptrRaw-class}} object 
#' @param plotShape if true plot each reference peak and the average peak 
#' (the peak shape)
#' @return A list of two vectors which are the reference peak normalized tof 
#' and intensity.
#' @keywords internal
determinePeakShape <- function(raw, plotShape = FALSE){
    # mass used for calibration
    massRef <- getCalibrationInfo(raw)$calibMassRef
    massRef.o <- massRef[order(massRef)]
    sp <- rowSums(getRawInfo(raw)$rawM)/ncol(getRawInfo(raw)$rawM)
    mz <- getRawInfo(raw)$mz
    # spectrum around calibration masses interval <- raw@calibSpectr
    interval <- lapply(massRef.o, function(x) {
        th <- 0.4
        # tof<-which(mz<= x+ th & mz>=x - th)
        mzx <- mz[mz <= x + th & mz >= x - th]
        sp <- sp[mz <= x + th & mz >= x - th]
        sp <- sp - snipBase(sp)
        list(mz = mzx, signal = sp)
    })
    
    # find center and width peak for each mass
    mz_centre <- vapply(interval, function(x) {
        sp <- x$signal
        mz <- x$mz
        localMax <- LocalMaximaSG(sp = sp, minPeakHeight = max(sp) * 0.1)
        if(length(localMax)==1){
            interpol <- stats::spline(mz[(localMax - 4):(localMax + 4)], sp[(localMax - 
                                                                                 4):(localMax + 4)], n = 1000)
            sg <- signal::sgolayfilt(interpol$y, n = 501, p = 3)  #n/
            center <- interpol$x[which.max(sg)]
            return(center)  
        } else return(0)
        
    }, FUN.VALUE = 1.1)
    interval<-interval[mz_centre!=0]
    massRef<-massRef[mz_centre!=0]
    n.mass <- length(massRef)
    mz_centre<-mz_centre[mz_centre!=0]
    deltaMz <- vapply(interval, function(x) width(tof = x$mz, peak = x$signal)$delta, 
                      FUN.VALUE = 1.1)
    deltaTof <- vapply(interval, function(x) width(tof = mzToTof(x$mz, calibCoef = getCalibrationInfo(raw)$calibCoef[[1]]), 
        peak = x$signal)$delta, FUN.VALUE = 1.1)
    tof_centre <- vapply(interval, function(x) {
        sp <- x$signal
        mz <- mzToTof(x$mz, calibCoef = getCalibrationInfo(raw)$calibCoef[[1]])
        localMax <- LocalMaximaSG(sp = sp, minPeakHeight = max(sp) * 0.2)
        interpol <- stats::spline(mz[(localMax - 4):(localMax + 4)], sp[(localMax - 
            4):(localMax + 4)], n = 1000)
        sg <- signal::sgolayfilt(interpol$y, n = 501, p = 3)  #n/
        center <- interpol$x[which.max(sg)]
        return(center)
    }, FUN.VALUE = 1.1)
    # normalized interval
    intervalnMz <- list()
    for (j in seq_len(n.mass)) {
        intervalnMz[[j]] <- (interval[[j]]$mz - mz_centre[[j]])/deltaMz[j]
    }
    intervalnTof <- list()
    for (j in seq_len(n.mass)) {
        intervalnTof[[j]] <- (mzToTof(interval[[j]]$mz, getCalibrationInfo(raw)$calibCoef[[1]]) - tof_centre[[j]])/deltaTof[j]
    }
    peakShape <- lapply(list(intervalnMz, intervalnTof), function(interval.n) {
        # reference interval : shorter
        interval.n.length <- vapply(interval.n, length, 1)
        interval.n.borne <- vapply(interval.n, range, c(1, 2))
        index.inter.longest <- which.max(interval.n.length)
        interval.n.longest <- interval.n[[index.inter.longest]]
        interval.ref.borne <- interval.n.borne[, which.min(apply(abs(interval.n.borne), 
            2, min))]
        interval.ref <- interval.n.longest[interval.ref.borne[1] < interval.n.longest & 
            interval.n.longest < interval.ref.borne[2]]
        peak <- matrix(0, ncol = n.mass, nrow = length(interval.ref))
        peak[, index.inter.longest] <- interval[[index.inter.longest]]$signal[interval.ref.borne[1] < 
            interval.n.longest & interval.n.longest < interval.ref.borne[2]]
        peak[, index.inter.longest] <- peak[, index.inter.longest]/max(peak[, index.inter.longest])
        # interpolation from the ref interval
        for (i in seq_len(n.mass)[-index.inter.longest]) {
            interpolation <- stats::spline(interval.n[[i]], interval[[i]]$signal, 
                xout = interval.ref)
            # baseline correction interpolation$y <- interpolation$y -
            # snipBase(interpolation$y,widthI = 4)
            # normalization
            indexPeakRef <- which(massRef == massRef[i])
            peak[, i] <- interpolation$y/stats::spline(interval.ref, interpolation$y, 
                xout = 0)$y
        }
        # average of cumulated signal
        peak_ref <- rowSums(peak)/n.mass
        return(list(intervalref = interval.ref, peakShape = peak_ref))
    })
    return(list(peakShapetof = list(tofRef = peakShape[[2]]$intervalref, peakRef = peakShape[[2]]$peakShape), 
        peakShapemz = list(tofRef = peakShape[[1]]$intervalref, peakRef = peakShape[[1]]$peakShape)))
}
# Baseline Correction --------------------------------------------
baselineEstimation <- function(sp, d = 2, eps = 1e-06) {
    x <- seq(1, length(sp))
    reg0 <- stats::lm(sp ~ poly(x, d))
    a0 <- reg0$coefficients
    yhat <- stats::predict(reg0)
    y <- (sp < yhat) * sp + yhat * (sp >= yhat)
    reg1 <- stats::lm(y ~ stats::poly(x, d))
    a1 <- reg1$coefficients
    yhat <- stats::predict(reg1)
    y <- (sp < yhat | yhat < min(sp)) * sp + yhat * (sp >= yhat & yhat >= min(sp))
    reg2 <- stats::lm(y ~ stats::poly(x, d))
    a2 <- reg2$coefficients
    c = 1
    while (sqrt(sum((a1 - a2)^2)) > (mean(sp) * eps) & c < 200) {
        a1 <- a2
        yhat <- stats::predict(reg2)
        y <- (sp < yhat) * sp + yhat * (sp >= yhat)
        reg2 <- stats::lm(y ~ stats::poly(x, d))
        a2 <- reg2$coefficients
        c = c + 1
    }
    return(yhat)
}

baselineEstimation2D <- function(rawM, dt = 1, dm = 2, p = 0.1, smoothPenalty = 0.1, 
    quantil = 0.01) {

    m <- splines::bs(seq(1, nrow(rawM)), df = dt, intercept = TRUE)
    t <- splines::bs(seq(1, ncol(rawM)), df = dm, intercept = TRUE)
    Y <- c(rawM)
    X <- t %x% m
    # rq<-quantreg::rq(Y ~X,tau = quantil)
    # bl<-matrix(rq$fitted.values,ncol=ncol(rawM),nrow=nrow(rawM)) first reg
    Dm <- diff(diag(ncol(m)), differences = 2)  #root square of the penality 
    # matrix of order 2
    Dt <- diff(diag(ncol(t)), differences = 2)
    n <- length(Y)
    Xa <- rbind(X, sqrt(smoothPenalty) * (diag(rep(1, ncol(t))) %x% Dm))
    Y[(n + 1):(n + (nrow(Dm) * ncol(t)))] <- 0
    w <- rep(1, length(Y))
    lm <- lm(Y ~ Xa - 1, weights = w)  # fit and return the penalized regression spline
    a1 <- lm$coefficients
    Z <- X %*% a1
    # w[which(Y[seq_len(n)]>Z & Z > min(Y[seq_len(n)]))]<-p w[which(Y[seq_len(n)]<=Z
    # | Z <= min(Y[seq_len(n)]))]<-1-p
    w[which(Y[seq_len(n)] > Z & Z > min(Y[seq_len(n)]))] <- 0
    w[which(Y[seq_len(n)] <= Z | Z <= min(Y[seq_len(n)]))] <- exp((Y[seq_len(n)] - 
        Z)[which(Y[seq_len(n)] <= Z | Z <= min(Y[seq_len(n)]))]/mean(abs((Y[seq_len(n)] - 
        Z)[Y[seq_len(n)] - Z < 0])) * c)
    # second reg
    reg2 <- stats::lm(Y ~ Xa - 1, weights = w)
    a2 <- reg2$coefficients
    c = 1
    # #(max(rawM)-min(rawM))*0.001 while( sqrt(mean((a1-a2)^2)) > 1e-10 & c < 20){
    while (mean(abs((Y[seq_len(n)] - Z)[Y[seq_len(n)] - Z < 0])) > 0.001 * mean(abs(c(rawM)))) {
        # while( sqrt(mean((bl1D-rowSums(matrix(Z2,ncol=ncol(rawM),nrow=nrow(rawM))))^2))
        # > sqrt(mean(bl1D)^2)*0.0001 & c < 100){
        a1 <- a2
        Z <- X %*% a1
        # w[which(Y[seq_len(n)]>Z & Z > min(Y[seq_len(n)]))]<-p w[which(Y[seq_len(n)]<=Z
        # | Z <= min(Y[seq_len(n)]))]<-1-p
        w[which(Y[seq_len(n)] > Z & Z > min(Y[seq_len(n)]))] <- 0
        w[which(Y[seq_len(n)] <= Z | Z <= min(Y[seq_len(n)]))] <- exp((Y[seq_len(n)] - 
            Z)[which(Y[seq_len(n)] <= Z | Z <= min(Y[seq_len(n)]))]/mean(abs((Y[seq_len(n)] - 
            Z)[Y[seq_len(n)] - Z < 0])) * c)
        reg2 <- stats::lm(Y ~ Xa - 1, weights = w)
        a2 <- reg2$coefficients
        c = c + 1
    }
    Z <- X %*% a2
    # Y <- (Y<Z | Z < min(Y))*Y + Z*(Y>=Z & Z >= min(Y) )
    bl <- matrix(Z, ncol = ncol(rawM), nrow = nrow(rawM))
    return(bl)
}
#'Baseline estimation
#'
#'@param sp an array with spectrum values
#'@param widthI width of interval
#'@param iteI number of iteration
#'
#'@return baseline estimation of the spectrum
#'@keywords internal
snipBase <- function(sp, widthI = 11, iteI = 5) {
    runave <- function(x, widI) {
        z <- x
        smoVi <- seq(widI + 1, length(x) - widI)
        z[smoVi] <- (z[smoVi - widI] + z[smoVi + widI])/2
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
#' @keywords internal
TransmissionCurve <- function(x, y) {
    model <- stats::lm(y ~ I(x^2) + x + I(x^0.5))
    curve <- stats::predict(model, new_data = data.frame(x = seq(1, utils::tail(y, 
        1))))
    extra <- Hmisc::approxExtrap(x, curve, xout = seq(1, 1000), method = "linear", 
        n = 50)
    return(extra)
}
# PPb concentration
ppbConvert <- function(peakList, transmission, U, Td, pd, k = 2 * 10^-9) {
    Tr_primary <- transmission[2, 1]
    x <- transmission[1, ][transmission[1, ] != 0]
    y <- transmission[2, ][transmission[1, ] != 0]
    TrCurve <- TransmissionCurve(x, y)
    TrCurve$y <- vapply(TrCurve$y, function(x) max(0.001, x), 0.1)
    ppb <- peakList[, "quanti"] * 10^9 * mean(U) * 2.8 * 22400 * 1013^2 * (273.15 + 
        mean(Td))^2 * Tr_primary/(k * 9.2^2 * mean(pd)^2 * 6022 * 10^23 * 273.15^2 * 
        TrCurve$y[peakList[, "Mz"]])
    ppb * 1000
}
## Other ----
#' Calculate the FWHM (Full Width at Half Maximum) in raw data
#'
#' @param tof A vector of tof interval
#' @param peak A vector of peak Intensity
#' @param fracMaxTIC the fraction of the maximum intenisty to compute the width
#' @return the delta FWHM in tof 
#' @keywords internal
width <- function(tof, peak, fracMaxTIC = 0.5) {
    hm <- max(peak) * fracMaxTIC
    sup1 <- findEqualGreaterM(peak[seq_len(which.max(peak))], hm)
    sup2 <- unname(which.max(peak)) + FindEqualLess(peak[(which.max(peak) + 1):length(peak)], 
        hm)
    # equation intrepolation linéaire entre sup et sup-1 = hm
    delta_tof <- unname(vapply(c(sup1, sup2), function(x) {
        (hm * (tof[x] - tof[x - 1]) - (tof[x] * peak[x - 1] - tof[x - 1] * peak[x]))/(peak[x] - 
            peak[x - 1])
    }, FUN.VALUE = 0.1))
    delta <- unname(abs(delta_tof[2] - delta_tof[1]))
    list(delta = delta, delta_tof = delta_tof - 1)
}
FindEqualLess <- function(x, values) {
    idx <- 1
    while (x[idx] > values) idx = idx + 1
    idx
}
