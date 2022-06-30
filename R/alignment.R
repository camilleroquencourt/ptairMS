utils::globalVariables(c("Mz", "quanti", "group", "background", "fracExp", 
                         "pValGreater",
                         "pValLess"))

#' Alignment with kernel gaussian density
#'
#' @param peakTab table with comlumn : mass, quantification, and groups 
#' number to aligned
#' @param ppmGroup width of sub group created beafore density estimation in ppm
#' @param dmzGroup width of sub group created beafore density estimation in Da
#' @return A list containing groups formed by alignment.
#' @keywords internal
align <- function(peakTab, ppmGroup = 70, dmzGroup = 0.001) {
    ### faire un test si c'est null renvoyer nul
    porder <- order(peakTab[, "Mz"])
    peaksL <- peakTab[porder, , drop = FALSE]
    
    # delete the unit in col names
    legList <- colnames(peaksL)
    quantiName <- legList[grep("quanti", legList)]
    legList[grep("quanti", legList)] <- "quanti"
    bgName <- legList[grep("background", legList)]
    legList[grep("background", legList)] <- "background"
    
    colnames(peaksL) <- legList
    
    # Get the number of different values in group
    n.group <- nlevels(as.factor(peaksL[, "group"]))
    if (n.group == 0) 
        return(list(peakTab))  #there is no peak 
    
    # Cuts the mass axis into 1 Da size intervals around the nominal mass
    
    nPoints <- 1024
    mzrange <- range(peaksL[, "Mz"])
    inter = 1
    massInter <- seq(floor(mzrange[1]) - 0.5, ceiling(mzrange[2]) + 0.5, inter)
    
    # Find the position of peaksL in massInter
    masspos <- findEqualGreaterM(peaksL[, "Mz"], massInter)
    
    # Initialize variables The list who will be return, contains each group of same
    # masses
    groupPeakList <- list()
    numGroupPeak <- 0  # the number of group formed
    
    # Loop on each group
    index_group <- which(diff(masspos, lag = 1) >= 1)
    for (i in index_group) {
        
        start <- masspos[i]  # Get start and end of the sub group
        end <- masspos[i + 1] - 1
        # end <- masspos[i + 2] - 1
        
        subpeakl <- peaksL[start:end, , drop = FALSE]  # Take submatrix
        
        # calculated the density Determining the base of the signal to skip for a
        # resolution.
            bw <- max(massInter[i] * ppmGroup/(10^6), dmzGroup)
            den <- stats::density(subpeakl[, "Mz"], bw, from = massInter[i] - bw, 
                                    to = massInter[i + 1] + bw, n = 1024)
        
        

        maxy <- max(den$y)
        den$y[which(den$y <= bw * maxy)] <- 0
        plim <- list(linflex = -1, rinflex = 0, state = 0)
        
        
        # While we detect peak in the density
        repeat {
            plim <- findLimDensity(den$y, plim$rinflex + 1, plim$state)
            if (plim$linflex == plim$rinflex) 
                break  ###In this case the last peak is found.
            
            subpeakl.group <- makeSubGroup(subpeakl, den, plim)
            if (nrow(subpeakl.group) != 0) {
                numGroupPeak <- numGroupPeak + 1
                groupPeakList[[numGroupPeak]] <- subpeakl.group
            }
            
            
        }  # end of repeat
    }  # end of error
    return(groupPeakList)
    
}

#' Use in align function. return a peak group thanks to kernel gaussian density 
#' in a peak matrix.
#' @param subpeakl a matrix with mz, ppb, background and group column. 
#' @param den the kernel gaussian density estimated on subpeakl
#' @param plim the limit of a peak in the density of the group who will pe fromed
#' @return the sub peakgroup  
#' @keywords internal
makeSubGroup <- function(subpeakl, den, plim) {
    
    selectedPGroup <- which(subpeakl[, "Mz"] >= den$x[plim$linflex] & subpeakl[, 
        "Mz"] <= den$x[plim$rinflex])
    
    subpeakl.group <- subpeakl[selectedPGroup, , drop = FALSE]
    
    ## if two variables of same level group are confused we keep the more intense one
    ## for mz center and sum theme
    if ((length(unique(subpeakl.group[, "group"])) < dim(subpeakl.group)[1])) {
        subpeakl.group <- data.table::as.data.table(subpeakl.group)
        if ("background" %in% colnames(subpeakl.group)) {
            if ("pValGreater" %in% colnames(subpeakl.group)) {
                subpeakl.group <- as.matrix(subpeakl.group[, list(Mz = Mz[which.max(quanti)], 
                  quanti = sum(quanti, na.rm = TRUE), background = sum(background, 
                    na.rm = TRUE), pValGreater = pValGreater[which.max(quanti)], 
                  pValLess = pValLess[which.max(quanti)]), by = group])
            } else {
                subpeakl.group <- as.matrix(subpeakl.group[, list(Mz = Mz[which.max(quanti)], 
                  quanti = sum(quanti, na.rm = TRUE), background = sum(background, 
                    na.rm = TRUE)), by = group])
            }
            
        } else subpeakl.group <- as.matrix(subpeakl.group[, list(Mz = Mz[which.max(quanti)], 
            quanti = sum(quanti, na.rm = TRUE)), by = group])
        
        
    }
    return(subpeakl.group)
}


#' aggregate peakgroup for align function  
#' @param subGroupPeak teh group tp aggregate
#' @param n.exp number of expiration done in the file
#' @return a matrix with the median of mz, mean of ppb, ppb in background, 
#'  and percentage of expiration where the peak is detected
#'  @keywords internal
aggregate <- function(subGroupPeak, n.exp) {
    if (nrow(subGroupPeak) == 0) 
        return(subGroupPeak)  #empty data frame
    subGroupPeak <- data.table::data.table(subGroupPeak)
    subGroupPeak <- subGroupPeak[, list(Mz = mean(Mz), 
                                        quanti = mean(quanti[group != 0]), 
                                        background = mean(quanti[group == 0]), 
                                        fracExp = round(sum(group !=0)/n.exp, 1))]
    if (is.na(subGroupPeak[, quanti])) 
        return(NULL)
    
    return(subGroupPeak)
}


aggregateTemporalFile <- function(time, indTimeLim, matPeak, funAggreg, dbl = 4) {
    # agregation of expirations and bg
    indLim <- indTimeLim$exp
    indBg <- indTimeLim$backGround
    
    XIC <- as.matrix(matPeak[, (ncol(matPeak) - length(time) + 1):ncol(matPeak), 
                             drop = FALSE])
    
    if (!is.null(indBg)) {

        ## baseline Corrected
        XICbl <- t(apply(XIC, 1, function(x) {
            lm <- stats::lm(cps ~ stats::poly(point, d = dbl), 
                            data = data.frame(cps = x[indBg],
                                              point = indBg))
            bl <- stats::predict(lm, newdata = data.frame(point = seq_along(time)))
            x <- x - bl
        }))
        
        
        for(j in seq_len(nrow(XIC))){
            x<-XIC[j,]
            lm <- stats::lm(cps ~ stats::poly(point, d = dbl), 
                            data = data.frame(cps = x[indBg],
                                              point = indBg))
            bl <- stats::predict(lm, newdata = data.frame(point = seq_along(time)))
            x <- x - bl
        }
        
        
        # test significativité comparaison de moyenne normal
        indExp <- Reduce(c, apply(indLim, 2, function(x) seq(x[1], x[2])))
        
        pValGreater <- apply(XIC, 1, 
                                function(x) stats::t.test(x[indExp],
                                                            x[indBg],
                                                            alternative = "greater")$p.value)
        pValLess <- apply(XIC, 1, function(x) stats::t.test(x[indExp], x[indBg], 
            alternative = "less")$p.value)
        pValGreater <- stats::p.adjust(pValGreater, method = "fdr")
        pValLess <- stats::p.adjust(pValLess, method = "fdr")
        
        ## aggregate
        bgIn <- time[indBg]
        exp <- time[indExp]
        colnames(XIC) <- as.character(time)
        quantiExp <- apply(XIC, 1, 
                            function(x) funAggreg(x[colnames(XIC) %in% exp]))
        quantiBg <- apply(XIC, 1, 
                            function(x) funAggreg(x[colnames(XIC) %in% bgIn]))
        
        colnames(XICbl) <- as.character(time)
        quantiExpBl <- apply(XICbl, 1, 
                                function(x) abs(funAggreg(x[colnames(XIC) %in% 
                                                                exp])))
        
        matPeakAg <- cbind(matPeak[, 1], quanti_cps = quantiExp, 
                            background_cps = quantiBg, 
            diffAbs_cps = quantiExpBl, pValLess, pValGreater)
    } else {
        quantiExp <- apply(XIC, 1, function(x) mean(x))
        matPeakAg <- cbind(matPeak[, 1], quanti_cps = quantiExp, 
                            background_cps = NA, 
            diffAbs_cps = NA)
    }
    
    return(matPeakAg)
}



### Alignsamples method -----

#' Alignment between samples
#' 
#' \code{AlignSamples} performs alignment between samples (i.e. the matching of 
#' variables between the peak lists within the \code{ptrSet} object) by using a 
#' kernel gaussian density (Delabriere et al, 2017).
#' This function returns an \code{\link[Biobase]{ExpressionSet}}, which contains 
#' the matrix of peak intensities, the sample metadata (borrowed from the
#' input ptrSet) and the variable metadata which contains the peak intensities in 
#' the background. 
#' Two filters may be applied to:
#' \itemize{
#' \item  keep only variables with a significant higher intensity in the 
#' expirations compared to the background (i.e., a p-value less than 
#' \code{pValGreaterThres}) for at least \code{fracExp} % of the samples
#' \item keep only variables which are detected in more 
#' than \code{fracGroup} percent of the samples (or \code{group})
#' }
#' If you do not want to apply those filters, set \code{fracGroup} to 0 and 
#'\code{pValGreaterThres} to 1.
#' @param X ptrSet already processed by the \code{\link[ptairMS]{detectPeak}} 
#' function
#' @param ppmGroup ppm maximal width for an mz group 
#' @param fracGroup only variables present in \code{fracGroup} 
#' percent of at least one \code{group} will be kept (if 0 the filter is 
#' not applied)
#' @param group character: sampleMetadata data column name. If \code{NULL},
#' variables not present in \code{fracGroup} percent of samples will be deleted. 
#' Else, variables not present in \code{fracGroup} percent in in at least one 
#' group group will 
#' be removed. 
#' @param pValGreaterThres threshold of the p-value for the unilateral
#' testing that quantification (in cps) of expiration points are greater than 
#' the intensities in the background. 
#' @param pValLessThres threshold of the p-value for the unilateral
#' testing that quantification (in cps) of expiration points are less than the 
#' intensities of the background.
#' @param fracExp fraction of samples which must have their p-value less than 
#' \code{pValGreaterThres} and \code{pValLessThres}
#' @param quantiUnit ppb, ncps or cps
#' @param bgCorrected logical: should the peak table contain the background 
#' corrected values?
#' @param dmzGroup minimum mz width to be used for grouping the features 
#' (required for low masses)
#' @return an \code{\link[Biobase]{ExpressionSet}} (Biobase object)
#' @examples 
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, 
#' setName="exhaledPtrset", mzCalibRef = c(21.022, 60.0525),
#' fracMaxTIC = 0.7, saveDir = NULL )
#' exhaledPtrset<-detectPeak(exhaledPtrset,mzNominal=c(21,60,79))
#' eset <- alignSamples(exhaledPtrset,pValGreaterThres=0.05)
#' Biobase::exprs(eset)
#' Biobase::fData(eset)
#' Biobase::pData(eset)
#' @references Delabriere et al., 2017
#' @rdname alignSamples
#' @import data.table
#'@export

setMethod(f = "alignSamples", signature = "ptrSet", 
          function(X, ppmGroup = 70, fracGroup = 0.8,
                    group = NULL, fracExp = 0.3, pValGreaterThres = 0.001,
                    pValLessThres = 0, quantiUnit = c("ppb", "ncps", "cps")[1], 
                    bgCorrected = TRUE, dmzGroup = 0.001) {
    
    sampleMetadata <- getSampleMetadata(X)
    peakList <- getPeakList(X)
    eSet <- alignSamplesFunc(peakList = peakList, 
                                sampleMetadata = sampleMetadata,
                                ppmGroup = ppmGroup, fracGroup = fracGroup, 
                                group = group, 
                                pValGreaterThres = pValGreaterThres, 
        pValLessThres = pValLessThres, fracExp = fracExp, dmzGroup = dmzGroup, 
        quantiUnit = quantiUnit, 
        bgCorrected = bgCorrected)
    
    return(eSet)
})


alignSamplesFunc <- function(peakList, sampleMetadata, ppmGroup = 100, 
                                fracGroup = 0.8,
                                group = NULL, dmzGroup = 0.001, 
                                pValGreaterThres = 0.005, pValLessThres = 0, 
                                fracExp = 0.3, 
                                quantiUnit = c("ppb", "ncps", "cps")[1], 
                                bgCorrected = TRUE) {
    
    
    peakList <- lapply(peakList,
                       function(x) {

                           x<-as.data.table(Biobase::fData(x)[, -c(2,3, 4)])


                           })
    
    testquantiUnit <- Reduce(c, lapply(peakList, function(x) {
        x <- as.matrix(x)
        all(is.na(x[, paste0("quanti_", quantiUnit)]))
    }))
    
    
    if (any(testquantiUnit)) {
        if (quantiUnit == "ppb") {
            testquantiUnitncps <- Reduce(c, lapply(peakList, function(x) all(is.na(x[, 
                paste0("quanti_", "ncps")]))))
            if (all(testquantiUnitncps)) {
                quantiUnit <- "ncps"
            } else {
                quantiUnit <- "cps"
            }
        } else {
            quantiUnit <- "cps"
        }
    }
   
    
    .SD <- data.table::.SD
        
    peakList <- lapply(peakList, function(x) {
            backgroundCOL<-which(grepl("background",colnames(x)))
            if(all(!is.na(x[,.SD,.SD=colnames(x)[backgroundCOL]]))){
                
                diff <- which(grepl("diffAbs", colnames(x)))
                quanti <- which(grepl("quanti", colnames(x)))
                bg <- which(grepl("background", colnames(x)))
                other <- seq_len(ncol(x))[-c(quanti, bg, diff)]
                
                if (bgCorrected) {
                    quanti <- diff[which(grepl(quantiUnit, colnames(x)[diff]))]
                } else quanti <- quanti[which(grepl(quantiUnit, colnames(x)[quanti]))]
                colnames(x)[quanti] <- paste("quanti", quantiUnit, sep = "_")
                bg <- bg[which(grepl(quantiUnit, colnames(x)[bg]))]
                
               
                x <- x[, .SD, .SD = colnames(x)[c(other, quanti, bg)]]
                
                
                x 
            }else {
                quanti <- which(grepl("quanti", colnames(x)))
                other <- seq_len(ncol(x))[-c(quanti,backgroundCOL)]
                quanti <- quanti[which(grepl(quantiUnit, colnames(x)[quanti]))]
                x <- x[, .SD, .SD = colnames(x)[c(other, quanti)]]
                x
            }
            
        })
    
    
    ## add column group with Samples group number
    mat <- NULL
    peakList <- lapply(peakList, function(x) as.matrix(x))
    
    for (i in seq_along(peakList)) {
        mat_temp <- cbind(peakList[[i]], group = i)
        mat <- rbind(mat, mat_temp)
    }
    
    
    # make subgroup of peak
    groupList <- align(peakTab = as.matrix(mat), ppmGroup, dmzGroup)
    
    # Get number of samples
    nSample <- length(peakList)
    
    # aggregate
    groupMat <- do.call(rbind, lapply(groupList, function(y) {
        y <- data.table::as.data.table(y)
        y <- y[, list(Mz = stats::median(Mz), quanti = paste(quanti, collapse = "_"), 
            background = if ("background" %in% colnames(y)) paste(background, collapse = "_"), 
            Samples = paste(group, collapse = "_"), nSamples = round(length(group)/nSample, 
                1))]
        y
    }))
    
    # formatting the final matrix
    mat.final.Exp <- apply(groupMat[, c("quanti", "Samples")], MARGIN = 1, function(x) {
        x<-as.matrix(x)
        output <- rep(NA, nSample)
        ch.Area <- x[1]
        ch.Area.split <- strsplit(ch.Area, split = "_")[[1]]
        vec.Area <- as.numeric(ch.Area.split)
        ch.samples <- x[2]
        vec.samples <- as.numeric(strsplit(ch.samples, split = "_")[[1]])
        output[vec.samples] <- vec.Area
        output
    })
    
    # matrix function for the case of one sample and mat.final is an array
    X <- t(matrix(mat.final.Exp, nrow = nSample))
    rownames(X) <- round(groupMat[, Mz], 4)
    colnames(X) <- names(peakList)
    
    filtered <- fliterEset(X, sampleMetadata, groupMat, groupList, peakList, group, 
        fracGroup, pValGreaterThres, pValLessThres, fracExp)
    
    if (is.null(filtered)) 
        return(NULL)
    
    X <- filtered$X
    Xbg <- filtered$Xbg
    
    if (quantiUnit == "ppb") {
        Xbg <- Xbg[(round(as.numeric(rownames(X))) >= 21), , drop = FALSE]
        X <- X[round(as.numeric(rownames(X))) >= 21, , drop = FALSE]
    }
    
    # adding the ion_mass as the first column in the fData
    Xbg <- cbind.data.frame(ion_mass = as.numeric(rownames(X)), as.data.frame(Xbg, 
        stringsAsFactors = FALSE), stringsAsFactors = FALSE)
    
    order <- order(as.numeric(rownames(X)))
    X <- X[order, , drop = FALSE]
    Xbg <- Xbg[order, , drop = FALSE]
    rownames(Xbg) <- rownames(X)
    message(nrow(X), " peaks aligned")
    featureData <- Biobase::AnnotatedDataFrame(as.data.frame(Xbg))
    sampleMetadata <- as.data.frame(sampleMetadata[colnames(X), , drop = FALSE])
    return(Biobase::ExpressionSet(assayData = X, 
                                    phenoData = Biobase::AnnotatedDataFrame(sampleMetadata), 
        featureData = featureData, annotation = quantiUnit))
}


fliterEset <- function(X, sampleMetadata, groupMat, groupList, peakList, group, 
                        fracGroup, 
    pValGreaterThres, pValLessThres, fracExp) {
    
    nSample <- ncol(X)
    
    # filter on samples frenquency
    if (!is.null(group)) {
        groupFac <- as.factor(sampleMetadata[, group])
        test <- apply(X, 1, function(x) vapply(levels(groupFac), 
                                                function(y) sum(!is.na(x[groupFac == 
            y]))/length(x[groupFac == y]) >= fracGroup, TRUE))
        keepVar <- apply(test, 2, any)
        
    } else {
        keepVar <- apply(X, 1, function(y) sum(!is.na(y))/length(y) >= fracGroup)
        
    }
    
    X <- X[keepVar, , drop = FALSE]
    
    if (nrow(X) == 0) {
        warning("peakList is empty, there is no peak presents in 
                    more thant ", 
            fracGroup, "% of samples")
        return(NULL)
    }
    
    if ("background" %in% colnames(groupMat)) {
        mat.final.bg <- apply(groupMat[, c("background", "Samples")], MARGIN = 1, 
            function(x) {
                output <- rep(NA, nSample)
                ch.Area <- x[1]
                ch.Area.split <- strsplit(ch.Area, split = "_")[[1]]
                vec.Area <- as.numeric(ch.Area.split)
                
                ch.samples <- x[2]
                vec.samples <- as.numeric(strsplit(ch.samples, split = "_")[[1]])
                
                output[vec.samples] <- vec.Area
                output
            })
        Xbg <- t(matrix(mat.final.bg, nrow = nSample))
        Xbg <- Xbg[keepVar, , drop = FALSE]
        rownames(Xbg) <- rownames(X)
        colnames(Xbg) <- paste("Background -", names(peakList))
        
        
        if (pValGreaterThres != 1 & "pValGreater" %in% colnames(peakList[[1]])) {
            groupMatpVal <- do.call(rbind, lapply(groupList, function(y) {
                y <- data.table::as.data.table(y)
                y <- y[, list(Mz = stats::median(Mz), pValGreater = paste(pValGreater, 
                  collapse = "_"), pValLess = paste(pValLess, collapse = "_"), 
                    Samples = paste(group, 
                  collapse = "_"), nSamples = round(length(group)/nSample, 1))]
                y
            }))
            
            matPvalGreater <- apply(groupMatpVal[, c("pValGreater", "Samples")], 
                MARGIN = 1, function(x) {
                  output <- rep(NA, nSample)
                  ch.Area <- x[1]
                  ch.Area.split <- strsplit(ch.Area, split = "_")[[1]]
                  vec.Area <- as.numeric(ch.Area.split)
                  ch.samples <- x[2]
                  vec.samples <- as.numeric(strsplit(ch.samples, split = "_")[[1]])
                  output[vec.samples] <- vec.Area
                  output
                })
            
            
            nFile <- length(peakList)
            matPvalGreater <- t(matrix(matPvalGreater, nrow = nSample))
            matPvalGreater <- matPvalGreater[keepVar, , drop = FALSE]
            if (!is.null(fracExp)) {
                Keep1 <- which(apply(matPvalGreater, 1, 
                                        function(x) sum(x < pValGreaterThres, 
                  na.rm = TRUE)/sum(!is.na(x)) >= fracExp))
            } else {
                Keep1 <- which(matPvalGreater < pValGreaterThres)
            }
            
            
            matPvalLess <- apply(groupMatpVal[, c("pValLess", "Samples")], MARGIN = 1, 
                function(x) {
                  output <- rep(NA, nSample)
                  ch.Area <- x[1]
                  ch.Area.split <- strsplit(ch.Area, split = "_")[[1]]
                  vec.Area <- as.numeric(ch.Area.split)
                  
                  ch.samples <- x[2]
                  vec.samples <- as.numeric(strsplit(ch.samples, split = "_")[[1]])
                  
                  output[vec.samples] <- vec.Area
                  output
                })
            
            matPvalLess <- t(matrix(matPvalLess, nrow = nSample))
            matPvalLess <- matPvalLess[keepVar, , drop = FALSE]
            
            if (!is.null(fracExp)) {
                Keep2 <- which(apply(matPvalGreater, 1, 
                                        function(x) sum(x < pValLessThres, 
                  na.rm = TRUE)/sum(!is.na(x)) >= fracExp))
            } else {
                Keep2 <- which(matPvalGreater < pValLessThres)
            }
            
            
            if (length(union(Keep1, Keep2)) != 0) {
                if (!is.null(fracExp)) {
                  X <- X[union(Keep1, Keep2), , drop = FALSE]
                  Xbg <- Xbg[union(Keep1, Keep2), , drop = FALSE]
                } else {
                  X[!seq(1, length(X)) %in% union(Keep1, Keep2)] <- NA
                  
                  # filter on samples frenquency
                  if (!is.null(group)) {
                    groupFac <- as.factor(sampleMetadata[, group])
                    test <- apply(X, 1, 
                                    function(x) vapply(levels(groupFac), 
                                                        function(y) sum(!is.na(x[groupFac == 
                      y]))/length(x[groupFac == y]) > fracGroup, TRUE))
                    keepVar <- apply(test, 2, any)
                    
                  } else {
                    keepVar <- apply(X, 1, function(y) sum(!is.na(y))/length(y) > 
                      fracGroup)
                    
                  }
                  
                  X <- X[keepVar, , drop = FALSE]
                  Xbg <- Xbg[keepVar, , drop = FALSE]
                  if (nrow(X) == 0) {
                    warning("peakList is empty, there is no peak presents in 
                            more thant ", 
                      fracGroup, "% of samples")
                    return(NULL)
                  }
                }
                
            } else {
                warning("peakList is empty")
                return(NULL)
            }
            
        }
        
    } else Xbg <- data.frame(row.names = rownames(X))
    
    return(list(X = X, Xbg = Xbg))
}

imputeRaw <- function(X, raw, timeLimit,quanti="ppb") {
    
    missingValues <- which(is.na(X), arr.ind = TRUE)
    indexFilesMissingValues <- unique(missingValues[, "col"])
    namesFilesMissingValues <- colnames(X)[indexFilesMissingValues]
    
    fctFit <- raw@fctFit
    l.shape <-raw@peakShape
    
    primaryIon <- raw@primaryIon
    filesFullName <- raw@name
    
    if (methods::is(filesFullName, "expression")) 
        filesFullName <- eval(filesFullName)
    

    filesFullName.j <- filesFullName
    
    # peak list raw
    peakList<-raw@peakList
    peakListRaw.j <- Biobase::fData(peakList)
    quantiImpute <- list()
    mzMissing <- as.numeric(rownames(missingValues))
    
    
    # open mz Axis
    mzAxis <- raw@mz 
    indexMzList <- lapply(unique(round(mzMissing)), function(m) which(m - 0.6 < mzAxis & 
                                                                          mzAxis < m + 0.6))
    names(indexMzList) <- unique(round(mzMissing))
    indexMz <- Reduce(union, indexMzList)
    
    # get index time
    indexLim <- timeLimit$exp
    indexTime <- Reduce(c, apply(indexLim, 2, function(x) seq(x["start"], x["end"])))
    
    nbExp <- ncol(indexLim)
    
    # open raw data
    #raw <- rhdf5::h5read(filesFullName.j, name = "/FullSpectra/TofData", 
    #                     index = list(indexMz, NULL, NULL, NULL))
    
    rawMn <- raw@rawM 
    # * 0.2 ns / 2.9 (single ion signal) if convert to cps
    
    # information for ppb convertion
    reaction <- try(reaction <- rhdf5::h5read(filesFullName.j, "AddTraces/PTR-Reaction"))
    transmission <- try(rhdf5::h5read(filesFullName.j, "PTR-Transmission"))
    
 
    
    for (m in unique(round(mzMissing))) {
        # exact missing mz
        mz <- mzMissing[round(mzMissing) == m]
        
        # mzAxis around m
        mzAxis.m <- mzAxis[which(m - 0.6 < mzAxis & 
                                     mzAxis < m + 0.6)]
        
        indexExp <- Reduce(c, apply(indexLim, 2, function(x) seq(x[1], x[2])))
        length.exp <- length(indexExp)
        
        quantiMat <- matrix(0, ncol = length(mz), nrow = nbExp)
        
        spectrum <- rowSums(rawMn[ which(m - 0.6 < mzAxis & 
                                             mzAxis < m + 0.6), 
                                  indexExp])/(length.exp * (diff(raw@time[c(1,2)])))
        
        spectrum <- spectrum - snipBase(spectrum)
        
        # substract fitted peak also find
        peakAlsoDetected <- peakListRaw.j[round(peakListRaw.j$Mz) == m, ]
        if (nrow(peakAlsoDetected) != 0) {
            # cumFitPeak
            
            fitPeaks <- apply(peakAlsoDetected, 1, 
                              function(x) eval(parse(text = fctFit))(x["Mz"], x["parameterPeak.delta1"], x["parameterPeak.delta2"], x["parameterPeak.heigh"], 
                                                                     x = mzAxis.m, l.shape = l.shape))
            if (nrow(peakAlsoDetected) > 1) 
                cumFitPeak <- rowSums(fitPeaks) else cumFitPeak <- c(fitPeaks)
                spectrum <- spectrum - cumFitPeak
        }
        
        # fit on the missing values
        resolutionrange<-raw@resolution
        resolution_upper <- resolutionrange[3]
        resolution_mean <- resolutionrange[2]
        resolution_lower <- resolutionrange[1]
        
        n.peak <- length(mz)
        delta <- rep(m/resolution_mean, 2 * n.peak)
        h <- vapply(mz, function(m) max(max(spectrum[which(abs(mzAxis.m - m) < (m * 
                                                                                    50/10^6))]), 1), 0)
        
        
        initMz <- matrix(c(mz, delta, h), nrow = n.peak)
        colnames(initMz) <- c("m", "delta1", "delta2", "h")
        
        
        lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), rep(0.1, 
                                                                                    n.peak)), ncol = 4) - matrix(c(initMz[, "m"]/(resolution_mean * 100), 
                                                                                                                   -initMz[, "m"]/(resolution_lower * 2), -initMz[, "m"]/(resolution_lower * 
                                                                                                                                                                              2), rep(0, n.peak)), ncol = 4)))
        
        upper.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), rep(Inf, 
                                                                                    n.peak)), ncol = 4) + matrix(c(initMz[, "m"]/(resolution_mean * 100), 
                                                                                                                   initMz[, "m"]/(resolution_upper * 2), initMz[, "m"]/(resolution_upper * 
                                                                                                                                                                            2), rep(0, n.peak)), ncol = 4)))
        
        
        fit <- fitPeak(initMz = initMz, sp = spectrum, mz.i = mzAxis.m, lower.cons, 
                       upper.cons, funcName = fctFit, l.shape = l.shape)
        
        fit.peak <- fit$fit.peak
        par_estimated <- fit$par_estimated
        
        
        quanti.m <- apply(par_estimated, 2, function(x) {
            th <- 10 * 0.5 * (log(sqrt(2) + 1)/x[2] + log(sqrt(2) + 1)/x[3])
            mz.x <- mzAxis.m[x[1] - th < mz & mz < x[1] + th]
            sum(eval(parse(text = fctFit))(x[1], x[2], x[3], x[4], mz.x, l.shape), 
                na.rm = TRUE)
        })
        
        list_peak <- cbind(Mz = mz, quanti = quanti.m/(primaryIon * 
                                                           488))
        
        # convert to ppb or ncps if there is reaction ans transmission information
        if (quanti == "ppb") {
            U <- c(reaction$TwData[1, , ])
            Td <- c(reaction$TwData[3, , ])
            pd <- c(reaction$TwData[2, , ])
            quanti.m <- ppbConvert(peakList = list_peak, transmission = transmission$Data, 
                                   U = U[indexExp], Td = Td[indexExp], pd = pd[indexExp])
            
        }
        if (quanti == "ncps") {
            # normalize by primary ions
            quanti.m <- quanti.m/(primaryIon * 488)
        }
        
       
            X[as.character(mz),] <- quanti.m
        
        
        
    }
    return(X)
}

imputeFunc <- function(file, missingValues, eSet, ptrSet) {
    
    fctFit <- getParameters(ptrSet)$detectPeakParam$fctFit
    if (is.null(fctFit)) 
        fctFit <- getPeaksInfo(ptrSet)$fctFit[[file]]
    
    l.shape <- getPeaksInfo(ptrSet)$peakShape[[file]]
    primaryIon <- getPTRInfo(ptrSet)$primaryIon
    filesFullName <- getParameters(ptrSet)$listFile
    
    if (methods::is(filesFullName, "expression")) 
        filesFullName <- eval(filesFullName)
    
    j <- which(file == colnames(Biobase::exprs(eSet)))
    filesFullName.j <- filesFullName[which(basename(filesFullName) == file)]
    
    
    mzMissing <- as.numeric(rownames(missingValues[missingValues[, "col"] == j, , 
        drop = FALSE]))
    
    # open mz Axis
    mzAxis <- rhdf5::h5read(filesFullName.j, name = "FullSpectra/MassAxis")
    indexMzList <- lapply(unique(round(mzMissing)), function(m) which(m - 0.6 < mzAxis & 
        mzAxis < m + 0.6))
    names(indexMzList) <- unique(round(mzMissing))
    indexMz <- Reduce(union, indexMzList)
    
    # get index time
    indexLim <- getTimeInfo(ptrSet)$timeLimit[[file]]$exp
    indexTime <- Reduce(c, apply(indexLim, 2, function(x) seq(x["start"], x["end"])))
    nbExp <- ncol(indexLim)
    
    # open raw data
    raw <- rhdf5::h5read(filesFullName.j, name = "/FullSpectra/TofData", 
                            index = list(indexMz, NULL, NULL, NULL))
    
    rawMn <- matrix(raw, nrow = dim(raw)[1], ncol = prod(utils::tail(dim(raw), 2)))
    # * 0.2 ns / 2.9 (single ion signal) if convert to cps
    
    # information for ppb convertion
    reaction <- try(reaction <- rhdf5::h5read(filesFullName.j, "AddTraces/PTR-Reaction"))
    transmission <- try(rhdf5::h5read(filesFullName.j, "PTR-Transmission"))
    
    # calibrate mass axis
    FirstcalibCoef <- try(rhdf5::h5read(filesFullName.j, "FullSpectra/MassCalibration", 
        index = list(NULL, 1)))
    attributCalib <- try(rhdf5::h5readAttributes(filesFullName.j, "/FullSpectra"))
    if (!is.null(attr(FirstcalibCoef, "condition")) & is.null(attr(attributCalib, "condition"))) {
        FirstcalibCoef <- matrix(c(attributCalib$`MassCalibration a`, attributCalib$`MassCalibration b`), 
                            ncol = 1)
    }
    tof <- sqrt(mzAxis) * FirstcalibCoef[1, 1] + FirstcalibCoef[2, 1]
    coefCalib <- getCalibrationInfo(ptrSet)$coefCalib[[file]][[1]]
    mzAxis <- ((tof - coefCalib["b", ])/coefCalib["a", ])^2
    
    # peak list raw
    peakList<-getPeakList(ptrSet)
    peakListRaw.j <- Biobase::fData(peakList[[file]])
    quantiImpute <- list()
    
    for (m in unique(round(mzMissing))) {
        # exact missing mz
        mz <- mzMissing[round(mzMissing) == m]
        
        # mzAxis around m
        mzAxis.m <- mzAxis[indexMzList[[as.character(m)]]]
        
        indexExp <- Reduce(c, apply(indexLim, 2, function(x) seq(x[1], x[2])))
        length.exp <- length(indexExp)
        
        quantiMat <- matrix(0, ncol = length(mz), nrow = nbExp)
        
        spectrum <- rowSums(rawMn[which(indexMz %in% indexMzList[[as.character(m)]]), 
            indexExp])/(length.exp * (diff(as.numeric(names(getTimeInfo(ptrSet)$TIC[[file]]))[c(1, 
            2)])))
        
        spectrum <- spectrum - snipBase(spectrum)
        
        # substract fitted peak also find
        peakAlsoDetected <- peakListRaw.j[round(peakListRaw.j$Mz) == m, ]
        if (nrow(peakAlsoDetected) != 0) {
            # cumFitPeak
            
            fitPeaks <- apply(peakAlsoDetected, 1, 
                                function(x) eval(parse(text = fctFit))(x["Mz"], x["parameterPeak.delta1"], x["parameterPeak.delta2"], x["parameterPeak.height"], 
                x = mzAxis.m, l.shape = l.shape))
            if (nrow(peakAlsoDetected) > 1) 
                cumFitPeak <- rowSums(fitPeaks) else cumFitPeak <- c(fitPeaks)
            spectrum <- spectrum - cumFitPeak
        }
        
        # fit on the missing values
        resolution_upper <- 8000
        resolution_mean <- 5000
        resolution_lower <- 3000
        
        n.peak <- length(mz)
        delta <- rep(m/resolution_mean, 2 * n.peak)
        h <- vapply(mz, function(m) max(max(spectrum[which(abs(mzAxis.m - m) < (m * 
            50/10^6))]), 1), 0)
        
        
        initMz <- matrix(c(mz, delta, h), nrow = n.peak)
        colnames(initMz) <- c("m", "delta1", "delta2", "h")
        
        
        lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), rep(0.1, 
            n.peak)), ncol = 4) - matrix(c(initMz[, "m"]/(resolution_mean * 100), 
            -initMz[, "m"]/(resolution_lower * 2), -initMz[, "m"]/(resolution_lower * 
                2), rep(0, n.peak)), ncol = 4)))
        
        upper.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), rep(Inf, 
            n.peak)), ncol = 4) + matrix(c(initMz[, "m"]/(resolution_mean * 100), 
            initMz[, "m"]/(resolution_upper * 2), initMz[, "m"]/(resolution_upper * 
                2), rep(0, n.peak)), ncol = 4)))
        
        
        fit <- fitPeak(initMz = initMz, sp = spectrum, mz.i = mzAxis.m, lower.cons, 
            upper.cons, funcName = fctFit, l.shape = l.shape)
        
        fit.peak <- fit$fit.peak
        par_estimated <- fit$par_estimated
        
        quanti.m <- apply(par_estimated, 2, function(x) {
            th <- 10 * 0.5 * (log(sqrt(2) + 1)/x[2] + log(sqrt(2) + 1)/x[3])
            mz.x <- mzAxis.m[x[1] - th < mz & mz < x[1] + th]
            sum(eval(parse(text = fctFit))(x[1], x[2], x[3], x[4], mz.x, l.shape), 
                na.rm = TRUE)
        })
        
        list_peak <- cbind(Mz = mz, quanti = quanti.m/(primaryIon[[file]]$primaryIon * 
            488))
        
        # convert to ppb or ncps if there is reaction ans transmission information
        if (Biobase::annotation(eSet) == "ppb") {
            U <- c(reaction$TwData[1, , ])
            Td <- c(reaction$TwData[3, , ])
            pd <- c(reaction$TwData[2, , ])
            quanti.m <- ppbConvert(peakList = list_peak, transmission = transmission$Data, 
                U = U[indexExp], Td = Td[indexExp], pd = pd[indexExp])
            
        }
        if (Biobase::annotation(eSet) == "ncps") {
            # normalize by primary ions
            quanti.m <- quanti.m/(primaryIon[[basename(file)]]$primaryIon * 488)
        }
        
        for (k in seq_along(mz)) {
            quantiImpute[[as.character(mz)[k]]] <- quanti.m[k]
        }
        
        
    }
    message(basename(file), " done")
    return(quantiImpute)
}

#' Imputes the missing values
#'
#' Imputes missing values by returning back to the raw data and fitting the 
#' peak shape function on the noise (or on the residuals signals if peaks have 
#' already been detected). 
#' @param eSet ExpressionSet returned by the \code{\link[ptairMS]{alignSamples}} 
#' function 
#' @param ptrSet \code{\link[ptairMS]{ptrSet-class}} object processed by the 
#' \code{\link[ptairMS]{detectPeak}} function
#' @param parallelize boolean. If \code{TRUE}, the loop over all files will be 
#' parallelized
#' @param nbCores number of clusters to use for parallel computation.
#' @return the same ExpressionSet as in input, with the imputed missing values 
#' in the assayData slot
#' @export 
#' @examples
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, 
#' setName="exhaledPtrset", mzCalibRef = c(21.022, 60.0525),
#' fracMaxTIC = 0.7, saveDir = NULL )
#' exhaledPtrset<-detectPeak(exhaledPtrset,mz=c(21,52))
#' eSet <- alignSamples(exhaledPtrset,pValGreaterThres=0.05,fracGroup=0)
#' Biobase::exprs(eSet)
#' eSet <- impute(eSet,exhaledPtrset)
#' Biobase::exprs(eSet)

impute <- function(eSet, ptrSet, parallelize = FALSE, nbCores = 2) {
    
    # get peak list peakList <- getPeakList(ptrSet)$aligned
    
    # get index of missing values
    missingValues <- which(is.na(Biobase::exprs(eSet)), arr.ind = TRUE)
    indexFilesMissingValues <- unique(missingValues[, "col"])
    namesFilesMissingValues <- colnames(Biobase::exprs(eSet))[indexFilesMissingValues]
    
    # get files full names in ptrSet object
    FUN <- function(file) {
        test <- try(imputeFunc(file = file, missingValues = missingValues, 
                                eSet = eSet, 
            ptrSet = ptrSet))
        if (!is.null(attr(test, "condition"))) {
            return(list(NULL))
        } else return(test)
        
    }
    
    if (parallelize) {
        cl <- parallel::makeCluster(nbCores)
        doParallel::registerDoParallel(cl)
        `%dopar%` <- foreach::`%dopar%`
        quantiMissing <- foreach::foreach(file = namesFilesMissingValues, 
                                            .packages = c("data.table")) %dopar% 
            {
                FUN(file)
            }
        parallel::stopCluster(cl)
    } else quantiMissing <- lapply(namesFilesMissingValues, FUN)
    
    names(quantiMissing) <- namesFilesMissingValues
    
    # add to peak table
    for (file in names(quantiMissing)) {
        for (mz in names(quantiMissing[[file]])) {
            Biobase::exprs(eSet)[mz, file] <- quantiMissing[[file]][[mz]]
        }
    }
    return(eSet)
}

#' Impute missing values on an matrix set from an ptrSet
#'
#' Imputing missing values by returning back to the raw data and fitting the 
#' peak shape function on the noise / residuals
#' @param X the peak table matrix with missing values
#' @param ptrSet processed by detectPeak function
#' @param quantiUnit the unit of the quantities in the matrix \code{X} (ppb, 
#' cps or ncps)
#' @return the same matrix as in input, with missing values imputing
#' @export 
#' @examples
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, 
#' setName="exhaledPtrset", mzCalibRef = c(21.022, 60.0525),
#' fracMaxTIC = 0.7, saveDir = NULL )
#' exhaledPtrset<-detectPeak(exhaledPtrset,mz=c(21,52))
#' eSet <- alignSamples(exhaledPtrset,pValGreaterThres=0.05,fracGroup=0)
#' X <-Biobase::exprs(eSet)
#' X <- imputeMat(X,exhaledPtrset,quantiUnit='ppb')
#' plotFeatures(exhaledPtrset,mz = 52.047,typePlot = "ggplot",colorBy = "subfolder")
imputeMat <- function(X, ptrSet, quantiUnit) {
    # peak func(ion)
    fctFit <- getPeaksInfo(ptrSet)$fctFit
    
    # get index of missing values
    missingValues <- which(is.na(X), arr.ind = TRUE)
    indexFilesMissingValues <- unique(missingValues[, "col"])
    namesFilesMissingValues <- colnames(X)[indexFilesMissingValues]
    
    # get primry ion quantity
    primaryIon <- getPTRInfo(ptrSet)$primaryIon
    
    # get files full names in ptrSet object
    filesFullName <- getParameters(ptrSet)$listFile
    if (methods::is(filesFullName, "expression")) 
        filesFullName <- eval(filesFullName)
    
    for (file in namesFilesMissingValues) {
        j <- which(file == colnames(X))
        filesFullName.j <- filesFullName[which(basename(filesFullName) == file)]
        
        mzMissing <- as.numeric(rownames(missingValues[missingValues[, "col"] == 
            j, , drop = FALSE]))
        
        # open mz Axis
        mzAxis <- rhdf5::h5read(filesFullName.j, name = "FullSpectra/MassAxis")
        indexMzList <- lapply(unique(round(mzMissing)), function(m) which(m - 0.6 < 
            mzAxis & mzAxis < m + 0.6))
        
        
        names(indexMzList) <- unique(round(mzMissing))
        indexMz <- Reduce(union, indexMzList)
        
        # get index time
        indexLim <- getTimeInfo(ptrSet)$timeLimit[[file]]$exp
        indexTime <- Reduce(c, apply(indexLim, 2, function(x) seq(x["start"], 
                                                                  x["end"])))
        nbExp <- ncol(indexLim)
        
        # open raw data
        raw <- rhdf5::h5read(filesFullName.j, name = "/FullSpectra/TofData", 
                                index = list(indexMz, 
            NULL, NULL, NULL))
        
        # * 0.2 ns / 2.9 (single ion signal) if convert to cps
        rawMn <- matrix(raw, nrow = dim(raw)[1], ncol = prod(utils::tail(dim(raw), 
            2)))
        
        # information for ppb convertion
        reaction <- try(reaction <- rhdf5::h5read(filesFullName.j, "AddTraces/PTR-Reaction"))
        transmission <- try(rhdf5::h5read(filesFullName.j, "PTR-Transmission"))
        
        # calibrate mass axis
        FirstcalibCoef <- rhdf5::h5read(filesFullName.j, "FullSpectra/MassCalibration", 
            index = list(NULL, 1))
        tof <- sqrt(mzAxis) * FirstcalibCoef[1, 1] + FirstcalibCoef[2, 1]
        coefCalib <- getCalibrationInfo(ptrSet)$coefCalib[[file]]
        mzAxis <- ((tof - coefCalib[[1]]["b", ])/coefCalib[[1]]["a", ])^2
        
        # peak list raw
        peakList<-getPeakList(ptrSet)
        peakListRaw.j <- Biobase::fData(peakList[[file]])
        
        for (m in unique(round(mzMissing))) {
            # exact missing mz
            mz <- mzMissing[round(mzMissing) == m]
            
            # mzAxis around m
            mzAxis.m <- mzAxis[indexMzList[[as.character(m)]]]
            
            indexExp <- Reduce(c, apply(indexLim, 2, function(x) seq(x[1], x[2])))
            length.exp <- length(indexExp)
            quantiMat <- matrix(0, ncol = length(mz), nrow = nbExp)
            
            spectrum <- rowSums(rawMn[which(indexMz %in% indexMzList[[as.character(m)]]), 
                indexExp])/(length.exp * (diff(as.numeric(names(getTimeInfo(ptrSet)$TIC[[file]]))[c(1, 
                2)])))
            
            spectrum <- spectrum - snipBase(spectrum)
            
            # substract fitted peak also find
            peakAlsoDetected <- peakListRaw.j[round(peakListRaw.j$Mz) == m, ]
            if (nrow(peakAlsoDetected) != 0) {
                # cumFitPeak
                funEVAL<-function(x) eval(parse(text = fctFit[[file]]))(x["Mz"], 
                                                                        x["parameterPeak.delta1"], 
                                                                        x["parameterPeak.delta2"], 
                                                                        x["parameterPeak.height"], 
                                                                        x = mzAxis.m, 
                                                                        l.shape = getPeaksInfo(ptrSet)$peakShape[[file]])
                fitPeaks <- apply(peakAlsoDetected, 1, funEVAL )
                if(all(is.na(fitPeaks))){
                    resolution_upper <- ceiling (max(ptrSet@resolution[[file]])/1000)*1000 +1000
                    resolution_mean <- round(mean(ptrSet@resolution[[file]])/100)*100
                    resolution_lower <- floor(min(ptrSet@resolution[[file]])/1000)*1000-500
                    
                    
                    
                    n.peak <- nrow(peakAlsoDetected)
                    mass<-peakAlsoDetected$Mz
                    delta <- rep(m/resolution_mean, 2 * n.peak)
                    h <- vapply(mass, function(m) max(max(spectrum[which(abs(mzAxis.m - m) < 
                                                                           (m * 50/10^6))]), 1), 0)
                    
                    
                    initMz <- matrix(c(mass, delta/2, h), nrow = n.peak)
                    colnames(initMz) <- c("m", "delta1", "delta2", "h")
                    
                    lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), 
                                                        rep(0.1, n.peak)), ncol = 4) - 
                                          matrix(c(initMz[, "m"]/(resolution_mean * 100), 
                                                   -initMz[, "m"]/(resolution_upper*2), 
                                                   -initMz[, "m"]/(resolution_upper*2), 
                                                   rep(0, n.peak)), 
                                                 ncol = 4)))
                    
                    upper.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), 
                                                        rep(Inf, n.peak)), ncol = 4) + 
                                          matrix(c(initMz[, "m"]/(resolution_mean *100), 
                                                   initMz[, "m"]/(resolution_lower*2), 
                                                   initMz[, "m"]/(resolution_lower*2), 
                                                   rep(0, n.peak)), ncol = 4)))
                    fit <- fitPeak(initMz, spectrum, mzAxis.m, lower.cons, upper.cons, fctFit[[file]], 
                                   l.shape = getPeaksInfo(ptrSet)$peakShape[[file]])
                    
                    fitPeaks <- fit$fit.peak
                    
                }
               
                if (length(dim(fitPeaks)) > 1 ) 
                    cumFitPeak <- rowSums(fitPeaks) else cumFitPeak <- c(fitPeaks)
                spectrum <- spectrum - cumFitPeak
                
                spectrum[spectrum<0]<-0
            }
            
            # fit on the missing values
            resolution_upper <- ceiling (max(ptrSet@resolution[[file]])/1000)*1000
            resolution_mean <- round(mean(ptrSet@resolution[[file]])/100)*100
            resolution_lower <- floor(min(ptrSet@resolution[[file]])/1000)*1000
            
            
            n.peak <- length(mz)
            delta <- rep(m/resolution_mean, 2 * n.peak)
            h <- vapply(mz, function(m) max(max(spectrum[which(abs(mzAxis.m - m) < 
                (m * 50/10^6))]), 1), 0)
            
            
            initMz <- matrix(c(mz, delta/2, h), nrow = n.peak)
            colnames(initMz) <- c("m", "delta1", "delta2", "h")
            
            lower.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), 
                rep(0.1, n.peak)), ncol = 4) - 
                    matrix(c(initMz[, "m"]/(resolution_mean * 100), 
                -initMz[, "m"]/(resolution_upper*2), 
                -initMz[, "m"]/(resolution_upper*2), 
                rep(0, n.peak)), 
                ncol = 4)))
            
            upper.cons <- c(t(initMz * matrix(c(rep(1, n.peak), rep(0, n.peak * 2), 
                rep(Inf, n.peak)), ncol = 4) + 
                    matrix(c(initMz[, "m"]/(resolution_mean *100), 
                             initMz[, "m"]/(resolution_lower*2), 
                             initMz[, "m"]/(resolution_lower*2), 
                rep(0, n.peak)), ncol = 4)))
            
            
            fit <- fitPeak(initMz, spectrum, mzAxis.m, lower.cons, upper.cons, fctFit[[file]], 
                l.shape = getPeaksInfo(ptrSet)$peakShape[[file]])
            
            fit.peak <- fit$fit.peak
            par_estimated <- fit$par_estimated
            
            quanti.m <- apply(par_estimated, 2, function(x) {
                th <- 10 * 0.5 * (log(sqrt(2) + 1)/x[2] + log(sqrt(2) + 1)/x[3])
                mz.x <- mzAxis.m[x[1] - th < mz & mz < x[1] + th]
                sum(eval(parse(text = fctFit))(x[1], x[2], x[3], x[4], mz.x, 
                                                l.shape = getPeaksInfo(ptrSet)$peakShape[[file]]), 
                  na.rm = TRUE)
            })
            
            list_peak <- cbind(Mz = mz, quanti = quanti.m/(primaryIon[[file]]$primaryIon * 
                488))
            
            # convert to ppb or ncps if there is reaction ans transmission information
            if (quantiUnit == "ppb") {
                U <- c(reaction$TwData[1, , ])
                Td <- c(reaction$TwData[3, , ])
                pd <- c(reaction$TwData[2, , ])
                quanti.m <- ppbConvert(peakList = list_peak, 
                                        transmission = transmission$Data, 
                  U = U[indexExp], Td = Td[indexExp], pd = pd[indexExp])
                
            }
            if (quantiUnit == "ncps") {
                # normalize by primary ions
                quanti.m <- quanti.m/(primaryIon[[basename(file)]]$primaryIon * 488)
            }
            
            
            
            # add to peak table
            X[as.character(mz), file] <- quanti.m
        }
        message(basename(file), " done")
    }  #end for file
    return(X)
}

