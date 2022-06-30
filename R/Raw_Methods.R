## Calibration ----
#' Calibrates the mass axis with references masses
#' 
#' To convert Time Of Flight (TOF) axis to mass axis, we use the formula:
#' mz = ((tof-b)/a )^2 (Muller et al. 2013) To estimate those 
#' parameters, references peaks with accurate know  masses and without 
#' overlapping peak are needed. The best is that the references masses covers a 
#' maximum of the mass range.
#' @param x a \code{\link[ptairMS]{ptrRaw-class}} or 
#' \code{\link[ptairMS]{ptrSet-class}} object
#' @param mzCalibRef Vector of accurate mass values of intensive peaks and 
#' 'unique' in a nominal mass interval (without overlapping)
#' @param calibrationPeriod in second, coefficient calibration are estimated for 
#' each sum spectrum of \code{calibrationPeriod} seconds
#' @param tol the maximum error tolerated in ppm. If more than \code{tol} 
#' warnings. 
#' @return the same ptrRaw or ptrSet as in input, with the following modified 
#' element:
#' \itemize{
#' \item mz: the new mz axis calibrated
#' \item rawM: same raw matrix with the new mz axis in rownames
#' \item calibMassRef: reference masses used for the calibration
#' \item calibMzToTof and calibTofToMz: function to convert TOF to mz
#' \item calibError: the calibration error to the reference masses in ppm
#' \item calibrationIndex: index time of each calibration period
#' } 
#' @examples 
#' 
#' ### ptrRaw object 
#' 
#' library(ptairData)
#' filePath <- system.file('extdata/exhaledAir/ind1', 'ind1-1.h5', 
#' package = 'ptairData')
#' raw <- readRaw(filePath, calib = FALSE)
#' rawCalibrated <- calibration(raw)
#' @rdname calibration
#' @export 
methods::setMethod(f = "calibration", signature = "ptrRaw", function(x, mzCalibRef = c(21.022, 
    29.013424, 41.03858, 59.049141, 75.04406, 203.943, 330.8495), calibrationPeriod = 60, 
    tol = 70) {
    object <- x
    # get mz axis and average spectrum
    mz <- getRawInfo(object)$mz
    time <- c(getRawInfo(object)$time)
    sp <- rowSums(getRawInfo(object)$rawM)/(dim(getRawInfo(object)$rawM)[2] * (time[3] - time[2]))
    width.window <- 0.4
    # check if mzCalibRef are in mz
    mzCalibRef<-mzCalibRef[mzCalibRef>0]
    outMz <- which(vapply(mzCalibRef, function(x) !any(round(x) - width.window < 
        mz & mz < round(x) + width.window), FUN.VALUE = TRUE))
    if (length(outMz) != 0) {
        message(paste(mzCalibRef[outMz], collapse = " "), "excluded, not contains in the mass axis \n")
        mzCalibRef <- mzCalibRef[-outMz]
    }
    # test if there is a only one peak on the TIS
    nLocalMax <- vapply(mzCalibRef, function(x) {
        spx <- sp[x - 0.4 < mz & mz < x + 0.4]
        length(LocalMaximaSG(sp = spx, minPeakHeight = 0.2 * max(spx)))
    }, FUN.VALUE = 0)
    badMass <- which(nLocalMax != 1)
    if (length(badMass) != 0) {
        mzCalibRef <- mzCalibRef[-badMass]
    }
    object@calibMassRef <- mzCalibRef
    # determine average peak shape on calibration masses
    peakShape <- determinePeakShape(raw = object)$peakShapetof
    object<-setPeakShape(object,peakShape) ####
    
    # performs calibration every steps second
    calib_List <- list(NULL)
    indexTime <- round(diff(time)[1], 3)
    nbIndex <- round(calibrationPeriod/indexTime)
    for (i in seq(0, (floor(length(time)/nbIndex) - 1))) {
        if (i >= 0) {
            index <- seq(from = i * nbIndex + 1, to = min((i * nbIndex + nbIndex), 
                length(getRawInfo(object)$time)))
            if (i == (floor(length(time)/nbIndex) - 1) & utils::tail(index, 1) < 
                length(time)) 
                index <- c(index, seq(utils::tail(index, 1) + 1, length(getRawInfo(object)$time)))
            sp.i <- rowSums(getRawInfo(object)$rawM[, index])
            calib_List[[i + 1]] <- c(calibrationFun(sp.i, mz, mzCalibRef, 
                                                    calibCoef = getCalibrationInfo(object)$calibCoef[[1]], 
                peakShape, tol), list(index = index))
        }
    }
    # use the mz axis of the first calibration
    mzVnbis <- calib_List[[1]]$mzVnbis
    object@mz <- mzVnbis
    rownames(object@rawM) <- mzVnbis
    # update object
    matrixEror <- Reduce(rbind, lapply(calib_List, function(x) x$error))
    if (length(calib_List) > 1) 
        matrixEror <- apply(matrixEror, 2, mean)
    
    calibrationInfo<-list(indexTimeCalib=lapply(calib_List, function(x) x$index),
                          calibSpectr=lapply(calib_List, function(x) x$calibSpectr),
                          calibError=matrixEror,
                          calibCoef=lapply(calib_List, function(x) x$coefs))
    object<-setCalibration(object,calibrationInfo)
    
    
    return(object)
})

#'calibration function
#'
#'Performs calibration on sp with mzCalibRef reference masses and mzToTofFunc 
#'as previous 
#'calibration function
#' @param sp spectrum
#' @param mz mass axis 
#' @param mzCalibRef masses of know reference peaks
#' @param calibCoef coeficient of the previous calibration
#' @param peakShape a list with reference axis and a reference peak shape 
#' centered in zero
#' @param tol maximum error tolarated in ppm
#' @return list 
#' @keywords internal
calibrationFun <- function(sp, mz, mzCalibRef, calibCoef, peakShape, tol) {

    width.window <- 0.4
    mzToTofFunc <- function(mz) mzToTof(mz, calibCoef)
    # calculate tof axis
    if (is.null(mzToTofFunc(1))) 
        tof <- seq(0, length(mz) - 1) else tof <- mzToTofFunc(mz)
    # spectrum of mass calib
    spTronc <- lapply(mzCalibRef, function(m) {
        index <- which((mz <= m + width.window) & (mz >= m - width.window))
        tof <- tof[index]
        noise <- sp[which(((m - 1 + width.window) < mz) & (mz < m - width.window))]
        mz <- mz[index]
        signal <- sp[index]
        return(list(signal = signal, mz = mz, tof = tof, noise = noise))
    })
    # test if there is a only one peak nLocalMax <-vapply(spTronc, function(x) {
    # length( LocalMaximaSG( sp = x$signal, minPeakHeight = 0.2*max(x$signal)) ) },
    # FUN.VALUE = 0 ) badMass <- which(nLocalMax !=1 ) if(length(badMass)!=0){
    # mzCalibRef <- mzCalibRef[-badMass] spTronc <- spTronc[-badMass] }
    if (length(mzCalibRef) < 2) 
        stop("To few references masses for calibration")
    # calculate the tof of reference masses with peak shape
    tofMax <- vapply(spTronc, function(x) {
        tofrange <- x$tof
        sp <- x$signal
        sp <- sp - snipBase(sp)
        acf <- stats::acf(x$noise, lag.max = 1, plot = FALSE)[1]$acf
        localMax <- LocalMaximaSG(sp = sp, minPeakHeight = max(sp) * 0.2, noiseacf = min(acf, 
            0.3))
        tcenter <- tofrange[localMax[which.max(sp[localMax])]]
        initTof <- matrix(c(tcenter, tcenter/10000, max(sp)), nrow = 1)
        largerfit <- initTof[2]
        fit <- fit_averagePeak(initTof, l.shape = peakShape, sp = sp[tcenter - largerfit < 
            tofrange & tofrange < tcenter + largerfit], bin = tofrange[tcenter - 
            largerfit < tofrange & tofrange < tcenter + largerfit], lower.cons = NULL, 
            upper.cons = NULL)
        tofMax <- fit$par_estimated[1, ]
        return(tofMax)
    }, FUN.VALUE = 0.1)
    # re estimated calibration coefficient with reference masses
    regression <- stats::nls(rep(1, length(mzCalibRef)) ~ I(((tofMax - b)/a)^2/mzCalibRef), 
        start = list(a = calibCoef["a", ], b = calibCoef["b", ]), algorithm = "port")
    coefs <- stats::coefficients(regression)
    coefs <- as.matrix(coefs)
    # the new mass axis calibrated
    mzVnbis <- tofToMz(tof, coefs)
    # the new position of reference masses
    mzRefRaw <- tofToMz(tofMax, coefs)
    # error in ppm
    error <- abs(mzRefRaw - mzCalibRef) * 10^6/mzCalibRef
    names(error) <- mzCalibRef
    # spectrum of calibration masses with new mass axis
    calibSpectr <- lapply(mzCalibRef, function(m) {
        index <- which((mzVnbis <= m + width.window) & (mzVnbis >= m - width.window))
        mz <- mzVnbis[index]
        signal <- sp[index]
        return(list(signal = signal, mz = mz))
    })
    if (any(abs(error) > tol)) 
        message("error greater than ", tol, "\n")
    return(list(mzVnbis = mzVnbis, mzCalibRef = mzCalibRef, calibSpectr = calibSpectr, 
        error = error, coefs = coefs))
}
alignCalibrationPeak <- function(calibSpectr, calibMassRef, ntimes) {
    lapply(seq(1, length(calibMassRef)), function(m) {
        mz1 <- calibSpectr[[1]][[m]]$mz
        signal <- Reduce(rbind, lapply(calibSpectr, function(x) stats::spline(x = x[[m]]$mz, 
            y = x[[m]]$signal, xout = mz1)$y))
        if (!is.null(nrow(signal))) {
            signal <- apply(signal, 2, function(x) sum(x)/ntimes)
        } else signal <- signal/ntimes
        list(mz = mz1, signal = signal)
    })
}
tofToMz <- function(tof, calibCoef) {
    ((tof - calibCoef["b", ])/calibCoef["a", ])^2
}
mzToTof <- function(m, calibCoef) {
    sqrt(m) * calibCoef["a", ] + calibCoef["b", ]
}
## plotRaw ----
#' @rdname plotRaw
#' @export
methods::setMethod(f = "plotRaw", signature = "ptrRaw", function(object, mzRange, 
    timeRange = c(NA, NA), type = c("classical", "plotly")[1], ppm = 2000, palette = c("heat", 
        "revHeat", "grey", "revGrey", "ramp")[1], showVocDB = TRUE, 
    figure.pdf = "interactive", ...) {
    mzVn <- getRawInfo(object)$mz
    timeVn <- getRawInfo(object)$time
    timeiNTER <- (timeVn[3] - timeVn[2])
    rawMN <- getRawInfo(object)$rawM
    if (length(mzRange) == 1) 
        mzRange <- mzRange + c(-1, 1) * mzRange * ppm * 1e-06
    mzRange[1] <- max(min(mzVn), mzRange[1], na.rm = TRUE)
    mzRange[2] <- min(max(mzVn), mzRange[2], na.rm = TRUE)
    timeRange[1] <- max(min(timeVn), timeRange[1], na.rm = TRUE)
    timeRange[2] <- min(max(timeVn), timeRange[2], na.rm = TRUE)
    mzVi <- which(mzVn > mzRange[1] & mzVn < mzRange[2])
    timeVi <- which(timeVn > timeRange[1] & timeVn < timeRange[2])
    rawSubMN <- rawMN[mzVi, timeVi]
    mzVn <- mzVn[mzVi]
    timeVn <- timeVn[timeVi]
    if (showVocDB) {
        vocdbDF <- .loadVocDB()
        vocdb_sel.vl <- vocdbDF[, "ion_mass"] >= mzRange[1] & vocdbDF[, "ion_mass"] <= 
            mzRange[2]
        if (sum(vocdb_sel.vl)) {
            vocdbDF <- vocdbDF[vocdb_sel.vl, , drop = FALSE]
        } else vocdbDF <- NULL
    } else vocdbDF <- NULL
    if (figure.pdf != "interactive") {
        if (type == "plotly") 
            stop("'plotly display is only available in the 'interactive' 
                     mode currently.", 
                call. = FALSE)
        filenameSplitVc <- unlist(strsplit(basename(figure.pdf), ".", fixed = TRUE))
        extC <- utils::tail(filenameSplitVc, 1)
        if (extC == "pdf") {
            grDevices::pdf(figure.pdf)
        } else stop("The extension of the 'figure.pdf' filename argument 
                     should be 'pdf'", 
            call. = FALSE)
    }
    switch(type, classical = {
        imageMN <- t(rawSubMN)[, seq_len(nrow(rawSubMN)), drop = FALSE]
        rownames(imageMN) <- round(as.numeric(rownames(imageMN)))
        colnames(imageMN) <- round(as.numeric(colnames(imageMN)), 4)
        paletteVc <- .palette(palette = palette)
        currentParLs <- graphics::par()
        for (parC in c("cin", "cra", "csi", "cxy", "din", "page")) currentParLs[[parC]] <- NULL
        marLs <- list(chr = c(0.6, 4.1, 1.1, 0), sca = c(0.6, 0.6, 1.1, 7.1), ima = c(4.1, 
            4.1, 0, 0), spe = c(4.1, 0.6, 0, 1.1))
        graphics::par(font = 2, font.axis = 2, font.lab = 2, pch = 18)
        graphics::layout(matrix(c(1, 2, 3, 4), byrow = TRUE, nrow = 2), heights = c(2, 
            5), widths = c(5, 2))
        ## chr: Current
        graphics::par(mar = marLs[["chr"]])
        chromVn <- apply(rawSubMN, 2, function(intVn) mean(intVn, na.rm = TRUE))
        plot(as.numeric(names(chromVn)), chromVn, cex = 0.7, pch = 16, xlab = "", 
            ylab = "", xaxs = "i", xaxt = "n", yaxs = "i")
        graphics::mtext("Mean of intensity", cex = 0.8, side = 2, line = 2.5)
        ## sca: Color scale
        graphics::par(mar = marLs[["sca"]])
        .drawScale(imageMN = imageMN, paletteVc = paletteVc)
        ## ima: Image
        graphics::par(mar = marLs[["ima"]])
        .drawImage(imageMN = imageMN, paletteVc = paletteVc)
        if (showVocDB && !is.null(vocdbDF)) {
            mzImaVn <- as.numeric(colnames(imageMN))
            graphics::abline(h = vapply(vocdbDF[, "ion_mass"], function(mzN) (mzN - 
                min(mzImaVn))/diff(range(mzImaVn)) * ncol(imageMN) + par("usr")[1], 
                FUN.VALUE = 1.1), lty = "dotted")
        }
        ## spe: Spectrum
        graphics::par(mar = marLs[["spe"]])
        specVn <- apply(rawSubMN, 1, function(intVn) mean(intVn, na.rm = TRUE)/timeiNTER)
        plot(specVn, as.numeric(names(specVn)), cex = 0.7, pch = 16, xlab = "", ylab = "", 
            xaxs = "i", yaxs = "i", yaxt = "n")
        graphics::mtext("Count Per Second (cps)", cex = 0.8, side = 1, line = 2.5)
        if (showVocDB && !is.null(vocdbDF)) {
            graphics::abline(h = vocdbDF[, "ion_mass"], lty = "dotted")
        }
        graphics::par(currentParLs)
        if (figure.pdf != "interactive") dev.off()
    }, plotly = {
        plotlyBuild <- function() {
            p <- plotly::subplot(plotly::plot_ly(x = as.numeric(colnames(rawSubMN)), 
                y = colSums(rawSubMN), type = "scatter", mode = "lines"), plotly::plotly_empty(), 
                plotly::plot_ly(z = rawSubMN, x = as.numeric(colnames(rawSubMN)), 
                  y = as.numeric(rownames(rawSubMN)), type = "heatmap"), plotly::plot_ly(x = rowSums(rawSubMN), 
                  y = as.numeric(rownames(rawSubMN)), type = "scatter", mode = "lines"), 
                nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), shareX = TRUE, 
                shareY = TRUE, titleX = FALSE, titleY = FALSE)
            return(p)
        }
        p <- suppressWarnings(suppressMessages(plotlyBuild()))
        p <- plotly::layout(p, showlegend = FALSE)
        print(plotly::ggplotly(p))
    })
    if (showVocDB & !is.null(vocdbDF)) {
        vocdbDF <- vocdbDF[nrow(vocdbDF):1, , drop = FALSE]
        print(vocdbDF[, c("ion_mass", "ion_formula", "name_iupac"), drop = FALSE])
    }
    return(invisible(list(rawsubM = rawSubMN, vocsubDB = vocdbDF)))
})


.palette <- function(palette) {
    switch(palette, heat = {
        return(rev(grDevices::rainbow(ceiling(256 * 1.5))[seq(1, 256)]))
    }, revHeat = {
        return(grDevices::rainbow(ceiling(256 * 1.5))[seq(1, 256)])
    }, grey = {
        return(grDevices::grey((0:255)/256))
    }, revGrey = {
        return(rev(grDevices::grey((0:255)/256)))
    }, ramp = {
        return((grDevices::colorRampPalette(c("blue", "orange", "red"), space = "rgb"))(256)[seq(1, 
            256)])
    })
}

.drawScale <- function(imageMN, paletteVc) {
    ylimVn <- c(0, 256)
    ybottomVn <- 0:255
    ytopVn <- seq_len(256)
    plot(x = 0, y = 0, font.axis = 2, font.lab = 2, type = "n", xlim = c(0, 1), ylim = ylimVn, 
        xlab = "", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
    graphics::rect(xleft = 0, ybottom = ybottomVn, xright = 1, ytop = ytopVn, col = paletteVc, 
        border = NA)
    graphics::axis(at = .prettyAxis(range(imageMN, na.rm = TRUE), 256)$atVn, font = 2, 
        font.axis = 2, labels = .prettyAxis(range(imageMN, na.rm = TRUE), 256)$labelVn, 
        las = 1, lwd = 2, lwd.ticks = 2, side = 4, xpd = TRUE)
    graphics::arrows(graphics::par("usr")[2], graphics::par("usr")[4], graphics::par("usr")[2], 
        graphics::par("usr")[3], code = 0, lwd = 2, xpd = TRUE)
    graphics::arrows(graphics::par("usr")[1], graphics::par("usr")[4], graphics::par("usr")[1], 
        graphics::par("usr")[3], code = 0, lwd = 2, xpd = TRUE)
    graphics::box(lwd = 2)
}

.prettyAxis <- function(axisValuesVn, opLengthN) {
    if (NA %in% axisValuesVn) {
        warning("NA in axisValuesVn")
        axisValuesVn <- as.vector(stats::na.omit(axisValuesVn))
    }
    if (opLengthN < length(axisValuesVn)) 
        stop("The length of in vector must be inferior to the length of the length 
         parameter.")
    if (length(axisValuesVn) < opLengthN) {
        axisValuesVn <- seq(from = min(axisValuesVn), to = max(axisValuesVn), length.out = opLengthN)
    }
    prettyAxisValues <- pretty(axisValuesVn)
    prettyLabelsVn <- prettyAtVn <- c()
    for (n in seq_along(prettyAxisValues)) if (min(axisValuesVn) < prettyAxisValues[n] && 
        prettyAxisValues[n] < max(axisValuesVn)) {
        prettyLabelsVn <- c(prettyLabelsVn, prettyAxisValues[n])
        prettyAtVn <- c(prettyAtVn, which(abs(axisValuesVn - prettyAxisValues[n]) == 
            min(abs(axisValuesVn - prettyAxisValues[n])))[1])
    }
    prettyAxisLs <- list(atVn = prettyAtVn, labelVn = prettyLabelsVn)
    return(prettyAxisLs)
}

.drawImage <- function(imageMN, paletteVc) {
    graphics::image(x = seq_len(nrow(imageMN)), y = seq_len(ncol(imageMN)), z = imageMN, 
        col = paletteVc, font.axis = 2, font.lab = 2, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "")
    timeVn <- as.numeric(rownames(imageMN))
    prettyTimeVn <- pretty(timeVn)
    prettyTimeVn <- prettyTimeVn[min(timeVn) <= prettyTimeVn & prettyTimeVn <= max(timeVn)]
    prettyTimeVi <- vapply(seq_along(prettyTimeVn), function(k) {
        which(abs(timeVn - prettyTimeVn[k]) == min(abs(timeVn - prettyTimeVn[k])))[1]
    }, FUN.VALUE = 1)
    graphics::axis(side = 1, at = prettyTimeVi, font = 2, labels = prettyTimeVn)
    mzVn <- as.numeric(colnames(imageMN))
    prettyMzVn <- pretty(mzVn)
    prettyMzVn <- prettyMzVn[min(mzVn) <= prettyMzVn & prettyMzVn <= max(mzVn)]
    prettyMzVi <- vapply(seq_along(prettyMzVn), function(k) {
        which(abs(mzVn - prettyMzVn[k]) == min(abs(mzVn - prettyMzVn[k])))[1]
    }, FUN.VALUE = 1)
    graphics::axis(side = 2, at = prettyMzVi, font = 2, labels = prettyMzVn)
    ## xlab
    graphics::mtext(line = 2.5, side = 2, text = "m/z", cex = 0.8)
    ## ylab
    graphics::mtext(line = 2.5, side = 1, text = "time (s)", cex = 0.8)
    ## border
    graphics::box(lwd = 2)
}
## plotCalib ----
#' @rdname plotCalib
#' @export
methods::setMethod(f = "plotCalib", signature = "ptrRaw", function(object, ppm = 2000, 
    ...) {
    raw <- object
    # get mass and specter
    mzCalibRef <- getCalibrationInfo(raw)$calibMassRef
    calibSpectr <- alignCalibrationPeak(getCalibrationInfo(raw)$calibSpectr, mzCalibRef, length(getRawInfo(raw)$time))
    error <- getCalibrationInfo(raw)$calibError
    # plot in a window of 2000 ppm
    nb_plot <- length(mzCalibRef)
    nb_row <- ceiling(sqrt(nb_plot))
    graphics::par(oma = c(0, 0, 3, 0))
    graphics::layout(matrix(seq(1, (nb_row^2)), nrow = nb_row, byrow = TRUE))
    # loop over masses
    for (i in seq_along(mzCalibRef)) {
        m <- mzCalibRef[i]
        th <- m * (ppm/2)/10^6
        mz <- calibSpectr[[i]]$mz
        index <- which(m - th < mz & mz < m + th)
        x <- mz[index]
        y <- calibSpectr[[i]]$signal[index]
        plot(x, y, type = "l", lwd = 2, ylab = "intenisty", xlab = "mz", main = c(m, 
            paste("error:", round(error[i], 2), "ppm")))
        graphics::abline(v = m, col = "red", lwd = 2)
    }
    title(main = getFileNames(raw), outer = TRUE, line = 0.5, cex.main = 2)
    graphics::par(oma = c(0, 0, 0, 0))
    graphics::layout(matrix(1))
})
## plotTIC----
#' @param fracMaxTIC Percentage (between 0 and 1) of the maximum of the Total 
#' Ion Current (TIC) amplitude with baseline removal. We will analyze only the 
#' part of the spectrum where the TIC intensity is higher than 
#' `fracMaxTIC * max(TIC) `. If you want to analyze the entire spectrum, set 
#' this parameter to 0. 
#' @rdname plotTIC
#' @export
#' @examples 
#' ### ptrRaw object
#' 
#' library(ptairData)
#' filePath <- system.file('extdata/exhaledAir/ind1', 'ind1-1.h5', 
#' package = 'ptairData')
#' raw <- readRaw(filePath)
#' p <- plotTIC(raw)
#' p
methods::setMethod(f = "plotTIC", signature = "ptrRaw", function(object, type, baselineRm, 
    showLimits, fracMaxTIC = 0.8, ...) {
    # get the TIC and time limit
    TIC <- colSums(getRawInfo(object)$rawM)
    # remove baseline
    if (baselineRm) {
        bl <- try(snipBase(TIC))
        if (is.null(attr(bl, "condition"))) 
            TIC <- TIC - bl else TIC <- TIC - TIC[1]
    }
    plot <- ggplot2::qplot(x = getRawInfo(object)$time, y = TIC, xlab = "time", ylab = "intensity", 
        main = paste("TIC of", getFileNames(object)))
    if (showLimits) {
        # calculate timeLimit
        indLim <- timeLimits(object, fracMaxTIC = fracMaxTIC, plotDel = FALSE)$exp
        plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = getRawInfo(object)$time[c(indLim)], 
            color = "time limits")) + ggplot2::scale_fill_manual("Legend")
    }
    plot <- plot + ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold"), 
        axis.title = ggplot2::element_text(size = 16), axis.text = ggplot2::element_text(size = 14), 
        legend.text = ggplot2::element_text(size = 14), legend.title = ggplot2::element_text(size = 16)) + 
        ggplot2::theme_classic()
    switch(type, ggplot = return(plot), plotly = return(plotly::ggplotly(plot)))
}  #end function
)
## timeLimit ----
#' Calculates time limits on the breath tracer
#' 
#' This function derives limits on the breath tracer indicated, where the 
#' intensity is greater than \code{fracMaxTIC*max(tracer)}. By setting
#' \code{fracMaxTIC} close to 1, the size of the limits will be restricted.
#' This function also determine the index corresponding to the background, where 
#' variation between two successive point can be control with \code{derivThreshold} 
#' parameter.
#' 
#' @param object a ptrRaw or ptrSet object
#' @param fracMaxTIC between 0 and 1. Percentage of the maximum of the tracer 
#' amplitude with baseline removal. If you want a finer limitation, increase 
#' \code{fracMaxTIC}, indeed decrease
#' @param fracMaxTICBg same as fracMaxTIC but for background detection (lower than 
#' fracMaxTIC*max(TIC))
#' @param derivThresholdExp the threshold of the difference between two
#' successive points of the expiration
#' @param derivThresholdBg the threshold of the difference between two
#' successive points of the background
#' @param mzBreathTracer NULL or a integer. Correspond to a nominal masses of 
#' Extract Ion Current (EIC)
#'  whose limits you want to compute. If NULL, the limits are calculated on the 
#'  Total Ion Current (TIC).
#' @param minPoints minimum duration of an expiration (in index).
#' @param degreeBaseline the degree of polynomial baseline function
#' @param baseline logical, should the trace be baseline corrected?
#' @param redefineKnots logical, should the knot location must be redefined with 
#' the new times limits ?
#' @param plotDel boolean. If TRUE, the trace is plotted with limits and 
#' threshold.
#' @return a list with expiration limits (a matrix of index, where each column 
#' correspond to one expiration, the first row it is the beginning and the second
#' the end, or NA if no limits are detected) and index of the background.
#' @rdname timeLimits
#' @examples
#' 
#' ## ptrRaw object
#' 
#' library(ptairData)
#' filePath <- system.file('extdata/exhaledAir/ind1', 'ind1-1.h5',
#' package = 'ptairData')
#' raw <- readRaw(filePath)
#' 
#' timLim <- timeLimits(raw, fracMaxTIC=0.9, plotDel=TRUE)
#' timLim_acetone <- timeLimits(raw, fracMaxTIC=0.5, mzBreathTracer = 59,
#' plotDel=TRUE)
#'@export
methods::setMethod(f = "timeLimits", signature = "ptrRaw", function(object, fracMaxTIC = 0.6, 
    fracMaxTICBg = 0.2, derivThresholdExp = 0.5, derivThresholdBg = 0.05, mzBreathTracer = NULL, 
    minPoints = 2, degreeBaseline = 1, baseline = TRUE, redefineKnots = TRUE, plotDel = FALSE) {
    rawM <- getRawInfo(object)$rawM
    mz <- getRawInfo(object)$mz
    if (is.null(dim(rawM))) 
        stop("rawM must be a matrix")
    if (fracMaxTIC < 0 || fracMaxTIC > 1) 
        stop("fracMaxTIC must be between 
                                                     0 and 1")
    if (is.null(mzBreathTracer)) {
        TIC <- colSums(rawM)
    } else {
        index <- lapply(mzBreathTracer, function(x) {
            th <- 350 * x/10^6
            which(x - th < mz & mz < x + th)
        })
        TIC <- colSums(rawM[unlist(index), ])
    }
    indLim <- timeLimitFun(TIC, fracMaxTIC, fracMaxTICBg, derivThresholdExp, derivThresholdBg, 
        mzBreathTracer, minPoints, degreeBaseline, baseline, plotDel)
    return(indLim)
})


timeLimitFun <- function(TIC, fracMaxTIC = 0.5, fracMaxTICBg = 0.5, derivThresholdExp = 0.5, 
    derivThresholdBg = 0.01, mzBreathTracer = NULL, minPoints = 3, degreeBaseline = 1, 
    baseline = TRUE, plotDel = FALSE) {
    if (fracMaxTIC == 0) {
        return(list(exp = matrix(c(1, length(TIC)), ncol = 1, nrow = 2, dimnames = list(c("start", 
            "end"))), backGround = NULL))
    }
    ## baseline corretion
    if (baseline) 
        bl <- try(baselineEstimation(TIC, d = degreeBaseline)) else bl <- 0
    if (is.null(attr(bl, "condition"))) 
        TIC.blrm <- TIC - bl else TIC.blrm <- TIC
    threshold <- (max(TIC.blrm) - min(TIC.blrm)) * fracMaxTIC
    thresholdBg <- (max(TIC.blrm) - min(TIC.blrm)) * fracMaxTICBg
    if (!baseline) {
        threshold <- threshold + min(TIC.blrm)
        thresholdBg <- thresholdBg + min(TIC.blrm)
    }
    ## delimitation
    ind.Exp <- which(TIC.blrm > (thresholdBg))
    dTIC <- diff(TIC.blrm)/max(TIC.blrm)
    bool <- abs(dTIC) < derivThresholdBg
    indBg <- seq(2, length(TIC))[bool]
    if (bool[1]) 
        indBg <- c(1, indBg)
    if (length(which(indBg %in% ind.Exp)) > 0) 
        indBg <- indBg[-which(indBg %in% ind.Exp)]
    if (length(indBg) == 0) {
        warning("no background detect")
        indBg <- NULL
    }
    ## delimitation
    hat <- which(TIC.blrm > (threshold))
    if (length(hat) == 0) {
        message("no limits detected")
        plot(TIC, type = "l", xlab = "Time (s)", ylab = "intensity", cex.lab = 1.5, 
            main = "Time limit")
        return(NA)
    }
    hat_end <- c(hat[which(diff(hat) != 1)], utils::tail(hat, 1))
    hat_begin <- c(hat[1], hat[which(diff(hat) != 1) + 1])
    hat_lim <- unname(rbind(hat_begin, hat_end))
    row.names(hat_lim) <- c("start", "end")
    hat_lim <- hat_lim[, hat_lim["end", ] - hat_lim["start", ] >= minPoints, drop = FALSE]
    for (j in seq_len(ncol(hat_lim))) {
        index <- seq(hat_lim[1, j], hat_lim[2, j])
        dTIC <- diff(TIC.blrm[index])/max(TIC.blrm[index])
        indexNEW <- index[abs(dTIC) < derivThresholdExp]
        hat_lim[, j] <- range(indexNEW)
    }
    hat_lim <- hat_lim[, hat_lim["end", ] - hat_lim["start", ] >= minPoints, drop = FALSE]
    inExp <- Reduce(c, apply(hat_lim, 2, function(x) seq(x[1], x[2])))
    if (plotDel) {
        if (is.null(mzBreathTracer)) 
            subtitle <- "TIC" else subtitle <- paste("EIC mz:", paste(round(mzBreathTracer, 2), collapse = "-"))
        plot(seq_along(TIC.blrm), TIC.blrm, type = "l", xlab = "Time (s)", ylab = "intensity", 
            cex.lab = 1.5, main = paste("Time limit", subtitle), ylim = c(min(TIC.blrm) - 
                0.2 * (max(TIC.blrm) - min(TIC.blrm)), max(TIC.blrm)), lwd = 2)
        graphics::points(seq_along(TIC.blrm)[inExp], TIC.blrm[inExp], col = "red", 
            pch = 19)
        graphics::abline(h = threshold, lty = 2)
        graphics::points(seq_along(TIC.blrm)[indBg], TIC.blrm[indBg], col = "blue", 
            pch = 19)
        graphics::legend("bottomleft", legend = c("threshold", "Exp", "Backgound"), 
            col = c("black", "red", "blue"), lty = c(2, NA, NA), pch = c(NA, 19, 
                19), horiz = TRUE)
    }
    return(list(backGround = indBg, exp = hat_lim))
}
bakgroundDetect <- function(TIC, derivThreshold = 0.01, minPoints = 4, plotDel = FALSE) {
    bl <- try(baselineEstimation(TIC, d = 1))
    if (is.null(attr(bl, "condition"))) 
        TIC.blrm <- TIC - bl else TIC.blrm <- TIC
    threshold <- (max(TIC.blrm) - min(TIC.blrm)) * 0.5
    ## delimitation
    ind.Exp <- which(TIC.blrm > (threshold))
    dTIC <- diff(TIC.blrm)/max(TIC.blrm)
    ind.Bg <- which(abs(dTIC) < derivThreshold)
    ind.Bg <- ind.Bg[-which(ind.Bg %in% ind.Exp)]
    ind.Bg_end <- c(ind.Bg[which(diff(ind.Bg) != 1)], utils::tail(ind.Bg, 1))
    ind.Bg_begin <- c(ind.Bg[1], ind.Bg[which(diff(ind.Bg) != 1) + 1])
    limBg <- unname(rbind(ind.Bg_begin, ind.Bg_end))
    row.names(limBg) <- c("start", "end")
    limBg <- limBg[, limBg["end", ] - limBg["start", ] >= minPoints, drop = FALSE]
    if (ncol(limBg) == 0) 
        warning("no background detected")
    if (plotDel) {
        subtitle <- "TIC"
        plot(TIC, type = "l", xlab = "Time (s)", ylab = "intensity", cex.lab = 1.5, 
            main = paste("Background limit", subtitle), ylim = c(min(TIC) - 0.2 * 
                (max(TIC) - min(TIC)), max(TIC)), lwd = 2)
        graphics::abline(v = c(limBg), col = "red", lty = 2, lwd = 2)
        graphics::legend("bottomleft", legend = c("limit"), col = c("red"), lty = c(2), 
            horiz = TRUE)
    }
    return(limBg)
}
### defineknots ----
#' @param timeLimit index time of the expiration limits and background 
#' @return numeric vector of knots 
#' @rdname defineKnots
#' @export
methods::setMethod(f = "defineKnots", signature = "ptrRaw", function(object, knotsPeriod = 3, 
    method = c("expiration", "uniform")[1], knotsList = NULL, timeLimit = list(NULL)) {
    if (knotsPeriod == 0) {
        knots <- list(NULL)
    } else {
        if (!is.null(knotsList)) {
            t <- getRawInfo(object)$time
            if (knotsList[1] >= t[1] & utils::tail(knotsList, 1) <= utils::tail(t, 
                1)) 
                stop("knots are not contained in the time axis \n")
        } else {
            background <- timeLimit$backGround
            t <- getRawInfo(object)$time
            knots <- try(defineKnotsFunc(t, background, knotsPeriod, method))
        }
    }
    return(knots)
})
defineKnotsFunc <- function(t, background, knotsPeriod, method, file = NULL) {
    # knot equally spaced
    duration <- utils::tail(t, 1)
    knots <- seq(t[1], utils::tail(t, 1), length.out = duration/knotsPeriod)
    if (method == "expiration" & !is.null(background)) {
        # reduce number of knots in bakcgrounds periods
        knots <- defineKnotsExpiration(t, background, knots)
    }
    test <- vapply(seq_along(knots[-1]), function(i) any(knots[i] < t & t < knots[i + 
        1]), FUN.VALUE = TRUE)
    if (!all(test)) {
        warning(file, ",knotsPeriod is to short, knots set to NULL")
        knots <- NULL
    }
    if (length(knots) > 100) 
        warning(file, " K= ", length(knots), " we suggest to set a highter knots 
                                      frequency \n")
    return(knots)
}
defineKnotsExpiration <- function(t, background, knots) {
    periods <- c(0, which(diff(background) > 1), length(background))
    # background periods
    bg <- lapply(seq(1, length(periods) - 1), function(i) {
        (t[background[periods[i] + 1]]):t[background[periods[i + 1]]]
    })
    # keep only first, middle and last point of each background periods, if their
    # exceed three points
    newknot <- list(NULL)
    for (j in seq_along(bg)) {
        k <- knots[bg[[j]][1] <= knots & knots <= utils::tail(bg[[j]], 1)]
        # knots in the background period
        if (length(k) > 3) {
            k <- c(k[1], mean(k), utils::tail(k, 1))
            # select first, middle and last point of the background periods
            knots <- knots[-which(bg[[j]][1] < knots & knots < utils::tail(bg[[j]], 
                1))]
            # delete knots
        }
        newknot[[j]] <- k
    }
    knots <- sort(unique(c(knots, Reduce(c, newknot))))  # add new knots
    return(knots)
}
## PeakList ----
#' Detection and quantification of peaks on a sum spectrum. 
#'
#' @param raw \code{\link[ptairMS]{ptrRaw-class}} object 
#' @param mzNominal the vector of nominal mass where peaks will be detected 
#' @param ppm the minimum distance between two peeks in ppm 
#' @param resolutionRange vector with resolution min, resolution Mean, and 
#' resolution max of the PTR
#' @param minIntensity the minimum intenisty for peaks detection. The final 
#' threshold for peak detection will be : max ( \code{minPeakDetect} , 
#' thresholdNoise ). The thresholdNoise correspond to 
#' max(\code{thNoiseRate} * max( noise around the nominal mass), \code{minIntensityRate} * 
#'  max( intenisty in the nominal mass). The noise around the nominal mass correspond : 
#'  \code{[m-windowSize-0.2,m-windowSize]U[m+windowSize,m+WindowSize+0.2]}.
#' @param fctFit the function for the quantification of Peak, should be average 
#' or Sech2
#' @param peakShape a list with reference axis and a reference peak shape 
#' centered in zero
#' @param maxIter maximum iteration of residual analysis
#' @param R2min R2 minimum to stop the iterative residual analysis
#' @param autocorNoiseMax the autocorelation threshold for Optimal windows 
#' Savitzky Golay 
#' filter in \code{OptimalWindowSG} ptairMS function. See \code{?OptimalWindowSG}
#' @param plotFinal boolean. If TRUE, plot the spectrum for all nominal masses, 
#' with the final fitted peaks
#' @param plotAll boolean. Tf TRUE, plot all step to get the final fitted peaks
#' @param thNoiseRate The rate which multiplies the max noise intensity
#' @param minIntensityRate The rate which multiplies the max signal intensity
#' @param countFacFWHM integer. We will sum the fitted peaks on a window's size 
#' of countFacFWHM * FWHM, centered in the mass peak center.
#' @param daSeparation the minimum distance between two peeks in Da for nominal 
#' mass < 17.
#' @param d the degree for the \code{Savitzky Golay} filtrer
#' @param windowSize peaks will be detected only around  m - windowSize ;
#'  m + windowSize, for all 
#' m in \code{mzNominal}
#' @return a list containing: \itemize{
#' \item peak: a data.frame, with for all peak detected: the mass center, the 
#' intensity count in cps, the peak width (delta_mz), correspond to the Full Width Half 
#' Maximum (FWHM),the resolution m/delta_m, the other parameters values estimated 
#' of \code{fitFunc}.
#' \item warnings: warnings generated by the peak detection algorithm per nominal masses
#' \item infoPlot: elements needed to plot the fitted peak per nominal masses
#' }
#' @examples 
#' library(ptairData)
#' filePath <- system.file('extdata/exhaledAir/ind1', 'ind1-1.h5', 
#' package = 'ptairData')
#' file <- readRaw(filePath)
#' 
#' peakList <- PeakList(file, mzNominal = c(21,63))
#' peakList$peak
#'@rdname PeakList 
#'@export
methods::setMethod(f = "PeakList", signature = "ptrRaw", 
                   function(raw, 
                            mzNominal = unique(round(getRawInfo(raw)$mz)), 
    ppm = 130, resolutionRange = c(300, 5000, 8000), minIntensity = 5, fctFit = c("sech2", 
        "averagePeak")[1], peakShape = NULL, maxIter = 3, R2min = 0.995, autocorNoiseMax = 0.3, 
    plotFinal = FALSE, plotAll = FALSE, thNoiseRate = 1.1, minIntensityRate = 0.01, 
    countFacFWHM = 10, daSeparation = 0.005, d = 3, windowSize = 0.4) {
    # get raw element
    sp <- rowSums(getRawInfo(raw)$rawM)/(ncol(getRawInfo(raw)$rawM) * (getRawInfo(raw)$time[3] - getRawInfo(raw)$time[2]))  # average spectrum 
    mz <- getRawInfo(raw)$mz  # mass axis
    mzCalibRef <- getCalibrationInfo(raw)$calibMassRef
    calibCoef <- getCalibrationInfo(raw)$calibCoef[[1]]
    if ( (fctFit == "averagePeak") & is.null(peakShape)) {
        peakShape <- determinePeakShape(raw)$peakShapemz
        raw<-setPeakShape(raw,peakShape)

    }
    prePeaklist <- lapply(mzNominal, function(m) {
        peakListNominalMass(i = m, mz = mz, sp = sp, ppmPeakMinSep = ppm, calibCoef = calibCoef, 
            resolutionRange = resolutionRange, minPeakDetect = minIntensity, fitFunc = fctFit, 
            maxIter = maxIter, R2min = R2min, autocorNoiseMax = autocorNoiseMax, 
            plotFinal, plotAll, thNoiseRate, minIntensityRate, countFacFWHM, daSeparation, 
            d, windowSize, peakShape)
    })
    peaklist <- do.call(rbind, lapply(prePeaklist, function(x) x[[1]]))
    warning <- do.call(rbind, lapply(prePeaklist, function(x) x[[2]]))
    infoPlot <- do.call(c, lapply(prePeaklist, function(x) x[[3]]))
    baseline <- do.call(c, lapply(prePeaklist, function(x) x[[4]]))
    return(list(peak = peaklist, warning = warning, infoPlot = infoPlot, baseline = baseline))
})
## TODO peakLIst method dor spectrum array

## detectpeak----
#' @param knots numeric vector corresponding to the knot values, which used for 
#' the two dimensional regression for each file. Should be provided 
#' by \code{\link[ptairMS]{defineKnots}} function
#' @param timeLimit index time of the expiration limits and background. 
#' Should be provided by \code{\link[ptairMS]{timeLimits}} function
#' @param mzPrimaryIon the exact mass of the primary ion isotope
#' @rdname detectPeak
#' @examples
#' 
#' ## For ptrRaw object
#' library(ptairData)
#' filePath <- system.file('extdata/exhaledAir/ind1', 'ind1-1.h5', 
#' package = 'ptairData')
#' raw <- readRaw(filePath,mzCalibRef=c(21.022,59.049))
#' timeLimit<-timeLimits(raw,fracMaxTIC=0.7)
#' knots<-defineKnots(object = raw,timeLimit=timeLimit)
#' raw <- detectPeak(raw, timeLimit=timeLimit, mzNominal = c(21,59),
#' smoothPenalty=0,knots=knots)
#' @export
methods::setMethod(f = "detectPeak", signature = "ptrRaw", function(x, ppm = 130, 
    minIntensity = 10, minIntensityRate = 0.01, mzNominal = NULL, resolutionRange = NULL, 
    fctFit = NULL, smoothPenalty = NULL, timeLimit, knots = NULL, mzPrimaryIon = 21.022, 
    ...) {
    
    raw <- x
    
    # get infomration
    massCalib <- getCalibrationInfo(raw)$calibMassRef
    
    # resolution
    if (is.null(resolutionRange)) {
        calibSpectr <- alignCalibrationPeak(getCalibrationInfo(raw)$calibSpectr, calibMassRef = massCalib, 
            ntimes = length(getRawInfo(raw)$time))
        resolutionEstimated <- estimateResol(getCalibrationInfo(raw)$calibMassRef, calibSpectr)
        resolutionRange <- c(floor(min(resolutionEstimated)/1000) * 1000, round(mean(resolutionEstimated)/1000) * 
            1000, ceiling(max(resolutionEstimated)/1000) * 1000)
    }
    
    
    
    # peakShape
    peakShape <- determinePeakShape(raw)$peakShapemz
    
    # check best fit
    sech2 <- mean(PeakList(raw, mzNominal = getCalibrationInfo(raw)$calibMassRef, 
                           fctFit = "sech2", 
                           maxIter = 1, ppm = 500, minIntensityRate = 0.2, 
                           windowSize = 0.2, resolutionRange = resolutionRange, 
                           peakShape = peakShape)$peak$R2)
    
    averagePeak <- mean(PeakList(raw, mzNominal = getCalibrationInfo(raw)$calibMassRef, fctFit = "averagePeak", 
                                 resolutionRange = resolutionRange, maxIter = 1, peakShape = peakShape, 
                                 ppm = 500, minIntensityRate = 0.2, windowSize = 0.2)$peak$R2)
    asymGauss <- mean(PeakList(raw, mzNominal = getCalibrationInfo(raw)$calibMassRef, fctFit = "asymGauss", 
                               resolutionRange = resolutionRange, maxIter = 1, peakShape = peakShape, 
                               ppm = 500, minIntensityRate = 0.2, windowSize = 0.2)$peak$R2)
    
    fctFit <- c("sech2", "averagePeak", "asymGauss")[which.max(c(sech2,averagePeak, asymGauss))]
    
    # primary ion
    p <- PeakList(raw, mzNominal = round(21.022), ppm = 700, minIntensity = 50, maxIter = 1,
                  fctFit = fctFit,resolutionRange =resolutionRange )
    raw@primaryIon<-p$peak$quanti_cps
    primaryIon <- list(primaryIon = p$peak$quanti_cps)
    # knot
    peakLists <- processFileTemporal(fullNamefile = raw, massCalib = massCalib, primaryIon = primaryIon, 
        indTimeLim = timeLimit, mzNominal = mzNominal, ppm = ppm, resolutionRange = resolutionRange, 
        minIntensity = minIntensity, knots = knots, smoothPenalty = smoothPenalty, 
        fctFit = fctFit, minIntensityRate = minIntensityRate, peakShape = peakShape, 
        ...)
    x <- peakLists
    infoPeak <- grep("parameter", colnames(x$raw))
    colnames(x$raw)[infoPeak] <- c("parameterPeak.delta1", "parameterPeak.delta2", 
        "parameterPeak.heigh")
    assayMatrix <- as.matrix(x$raw)[, -c(1, 2, 3, 4, infoPeak), drop = FALSE]
    rownames(assayMatrix) <- round(x$raw$Mz, 4)
    featuresMatrix <- data.frame(cbind((as.matrix(x$raw)[, c(1, 2, 3, 4, infoPeak), 
        drop = FALSE]), (x$aligned[, -1])))
    rownames(featuresMatrix) <- rownames(assayMatrix)
    peakLists <- Biobase::ExpressionSet(assayData = assayMatrix, featureData = Biobase::AnnotatedDataFrame(featuresMatrix))
    
    raw@peakList<-peakLists
    raw@fctFit<-fctFit
    raw@resolution<-resolutionRange
    
    return(raw)
})


estimateResol <- function(calibMassRef, calibSpectr) {
    m <- calibMassRef
    delta <- vapply(calibSpectr, function(x) {
        mzM <- x$mz
        spM <- x$signal
        # half maximum
        hm <- max(spM)/2
        # limits
        lim1 <- findEqualGreaterM(spM[seq_len(which.max(spM))], hm)
        lim2 <- unname(which.max(spM)) + FindEqualLess(spM[(which.max(spM) + 1):length(spM)], 
            hm)
        # equation : (intrepolation linéaire entre lim et (lim-1)) = hm
        deltaBorne <- vapply(c(lim1, lim2), function(x) (hm * (mzM[x] - mzM[x - 1]) - 
            (mzM[x] * spM[x - 1] - mzM[x - 1] * spM[x]))/(spM[x] - spM[x - 1]), FUN.VALUE = 1.1)
        delta <- diff(deltaBorne)
        return(delta)
    }, FUN.VALUE = 0)
    resol <- m/delta
    names(resol) <- m
    return(resol)
}
## show ----
#' show a ptrRaw object 
#' 
#' It indicates the files, the mz range, time acquisition range, and calibration error.
#' @param object a ptrRaw object
#' @return nothing
#' @export 
methods::setMethod("show", "ptrRaw", function(object) {
    cat(getFileNames(object), "\n")
    cat("  mz range: ", paste(round(range(getRawInfo(object)$mz), 2), collapse = " - "), "\n")
    cat("  time range: ", paste(round(range(getRawInfo(object)$time), 2), collapse = " - "), 
        "\n")
    cat("  Calibration error in ppm: \n")
    if (length(getCalibrationInfo(object)$calibError) == 1) {
        cat("     No calibration performs")
    } else {
        print(round(getCalibrationInfo(object)$calibError, 2))
    }
})
#### dead time correction -----
#'Dead time correction on raw data
#'@param raw ptrRaw object
#'@param ve extending dead time
#'@param vne non extending dead time
#'@param r number of extraction
#'@param threshold only bin of intensity more then threshold*r which be corrected
#'@return a ptrRaw object with the raw matrix corrected
#'@keywords internal
deadTimeCorr <- function(raw, ve, vne, r, threshold = 0.1) {
    rawM <- getRawInfo(raw)$rawM
    index <- which(rawM > threshold * r)
    for (j in index) {
        rawM[j] <- -r * log(1 - (rawM[j]/r * (1 - sum(raw[(j - vne):(j - 1 - ve)])/r)^-1 * 
            exp(sum(rawM[(j - ve):(j - 1)])/r)))
    }
    raw@rawM <- rawM
    return(raw)
}
