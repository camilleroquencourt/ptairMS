utils::globalVariables(c("error", "name", "out", "intervalRef", "signal", "signal0", 
    "signal1", "Mz"))
### plot ----
#' Plot a ptrSet object
#' 
#' plot a ptrSet object
#' @aliases plot.ptrSet plot,ptrSet-method
#' @param x a ptrSet object 
#' @param y not use
#' @param typePlot could be : \code{calibError}, \code{resolution},  
#' \code{peakShape}, or 
#' a empty character if you want all. 
#' @return plot 
#' @rdname plot
#' @export
#' @examples 
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' plot(exhaledPtrset )
#' plot(exhaledPtrset ,typePlot='calibError')
#' plot(exhaledPtrset ,typePlot='resolution')
#' plot(exhaledPtrset ,typePlot='peakShape')
methods::setMethod(f = "plot", signature = "ptrSet", function(x, y, typePlot = "") {
    # if(!is.null(saveFile)) { ext<- tools::file_ext(saveFile) if(ext == 'pdf')
    # pdf(saveFile) if(ext == 'png') png(saveFile) if(ext == 'jpeg') jpeg(saveFile) }
    if (!typePlot %in% c("", "calibError", "resolution", "peakShape")) 
        warning("typePolt should be calibError,resolution or peakShape")
    if (typePlot == "calibError") {
        return(plotCalibError(x))
    } else if (typePlot == "resolution") {
        return(plotResolution(x)$plot)
    } else if (typePlot == "peakShape") {
        return(plotPeakShape(x))
    } else {
        reso <- try(plotResolution(x))
        calib <- try(plotCalibError(x))
        peakShape <- try(plotPeakShape(x))
        reaction <- try(plotPtrReaction(x))
        left <- ggpubr::ggarrange(calib, reso$plot, reso$table, nrow = 3, 
                                  heights = c(0.5, 
            0.4, 0.1))
        right <- ggpubr::ggarrange(peakShape, reaction, nrow = 2)
        return(ggpubr::ggarrange(left, right, ncol = 2, align = "hv"))
    }
    # if(!is.null(saveFile)) dev.off()
})
# plot resolution boxplot for a ptrSet
plotResolution <- function(set) {
    mzCalibRef <- getParameters(set)$mzCalibRef
    resolution <- getPTRInfo(set)$resolution
    
    # complete missing masses
    index <- which(vapply(resolution, length, 1) != length(mzCalibRef))
    resolution[index] <- lapply(resolution[index], function(x) {
        missing <- which(!mzCalibRef %in% names(x))
        x_new <- c(x, rep(as.numeric(NA), length(missing)))
        names(x_new) <- c(names(x), mzCalibRef[missing])
        x_new <- x_new[order(as.numeric(names(x_new)))]
        x_new
    })
    # concatenate
    resolutiontemp <- lapply(resolution, function(x) {
        mat <- cbind(resolution = x, Mz = as.numeric(names(x)))
        rownames(mat) <- NULL
        mat
    })
    resolMat <- do.call(rbind, resolutiontemp)
    resolMat <- data.table::as.data.table(resolMat)
    `:=` <- data.table::`:=`
    resolMat[, `:=`(Mz, factor(Mz))]
    resolMat <- cbind(resolMat, name = rep(names(resolution), 
                                           each = nlevels(resolMat$Mz)))
    # identify outliers
    resolMat <- resolMat[!is.na(resolMat$resolution), ]
    resolMat <- resolMat[, list(resolution, name, out = ifelse(is_outlier(resolution), 
        resolution, as.numeric(NA))), by = Mz]
    # boxplot with outlier labels
    g <- ggplot2::ggplot(subset(resolMat, !is.na(resolution)), ggplot2::aes(y = resolution, 
        x = Mz)) + ggplot2::geom_boxplot() + ggplot2::geom_text(data = resolMat[!is.na(out), 
        ], ggplot2::aes(x = Mz, y = resolution, label = name), vjust = -0.5, size = 3) + 
        ggplot2::ggtitle(label = expression("Resolution m/" ~ Delta ~ "(m)")) + ggplot2::ylab("m/delta(m)") + 
        ggplot2::theme_classic()
    resolutionEstimated <- Reduce(c, getPTRInfo(set)$resolution)
    info <- gridExtra::tableGrob(t(data.frame(resolution = c(min = floor(min(resolutionEstimated)/1000) * 
        1000, mean = round(mean(resolutionEstimated)/1000) * 1000, 
        max = ceiling(max(resolutionEstimated)/1000) * 
        1000))), theme = gridExtra::ttheme_minimal(base_size = 12))
    return(list(plot = g, table = info))
}
#' @importFrom rlang .data    
plotCalibError <- function(set) {
    massCalib <- getCalibrationInfo(set)$mzCalibRef
    errorList <- getCalibrationInfo(set)$errorCalibPpm
    names(errorList) <- paste(seq(1, length(errorList)), names(errorList),
                              sep = "-")
    # get list files
    listFiles <- getParameters(set)$listFile
    if (methods::is(listFiles, "expression")) 
        listFiles <- eval(listFiles)
    # number of files
    nfiles <- length(massCalib)
    # calibration reference
    mzCalibRef <- getParameters(set)$mzCalibRef
    # count how many time masses where excluded
    count <- table(unlist(massCalib))
    exclud <- nfiles - count
    # TODO: affichier cette info sur le graph
    # complete missing masses
    index <- which(vapply(errorList, length, 1) != length(mzCalibRef))
    errorList[index] <- lapply(errorList[index], function(x) {
        missing <- which(!mzCalibRef %in% names(x))
        x_new <- c(x, rep(NA, length(missing)))
        names(x_new) <- c(names(x), mzCalibRef[missing])
        x_new <- x_new[order(as.numeric(names(x_new)))]
        x_new
    })
    # concatenate
    errorListemp <- lapply(errorList, function(x) {
        mat <- cbind(error = x, Mz = as.numeric(names(x)))
        rownames(mat) <- NULL
        mat
    })
    calibErrorMat <- do.call(rbind, errorListemp)
    calibErrorMat <- data.table::as.data.table(calibErrorMat)
    `:=` <- data.table::`:=`
    calibErrorMat[, `:=`(Mz, factor(Mz))]
    calibErrorMat <- cbind(calibErrorMat, name = rep(names(errorList), 
                                                     each = nlevels(calibErrorMat$Mz)))
    # identify outliers
    calibErrorMat <- calibErrorMat[!is.na(calibErrorMat$error), ]
    calibErrorMat <- calibErrorMat[, list(error, name, 
                                          out = ifelse(is_outlier(error), 
        error, as.numeric(NA))), by = Mz]
    # boxplot with labels
    p <- ggplot2::ggplot(subset(calibErrorMat, !is.na(error)), 
                         ggplot2::aes(y = .data$error,  x = .data$Mz)) + 
        ggplot2::geom_boxplot() + 
        ggplot2::geom_text(data = calibErrorMat[!is.na(out),], 
                           ggplot2::aes(x = .data$Mz, 
                                        y = .data$error, 
                                        label = .data$name), vjust = -0.5, size = 3.5) + 
        ggplot2::ggtitle("Calibration error") + ggplot2::ylab("ppm") + 
        ggplot2::theme_classic()
    return(p)
}
is_outlier <- function(x) {
    return(x < stats::quantile(x, 0.25, na.rm = TRUE) - 1.5 * stats::IQR(x, na.rm = TRUE) | 
        x > stats::quantile(x, 0.75, na.rm = TRUE) + 1.5 * stats::IQR(x, na.rm = TRUE))
}
#' plot the average peak shape of reference calibration masses for a ptrSet
#' @param set ptrSet object
#' @param showAverage boolean 
#' @return ggplot object
#' @importFrom scales hue_pal 
#' @keywords internal
plotPeakShape <- function(set, showAverage = FALSE) {
    mzRef <- getParameters(set)$mzCalibRef
    n.files <- length(getCalibrationInfo(set)$mzCalibRef)
    npoints <- 50
    interval.ref <- seq(-3, 3, length = npoints)  # +- 3 * peak width(=1) 
    # loop over file
    peak_ref <- array(NA, dim = c(length(mzRef), n.files, npoints))
    for (n.file in seq_len(n.files)) {
        massRef <- getCalibrationInfo(set)$mzCalibRef[[n.file]]
        interval <- getCalibrationInfo(set)$signalCalibRef[[n.file]]
        names(interval) <- massRef
        n.mass <- length(interval)
        delta <- vapply(interval, function(x) width(tof = x$mz, peak = x$signal)$delta, 
            FUN.VALUE = 1.1)
        t_centre <- vapply(massRef, function(x) {
            sp <- interval[[as.character(x)]]$signal
            mz <- interval[[as.character(x)]]$mz
            localMax <- LocalMaximaSG(sp = sp, minPeakHeight = max(sp) * 0.2)
            interpol <- stats::spline(mz[(localMax - 4):(localMax + 4)], sp[(localMax - 
                4):(localMax + 4)], n = 100)
            sg <- signal::sgolayfilt(interpol$y, n = 51, p = 3)  #n/2
            center <- interpol$x[which.max(sg)]
            return(center)
        }, FUN.VALUE = 1.1)
        # normalization of intreval
        interval.n <- list()
        for (j in seq_len(n.mass)) {
            interval.n[[j]] <- (interval[[j]]$mz - t_centre[j])/delta[j]
        }
        names(interval.n) <- massRef
        peak <- matrix(0, ncol = n.mass, nrow = length(interval.ref))
        for (i in seq_len(n.mass)) {
            # interpolation
            interpolation <- stats::spline(interval.n[[i]], interval[[i]]$signal, 
                xout = interval.ref)
            # baseline correction
            interpolation$y <- interpolation$y - snipBase(interpolation$y, widthI = 4)
            # normalization
            indexPeakRef <- which(massRef == massRef[i])
            peak_ref[indexPeakRef, n.file, ] <- interpolation$y/stats::spline(interval.ref, 
                interpolation$y, xout = 0)$y
            # cumsum(interpolation$y)/sum(interpolation$y)
        }
    }
    # aggreate
    # mean
    peaksMAT <- apply(peak_ref, 1, function(x) {
        if (all(is.na(x))) 
            return(rep(NA, 4 * npoints))
        peak <- apply(x, 2, function(x) mean(x, na.rm = TRUE))
        peak <- stats::spline(interval.ref, peak, n = 4 * npoints)
        return(peak$y)
    })
    interval.ref2 <- seq(-3, 3, length = 4 * npoints)
    peaks <- matrix(peaksMAT, ncol = 1)
    peakData <- data.frame(signal = peaks, Mz = as.factor(rep(round(mzRef, 3), 
                                                              each = 4 *npoints)), 
                           intervalRef = rep(interval.ref2, length(mzRef)))
    peakData <- peakData[!is.na(peakData$signal), ]
    # plot
    if (showAverage) {
        p <- ggplot2::ggplot() + 
            ggplot2::geom_line(data = peakData, ggplot2::aes(x = .data$intervalRef,
                                                             y = .data$signal, 
                                                             group = .data$Mz, 
                                                             colour = .data$Mz), 
                               size = 1) + 
            ggplot2::geom_line(data = data.frame(x = interval.ref2,
                                                 y = apply(peaksMAT, 
                                                           which.max(dim(peaksMAT)), 
                                                           mean)), 
                               ggplot2::aes(x = .data$x,y = .data$y, colour = "average"), 
                               size = 1.2) + 
            ggplot2::ggtitle("Average normalized peak shape of 
                             calibration peaks") + 
            ggplot2::xlab("Mz interval normalized") + 
            ggplot2::ylab("Intenisty normalized") + 
            ggplot2::scale_color_manual(values = c((scales::hue_pal())(length(levels(peakData$Mz))), 
                "black"))
    } else {
        p <- ggplot2::ggplot() + 
            ggplot2::geom_line(data = peakData, 
                               ggplot2::aes(x = .data$intervalRef, 
                                            y = .data$signal, group = .data$Mz, 
                                            colour = .data$Mz), size = 1) + 
            ggplot2::ggtitle("Average normalized peak shape of 
                             calibration peaks") + 
            ggplot2::xlab("Mz interval normalized") + 
            ggplot2::ylab("Intenisty normalized") + 
            ggplot2::scale_color_manual(values = c((scales::hue_pal())(length(levels(peakData$Mz)))))
    }
    if (n.files > 1) {
        # confidence interval
        peaksUp <- apply(peak_ref, 1, function(x) {
            if (all(is.na(x))) 
                return(rep(NA, 4 * npoints))
            peak <- apply(x, 2, function(y) mean(y, na.rm = TRUE) + 
                              stats::qnorm(0.975) * stats::sd(y, na.rm = TRUE)/
                              sqrt(length(y)))
            peak <- stats::spline(interval.ref, peak, n = 4 * npoints)
            return(peak$y)
        })
        peaksDown <- apply(peak_ref, 1, function(x) {
            if (all(is.na(x))) 
                return(rep(NA, 4 * npoints))
            peak <- apply(x, 2, function(y) mean(y, na.rm = TRUE) - stats::qnorm(0.975) * 
                stats::sd(y, na.rm = TRUE)/sqrt(length(y)))
            peak <- stats::spline(interval.ref, peak, n = 4 * npoints)
            return(peak$y)
        })
        peaksUp <- matrix(peaksUp, ncol = 1)
        peaksDown <- matrix(peaksDown, ncol = 1)
        peakData <- data.frame(signal0 = peaksDown, signal1 = peaksUp, 
                               Mz = as.factor(rep(round(mzRef, 3), each = 4 * npoints)), 
                               intervalRef = rep(interval.ref2, length(mzRef)))
        peakData <- peakData[!is.na(peakData$signal0), ]
        p <- p + ggplot2::geom_line(data = peakData, ggplot2::aes(x = .data$intervalRef, 
            y = .data$signal0, group = .data$Mz, colour = .data$Mz), size = 0.6, 
            linetype = "dashed") + ggplot2::geom_line(data = peakData, mapping = ggplot2::aes(x = .data$intervalRef, 
            y = .data$signal1, group = .data$Mz, colour = .data$Mz), size = 0.6, 
            linetype = "dashed")
    }
    return(p + ggplot2::theme_classic())
}
plotPeakShapeTof <- function(set) {
    mzRef <- getParameters(set)$mzCalibRef
    n.files <- length(getCalibrationInfo(set)$mzCalibRef)
    npoints <- 50
    interval.ref <- seq(-3, 3, length = npoints)  # +- 3 * peak width(=1) 
    # loop over file
    peak_ref <- array(NA, dim = c(length(mzRef), n.files, npoints))
    for (n.file in seq_len(n.files)) {
        massRef <- getCalibrationInfo(set)$mzCalibRef[[n.file]]
        interval <- alignCalibrationPeak(getCalibrationInfo(set)$signalCalibRef[[n.file]], massRef, 
                                         length(getTimeInfo(set)$TIC[[n.files]]))
        names(interval) <- massRef
        n.mass <- length(interval)
        coefs <- getCalibrationInfo(set)$coefCalib[[n.file]]
        # delta <- vapply(interval, function(x) width(tof =
        # sqrt(x$mz)*coefs['a',]+coefs['b',], peak = x$signal)$delta, FUN.VALUE = 1.1)
        delta <- vapply(interval, function(x) width(tof = x$mz, peak = x$signal)$delta, 
            FUN.VALUE = 1.1)
        t_centre <- vapply(massRef, function(x) {
            sp <- interval[[as.character(x)]]$signal
            mz <- interval[[as.character(x)]]$mz
            localMax <- LocalMaximaSG(sp = sp, minPeakHeight = max(sp) * 0.2)
            interpol <- stats::spline(mz[(localMax - 4):(localMax + 4)], sp[(localMax - 
                4):(localMax + 4)], n = 1000)
            sg <- signal::sgolayfilt(interpol$y, n = 501, p = 3)
            tofMax <- interpol$x[which.max(sg)]
            return(tofMax)
        }, FUN.VALUE = 1.1)
        # t_centre<-sqrt(t_centre)*coefs['a',]+coefs['b',]
        # normalization of intreval
        interval.n <- list()
        for (j in seq_len(n.mass)) {
            # interval.n[[j]]<-(sqrt(interval[[j]]$mz)*coefs['a',]+coefs['b',]-
            # t_centre[j])/delta[j]
            interval.n[[j]] <- (interval[[j]]$mz - t_centre[j])/delta[j]
        }
        names(interval.n) <- massRef
        peak <- matrix(0, ncol = n.mass, nrow = length(interval.ref))
        for (i in seq_len(n.mass)) {
            # interpolation
            interpolation <- stats::spline(interval.n[[i]], interval[[i]]$signal, 
                xout = interval.ref)
            # baseline correction
            interpolation$y <- interpolation$y - snipBase(interpolation$y, widthI = 4)
            # normalization
            indexPeakRef <- which(massRef == massRef[i])
            peak_ref[indexPeakRef, n.file, ] <- interpolation$y/stats::spline(interval.ref, 
                interpolation$y, xout = 0)$y
            # cumsum(interpolation$y)/sum(interpolation$y)
        }
    }
    # aggreate
    # mean
    peaks <- apply(peak_ref, 1, function(x) {
        if (all(is.na(x))) 
            return(rep(NA, 4 * npoints))
        peak <- apply(x, 2, function(x) mean(x, na.rm = TRUE))
        peak <- stats::spline(interval.ref, peak, n = 4 * npoints)
        return(peak$y)
    })
    interval.ref2 <- seq(-3, 3, length = 4 * npoints)
    peaks <- matrix(peaks, ncol = 1)
    peakData <- data.frame(signal = peaks, Mz = as.factor(rep(mzRef, each = 4 * npoints)), 
        intervalRef = rep(interval.ref2, length(mzRef)))
    peakData <- peakData[!is.na(peakData$signal), ]
    # plot
    p <- ggplot2::ggplot(data = peakData) + 
        ggplot2::geom_line(ggplot2::aes(x = intervalRef, y = signal, 
                                        group = Mz, colour = Mz), size = 1) + 
        ggplot2::ggtitle("Average normalized peak shape") + 
        ggplot2::xlab("Mz interval normalized") + 
        ggplot2::ylab("Intenisty normalized")
    if (n.files > 1) {
        # confidence interval
        peaksUp <- apply(peak_ref, 1, function(x) {
            if (all(is.na(x))) 
                return(rep(NA, 4 * npoints))
            peak <- apply(x, 2, function(y) mean(y, na.rm = TRUE) + 
                                stats::qnorm(0.975) * stats::sd(y, na.rm = TRUE)/
                                sqrt(length(y)))
            peak <- stats::spline(interval.ref, peak, n = 4 * npoints)
            return(peak$y)
        })
        peaksDown <- apply(peak_ref, 1, function(x) {
            if (all(is.na(x))) 
                return(rep(NA, 4 * npoints))
            peak <- apply(x, 2, function(y) mean(y, na.rm = TRUE) - 
                                stats::qnorm(0.975) * stats::sd(y, na.rm = TRUE)/
                                sqrt(length(y)))
            peak <- stats::spline(interval.ref, peak, n = 4 * npoints)
            return(peak$y)
        })
        peaksUp <- matrix(peaksUp, ncol = 1)
        peaksDown <- matrix(peaksDown, ncol = 1)
        peakData <- data.frame(signal0 = peaksDown, signal1 = peaksUp, 
                                Mz = as.factor(rep(mzRef,each = 4 * npoints)), 
                                intervalRef = rep(interval.ref2, length(mzRef)))
        peakData <- peakData[!is.na(peakData$signal0), ]
        p <- p + ggplot2::geom_line(data = peakData, ggplot2::aes(x = intervalRef, 
            y = signal0, group = Mz, colour = Mz), size = 0.6, 
            linetype = "dashed") + 
            ggplot2::geom_line(data = peakData, ggplot2::aes(x = intervalRef, 
                                                                y = signal1, 
                group = Mz, colour = Mz), size = 0.6, linetype = "dashed")
    }
    return(p)
}

plotPtrReaction <- function(pSet) {
    U <- Reduce(c, lapply(getPTRInfo(pSet)$prtReaction, function(x) {
        y <- c(x$TwData[1, , ])
        mean(y[y != 0])
    }))
    PD <- Reduce(c, lapply(getPTRInfo(pSet)$prtReaction, function(x) {
        y <- c(x$TwData[2, , ])
        mean(y[y != 0])
    }))
    TD <- Reduce(c, lapply(getPTRInfo(pSet)$prtReaction, function(x) {
        y <- c(x$TwData[3, , ])
        mean(y[y != 0])
    }))
    EN <- Reduce(c, lapply(getPTRInfo(pSet)$prtReaction, function(x) {
        y <- c(x$TwData[4, , ])
        mean(y[y != 0])
    }))
    primaryIon <- Reduce(c, lapply(getPTRInfo(pSet)$primaryIon, function(x) x$primaryIon))
    date <- Reduce(c, getDate(pSet))
    date <- vapply(date, function(x) chron::chron(dates. = strsplit(x, " ")[[1]][1], 
        times. = strsplit(x, " ")[[1]][2], format = c("d/m/y", "h:m:s")), 
        FUN.VALUE = chron::chron(1))
    
    Udrift <- ggplot2::ggplot() + 
        ggplot2::geom_point(mapping = ggplot2::aes(x = .data$date,y = .data$U), 
                            data = data.frame(date = as.Date(chron::as.dates(date)), 
                                                U = U)) + 
        ggplot2::ggtitle("Drift voltage") + ggplot2::ylab("V") + 
        ggplot2::theme_classic() + 
        ggplot2::theme(title = ggplot2::element_text(size = 9))
    
    Tdrift <- ggplot2::ggplot() + 
        ggplot2::geom_point(mapping = ggplot2::aes(x = .data$date,y = .data$TD), 
                            data = data.frame(date = as.Date(chron::as.dates(date)), TD = TD)) + 
        ggplot2::ggtitle("Drift temperature") + 
        ggplot2::ylab("degree") + 
        ggplot2::theme_classic() + 
        ggplot2::theme(title = ggplot2::element_text(size = 9))
    
    Pdrift <- ggplot2::ggplot() + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$date, 
        y = .data$PD), data = data.frame(date = as.Date(chron::as.dates(date)), PD = PD)) + 
        ggplot2::ggtitle("Drift pressure") + 
        ggplot2::ylab("mbar") + 
        ggplot2::theme_classic() + 
        ggplot2::theme(title = ggplot2::element_text(size = 9))
    
    primaryIonPlot <- ggplot2::ggplot() + 
        ggplot2::geom_point(mapping = ggplot2::aes(x = .data$date, y = .data$cps), 
        data = data.frame(date = as.Date(chron::as.dates(date)),
                          cps = primaryIon)) + 
        ggplot2::ggtitle("Primary ion isotope intensity") + 
        ggplot2::xlab("Date") + ggplot2::theme_classic() + 
        ggplot2::theme(title = ggplot2::element_text(size = 9))
    
    blank<-ggplot2::element_blank()
    reaction <- ggpubr::ggarrange(Udrift + 
                                      ggplot2::theme(axis.ticks.x =blank , 
                                                     axis.text.x = blank, 
                                                     axis.title.x = blank),
                                  Tdrift + ggplot2::theme(axis.ticks.x =blank, 
                                                          axis.text.x = blank,
                                                          axis.title.x = blank), 
                                  Pdrift + ggplot2::theme(axis.ticks.x = blank,
                                                          axis.text.x = blank, 
                                                          axis.title.x = blank), 
        primaryIonPlot, ncol = 1, heights = c(0.23, 0.23, 0.23, 0.31), align = "v")
    return(primaryIonPlot)
}

### plotFiles----

#' Plot the calibration peaks after calibration
#' 
#' @param object a ptrSet or ptrRaw object 
#' @param ppm the width of plot windows
#' @param pdfFile is different of \code{NULL}, the file path to save the 
#' plots in pdf
#' @param fileNames the name of the files in the ptrSet object to plot. 
#' If \code{NULL}, all files will be
#' plotted
#' @param ... not used
#' @return plot 
#' @export
#' @rdname plotCalib
#' @examples 
#' ## ptrSet
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' plotCalib(exhaledPtrset ,fileNames=getFileNames(exhaledPtrset )[1])
#'
#' ## ptrRaw 
#' filePath<-system.file('extdata/exhaledAir/ind1/ind1-1.h5', 
#' package = 'ptairData')
#' raw <- readRaw(filePath,mzCalibRef=c(21.022,59.049))
#' plotCalib(raw)
methods::setMethod(f = "plotCalib", signature = "ptrSet", function(object, ppm = 2000, 
    pdfFile = NULL, fileNames = NULL, ...) {
    fileNamesObject <- getParameters(object)$listFile
    if (methods::is(fileNamesObject, "expression")) 
        fileNamesObject <- eval(fileNamesObject)
    set <- object
    # get list files
    if (is.null(fileNames)) {
        fileNames <- basename(fileNamesObject)
    } else {
        # check if fileNames are in the ptrSet object put in basename
        fileNames <- basename(fileNames)
        if (methods::is(fileNames, "expression")) 
            fileNames <- eval(fileNames)
        test <- fileNames %in% basename(fileNamesObject)
        if (!all(test)) 
            stop("This file(s) names are not in 
                                   the directory: \n", 
                paste(!fileNames[test], collapse = "\n"))
    }
    # get information for plot
    massCalib <- getCalibrationInfo(set)$mzCalibRef[fileNames] 
    # spectrCalib<-lapply(fileNames, function(x) 
    # alignCalibrationPeak(getCalibrationInfo(set)$signalCalibRef[[x]], length(getTimeInfo(set)$TIC[[x]])))
    # names(spectrCalib)<-fileNames
    spectrCalib <- getCalibrationInfo(set)$signalCalibRef[fileNames]
    errorList <- getCalibrationInfo(set)$errorCalibPpm[fileNames]
    # get the extension of the file
    if (!is.null(pdfFile)) {
        # check is pdf extension
        file_ext <- tools::file_ext(pdfFile)
        if (file_ext != "pdf") {
            warning("the plot will be saved in pdf extension")
            pdfFile <- gsub(paste0(".", file_ext, "$"), ".pdf", pdfFile)
        }
        pdf(file = pdfFile, width = 7 * 1.8, height = 5 * 1.8)
    }
    # plot all masses use for calibration for each files loop over files
    for (i in seq_along(massCalib)) {
        error <- errorList[[i]]
        # plot in a window of 2000 ppm
        nb_plot <- length(massCalib[[i]])
        nb_col <- min(3, nb_plot)
        nb_row <- ceiling(nb_plot/nb_col)
        graphics::par(oma = c(0, 0, 3, 0))
        layout(matrix(seq(1, nb_row * nb_col), nrow = nb_row, ncol = nb_col, 
                        byrow = TRUE))
        for (j in seq_along(massCalib[[i]])) {
            m <- massCalib[[i]][j]
            th <- m * (ppm/2)/10^6
            mz <- spectrCalib[[i]][[j]]$mz
            index <- which(m - th < mz & mz < m + th)
            x <- mz[index]
            y <- spectrCalib[[i]][[j]]$signal[index]
            plot(x, y, type = "l", lwd = 2, ylab = "intenisty", 
                    xlab = "mz", main = c(m, paste("error:", round(error[j], 2),
                                                    "ppm")))
            abline(v = m, col = "red", lwd = 2)
        }  #end loop masses
        title(main = paste(i, "-", names(spectrCalib)[i]), outer = TRUE, 
            line = 0.5, 
            cex.main = 2)
    }  #end loop files
    if (!is.null(pdfFile)) 
        dev.off()
    graphics::par(oma = c(0, 0, 0, 0))
})

#' plot the Total Ion sptectrum (TIC) for one or several files.
#' @param object ptrSet or ptrRaw S4 object
#' @param type set 'plotly' to get an interactive plot, and 'ggplot' for 
#' classical plot.  
#' @param baselineRm logical. If \code{TRUE}, remove the baseline of the TIC
#' @param showLimits logical. If \code{TRUE}, add the time limits to the plot 
#' (obtain with the `fracMaxTIC` argument or `createPtrSet` function)
#' @return a plotly of ggplot2 object. 
#' @param pdfFile a absolute file path. A pdf will be generated with a plot 
#' for each file, caints TIC and 
#'time limits.
#' @param fileNames vector of character. The file names of the ptrSer that 
#' you want to plot. Could be in 
#' basename or fullname.
#' @param colorBy character. A name of the ptrSet's sampleMetaData column, 
#' to display with
#' the same color files of same attributes. 
#' @param normalizePrimariIon should the TIC be normalized by the primary ion
#' @param ... not used
#' @rdname plotTIC
#' @importFrom grDevices dev.off pdf
#' @export
#' @examples
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' plotTIC(exhaledPtrset ,type='ggplot')
methods::setMethod(f = "plotTIC", signature = "ptrSet", 
                    function(object, type, baselineRm,
                                showLimits, pdfFile = NULL, 
                                fileNames = NULL, colorBy = "rownames", 
                                normalizePrimariIon = FALSE, ...) {
    fileNamesObject <- getParameters(object)$listFile
    if (methods::is(fileNamesObject, "expression")) 
        fileNamesObject <- eval(fileNamesObject)
    set <- object
    # get list files
    if (is.null(fileNames)) {
        fileNames <- fileNamesObject
    } else {
        # check if fileNames are in the ptrset object
        test <- basename(fileNames) %in% basename(fileNamesObject)
        if (!all(test)) 
            stop("This file(s) are not in the directory: \n", 
                paste(basename(fileNames)[which(!test)], 
                collapse = "\n"))
        # put in full name
        if (any(dirname(fileNames) == ".")) {
            fileNames <- fileNamesObject[basename(fileNamesObject) %in% 
                                            basename(fileNames)]
        }
    }
    # get the TIC and time limit
    TIC <- getTimeInfo(set)$TIC[basename(fileNames)]
    indLim <- getTimeInfo(set)$timeLimit[basename(fileNames)]$exp
    # getSampleMetadata
    SMD <- getSampleMetadata(set)
    # test if colorBy is a colnames of samplemetadata
    if (colorBy != "rownames" & !(colorBy %in% colnames(SMD))) {
        message(colorBy, " is not a column of sample metadata")
        colorBy <- "rownames"
    }
    # remove baseline
    if (baselineRm) {
        TIC <- lapply(TIC, function(x) {
            bl <- try(baselineEstimation(x, d = 1))
            if (is.null(attr(bl, "condition"))) 
                TIC.blrm <- x - bl else TIC.blrm <- x - x[1]
            TIC.blrm
        })
    }
    if (!is.null(pdfFile)) {
        # check is pdf extension
        file_ext <- tools::file_ext(pdfFile)
        if (file_ext != "pdf") {
            warning("the plot will be saved in pdf extension")
            pdfFile <- gsub(paste0(".", file_ext, "$"), ".pdf", pdfFile)
        }
        pdf(pdfFile, width = 7 * 1.8, height = 5 * 1.8)
    }
    p <- ggplot2::ggplot()
    for (j in seq_along(TIC)) {
        ticPlot <- TIC[[j]]
        time <- as.numeric(names(ticPlot))
        if (normalizePrimariIon) 
            ticPlot <- ticPlot/getPTRInfo(set)$primaryIon[[j]]
        if (!is.null(pdfFile)) {
            plot <- ggplot2::qplot(x = time, y = ticPlot/(time[2] - time[1]), 
                                    xlab = "time", ylab = "Sum cps", 
                                    main = paste(j, names(TIC)[j], sep = " - "),
                                    size = 0.8)
            if (ncol(indLim[[j]]) > 0) 
                plot <- plot + 
                    ggplot2::geom_vline(ggplot2::aes(xintercept = as.numeric(names(TIC[[j]]))[c(indLim[[j]])],
                                                     color = "time limits")) + 
                    ggplot2::scale_fill_manual("Legend") +
                    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, 
                                                                      face = "bold"),
                                   axis.title = ggplot2::element_text(size = 16), 
                                   axis.text = ggplot2::element_text(size = 14),
                                   legend.text = ggplot2::element_text(size = 14), 
                                   legend.title = ggplot2::element_text(size = 16))
            print(plot)
        }
        if (colorBy == "rownames") {
            colour <- names(TIC)[j]
        } else {
            colour <- as.factor(SMD[names(TIC)[j], colorBy])
        }
        data <- data.frame(time = time, 
                            Sum_cps = ticPlot/(time[2] - time[1]), 
                            Legend = colour)
        p <- p + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$time, 
                                                            y = .data$Sum_cps,
                                                            color = .data$Legend), 
                                    data = data) + 
            ggplot2::geom_line(mapping = ggplot2::aes(x = .data$time, 
                                                     y = .data$Sum_cps, 
                                                     color = .data$Legend), 
                               data = data, size = 1)
    }
    
    if (!is.null(pdfFile)) 
        dev.off()
    
    if (showLimits) {
        if (colorBy == "rownames") {
            Legend <- names(TIC)
        } else {
            Legend <- as.factor(SMD[, colorBy])
        }
        
        
        limitdf <- lapply(Legend, 
                        function(j) {
                            t<-as.numeric(names(getTimeInfo(set)$TIC[[j]]))[c(getTimeInfo(set)$timeLimit[[j]]$exp)]
                            data.frame(x = t, Legend = rep(j, 2))}
            )
        limitdf <- Reduce(rbind, limitdf)
        p <- p + 
            ggplot2::geom_vline(mapping = ggplot2::aes(xintercept = .data$x, 
            color = .data$Legend), data = limitdf, size = 0.9)
    }
    # set title and legend
    p <- p + ggplot2::ggtitle(paste("TIC of", getParameters(set)$name)) + 
        ggplot2::theme_classic()
    switch(type, ggplot = return(p), plotly = return(plotly::ggplotly(p)))
}  #end function
)
#' Plot method for 'ptrRaw' objects
#'
#' Displays the image of the matrix of intensities, the TIC and the TIS,
#' for the selected m/z and time ranges
#' 
#' @param object An S4 object of class \code{ptrRaw-class} or \code{ptrSet}
#' @param mzRange Either a vector of 2 numerics indicating the m/z limits
#' or an integer giving a nominal m/z
#' @param timeRange Vector of 2 numerics giving the time limits
#' @param type Character: plot type; either 'classical' [default] or 'plotly'
#' @param ppm Integer: Half size of the m/z window when mzRange is set to a
#' nominal mass
#' @param palette Character: Color palette for the 'classical' plot; either 'heat'
#' [default], 'revHeat', 'grey', 'revGrey' or 'ramp'
#' @param showVocDB Logical: Should putative m/z annotations from the internal
#' package database be displayed (default is TRUE)
#' @param figure.pdf Character: Either 'interactive' [default], or the filename
#' of the figure to be saved (with the 'pdf' extension); only available for the
#' 'classical' display
#' @return Invisibly returns a list of the raw (sub)matrix 'rawsubM' and
#' the voc (sub)database 'vocsubDB'
#' @param fileNames vector of character. The file names of the ptrSer that you 
#' want to plot. Could be in 
#' basename or fullname.
#' @param ... not used
#' @rdname plotRaw
#' @examples 
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' ptairMS::plotRaw(exhaledPtrset ,mzRange=59,fileNames='ind1-1.h5')
#'
#' patientRaw <- ptairMS::readRaw(system.file('extdata/exhaledAir/ind1/ind1-1.h5',  
#' package = 'ptairData'),  mzCalibRef=c(21.022,59.049,75.05))
#' ptairMS::plotRaw(patientRaw, mzRange = 59)
#' ptairMS::plotRaw(patientRaw, mzRange = 59, type = 'plotly')
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline layout par title
#' @export 
methods::setMethod(f = "plotRaw", signature = "ptrSet", 
                   function(object, mzRange, timeRange = c(NA, NA), 
                            type = c("classical", "plotly")[1], ppm = 2000, 
                            palette = c("heat", "revHeat", "grey", "revGrey", 
                                        "ramp")[1], showVocDB = TRUE, 
                            figure.pdf = "interactive",
                            fileNames = NULL, ...) {
    set <- object
    fileNamesObject <- getParameters(set)$listFile
    if (methods::is(fileNamesObject, "expression")) 
        fileNamesObject <- eval(fileNamesObject)
    # get list files
    if (is.null(fileNames)) {
        fileNames <- fileNamesObject
    } else {
        # check if fileNames are in the ptrset object
        test <- basename(fileNames) %in% basename(fileNamesObject)
        if (!all(test)) 
            stop("This file(s) are not in the directory: \n", 
                 paste(basename(fileNames)[which(!test)], 
                collapse = "\n"))
        # put in full name
        if (any(dirname(fileNames) == ".")) {
            fileNames <- fileNamesObject[basename(fileNamesObject) %in% 
                                             basename(fileNames)]
        }
    }
    if (length(mzRange) == 1) 
        mzRange <- mzRange + c(-1, 1) * mzRange * ppm * 1e-06
    if (figure.pdf != "interactive") {
        if (type == "plotly") 
            stop("'plotly display is only available in the 'interactive' 
           mode currently.", 
                call. = FALSE)
        filenameSplitVc <- unlist(strsplit(basename(figure.pdf), ".", fixed = TRUE))
        extC <- utils::tail(filenameSplitVc, 1)
        if (extC == "pdf") {
            pdf(figure.pdf)
        } else stop("The extension of the 'figure.pdf' filename argument 
                    should be 'pdf'", 
            call. = FALSE)
    } else if (is.null(fileNames)) 
        graphics::par(ask = TRUE)
    for (file in fileNames) {
        # open a mz range
        mz <- rhdf5::h5read(file, name = "FullSpectra/MassAxis")
        mzRange[1] <- max(min(mz), mzRange[1], na.rm = TRUE)
        mzRange[2] <- min(max(mz), mzRange[2], na.rm = TRUE)
        # get the index of the aroud mzRange
        index <- which(mzRange[1] < mz & mz < mzRange[2])
        # object data
        object <- rhdf5::h5read(file, name = "/FullSpectra/TofData", 
                                index = list(index, NULL, NULL, NULL))
        
        # time axis
        timeVn <- as.numeric(names(getTimeInfo(set)$TIC[[basename(file)]]))
        timeiNTER <- (timeVn[3] - timeVn[2])
        timeRange.file <- rep(0, 0)
        timeRange.file[1] <- max(min(timeVn), timeRange[1], na.rm = TRUE)
        timeRange.file[2] <- min(max(timeVn), timeRange[2], na.rm = TRUE)
        indexTime <- which(timeVn >= timeRange.file[1] & 
                               timeVn <= timeRange.file[2])
        timeVn <- timeVn[indexTime]
        # calibrate mass axis
        
        FirstcalibCoef <- try(rhdf5::h5read(file, "FullSpectra/MassCalibration", index = list(NULL, 
                                                                                             1)))
        attributCalib <- try(rhdf5::h5readAttributes(file, "/FullSpectra"))
        if (!is.null(attr(FirstcalibCoef, "condition")) & is.null(attr(attributCalib, "condition"))) {
            FirstcalibCoef <- matrix(c(attributCalib$`MassCalibration a`, attributCalib$`MassCalibration b`), 
                                ncol = 1)
        }
        tof <- sqrt(mz) * FirstcalibCoef[1, 1] + FirstcalibCoef[2, 1]
        # tof<- seq(0,length(mz))
        coefCalib <- getCalibrationInfo(set)$coefCalib[[basename(file)]][[1]]
        mzVn <- ((tof - coefCalib["b", ])/coefCalib["a", ])^2
        mzVn <- mzVn[index]
        # formate object matix
        rawSubMN <- matrix(object, nrow = dim(object)[1], 
                           ncol = prod(utils::tail(dim(object), 2)))
        rawSubMN <- rawSubMN[, indexTime]
        dimnames(rawSubMN) <- list(mzVn, timeVn)
        if (showVocDB) {
            vocdbDF <- .loadVocDB()
            vocdb_sel.vl <- vocdbDF[, "ion_mass"] >= mzRange[1] & 
                vocdbDF[, "ion_mass"] <= mzRange[2]
            if (sum(vocdb_sel.vl)) {
                vocdbDF <- vocdbDF[vocdb_sel.vl, , drop = FALSE]
            } else vocdbDF <- NULL
        } else vocdbDF <- NULL
        switch(type, classical = {
            imageMN <- t(rawSubMN)[, seq_len(nrow(rawSubMN)), drop = FALSE]
            rownames(imageMN) <- round(as.numeric(rownames(imageMN)))
            colnames(imageMN) <- round(as.numeric(colnames(imageMN)), 4)
            paletteVc <- .palette(palette = palette)
            currentParLs <- graphics::par()
            for (parC in c("cin", "cra", "csi", "cxy", "din", "page")) currentParLs[[parC]] <- NULL
            marLs <- list(chr = c(0.6, 4.1, 4, 0), sca = c(0.6, 0.6, 4, 7.1), ima = c(4.1, 
                4.1, 0, 0), spe = c(4.1, 0.6, 0, 1.1))
            graphics::par(font = 2, font.axis = 2, font.lab = 2, pch = 18)
            layout(matrix(c(1, 2, 3, 4), byrow = TRUE, nrow = 2), heights = c(2, 
                5), widths = c(5, 2))
            ## chr: Chromatogram
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
                abline(h = vapply(vocdbDF[, "ion_mass"], 
                                function(mzN) (mzN - min(mzImaVn))/
                                    diff(range(mzImaVn)) * ncol(imageMN) +
                                    par("usr")[1], 
                                FUN.VALUE = 1.1), lty = "dotted")
            }
            ## spe: Spectrum
            graphics::par(mar = marLs[["spe"]])
            specVn <- apply(rawSubMN, 1, function(intVn) mean(intVn, na.rm = TRUE)/timeiNTER)
            plot(specVn, as.numeric(names(specVn)), cex = 0.7, pch = 16, xlab = "", 
                ylab = "", xaxs = "i", yaxs = "i", yaxt = "n")
            graphics::mtext("Count Per Second (cps)", cex = 0.8, side = 1, line = 2.5)
            if (showVocDB && !is.null(vocdbDF)) {
                abline(h = vocdbDF[, "ion_mass"], lty = "dotted")
            }
            title(main = basename(file), outer = TRUE, line = -2, cex.main = 2)
            graphics::par(currentParLs)
        }, plotly = {
            plotlyBuild <- function() {
                p <- plotly::subplot(plotly::plot_ly(x = as.numeric(colnames(rawSubMN)), 
                                                    y = colSums(rawSubMN), 
                                                    type = "scatter", 
                                                    mode = "lines"),
                                    plotly::plotly_empty(), 
                                    plotly::plot_ly(z = rawSubMN, 
                                                    x = as.numeric(colnames(rawSubMN)),
                                                    y = as.numeric(rownames(rawSubMN)), 
                                                    type = "heatmap"), 
                                    plotly::plot_ly(x = rowSums(rawSubMN),
                                                    y = as.numeric(rownames(rawSubMN)), 
                                                    type = "scatter", 
                                                    mode = "lines"),
                                    nrows = 2, heights = c(0.2, 0.8), 
                                    widths = c(0.8, 0.2), shareX = TRUE, 
                                    shareY = TRUE, titleX = FALSE,
                                    titleY = FALSE)
                `%>%` <- plotly::`%>%`
                p <- p %>% plotly::layout(title = basename(file))
                return(p)
            }
            p <- suppressWarnings(suppressMessages(plotlyBuild()))
            p <- plotly::layout(p, showlegend = FALSE)
            print(plotly::ggplotly(p))
        })
    }  #end loop file
    # reset parameter plot
    graphics::par(ask = FALSE)
    # close pdf
    if (figure.pdf != "interactive") 
        dev.off()
    # print voc ref
    if (!is.null(vocdbDF)) {
        vocdbDF <- vocdbDF[nrow(vocdbDF):1, , drop = FALSE]
        print(vocdbDF[, c("ion_mass", "ion_formula", "name_iupac"), drop = FALSE])
    }
})
#' @rdname plotFeatures
#' @export
methods::setMethod(f = "plotFeatures", signature = "ptrSet", 
                function(set, mz, 
                        typePlot,
                        addFeatureLine, ppm, pdfFile, 
                        fileNames, colorBy) {
    
                    fileNamesObject <- getParameters(set)$listFile
    if (methods::is(fileNamesObject, "expression")) 
        fileNamesObject <- eval(fileNamesObject)
    
    # get list files
    if (is.null(fileNames)) {
        fileNames <- fileNamesObject
    } else {
        # check if fileNames are in the ptrSet object put in basename
        test <- basename(fileNames) %in% basename(fileNamesObject)
        if (!all(test)) 
            stop("This file(s) names are not in 
                                     the directory: \n", 
                paste(!basename(fileNames)[test], collapse = "\n"))
        # Put in fullNames
        if (any(dirname(fileNames) == ".")) {
            fileNames <- fileNamesObject[which(basename(fileNamesObject) %in% fileNames)]
        }
    }
    
    # getSampleMetadata
    SMD <- getSampleMetadata(set)
    
    # test if colorBy is a colnames of samplemetadata
    if (colorBy != "rownames" & !(colorBy %in% colnames(SMD))) {
        message(colorBy, " is not a column of sample metadata")
        colorBy <- "rownames"
    }
    
    # get the extension file
    if (!is.null(pdfFile)) {
        # check is pdf extension
        file_ext <- tools::file_ext(pdfFile)
        if (file_ext != "pdf") {
            warning("the plot will be saved in pdf extension")
            pdfFile <- gsub(paste0(".", file_ext, "$"), ".pdf", pdfFile)
        }
        pdf(file = pdfFile, width = 7 * 1.8, height = 5 * 1.8)
    }
    
    # loop over file
    plotAll <- ggplot2::ggplot()
    listPlotFile <- list()
    j <- 0
    for (file in fileNames) {
        j <- j + 1
        # get the mass axis, must be the same for all files
        mzAxis <- rhdf5::h5read(file, name = "FullSpectra/MassAxis")
        # get the index of the mzAxis
        thLarge <- max(0.4, mz * (ppm/2)/10^6)
        indexMz <- which(mz - thLarge < mzAxis & mzAxis < mz + thLarge)
        indLim <- getTimeInfo(set)$timeLimit[[basename(file)]]$exp
        indLimBg <- getTimeInfo(set)$timeLimit[[basename(file)]]$backGround
        n.limit <- dim(indLim)[2]
        raw <- rhdf5::h5read(file, name = "/FullSpectra/TofData", 
                            index = list(indexMz, NULL, NULL, NULL))
        time <- c(rhdf5::h5read(file, name = "/TimingData/BufTimes"))
        index_zero <- which(time == 0)[-1]
        if (length(index_zero)) 
            time <- time[-index_zero]
        mzAxis.j <- mzAxis[indexMz]
        rawMn <- matrix(raw, nrow = dim(raw)[1], ncol = prod(utils::tail(dim(raw), 
            2)))
        # * 0.2 ns / 2.9 (single ion signal) if convert to cps
        
        FirstcalibCoef <- try(rhdf5::h5read(file, "FullSpectra/MassCalibration", index = list(NULL, 
                                                                                              1)))
        attributCalib <- try(rhdf5::h5readAttributes(file, "/FullSpectra"))
        if (!is.null(attr(FirstcalibCoef, "condition")) & is.null(attr(attributCalib, "condition"))) {
            FirstcalibCoef <- matrix(c(attributCalib$`MassCalibration a`, attributCalib$`MassCalibration b`), 
                                     ncol = 1)
        }
        
        tof <- sqrt(mzAxis.j) * FirstcalibCoef[1, 1] + FirstcalibCoef[2, 1]
        coefCalib <- getCalibrationInfo(set)$coefCalib[[basename(file)]][[1]]
        mzNew <- ((tof - coefCalib["b", ])/coefCalib["a", ])^2
        # get smaller windows
        th <- mz * (ppm/2)/10^6
        indexSub <- which(mz - th < mzNew & mzNew < mz + th)
        # plot mean spectrum for expriration/time periods
        plotFile <- ggplot2::ggplot()
        indexTimeVec <- NULL
        for (i in seq_len(n.limit)) {
            # get index of the time period
            indexTime <- seq(indLim["start", i], indLim["end", i])
            indexTimeVec <- c(indexTimeVec, indexTime)
            spectrum <- rowSums(rawMn[, indexTime])/(
                ncol(rawMn[, indexTime]) * (time[3] - time[2]))
            data <- data.frame(mz = mzNew[indexSub], cps = spectrum[indexSub], 
                                timePeriods = as.character(i))
            plotFile <- plotFile + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$mz, 
                                                                            y = .data$cps, 
                                                                            color = .data$timePeriods), 
                                                        data = data) + 
                ggplot2::geom_line(mapping = ggplot2::aes(x = .data$mz, 
                                                        y = .data$cps, 
                                                        color = .data$timePeriods), 
                                    data = data, size = 1)
        }
        # get the background
        if (!is.null(indLimBg)) {
            background <- rowSums(rawMn[, indLimBg])/
                (length(indLimBg) * (time[3] - time[2]))
            data = data.frame(mz = mzNew[indexSub], cps = background[indexSub], 
                            timePeriods = "background")
            
            # plot background to the file plot
            plotFile <- plotFile + 
                ggplot2::geom_point(mapping = ggplot2::aes(x = .data$mz,
                                                            y = .data$cps,
                                                            color = .data$timePeriods), 
                                    data = data) + 
                ggplot2::geom_line(mapping = ggplot2::aes(x = .data$mz,
                                                          y = .data$cps, 
                                                          color = .data$timePeriods), 
                                    data = data, linetype = "dashed")
        }
        
        plotFile <- plotFile + ggplot2::ggtitle(basename(file))
        listPlotFile[[j]] <- plotFile
        
        ## add summary line to plotAll get the average time periods spectrum
        spectrum <- rowSums(rawMn[, indexTimeVec])/(length(indexTimeVec) * (time[3] - 
            time[2]))
        
        # spline intrepolation for spectrum
        splineInterpol <- stats::splinefun(mzNew, spectrum)
        
        # plot spectrum
        if (colorBy == "rownames") {
            colour <- basename(file)
        } else {
            colour <- SMD[basename(file), colorBy]
        }
        
        data = data.frame(mz = mzNew[indexSub], cps = spectrum[indexSub], 
                        Legend = rep(colour,length(indexSub)))
        
        plotAll <- plotAll + 
            ggplot2::geom_point(ggplot2::aes(x = .data$mz, 
                                            y = .data$cps, 
                                            color = .data$Legend), 
                                data = data) + 
            ggplot2::stat_function(mapping = ggplot2::aes(x = .data$mz, 
                                                        color = .data$Legend), 
                                    data = data, fun = splineInterpol, 
                                    n = 1000, size = 1)
        
        # spline interpolation for background
        if (!is.null(indLimBg)) {
            background <- rowSums(rawMn[, indLimBg])/
                (length(indLimBg) * (time[3] - time[2]))
            splineInterpol <- stats::splinefun(mzNew, background)
            
            # plot background
            data = data.frame(mz = mzNew[indexSub], 
                                Legend = rep(colour, length(indexSub)))
            
            plotAll <- plotAll + 
                ggplot2::stat_function(mapping = ggplot2::aes(x = .data$mz,
                                                            color = .data$Legend), 
                                        data = data, 
                                        fun = splineInterpol, n = 1000, 
                                        linetype = "dashed", size = 1)
        }
        # plot vertical line line
        if (addFeatureLine) 
            plotAll <- plotAll + ggplot2::geom_vline(xintercept = mz)
    }  #END loop file
    plotAll <- plotAll + ggplot2::ggtitle(paste("Features", round(mz, 3), 
                                                "of", getParameters(set)$name), 
        subtitle = paste(annotateVOC(mz)[2:3], collapse = " ")) + 
        ggplot2::labs(color = colorBy) + 
        ggplot2::theme_classic()
    if (!is.null(pdfFile)) {
        for (j in seq_along(listPlotFile)) print(listPlotFile[[j]])
        dev.off()
    }
    switch(typePlot, ggplot = return(plotAll), 
            plotly = return(plotly::ggplotly(plotAll)))
})
## sampleMetadata ----
#' reset the default sampleMetadata, containing the suborders names and the 
#' acquisition dates as columns. 
#' @param ptrset a ptrser object
#' @return a data.frame 
#' @export
#' @examples 
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' SMD<- resetSampleMetadata(exhaledPtrset)
resetSampleMetadata <- function(ptrset) {
    dir <- getParameters(ptrset)$dir
    fileNamesObject <- getParameters(ptrset)$listFile
    if (methods::is(dir, "expression")) {
        dir <- eval(dir)
        fileNamesObject <- eval(fileNamesObject)
    }
    filesFullName <- list.files(dir, recursive = TRUE, 
                                pattern = "\\.h5$", full.names = TRUE)
    fileDir <- dirname(list.files(dir, 
                                    recursive = TRUE, pattern = "\\.h5$"))
    filesProcessed <- basename(fileNamesObject)
    # delete new files
    newFilesIndex <- which(!basename(filesFullName) %in% filesProcessed)
    if (length(newFilesIndex > 0)) {
        filesFullName <- filesFullName[-newFilesIndex]
        fileDir <- fileDir[-newFilesIndex]
    }
    fileName <- basename(filesFullName)
    # test if there is subfolder if there is no subfolder
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
        sampleMetadata <- data.frame(subfolder = group, 
                                    row.names = fileName, 
                                    stringsAsFactors = FALSE)
    }
    sampleMetadata$date <- Reduce(c, getDate(ptrset))
    return(sampleMetadata)
}
#' get sampleMetadata of a ptrSet
#' @param set a ptrSet object
#' @return a data.frame 
#' @export
#' @examples 
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL)
#' SMD<-getSampleMetadata(exhaledPtrset)
getSampleMetadata <- function(set) {
    if (!methods::is(set, "ptrSet")) 
        stop("set is not a ptrSet object")
    sampleMetadata <- set@sampleMetadata
    return(sampleMetadata)
}
#' set sampleMetadata in a ptrSet
#' 
#' Insert a samplemetada data.frame in a ptrSet object. The dataframe must 
#' have all file names in rownames. 
#' @param set a ptrSet object
#' @param sampleMetadata a data.frame with all file names of the ptrSet in 
#' row names
#' @return the ptrSet object in argument with the sampleMetadata modified
#' @export
#' @examples 
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' SMD<-getSampleMetadata(exhaledPtrset )
#' colnames(SMD)[1]<-'individual'
#' exhaledPtrset<-setSampleMetadata(exhaledPtrset ,SMD)
setSampleMetadata <- function(set, sampleMetadata) {
    # check if set is ptrSet
    if (!methods::is(set, "ptrSet")) 
        stop("set is not a ptrSet object")
    # check if row names contains all files
    files <- names(set@mzCalibRef)
    fileName <- basename(files)
    testFilesName <- fileName %in% row.names(sampleMetadata)
    if (!all(testFilesName)) {
        stop(fileName[!testFilesName], "not in sampleMetadata 
                          row names, please complete theme \n")
    }
    set@sampleMetadata <- sampleMetadata
    ## save
    if (!is.null(getParameters(set)$saveDir)) {
        changeName <- parse(text = paste0(getParameters(set)$name, "<- set "))
        eval(changeName)
        eval(parse(text = paste0("save(", getParameters(set)$name, ",file= 
                                         paste0( getParameters(set)$saveDir,'/', '", 
                                 getParameters(set)$name, ".RData '))")))
    }
    return(set)
}
#' export sampleMetadata
#' @param set a ptrSet object
#' @param saveFile a file path in tsv extension where the data.frame will be 
#' exported
#' @return nothing
#' @export
#' @examples 
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' saveFile<-file.path(getwd(),'sampleMetadata.tsv')
#' #exportSampleMetada(exhaledPtrset ,saveFile)
exportSampleMetada <- function(set, saveFile) {
    if (!methods::is(set, "ptrSet")) 
        stop("set is not a ptrSet object")
    if (tools::file_ext(saveFile) != "tsv") 
        stop("saveFile must be in 
                                                        .tsv extension")
    # get sampleMetadata
    sampleMetadata <- getSampleMetadata(set)
    # write the data.frame in a tsv file
    utils::write.table(sampleMetadata, file = saveFile, sep = "\t", 
                        col.names = NA)
}
#' import a sampleMetadata from a tsv file to a ptrSet object 
#' @param set a ptrSet object
#' @param file a tsv file contains the sample metada to import, with all 
#' file names in row name 
#' (the fist column on th excel).
#' @return a ptrSet with th enew sample Metadata
#' @export
#' @examples 
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' saveFile<-file.path(getwd(),'sampleMetadata.tsv')
#' #exportSampleMetada(exhaledPtrset ,saveFile)
#' #exhaledPtrset<-importSampleMetadata(exhaledPtrset ,saveFile)
importSampleMetadata <- function(set, file) {
    if (!methods::is(set, "ptrSet")) 
        stop("set is not a ptrSet object")
    sampleMetadata <- try(utils::read.table(file = file, sep = "\t", 
                                            header = TRUE, 
        row.names = 1, quote = ""))
    # check if row names contains all files
    dir <- getParameters(set)$dir
    if (methods::is(dir, "expression")) 
        dir <- eval(dir)
    files <- list.files(dir, recursive = TRUE, pattern = "\\.h5$")
    fileName <- basename(files)
    testFilesName <- fileName %in% row.names(sampleMetadata)
    if (!all(testFilesName)) {
        stop(fileName[!testFilesName], "not in sampleMetadata 
                            row names, please complete theme \n")
    }
    set<- setSampleMetadata(set = set,sampleMetadata = sampleMetadata)
    return(set)
}
## time limits ----
#' @rdname timeLimits
#' @export
methods::setMethod(f = "timeLimits", signature = "ptrSet", 
                    function(object, fracMaxTIC = 0.5,fracMaxTICBg = 0.2, 
                            derivThresholdExp = 1, derivThresholdBg = 0.05,
                            minPoints = 2, 
                            degreeBaseline = 1, baseline = TRUE, 
                            redefineKnots = TRUE, plotDel = FALSE) {
    fileNames <- getParameters(object)$listFile
    if (methods::is(fileNames, "expression")) 
        fileNames <- eval(fileNames)
    fileNames <- basename(fileNames)
    for (file in fileNames) {
        TIC <- getTimeInfo(object)$breathTracer[[file]]
        indLim <- timeLimitFun(TIC, fracMaxTIC = fracMaxTIC, 
                                fracMaxTICBg = fracMaxTICBg,
                                derivThresholdExp = derivThresholdExp, 
                                derivThresholdBg = derivThresholdBg,
                                degreeBaseline = degreeBaseline, 
                                baseline = baseline, 
                                minPoints = minPoints, 
                               plotDel = plotDel)
        
        
        object<-setTimeLimits(object,newTimeLimites =indLim,index = file )
    }
    
    paramterTimeLimit <- list(fracMaxTIC = fracMaxTIC, 
                            fracMaxTICBg = fracMaxTICBg, 
                            derivThresholdExp = derivThresholdExp, 
                            derivThresholdBg = derivThresholdBg, 
                            degreeBaseline = degreeBaseline, 
                            minPoints = minPoints)
    
    parameter<-getParameters(object)
   parameter$timeLimit <- paramterTimeLimit
   object<-setParameters(object,parameter) 
   
   
    if (redefineKnots) 
        object <- defineKnots(object, knotsPeriod = getParameters(object)$knotsPeriod)
    saveDir <- getParameters(object)$saveDir
    objName <- getParameters(object)$name
    
    if (!is.null(saveDir)) {
        if (!is.null(objName)) {
            changeName <- parse(text = paste0(objName, "<- object "))
            eval(changeName)
            eval(parse(text = paste0("save(", objName, ",file= paste0( 
                                           saveDir,'/', '", 
                objName, ".RData '))")))
        } else save(object, file = paste0(saveDir, "/ptrSet.RData"))
    }
    return(object)
})

## define knots-----
#' Define the knots location
#' 
#' \code{defineKnots} function determine the knots location for a ptrSet or 
#' ptrRaw object.
#' There is three possibilities :  
#' \itemize{
#' \item \code{method = expiration} in the expiration periods, a knots is 
#' placed every \code{knotsPeriod} seconds, and 1 knots in the middle of two 
#' expiration, one at begin and at the end
#' \item \code{method= uniforme}, the knots are placed uniformly every 
#' \code{knotsPeriod} time points
#' \item give in \code{knotsList} a list of knot, with all base name file in 
#' name of the list element.  All file must be informed. The knots location 
#' must be contained in the time axis
#' }
#' @param  object ptrSet object
#' @param knotsPeriod the period in second (times scale) between two knots for 
#' the two dimensional modelization
#' @param method expiration or uniform
#' @param knotsList a list of knot location for each files, with all base name 
#' file in name of the list element
#' @return a list with numeric vector of knots for each file
#' @examples 
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#'
#' #### placed knots every 2 times points
#' exhaledPtrset <- defineKnots(exhaledPtrset ,knotsPeriod=2,method='uniform')
#' 
#' #### placed knots every 3 times points in the expiration (default)
#' exhaledPtrset <- defineKnots(exhaledPtrset ,knotsPeriod=3,method='expiration')
#' @rdname defineKnots
#' @export
methods::setMethod(f = "defineKnots", signature = "ptrSet", 
                    function(object, knotsPeriod = 3, 
                            method = c("expiration", "uniform")[1], 
                            knotsList = NULL) {
    if (knotsPeriod == 0) {
        knots <- list(NULL)
    } else {
        if (!is.null(knotsList)) {
            # check if all file are set
            if (any(!names(knotsList) %in% names(getTimeInfo(object)$TIC))) 
                stop(names(knotsList)[!names(knotsList) %in% names(getTimeInfo(object)$TIC)], 
                  " missing")
            # check if interior knots are in the times axis
            test <- apply(names(knotsList), function(file) {
                knots <- knotsList[[file]]
                t <- as.numeric(names(getTimeInfo(object)$TIC[[file]]))
                return(knots[1] >= t[1] & utils::tail(knots, 1) <= utils::tail(t, 
                  1))
            })
            if (any(!test)) 
                stop("knots are not contained in 
                                                   the time axis for file ", 
                  names(knotsList)[which(!test)], "\n")
            knots <- knotsList
        } else {
            knots <- lapply(names(getTimeInfo(object)$TIC), function(file) {
                background <- getTimeInfo(object)$timeLimit[[file]]$backGround
                t <- as.numeric(names(getTimeInfo(object)$TIC[[file]]))
                res <- try(defineKnotsFunc(t, background, knotsPeriod, method, file))
                if (!is.null(attr(res, "condition"))) 
                  res <- NULL
                res
            })
        }
        names(knots) <- names(getTimeInfo(object)$TIC)
    }
        
                                 
    object <- setKnots(object,knots)
    parameter<-getParameters(object)
    saveDir <- parameter$saveDir
    objName <- parameter$name
    parameter$knotsPeriod <- knotsPeriod
    object<-setParameters(object,parameter)    
   
    if (!is.null(saveDir)) {
        if (!is.null(objName)) {
            changeName <- parse(text = paste0(objName, "<- object "))
            eval(changeName)
            eval(parse(text = paste0("save(", objName, ",file= 
                                                    paste0( saveDir,'/', '", 
                objName, ".RData '))")))
        } else save(object, file = paste0(saveDir, "/ptrSet.RData"))
    }
    return(object)
})

## calibration----
#' @rdname calibration
#' @export 
methods::setMethod(f = "calibration", signature = "ptrSet", 
                   function(x, mzCalibRef = c(21.022, 
    29.013424, 41.03858, 75.04406, 203.943, 330.8495), 
    calibrationPeriod = 60, tol = 70) {
                       
    fileNames <- getParameters(x)$listFile
    
    if (methods::is(fileNames, "expression")) 
        fileNames <- eval(fileNames)
    
    for (file in fileNames) {
        raw <- readRaw(file, mzCalibRef = mzCalibRef, 
                        calibrationPeriod = calibrationPeriod)
        calibSpectr <- alignCalibrationPeak(getCalibrationInfo(raw)$calibSpectr, 
                                            getCalibrationInfo(raw)$calibMassRef, 
                                            length(getRawInfo(raw)$time))
        resolutionEstimated <- estimateResol(calibMassRef = getCalibrationInfo(raw)$calibMassRef, 
                                             calibSpectr = calibSpectr)
        calibList<-list(coefCalib=getCalibrationInfo(raw)$calibCoef,
                        mzCalibRef=getCalibrationInfo(raw)$calibMassRef,
                        errorCalibPpm=getCalibrationInfo(raw)$calibError,
                        signalCalibRef=calibSpectr,
                        resolution=resolutionEstimated,
                        indexTimeCalib=raw@indexTimeCalib)
        
        x<-setCalibration(x,calibList,index=basename(file))
    }
    
    parameter<-getParameters(x)
    parameter$mzCalibRef <- mzCalibRef
    x<-setParameters(x,parameter)
    saveDir <- parameter$saveDir
    objName <- parameter$name
    if (!is.null(saveDir)) {
        if (!is.null(objName)) {
            changeName <- parse(text = paste0(objName, "<- x "))
            eval(changeName)
            eval(parse(text = paste0("save(", objName, ",file= paste0( 
                                           saveDir,'/', '", 
                objName, ".RData '))")))
        } else save(x, file = paste0(saveDir, "/ptrSet.RData"))
    }
    return(x)
})

## other ----
#' get the files directory of a ptrSet
#' @param ptrSet ptrSet object 
#' @return the directory in absolute path as character
#' @examples 
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' getDirectory(exhaledPtrset )
#' @export
getDirectory <- function(ptrSet) {
    dir <- getParameters(ptrSet)$dir
    if (methods::is(dir, "expression")) 
        dir <- eval(dir)
    return(dir)
}
#' remove the peakList of an ptrSet object 
#' 
#' This function is useful when you want to change the parameters of the detect 
#' peak function. First delete the peakList with \code{rmPeakList}, and apply 
#' \code{detectPeak}with the new parameters.
#' @param object ptrSet object 
#' @return a ptrSet
#' @export
#' @examples 
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' exhaledPtrset <-rmPeakList(exhaledPtrset )
rmPeakList <- function(object) {
    object@peakList <- list()
    parameter<-getParameters(object)
    parameter$detectPeakParam <- NULL
    object<-setParameters(object,parameter)
    return(object)
}

#' show a ptrSet object 
#' 
#' It indicates the directory, the numbre of files that contain the directory at 
#' the moment, and the number of processed files.
#' The two numbers are different, use \code{updatePtrSet} function.
#' @param object a ptrSet object
#' @return nothing
#' @export 
methods::setMethod("show", "ptrSet", function(object) {
    dir <- getParameters(object)$dir
    fileCheck <- getTimeInfo(object)$TIC
    if (methods::is(dir, "expression")) {
        dir <- eval(dir)
    }
    nFiles <- length(list.files(dir, recursive = TRUE, pattern = "\\.h5$"))
    nFilesCheck <- length(fileCheck)
    nFilesProcess <- length( getPeakList(object) )
    cat("ptrSet object :", getParameters(object)$name, "\n")
    cat("directory:", dir, "\n")
    cat("   ", nFiles, "files contained in the directory \n")
    cat("   ", nFilesCheck, "files checked", "\n")
    cat("   ", nFilesProcess, "files processed by detectPeak", "\n")
    if (!is.null(getParameters(object)$saveDir)) 
        cat("object save in ", getParameters(object)$saveDir)
})
#' get the peak list of a ptrSet object 
#' @param set ptrSet object 
#' @return  a list of expressionSet, where:
#'\itemize{
#'\item assay Data contains the matrix with m/z peak center in row names, and the 
#'quantification in cps at each time point
#'\item feature Data the matrix with  m/z peak center in row names, and the 
#'following columns: 
#'   \itemize{
#'   \item quanti_unit: the mean of the quantification over the expiration/headspace 
#'time limits defined 
#'    \item background_unit: the mean of the quantification over the background time 
#'limits defined 
#'   \item diffAbs_unit: the mean of the quantification over the expiration/headspace
#'time limits defined after subtracting the baseline estimated from the background 
#'points defined 
#'    \item pValLess/ pValGreater: the p-value of the unilateral t-test, who test 
#'that quantification (in cps) of expiration points are less/greater than the 
#'intensity of the background. 
#'    \item lower/upper: integration boundaries
#'    \item parameter peak: estimated peak parameter 
#'   }
#'} 
#' @examples 
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' exhaledPtrset<-detectPeak(exhaledPtrset,mz=c(21,59))
#' peakList<- getPeakList(exhaledPtrset )
#' X<-Biobase::exprs(peakList[[1]])
#' Y<- Biobase::fData(peakList[[1]])
#' head(Y)
#' @export
getPeakList <- function(set) {
    return(set@peakList)
}
#' get the file names containing in the directory of a ptrSet or ptrRaw
#' @param object ptrSet object 
#' @param fullNames logical: if \code{TRUE}, it return the the directory path is 
#' prepended to the file names.
#' @return a vector of character that contains all file names
#' @examples 
#' library(ptairData)
#' dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")
#' exhaledPtrset <- createPtrSet(dir=dirRaw, setName="exhaledPtrset", 
#' mzCalibRef = c(21.022, 60.0525), fracMaxTIC = 0.7, saveDir = NULL )
#' getFileNames(exhaledPtrset )
#' @rdname getFileNames
#' @export
methods::setMethod("getFileNames", signature = "ptrSet", 
                   function(object, fullNames) {
                       
    fileFullNames <- getParameters(object)$listFile
    if (methods::is(fileFullNames, "expression")) 
        fileFullNames <- eval(fileFullNames)
    if (fullNames) 
        return(fileFullNames) else return(basename(fileFullNames))
})



#' get the file names containing in the directory of a ptrSet
#' @rdname getFileNames
#' @export
methods::setMethod("getFileNames", signature = "ptrRaw", 
                   function(object, fullNames) {
                       return(object@name)
                   })


methods::setMethod("getParameters", signature = "ptrSet", 
                   function(object) {
                       
                      return(object@parameter)
                   })

methods::setMethod("setParameters", signature = "ptrSet", 
                   function(object,newparam) {
                       object@parameter<-newparam
                       return(object)
                   })


methods::setMethod("getDate", signature = "ptrSet", 
                   function(object) {
                       
                       return(object@date)
                   })

methods::setMethod("getDate", signature = "ptrRaw", 
                   function(object) {
                       
                       return(object@date)
                   })



methods::setMethod("getCalibrationInfo", signature = "ptrSet", 
                   function(object) {
                       return(list(mzCalibRef=object@mzCalibRef,
                                   signalCalibRef = object@signalCalibRef,
                                   errorCalibPpm =object@errorCalibPpm,
                                   coefCalib=object@coefCalib,
                                   indexTimeCalib=object@indexTimeCalib)
                                   )
                   })

methods::setMethod("getCalibrationInfo", signature = "ptrRaw", 
                   function(object) {
                       return(list(
                                    calibMzToTof = object@calibMzToTof,
                                    calibToftoMz = object@calibToftoMz,
                                    calibCoef=object@calibCoef,
                                    calibMassRef = object@calibMassRef,
                                    calibError=object@calibError,
                                    calibSpectr = object@calibSpectr,
                                    indexTimeCalib=object@indexTimeCalib)
                       )
                   })


methods::setMethod("setCalibration", signature = "ptrRaw", 
                   function(object,newCalibrationInfo) {
    object@indexTimeCalib <- newCalibrationInfo$indexTimeCalib
    object@calibSpectr <- newCalibrationInfo$calibSpectr
    object@calibError <- newCalibrationInfo$calibError
    object@calibCoef <-newCalibrationInfo$calibCoef
    return(object)
                   })

methods::setMethod("setCalibration", signature = "ptrSet", 
                   function(object,newCalibrationInfo,index=NULL) {
                        if(is.null(index)) index<- seq_along(newCalibrationInfo)
                        if( length(index == 1) ){
                                object@coefCalib[[index]] <- newCalibrationInfo$coefCalib
                                object@mzCalibRef[[index]] <- newCalibrationInfo$mzCalibRef
                                object@errorCalibPpm[[index]] <- newCalibrationInfo$errorCalibPpm
                                object@signalCalibRef[[index]] <- newCalibrationInfo$signalCalibRef
                                object@resolution[[index]] <- newCalibrationInfo$resolution
                                object@indexTimeCalib[[index]] <-newCalibrationInfo$indexTimeCalib
                                
                        } else{
                            object@coefCalib[index] <- newCalibrationInfo$coefCalib
                            object@mzCalibRef[index] <- newCalibrationInfo$mzCalibRef
                            object@errorCalibPpm[index] <- newCalibrationInfo$errorCalibPpm
                            object@signalCalibRef[index] <- newCalibrationInfo$signalCalibRef
                            object@resolution[index] <- newCalibrationInfo$resolution
                            object@indexTimeCalib[index] <-newCalibrationInfo$indexTimeCalib
                        }
                        
                       return(object)
                   })





methods::setMethod("getPeaksInfo", signature = "ptrSet", 
                   function(object) {
                       return(list(knots=object@knots,
                       fctFit=object@fctFit,
                       peakShape=object@peakShape))
                   })

methods::setMethod("getPeaksInfo", signature = "ptrRaw", 
                   function(object) {
                       return(list(
                                   peakShape=object@peakShape))
                   })

methods::setMethod("getTimeInfo", signature = "ptrSet", 
                   function(object) {
                       return(list(TIC =object@TIC,
                       breathTracer=object@breathTracer,
                       timeLimit=object@timeLimit))
                   })


methods::setMethod("setTimeLimits", signature = "ptrSet",
                   function(object,newTimeLimites,index=NULL) {
                       if(is.null(index)){
                           object@timeLimit <-  newTimeLimites
                       } else if(length(index) == 1 ) {
                           object@timeLimit[[index]] <-  newTimeLimites 
                           } else object@timeLimit[index] <-  newTimeLimites 
                        return(object)
                   })

methods::setMethod("setKnots", signature = "ptrSet",
                   function(object,newKnots,index=NULL) {
                       if(is.null(index)){
                           object@knots <-  newKnots
                       } else if(length(index) == 1 ) {
                           object@knots[[index]] <-  newKnots 
                       } else object@knots[index] <-  newKnots 
                       return(object)
                   })



methods::setMethod("getPTRInfo", signature = "ptrSet", 
                   function(object) {
                       
                       return(list(primaryIon=object@primaryIon,
                       ptrTransmisison = object@ptrTransmisison,
                       prtReaction = object@prtReaction,
                       resolution = object@resolution))
                   })

methods::setMethod("getPTRInfo", signature = "ptrRaw", 
                   function(object) {
                       
                       return(list(
                                   ptrTransmisison = object@ptrTransmisison,
                                   prtReaction = object@prtReaction))
                   })


methods::setMethod("getRawInfo", signature = "ptrRaw", 
                   function(object) {
                       return(list(
                           rawM= object@rawM,
                           mz =  object@mz,
                           time = object@time))
                   })


methods::setMethod("setPeakShape", signature = "ptrRaw", 
                   function(object,peakShape) {
                      object@peakShape<-peakShape
                      object
                   })




methods::setMethod("deleteFilePtrSet", signature = "ptrSet", 
                   function(object,deletedFiles) {
   
    object@mzCalibRef[deletedFiles] <- NULL
    object@signalCalibRef[deletedFiles] <- NULL
    object@errorCalibPpm[deletedFiles] <- NULL
    object@coefCalib[deletedFiles] <- NULL
    object@indexTimeCalib[deletedFiles] <- NULL
    object@knots[deletedFiles] <- NULL
    object@fctFit[deletedFiles] <- NULL
    object@peakShape[deletedFiles] <- NULL
    object@primaryIon[deletedFiles] <- NULL
    object@resolution[deletedFiles] <- NULL
    object@TIC[deletedFiles] <- NULL
    object@timeLimit[deletedFiles] <- NULL
    object@peakList[deletedFiles] <- NULL
    object@breathTracer[deletedFiles] <- NULL
    object@ptrTransmisison[deletedFiles] <- NULL
    object@prtReaction[deletedFiles] <- NULL
    object@date[deletedFiles] <- NULL

    object
})



methods::setMethod("mergedPtrSet", signature = "ptrSet", 
                   function(object,ptrSetNewfile,orderedFile) {

    object@mzCalibRef <- c(getCalibrationInfo(object)$mzCalibRef, ptrSetNewfile$mzCalibRefList)[basename(orderedFile)]
    object@signalCalibRef <- c(getCalibrationInfo(object)$signalCalibRef, ptrSetNewfile$signalCalibRef)[basename(orderedFile)]
    object@errorCalibPpm <- c(getCalibrationInfo(object)$errorCalibPpm, ptrSetNewfile$errorCalibPpm)[basename(orderedFile)]
    object@coefCalib <- c(getCalibrationInfo(object)$coefCalib, ptrSetNewfile$coefCalibList)[basename(orderedFile)]
    object@resolution <- c(getPTRInfo(object)$resolution, ptrSetNewfile$resolution)[basename(orderedFile)]
    object@TIC <- c(getTimeInfo(object)$TIC, ptrSetNewfile$TIC)[basename(orderedFile)]
    object@timeLimit <- c(getTimeInfo(object)$timeLimit, ptrSetNewfile$timeLimit)[basename(orderedFile)]
    object@primaryIon <- c(getPTRInfo(object)$primaryIon, ptrSetNewfile$primaryIon)[basename(orderedFile)]
    object@breathTracer <- c(getTimeInfo(object)$breathTracer, ptrSetNewfile$breathTracer)[basename(orderedFile)]
    object@ptrTransmisison <- c(getPTRInfo(object)$ptrTransmisison, ptrSetNewfile$ptrTransmisison)[basename(orderedFile)]
    object@prtReaction <- c(getPTRInfo(object)$prtReaction, ptrSetNewfile$prtReaction)[basename(orderedFile)]
    object@date <- c(getDate(object), ptrSetNewfile$date)[basename(orderedFile)]
    object@peakShape <- c(object@peakShape, ptrSetNewfile$peakShape)[basename(orderedFile)]
    object@knots<-c(object@knots, ptrSetNewfile$knots)[basename(orderedFile)]
    object@indexTimeCalib<-c(object@indexTimeCalib, ptrSetNewfile$indexTimeCalib)[basename(orderedFile)]
    object@fctFit<-c(object@fctFit, ptrSetNewfile$fctFit)[basename(orderedFile)]
    
    
    return(object)
    })
