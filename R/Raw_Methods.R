## Calibration ----

#' Calibrate the mass axis with references masses
#' 
#' To convert Time Of Flight (TOF) axis to mass axis, different formula are proposes in the 
#' literature (average and al. 2013, Cappelin and al. 2010) mz = ((tof-b)/a )^2 and 
#' mz = a + b (tof) +c  (tof)^2. To estimate those parameters, references peaks with accurate know 
#' masses and without overlapping peak are needed. The best is that the references masses 
#' covers a maximum of the mass range.
#' @param x a prtRaw or ptrSet object
#' @param mzCalibRef Vector of accurate mass values of intensive peaks and 'unique' in a 
#' nominal mass interval (without overlapping)
#' @param tol the maximum error tolerated in ppm. If more than \code{tol} warnings. 
#' @return the same ptrRaw or ptrSet as in input, with the following modified element:
#' \itemize{
#' \item mz: the new mz axiscalibrated
#' \item ramM: same raw matrix with the new mz axis in rownames
#' \item calibMassRef: reference masses used for the calibration
#' \item calibMzToTof and calibTofToMz: function to convert TOF to mz
#' \item calibError: the calibration error to the reference masses in ppm
#' } 
#' @examples 
#' library(ptairData)
#' filePath <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
#' raw <- readRaw(filePath, calibTIS = FALSE)
#' rawCalibrated <- calibration(raw)
#' @rdname calibration
#' @export 
methods::setMethod(f = "calibration",
          signature = "ptrRaw", 
          function(x, mzCalibRef = c(21.022, 29.013424,41.03858,59.049141,75.04406, 
                                              203.943, 330.8495), tol=70){
            
            object <- x
            
            # get mz axis and average spectrum
            mz <- object@mz
            sp <- rowSums(object@rawM)/dim(object@rawM)[2]

            #performs calibration
            calib<-calibrationFun(sp,mz,mzCalibRef,object@calibMzToTof,tol)
           
            # update object object
            object@mz <- calib$mzVnbis
            rownames(object@rawM) <- calib$mzVnbis
            object@calibMassRef <- calib$mzCalibRef
            object@calibToftoMz <- calib$calib_formula
            object@calibMzToTof <- calib$calib_invformula
            object@calibSpectr <- calib$calibSpectr
            object@calibError <- calib$error
            object@calibCoef <- calib$coefs

            return(object)
          } )

#'calibration function
#'
#'Performs calibration on sp with mzCalibRef reference masses and mzToTofFunc as previous 
#'calibration function
#' @param sp spectrum
#' @param mz mass axis 
#' @param mzCalibRef masses of know reference peaks
#' @param mzToTofFunc function to convert mz to tof (previous calibration)
#' @param tol maximum error tolarated in ppm
#' @return list 
calibrationFun<-function(sp,mz,mzCalibRef,mzToTofFunc,tol){
  width.window<-0.4
  
  # check if mzCalibRef are in mz
  outMz <- which(vapply(mzCalibRef, 
                        function(x) !any( round(x) -width.window < mz & 
                                            mz < round(x) + width.window),
                        FUN.VALUE = TRUE ))
  if(length(outMz)!=0) {
    message(paste( paste(mzCalibRef[outMz],collapse = " "),
                   "excluded, not contains in the mass axis \n"))
    mzCalibRef <- mzCalibRef[-outMz]
  }
  
  # calculate tof axis
  if(is.null(mzToTofFunc(1))) tof<-seq(0,length(mz)-1) else tof <- mzToTofFunc(mz)
  
  #spectrum of mass calib 
  spTronc <- lapply(mzCalibRef,function(m) {
    index<- which((mz < m + width.window) & (mz > m - width.window))
    tof<-tof[index]
    mz <- mz[index]
    signal <- sp[index] 
    return(list(signal=signal,mz=mz,tof=tof))})
  
  # test if there is a only one peak
  nLocalMax <-vapply(spTronc, function(x) {
    length( LocalMaximaSG( sp = x$signal,
                           minPeakHeight = 0.1*max(x$signal)) )
    
  }, FUN.VALUE = 0 )
  
  
  badMass <- which(nLocalMax !=1 ) 
  if(length(badMass)!=0){
    
    mzCalibRef <- mzCalibRef[-badMass]
    spTronc <- spTronc[-badMass]
  }
  
  if( length(mzCalibRef) < 2 ) stop("To few references masses for calibration")
  
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
    
    
    #interpol<- stats::spline( tofrange, sp, n=1000*length(tofrange) )
    #tofMax<-interpol$x[ which.max(interpol$y)]
    return(tofMax)
  },FUN.VALUE = 0.1)
  
  # re estimated calibration coefficient with reference masses
  regression <- stats::nls( rep(1,length(mzCalibRef))  ~  I( ( (tofMax - b) / a) ^ 2 /mzCalibRef ),
                     start = list(a = 8838, b= -219 ), algorithm = "port")
  coefs <- stats::coefficients(regression)
  calib_formula <- function(tof) ((tof - coefs['b']) / coefs['a']) ^ 2
  calib_invformula <- function(m) sqrt(m)*coefs['a'] + coefs['b']
  
  # the new mass axis calibrated
  mzVnbis <- calib_formula(tof) 
  
  # the new position of reference masses
  mzRefRaw <- calib_formula(tofMax)
  
  # error in ppm
  error <- abs(mzRefRaw-mzCalibRef)*10^6/mzCalibRef
  names(error)<-mzCalibRef
  
  # spectrum of calibration masses with new mass axis
  calibSpectr <- lapply(mzCalibRef,function(m) {
    index<- which((mzVnbis < m + width.window) & (mzVnbis > m - width.window))
    mz <- mzVnbis[index]
    signal <- sp[index] 
    return(list(signal=signal,mz=mz))})
  
  if(any(abs(error) > tol)) message(paste("error greater than ",tol,"\n"))
  
  return(list(mzVnbis =mzVnbis,
              mzCalibRef = mzCalibRef,
              calib_formula = calib_formula,
              calib_invformula= calib_invformula,
              calibSpectr= calibSpectr,
              error= error,
              coefs= as.matrix(coefs)))
  
}
## plotRaw ----
#' @rdname plotRaw
#' @export
methods::setMethod(f = "plotRaw",
          signature = "ptrRaw",
          function(object,
                   mzRange ,
                   timeRange = c(NA, NA),
                   type = c("classical", "plotly")[1],
                   ppm = 2000,
                   palette = c("heat",
                               "revHeat",
                               "grey",
                               "revGrey",
                               "ramp")[1],
                   showVocDB = TRUE,
                   figure.pdf = "interactive",...) 
          {
            
            
            mzVn <- object@mz
            timeVn <- object@time
            rawMN <- object@rawM
            
            
            if (length(mzRange) == 1)
              mzRange <- mzRange + c(-1, 1) * mzRange * ppm * 1e-6
            
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
              
              vocdb_sel.vl <- vocdbDF[, "mass_Hplus"] >= mzRange[1] &
                vocdbDF[, "mass_Hplus"] <= mzRange[2]
              
              if (sum(vocdb_sel.vl)) {
                vocdbDF <- vocdbDF[vocdb_sel.vl, , drop = FALSE]
              } else
                vocdbDF <- NULL
              
            } else
              vocdbDF <- NULL
            
            
            if (figure.pdf != "interactive") {
              if (type == "plotly")
                stop("'plotly display is only available in the 'interactive' mode currently.",
                     call. = FALSE)
              filenameSplitVc <- unlist(strsplit(basename(figure.pdf), ".", fixed = TRUE))
              extC <- utils::tail(filenameSplitVc, 1)
              if (extC == "pdf") {
                grDevices::pdf(figure.pdf)
              } else
                stop("The extension of the 'figure.pdf' filename argument should be 'pdf'",
                     call. = FALSE)
            }
            
            switch(type,
                   
                   classical = {
                     
                     imageMN <- t(rawSubMN)[, seq_len(nrow(rawSubMN)), drop = FALSE]
                     rownames(imageMN) <- round(as.numeric(rownames(imageMN)))
                     colnames(imageMN) <- round(as.numeric(colnames(imageMN)), 4)
                     
                     paletteVc <- .palette(palette = palette)
                     
                     currentParLs <- graphics::par()
                     for (parC in c("cin", "cra", "csi", "cxy", "din", "page"))
                       currentParLs[[parC]] <- NULL   
                     
                     marLs <- list(chr = c(0.6, 4.1, 1.1, 0),
                                   sca = c(0.6, 0.6, 1.1, 7.1),
                                   ima = c(4.1, 4.1, 0, 0),
                                   spe = c(4.1, 0.6, 0, 1.1))
                     
                     graphics::par(font = 2,
                         font.axis = 2,
                         font.lab = 2,
                         pch = 18)
                     
                     graphics::layout(matrix(c(1, 2,
                                     3, 4),
                                   byrow = TRUE,
                                   nrow = 2),
                            heights = c(2, 5),
                            widths = c(5, 2))
                     
                     ## chr: Chromatogram
                     
                     graphics::par(mar = marLs[["chr"]])
                     
                     chromVn <- apply(rawSubMN, 2,
                                      function(intVn)
                                        mean(intVn, na.rm = TRUE))
                     plot(as.numeric(names(chromVn)),
                          chromVn,
                          cex = 0.7,
                          pch = 16,
                          xlab = "",
                          ylab = "",
                          xaxs = "i",
                          xaxt = "n",
                          yaxs = "i")
                     
                     graphics::mtext("Mean of intensity",
                           cex = 0.8,
                           side = 2,
                           line = 2.5)
                     
                     ## sca: Color scale
                     
                     graphics::par(mar = marLs[["sca"]])
                     
                     .drawScale(imageMN = imageMN,
                                paletteVc = paletteVc)
                     
                     ## ima: Image
                     
                     graphics::par(mar = marLs[["ima"]])
                     
                     .drawImage(imageMN = imageMN, paletteVc = paletteVc)
                     
                     if (showVocDB && !is.null(vocdbDF)) {
                       mzImaVn <- as.numeric(colnames(imageMN))
                       graphics::abline(h = vapply(vocdbDF[, "mass_Hplus"],
                                         function(mzN)
                                           (mzN - min(mzImaVn))/diff(range(mzImaVn)) * ncol(imageMN) + par("usr")[1],
                                         FUN.VALUE =1.1 ),
                              lty = "dotted")
                     }
                     
                     ## spe: Spectrum
                     
                     graphics::par(mar = marLs[["spe"]])
                     
                     specVn <- apply(rawSubMN, 1,
                                     function(intVn)
                                       mean(intVn, na.rm = TRUE))
                     
                     plot(specVn,
                          as.numeric(names(specVn)),
                          cex = 0.7,
                          pch = 16,
                          xlab = "",
                          ylab = "",
                          xaxs = "i",
                          yaxs = "i",
                          yaxt = "n")
                     
                     graphics::mtext("Mean of intensity",
                           cex = 0.8,
                           side = 1,
                           line = 2.5)
                     
                     if (showVocDB && !is.null(vocdbDF)) {
                       
                       graphics::abline(h = vocdbDF[, "mass_Hplus"], lty = "dotted")
                       
                     }
                     
                     graphics::par(currentParLs)
                     
                     if (figure.pdf != "interactive")
                       dev.off()
                     
                   },
                   
                   plotly = {
                     
                     plotlyBuild <- function() {
                       
                       p <- plotly::subplot(
                         
                         plotly::plot_ly(x = as.numeric(colnames(rawSubMN)),
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
                                         type = 'scatter', mode = 'lines'),
                         
                         nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), 
                         shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
                         
                       )
                       
                       return(p)
                       
                     }
                     
                     p <- suppressWarnings(suppressMessages(plotlyBuild()))
                     
                     p <- plotly::layout(p, showlegend = FALSE)
                     
                     print(plotly::ggplotly(p))
                     
                   })
            
            if (showVocDB & !is.null(vocdbDF)) {
              vocdbDF <- vocdbDF[nrow(vocdbDF):1, , drop = FALSE]
              print(vocdbDF[, c("mass_Hplus", "formula_Hplus", "name_iupac"),
                            drop = FALSE])
            }
            
            return(invisible(list(rawsubM = rawSubMN,
                                  vocsubDB = vocdbDF)))
            
          })


.palette <- function(palette) {
  
  switch(palette,
         heat = {return(rev(grDevices::rainbow(ceiling(256 * 1.5))[seq(1,256)]))},
         revHeat = {return(grDevices::rainbow(ceiling(256 * 1.5))[seq(1,256)])},
         grey = {return(grDevices::grey((0:255) / 256))},
         revGrey = {return(rev(grDevices::grey((0:255) / 256)))},
         ramp = {return(grDevices::colorRampPalette(c("blue", "orange", "red"),
                                                    space = "rgb")(256)[seq(1,256)])})
  
}

.drawScale <- function(imageMN,
                       paletteVc) {
  
  
  ylimVn <- c(0, 256)
  ybottomVn <- 0:255
  ytopVn <- seq_len(256)
  
  plot(x = 0,
       y = 0,
       font.axis = 2,
       font.lab = 2,
       type = "n",
       xlim = c(0, 1),
       ylim = ylimVn,
       xlab = "",
       ylab = "",
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n")
  
  graphics::rect(xleft = 0,
       ybottom = ybottomVn,
       xright = 1,
       ytop = ytopVn,
       col = paletteVc,
       border = NA)
  
  graphics::axis(at = .prettyAxis(range(imageMN, na.rm = TRUE), 256)$atVn,
       font = 2,
       font.axis = 2,
       labels = .prettyAxis(range(imageMN, na.rm = TRUE), 256)$labelVn,
       las = 1,
       lwd = 2,
       lwd.ticks = 2,
       side = 4,
       xpd = TRUE)
  
  graphics::arrows(graphics::par("usr")[2],
                   graphics::par("usr")[4],
                   graphics::par("usr")[2],
                   graphics::par("usr")[3],
                   code = 0,
                   lwd = 2,
                   xpd = TRUE)
  
  graphics::arrows(graphics::par("usr")[1],
                   graphics::par("usr")[4],
                   graphics::par("usr")[1],
                   graphics::par("usr")[3],
                   code = 0,
                   lwd = 2,
                   xpd = TRUE)
  
  graphics::box(lwd = 2)
  
  
}


.prettyAxis <- function(axisValuesVn,
                        opLengthN) {
  
  if (NA %in% axisValuesVn) {
    
    warning("NA in axisValuesVn")
    
    axisValuesVn <- as.vector(stats::na.omit(axisValuesVn))
    
  }
  
  if (opLengthN < length(axisValuesVn))
    stop("The length of in vector must be inferior to the length of the length parameter.")
  
  if (length(axisValuesVn) < opLengthN) {
    
    axisValuesVn <- seq(from = min(axisValuesVn), to = max(axisValuesVn), length.out = opLengthN)
    
  }
  
  prettyAxisValues <- pretty(axisValuesVn)
  
  prettyLabelsVn <- prettyAtVn <- c()
  
  for (n in seq_along(prettyAxisValues))
    if (min(axisValuesVn) < prettyAxisValues[n] && prettyAxisValues[n] < max(axisValuesVn)) {
      prettyLabelsVn <- c(prettyLabelsVn, prettyAxisValues[n])
      prettyAtVn <- c(prettyAtVn, which(abs(axisValuesVn - prettyAxisValues[n]) == min(abs(axisValuesVn - prettyAxisValues[n])))[1])
    }
  
  prettyAxisLs <- list(atVn = prettyAtVn,
                       labelVn = prettyLabelsVn)
  
  
  return(prettyAxisLs)
  
}


.drawImage <- function(imageMN,
                       paletteVc) {
  
  graphics::image(x = seq_len(nrow(imageMN)),
                  y = seq_len(ncol(imageMN)),
                  z = imageMN,
                  col = paletteVc,
                  font.axis = 2,
                  font.lab = 2,
                  xaxt = "n",
                  yaxt = "n",
                  xlab = "",
                  ylab = "")
  
  timeVn <- as.numeric(rownames(imageMN))
  
  prettyTimeVn <- pretty(timeVn)
  
  prettyTimeVn <- prettyTimeVn[min(timeVn) <= prettyTimeVn &
                                 prettyTimeVn <= max(timeVn)]
  
  prettyTimeVi <- vapply(seq_along(prettyTimeVn), function(k) {
    which(abs(timeVn - prettyTimeVn[k]) == min(abs(timeVn - prettyTimeVn[k])))[1]
  },FUN.VALUE = 1)
  
  graphics::axis(side = 1,
       at = prettyTimeVi,
       font = 2,
       labels = prettyTimeVn)
  
  mzVn <- as.numeric(colnames(imageMN))
  
  prettyMzVn <- pretty(mzVn)
  
  prettyMzVn <- prettyMzVn[min(mzVn) <= prettyMzVn &
                             prettyMzVn <= max(mzVn)]
  
  prettyMzVi <- vapply(seq_along(prettyMzVn), function(k) {
    which(abs(mzVn - prettyMzVn[k]) == min(abs(mzVn - prettyMzVn[k])))[1]
  },FUN.VALUE = 1)
  
  graphics::axis(side = 2,
       at = prettyMzVi,
       font = 2,
       labels = prettyMzVn)
  
  ## xlab
  
  graphics::mtext(line = 2.5,
        side = 2,
        text = "m/z",
        cex = 0.8)
  
  ## ylab
  
  graphics::mtext(line = 2.5,
        side = 1,
        text = "time (s)",
        cex = 0.8)
  
  ## border
  
  graphics::box(lwd = 2)
  
}


## plotCalib ----
#' @rdname plotCalib
#' @export
methods::setMethod(f="plotCalib",
          signature = "ptrRaw",
          function(object,ppm=2000,...){
  
        raw <- object
        # get mass and specter
        mzCalibRef<- raw@calibMassRef
        calibSpectr <- raw@calibSpectr
        error <- raw@calibError
  
        # plot in a window of 2000 ppm
        nb_plot<-length(mzCalibRef)
        nb_row <- ceiling(sqrt(nb_plot))
        graphics::par(oma = c(0, 0, 3, 0))
        graphics::layout(matrix(seq(1,(nb_row^2)),nrow=nb_row,byrow = TRUE))
        
        #loop over masses
        for (i in seq_along(mzCalibRef)){
          m <- mzCalibRef[i]
          th<- m*(ppm/2)/10^6
          mz <- calibSpectr[[i]]$mz
          index<- which(m-th < mz & mz < m+ th)
          x <- mz[index]
          y <- calibSpectr[[i]]$signal[index]
          plot( x , y , type="l",
          lwd=2, ylab="intenisty", xlab="mz", 
          main  = c(m,paste("error:", round(error[i],2),"ppm")))
          graphics::abline(v=m, col="red", lwd=2)
      }
        title(main=raw@name,outer = TRUE,line =0.5,cex.main=2)
        graphics::par(oma = c(0, 0, 0, 0))
        graphics::layout(matrix(1))
} )


##plotTIC----
#' @param fracMaxTIC Percentage (between 0 and 1) of the maximum of the Total Ion Chromatogram (TIC) 
#' amplitude with baseline removal. We will analyze only the part of the spectrum where 
#' the TIC intensity is higher than `fracMaxTIC * max(TIC) `. If you want to analyze the entire spectrum, 
#' set this parameter to 0. 
#' @rdname plotTIC
#' @export
#' @examples 
#' library(ptairData)
#' filePath <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
#' raw <- readRaw(filePath)
#' p <- plotTIC(raw)
#' p
#' 
methods::setMethod(f="plotTIC",
          signature = "ptrRaw",
          function(object, type, baselineRm, showLimits,fracMaxTIC=0.5,...){
            
            
            #get the TIC and time limit
            TIC<-colSums(object@rawM)
            
            
            #remove baseline
            if(baselineRm) {
              bl <- try(snipBase(TIC))
              if(is.null(attr(bl,"condition"))) TIC <-TIC - bl else TIC<- TIC - TIC[1]
  
            }
            
            
      
            plot<- ggplot2::qplot(x=object@time,y=TIC,
                                  xlab="time",ylab="intensity",main=paste("TIC of",object@name)) 
            if(showLimits){
              #calculate timeLimit 
              indLim <- timeLimits(object, fracMaxTIC = fracMaxTIC, plotDel = FALSE)
              plot<- plot +
                geom_vline(aes(xintercept = object@time[c(indLim)],
                               color="time limits")) + scale_fill_manual("Legend")
            }
             plot <- plot + ggplot2::theme( 
               plot.title = ggplot2::element_text(size=20, face="bold"),
               axis.title = ggplot2::element_text(size=16),
               axis.text = ggplot2::element_text(size=14),
               legend.text =  ggplot2::element_text(size=14),
               legend.title = ggplot2::element_text(size=16))
                
            switch (type,
                    ggplot = return(plot),
                    plotly = return(plotly::ggplotly(plot))
            )
          }#end function
)


## timeLimit ----
#' Calculate time limits on the Chromatogram
#' 
#' This function derives limits on the Total Ion Chromatogram TIC, where the intenisty is greater than \code{fracMaxTIC*max(TIC)}, 
#' where max(TIC)  is the maximum of teh TIC with baseline removal.
#' In this way,  the expiration limits, or headsapce analysis limits can be detected. So, by setting 
#' \code{fracMaxTIC} close to 1, the size of teh limits will be restricted.
#' 
#' @param object a ptrRaw or ptrSet object
#' @param fracMaxTIC between 0 and 1. Percentage of the maximum of the Chromatogram amplitude with baseline removal. 
#' If you want a finer limitation, increase \code{fracMaxTIC}, indeed decrease
#' @param traceMasses NULL or a integer. Correspond to a nominal masses of Extract Ion Chromatogram (EIC)
#'  whose limits you want to compute. If NULL, the limits are calculated on the Total Ion Chromatogram (TIC).
#' @param minPoints minimum duration of an expiration (in index).
#' @param plotDel boolean. If TRUE, the Chormatogram is ploted with limits and threshold.
#' @return a matrix of index, where each colomn correspond to one expriration, the first row 
#' it is the beginning and the seconde the end, or NA if no limits are detected.
#' @rdname timeLimits
#' @examples
#' library(ptairData)
#' filePath <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
#' raw <- readRaw(filePath)
#' 
#' ind_lim <- timeLimits(raw, fracMaxTIC=0.9, plotDel=TRUE)
#' ind_lim_acetone <- timeLimits(raw, fracMaxTIC=0.5, traceMasses = 59,plotDel=TRUE)
#'@export
methods::setMethod(f="timeLimits",
          signature = "ptrRaw",
          function(object,fracMaxTIC=0.5, traceMasses= NULL, minPoints = 2, plotDel=FALSE){
            
            rawM <-object@rawM
            mz <- object@mz
            
            if(is.null(dim(rawM))) stop("rawM must be a matrix")
            if(fracMaxTIC <0 || fracMaxTIC > 1) stop("fracMaxTIC must be between 0 and 1")
            if(fracMaxTIC==0) {
              return(matrix(c(1,ncol(rawM)),ncol=1,nrow=2,dimnames = list(c("start","end"))))
                     }
            if(is.null(traceMasses)) { TIC <- colSums(rawM)
            } else {
              index<- lapply(traceMasses, function(x) {
                which( x - 0.4 < mz & mz < x + 0.4)
              })
              TIC <- colSums(rawM[unlist(index),]) 
            }
            
            indLim<-timeLimitFun(TIC,fracMaxTIC, traceMasses, minPoints, plotDel)
            
            return(indLim)
          }
          )

timeLimitFun<-function(TIC,fracMaxTIC=0.5, traceMasses= NULL, minPoints = 3, plotDel=FALSE){
  
  ## baseline corretion
  bl <- try(baselineEstimation(TIC,d=1))
  if(is.null(attr(bl,"condition"))) TIC.blrm<-TIC - bl else TIC.blrm<-TIC
  threshold<-(max(TIC.blrm)-min(TIC.blrm))*fracMaxTIC
  
  ## delimitation
  hat <- which(TIC.blrm > (threshold) )
  if(length(hat)==0) {
    message("no limits detected")
    plot(TIC,type='l',xlab="Time (s)",ylab="intensity",cex.lab=1.5, main = "Time limit") 
    return(NA)
  }
  
  hat_end <- c(hat[which(diff(hat) !=1)],utils::tail(hat,1))
  hat_begin<-c(hat[1],hat[which(diff(hat) !=1)+1])
  hat_lim<-unname(rbind(hat_begin,hat_end))
  row.names(hat_lim)<-c("start","end")
  
  hat_lim<-hat_lim[, hat_lim["end",]- hat_lim["start",] >= minPoints,drop=FALSE]
  
  if(plotDel){
    if(is.null(traceMasses)) subtitle<-"TIC" 
    else subtitle <- paste("EIC mz:",paste(round(traceMasses,2),collapse = "-"))
    plot(TIC.blrm,type='l',xlab="Time (s)",ylab="intensity",cex.lab=1.5, main = paste("Time limit",subtitle), 
         ylim=c(min(TIC.blrm)-0.2*(max(TIC.blrm)-min(TIC.blrm)),max(TIC.blrm)),lwd=2)
    graphics::abline(v=c(hat_lim),col="red",lty=2,lwd=2)
    graphics::abline(h=threshold,lty=2)
    graphics::legend("bottomleft",legend=c("threshold", "limit"),col=c("black","red"),lty=c(2,2),horiz = TRUE)
  }
  
  return(hat_lim)
}

bakgroundDetect<-function(TIC,derivThreshold=0.01,  minPoints = 4, plotDel=FALSE){
  
  bl <- try( baselineEstimation(TIC,d=1))
  if(is.null(attr(bl,"condition"))) TIC.blrm<-TIC - bl else TIC.blrm<-TIC
  threshold<-(max(TIC.blrm)-min(TIC.blrm))*0.5
  
  ## delimitation
  ind.Exp <- which(TIC.blrm > (threshold) )

  dTIC <- diff(TIC.blrm)/max(TIC.blrm)
  ind.Bg<-which(abs(dTIC)<derivThreshold)
  ind.Bg<-ind.Bg[ -which(ind.Bg %in%ind.Exp) ]
  
  ind.Bg_end <- c(ind.Bg[which(diff(ind.Bg) !=1)],utils::tail(ind.Bg,1))
  ind.Bg_begin<-c(ind.Bg[1],ind.Bg[which(diff(ind.Bg) !=1)+1])
  limBg<-unname(rbind(ind.Bg_begin,ind.Bg_end))
  row.names(limBg)<-c("start","end")
  
  limBg<-limBg[, limBg["end",]- limBg["start",] >= minPoints,drop=FALSE]
  
  if(ncol(limBg)==0) warning("no background detected")
  
    if(plotDel){
    subtitle<-"TIC" 
    plot(TIC,type='l',xlab="Time (s)",ylab="intensity",cex.lab=1.5, main = paste("Background limit",subtitle), 
         ylim=c(min(TIC)-0.2*(max(TIC)-min(TIC)),max(TIC)),lwd=2)
    graphics::abline(v=c(limBg),col="red",lty=2,lwd=2)
    graphics::legend("bottomleft",legend=c( "limit"),col=c("red"),lty=c(2),horiz = TRUE)
  }
  return(limBg)
}

## PeakList ----
#' Detection and quantification of peaks on spectrum. 
#'
#' @param raw ptrRaw object 
#' @param mzNominal the vector of nominal mass where peaks will be detected 
#' @param ppm the minimum distance between two peeks in ppm 
#' @param minIntensity the minimum intenisty for peaks detection. The final threshold for peak detection
#' will be : max ( \code{minPeakDetect} , thresholdNoise ). The thresholdNoise correspond to
#'  max(\code{thNoiseRate} * max( noise around the nominal mass), \code{thIntensityRate} * 
#'  max( intenisty in the nominal mass). The noise around the nominal mass correspond : 
#'  \code{[m-windowSize-0.2,m-windowSize]U[m+windowSize,m+WindowSize+0.2]}.
#' @param fctFit the function for the quantification of Peak, should be average or Sech2
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
#' @return a list containning: \itemize{
#' \item peak: a data.frame, with for all peak detected: the mass center, the intensity count, the peak width (delta_mz),
#' correspond to the Full Width Half Maximum (FWHM),the resolution m/delta_m, the other 
#' parameters values estimated of \code{fitFunc}.
#' \item warnings: warnings generated by the peak detection algorithm per nominal masses
#' \item infoPlot: elements nedded to plot the fitted peak per nominal masses
#' }
#' @examples 
#' library(ptairData)
#' filePath <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
#' file <- readRaw(filePath)
#' 
#' peakList <- PeakList(file, mzNominal = c(21,63))
#' peakList$peak
#'@rdname PeakList 
#'@export
methods::setMethod(f="PeakList",
          signature = "ptrRaw",
          function(raw,
                   mzNominal = unique(round(raw@mz)), ppm = 130, 
                   minIntensity=5, fctFit=c("Sech2","average")[1], maxIter=2, autocorNoiseMax = 0.3,
                   plotFinal=FALSE, plotAll=FALSE, thNoiseRate=1.1, thIntensityRate = 0.01,
                   countFacFWHM=10, daSeparation=0.005, d=3, windowSize=0.4) {
  
  #get raw element 
  sp <- rowSums(raw@rawM)/(ncol(raw@rawM)*round(raw@time[2]-raw@time[1])) # average spectrum 
  mz <- raw@mz # mass axis
  mzCalibRef <- raw@calibMassRef 
  tofToMz <- raw@calibToftoMz 
  mzToTof <- raw@calibMzToTof
  
  if(fctFit=="average") l.shape<-determinePeakShape(sp,mz,massRef = mzCalibRef)
  
  prePeaklist <- lapply(mzNominal, function(m) peakListNominalMass(m,mz,sp,ppm, mzToTof,tofToMz,
                                                                    minIntensity, fctFit, maxIter, autocorNoiseMax ,
                                                                    plotFinal, plotAll, thNoiseRate, thIntensityRate ,
                                                                    countFacFWHM, daSeparation, d, windowSize) )
  
  peaklist<-do.call(rbind,lapply(prePeaklist, function(x) x[[1]]))
  warning <-do.call(rbind,lapply(prePeaklist, function(x) x[[2]]))
  infoPlot<- do.call(c,lapply(prePeaklist, function(x) x[[3]]))
  
  return(list(peak=peaklist, warning=warning,infoPlot=infoPlot))
})

##TODO peakLIst method dor spectrum array

##detectpeak----
#' @rdname detectPeak
#' @examples 
#' 
#' library(ptairData)
#' filePath <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
#' raw <- readRaw(filePath,mzCalibRef=c(21.022,59.049))
#' peakList <- detectPeak(raw, mzNominal = c(21,63))
#' peakList$aligned
#' @export
methods::setMethod(f="detectPeak",
          signature = "ptrRaw",
          function(x, 
                   mzNominal=NULL , ppm=130, ppmGroupBkg=50, fracGroup=0.8, minIntensity=10, 
                   fctFit=c("Sech2","average")[1],normalize=TRUE,fracMaxTIC=0.5,processFun=processFileSepExp,...)
          {
            raw<-x
            #get infomration
            massCalib<-raw@calibMassRef
            
            sp<-rowSums(raw@rawM)/ncol(raw@rawM)
            mz<-raw@mz
            fit<-peakListNominalMass(21,mz,sp,ppmPeakMinSep=500, mzToTof= raw@calibMzToTof,
                                     minPeakDetect=10, fitFunc="Sech2", maxIter=1, autocorNoiseMax=0.3 ,
                                     plotFinal=F, plotAll=F, thNoiseRate=1.1, thIntensityRate=0.01 ,
                                     countFacFWHM=10, daSeparation=0.1, d=3, windowSize=0.2 )
            primaryIon<-max(fit$peak[,"quanti"]) #max in case of false positif detected
            
            indTimeLim<- timeLimits(raw, fracMaxTIC = fracMaxTIC)
            
            peakLists<-processFun(raw,massCalib,primaryIon,indTimeLim, mzNominal,
                                                  ppm, ppmGroupBkg, fracGroup,minIntensity, 
                                                  fctFit,normalize)

            return(peakLists)
          } )


estimateResol<-function(calibMassRef,calibSpectr){
              m <- calibMassRef
              
              delta<-vapply(calibSpectr,function(x) {
                mzM <- x$mz
                spM <- x$signal
                #half maximum 
                hm <- max(spM)/2
                #limits 
                lim1 <- findEqualGreaterM(spM[seq_len(which.max(spM))],hm)
                lim2 <- unname(which.max(spM))+ FindEqualLess(spM[(which.max(spM)+1):length(spM)],hm)
                # equation : (intrepolation linéaire entre lim et (lim-1)) = hm
                deltaBorne<- vapply( c(lim1,lim2) ,
                                      function(x) (hm*(mzM[x]-mzM[x-1])-(mzM[x]*spM[x-1]-mzM[x-1]*spM[x]))/
                                           (spM[x]-spM[x-1]),FUN.VALUE = 1.1 )
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
methods::setMethod("show","ptrRaw",
          function(object){
            cat(object@name,"\n")
            cat("  mz range: ",paste(round(range(object@mz),2),collapse = " - "),"\n")
            cat("  time range: ",paste(round(range(object@time),2),collapse = " - "),"\n")
            cat("  Calibration error in ppm: \n")
            if(length(object@calibError)==1) {
              cat("     No calibration performs") 
              } else {
                print(round(object@calibError,2))
              }
            
          })

#### dead time correction -----
#'Dead time correction on raw data
#'@param raw ptrRaw object
#'@param ve exenting dead time
#'@param vne non extending dead time
#'@param r number of extraction
#'@param threshold only bin of intenisty more then threshold*r which be corrected
#'@return a ptrRaw object with the raw matrix corrected
deadTimeCorr<-function(raw,ve,vne,r,threshold=0.1){
  rawM <- raw@rawM  
  index<-which(rawM > threshold*r)
  for(j in index){
    rawM[j]<- -r*log( 1-( rawM[j]/r * 
                            (1-sum(raw[(j-vne):(j-1-ve)])/r)^-1*
                            exp(sum(rawM[(j-ve):(j-1)])/r)
                          )
                      )
  }
  raw@rawM<-rawM
  return(raw)
}
