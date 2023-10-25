## ----set, include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(fig.width=12, fig.height=8) 

## ---- eval=FALSE--------------------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("ptairMS")

## ----package, message = FALSE, warning = FALSE--------------------------------
library(ptairMS)
library(ptairData)

## ---- eval=FALSE--------------------------------------------------------------
#  ptairMS::RunShinnyApp()

## ----ptairData----------------------------------------------------------------
dirRaw <- system.file("extdata/exhaledAir", package = "ptairData")

## ----createPtrSet-------------------------------------------------------------
exhaledPtrset <- createPtrSet(dir=dirRaw,
                     setName="exhaledPtrset",
                     mzCalibRef = c(21.022, 60.0525),fracMaxTIC = 0.7,
                     saveDir = NULL )

## -----------------------------------------------------------------------------
exhaledPtrset

## ----getFileNames-------------------------------------------------------------
getFileNames(exhaledPtrset)

## ----plot---------------------------------------------------------------------
plot(exhaledPtrset)

## ----calib table--------------------------------------------------------------
calib_table<-read.csv(system.file("extdata", "reference_tables/calib_table.tsv", package = "ptairMS"),sep="\t")
knitr::kable(calib_table)

## ----calibration--------------------------------------------------------------
plotCalib(exhaledPtrset,fileNames=getFileNames(exhaledPtrset)[1])

## ----plotCalib----------------------------------------------------------------
exhaledPtrset <- calibration(exhaledPtrset, mzCalibRef =  c(21.022, 60.0525,75.04406))
plot(exhaledPtrset,type="calibError")

## ----shinny app, eval=FALSE---------------------------------------------------
#  exhaledPtrset <- changeTimeLimits(exhaledPtrset)

## ----timeLimits1--------------------------------------------------------------
samplePath <-getFileNames(exhaledPtrset,fullNames = TRUE)[1]
sampleRaw <- readRaw(samplePath, calib = FALSE)
expirationLimit <- timeLimits(sampleRaw,fracMaxTIC =  0.5,plotDel = TRUE, mzBreathTracer = 60.05)
expirationLimit <- timeLimits(sampleRaw,fracMaxTIC =  0.9,plotDel = TRUE,mzBreathTracer = NULL)

## ----plotTIC1-----------------------------------------------------------------
plotTIC(object = exhaledPtrset,baselineRm = TRUE,type = "ggplot")

## ----getSampleMetadata--------------------------------------------------------
getSampleMetadata(exhaledPtrset)

## ----setSampleMetadata--------------------------------------------------------
sampleMD <- getSampleMetadata(exhaledPtrset)
colnames(sampleMD)[1] <- "individual"  

exhaledPtrset <- setSampleMetadata(exhaledPtrset,sampleMD)
getSampleMetadata(exhaledPtrset)

## ----exportSampleMetada,eval=FALSE--------------------------------------------
#  exportSampleMetada(exhaledPtrset, saveFile = file.path(DirBacteria,"sampleMetadata.tsv"))
#  exhaledPtrset <- importSampleMetadata(exhaledPtrset, file = file.path(DirBacteria,"sampleMetadata.tsv"))

## ----plotRaw_ptrSet-----------------------------------------------------------
plotRaw(exhaledPtrset, mzRange = 59 , fileNames = getFileNames(exhaledPtrset)[1],showVocDB = TRUE)

## ----plotFeatures, message=FALSE, warning=FALSE-------------------------------
plotFeatures(exhaledPtrset,mz=59.049,type="ggplot",colorBy = "individual")

## ---- update------------------------------------------------------------------
exhaledPtrset <- updatePtrSet(exhaledPtrset)

