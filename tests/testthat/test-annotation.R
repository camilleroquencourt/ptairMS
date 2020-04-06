testthat::context("Test m/z annotation")

test_annotateVOC <- function(){
  library(ptairData)
  bacteria_dir.c <- system.file("extdata/mycobacteria",  package = "ptairData")
  bacteria.ptrset <- createPtrSet(bacteria_dir.c, setName = "bacteria",
                                  mzCalibRef = c(21.022,59.049), fracMaxTIC = 0.8, saveDir = NULL)
  bacteria.ptrset <- detectPeak(bacteria.ptrset, mzNominal = c(21,59,60))
  bacteria.eset <- alignSamples(bacteria.ptrset)
  
  # numeric vector
  annotateVectorDF <- ptairMS::annotateVOC(as.numeric(Biobase::featureNames(bacteria.eset)))
  testthat::expect_identical(annotateVectorDF["59.0492", "vocDB_name_iupac"],
                             "prop-2-en-1-ol, propan-2-one, propanal")
  
  # data.frame
  fdataDF <- Biobase::fData(bacteria.eset)
  annotateDataFrameDF <- annotateVOC(fdataDF)
  testthat::expect_identical(annotateDataFrameDF["59.0492", "vocDB_name_iupac"],
                             "prop-2-en-1-ol, propan-2-one, propanal")
  
  # ExpressionSet
  bacteria.eset <- annotateVOC(bacteria.eset)
  testthat::expect_identical(Biobase::fData(bacteria.eset)["59.0492", "vocDB_name_iupac"],
                             "prop-2-en-1-ol, propan-2-one, propanal")
  
}

test_isotope<-function(){
  library(ptairData)
  directory <- system.file("extdata/mycobacteria",  package = "ptairData")
  bacteria.ptrset <- createPtrSet(directory, setName = "bacteria",
  mzCalibRef = c(21.022,59.049))
  bacteria.ptrset <- detectPeak(bacteria.ptrset,mz=c(59,60))
  bacteria.eset <- alignSamples(bacteria.ptrset,fracGroup=1,bgThreshold = 0)
  bacteria.eset <-findIsotope(bacteria.eset)
  testthat::expect_equal( Biobase::fData(bacteria.eset)[1,"isotope"],row.names(Biobase::fData(bacteria.eset))[2])
}


test_that("annotateVOC function",test_annotateVOC())
test_that("findIsotope function",test_isotope())
  