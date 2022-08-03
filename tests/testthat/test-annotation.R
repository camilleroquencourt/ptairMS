testthat::context("Test m/z annotation")

test_annotateVOC_and_isotope <- function(){
  library(ptairData)
  bacteria_dir.c <- system.file("extdata/mycobacteria",  package = "ptairData")
  bacteria.ptrset <- createPtrSet(bacteria_dir.c, setName = "bacteria",
                                  mzCalibRef = c(21.022,59.049), fracMaxTIC = 0.8, saveDir = NULL)
  bacteria.ptrset <- detectPeak(bacteria.ptrset, mzNominal = c(21,59,60),resolutionRange = c(3000,5000,10000))
  bacteria.eset <-  alignSamples(bacteria.ptrset,pValGreaterThres = 0.05,quanti="ppb",fracExp = 1)
  
  # numeric vector
  annotateVectorDF <- ptairMS::annotateVOC(as.numeric(Biobase::featureNames(bacteria.eset)))
  testthat::expect_identical(annotateVectorDF[1, "vocDB_name_iupac"],
                             "prop-2-en-1-ol, propan-2-one, propanal")
  
  # data.frame
  fdataDF <- Biobase::fData(bacteria.eset)
  annotateDataFrameDF <- annotateVOC(fdataDF)
  testthat::expect_identical(annotateDataFrameDF[1, "vocDB_name_iupac"],
                             "prop-2-en-1-ol, propan-2-one, propanal")
  
  # ExpressionSet
  bacteria.eset <- annotateVOC(bacteria.eset)
  testthat::expect_identical(Biobase::fData(bacteria.eset)[1, "vocDB_name_iupac"],
                             "prop-2-en-1-ol, propan-2-one, propanal")
  
  bacteria.eset <-findIsotope(bacteria.eset)
  testthat::expect_equal( Biobase::fData(bacteria.eset)[1,"isotope"],row.names(Biobase::fData(bacteria.eset))[2])
  
}

test_that("annotateVOC function and findIsotope",test_annotateVOC_and_isotope())

  