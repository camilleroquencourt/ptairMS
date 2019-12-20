testthat::context("Test m/z annotation")

test_that("annotateVOC function", {
  
  bacteria_dir.c <- system.file("extdata/mycobacteria",  package = "ptairData")
  bacteria.ptrset <- ptairMS::createPtrSet(bacteria_dir.c, setName = "bacteria",
                                           mzCalibRef = c(21.022,59.049), fracMaxTIC = 0.8, saveDir = NULL)
  bacteria.ptrset <- detectPeak(bacteria.ptrset, mzNominal = c(21,59))
  bacteria.eset <- alignSamples(bacteria.ptrset)
  
  # numeric vector
  
  annotateVectorDF <- ptairMS::annotateVOC(as.numeric(Biobase::featureNames(bacteria.eset)))
  
  testthat::expect_identical(annotateVectorDF["59.0492", "vocDB_cas.name"],
                             "2-propen-1-ol, acetone, propanal")
  
  # data.frame
  
  fdataDF <- Biobase::fData(bacteria.eset)
  
  annotateDataFrameDF <- ptairMS::annotateVOC(fdataDF)
  
  testthat::expect_identical(annotateDataFrameDF["59.0492", "vocDB_cas.name"],
                             "2-propen-1-ol, acetone, propanal")
  
  # ExpressionSet
  
  bacteria.eset <- ptairMS::annotateVOC(bacteria.eset)
  
  testthat::expect_identical(Biobase::fData(bacteria.eset)["59.0492", "vocDB_cas.name"],
                             "2-propen-1-ol, acetone, propanal")
  
})
  