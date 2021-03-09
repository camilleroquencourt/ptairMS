testthat::context("Test Peak detection")

test_PeakList<-function(){
  ## error
  testthat::expect_error(PeakList(NULL))
  
  library(ptairData)
  filePath <-  system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
  file <- readRaw(filePath,mzCalibRef = c(21.022,59.049))
  
  peakList <- PeakList(file, mzNominal = 21,fctFit = "averagePeak")
  peakList2 <- PeakList(file, mzNominal = 21,fctFit = "sech2")
  peakList3 <- PeakList(file, mzNominal = 21,fctFit = "asymGauss")
  
  ## find isotope of primary ion (always present)
  testthat::expect_true(dim(peakList$peak)[1] >= 1)
  testthat::expect_true(dim(peakList2$peak)[1] >= 1)
  testthat::expect_true(dim(peakList3$peak)[1] >= 1)
  
  peakList <- PeakList(file, mzNominal = 63)
  testthat::expect_true(nrow(peakList$peak) == 3)
}

test_peakSet<-function(){
  
  library(ptairData)
  dir <- system.file("extdata/mycobacteria",  package = "ptairData")
  dirSet <- createPtrSet(dir = dir, setName = "test", mzCalibRef = c(21.022,59.049))
  ListFiles <- names(dirSet@TIC)
  peakLists <- detectPeak(dirSet, smoothPenalty = 0,
                          resolutionRange =  c(3000,5000,9000),
                          mzNominal = c(21,59))

  # type,
  testthat::expect_is(peakLists,'ptrSet')
  testthat::expect_true(all(sapply(peakLists@peakList,function(x) nrow(Biobase::exprs(x))>=2)))
  
  
  }

testthat::test_that("PeakList() work", test_PeakList())
testthat::test_that("peakSet() work", test_peakSet())

