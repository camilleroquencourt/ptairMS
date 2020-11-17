testthat::context("Test Peak detection")

test_PeakList<-function(){
  
  ## error
  testthat::expect_error(PeakList(NULL))
  
  library(ptairData)
  filePath <-  system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
  file <- readRaw(filePath,mzCalibRef = c(21.022,59.049))
  
  peakList <- PeakList(file, mzNominal = 21)

  ## find isotope of primary ion (always present)
  testthat::expect_true(dim(peakList$peak)[1] >= 1)
  
  peakList <- PeakList(file, mzNominal = 63)
  testthat::expect_true(nrow(peakList$peak) == 3)
}

test_peakSet<-function(){
  
  library(ptairData)
  dir <- system.file("extdata/mycobacteria",  package = "ptairData")
  dirSet <- createPtrSet(dir, setName = "test", mzCalibRef = c(21.022,59.049))
  ListFiles <- names(dirSet@TIC)
  peakLists <- detectPeak(dirSet, mzNominal  = c(21, 60),processFun=ptairMS:::processFileSepExp)
  peakLists2 <- detectPeak(dirSet, mzNominal  = c(21, 60),processFun=ptairMS:::processFileTemporal)

  testthat::expect_equal(peakLists@peakListAligned$Control1.h5$quanti_ncps,
                         peakLists@peakListAligned$Control1.h5$quanti_cps/(peakLists@primaryIon$Control1.h5$primaryIon*488))
  # type
  testthat::expect_is(peakLists,'ptrSet')
  
  # One matrix for each file
  testthat::expect_true( length(peakLists@peakListAligned) == length(ListFiles) )
  testthat::expect_true(all(sapply(peakLists@peakListAligned,function(x) nrow(x)==2)))
  testthat::expect_true(all(sapply(peakLists@peakListRaw,function(x) ncol(x)==9)))
  
  }

testthat::test_that("PeakList() work", test_PeakList())
testthat::test_that("peakSet() work", test_peakSet())

