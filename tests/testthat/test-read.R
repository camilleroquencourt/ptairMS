testthat::context("Test file reading.")

test_readRaw <- function() {

  library(ptairData)
  file <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")

  # No file
  testthat::expect_error(readRaw(NULL))
  testthat::expect_error(readRaw(NA_character_))
  testthat::expect_error(readRaw(character()))
  testthat::expect_error(readRaw(10L))
  testthat::expect_error(readRaw('foo'))

  # Ouput 
  ret <- readRaw(file,calib=FALSE)
  mz <- rhdf5::h5read(file,"/FullSpectra/MassAxis")
  testthat::expect_is(ret, 'ptrRaw')
  testthat::expect_equal(length(ret@mz),length(mz))
  
}

test_createSet <- function(){
  
  library(ptairData)
  dir <- system.file("extdata/mycobacteria", package = "ptairData")
  pSet <- createPtrSet(dir,setName = "test", mzCalibRef = c(21.022,59.049))
  testthat::expect_is(pSet,"ptrSet")
  testthat::expect_equal(length(pSet@mzCalibRef),6)
  testthat::expect_equal(length(pSet@signalCalibRef),6)
  testthat::expect_equal(length(pSet@errorCalibPpm),6)
  testthat::expect_equal(length(pSet@resolution),6)
  testthat::expect_equal(length(pSet@TIC),6)
  testthat::expect_equal(length(pSet@timeLimit),6)
  
  }


test_knots<- function(){
  library(ptairData)
  dir <- system.file("extdata/mycobacteria", package = "ptairData")
  
  pSet <- createPtrSet(dir,setName = "test", 
                       mzCalibRef = c(21.022,59.049))
  
  knot1<-pSet@knots
  pSet<-defineKnots(pSet)
  knot2<-pSet@knots
  testthat::expect_equal(knot1,knot2)
  testthat::expect_equal(knot2[[1]],
                         c(0.000000,  3.117667,  6.235335,  9.353002, 12.470669, 15.588337, 18.706004, 
                           21.823671, 24.941339, 28.059006, 31.176673, 34.294341,
                           37.412008, 40.529675, 43.647343, 46.765010, 49.882677, 53.000345))
  
  pSeterror<-suppressWarnings(defineKnots(pSet,knotsPeriod = 0.5))
  testthat::expect_null(pSeterror@knots$Control1.h5)
  
  pSet2<-defineKnots(pSet,method = "uniform",knotsPeriod = 5)
  testthat::expect_equal(pSet2@knots[[1]],c(0.000000,  5.300034, 10.600069, 15.900103, 21.200138,
                                            26.500172, 31.800207, 37.100241, 42.400276, 47.700310, 53.000345))
  pSet3<-timeLimits(pSet,fracMaxTIC = 0,redefineKnots = TRUE)
  testthat::expect_equal(pSet3@knots[[1]], c(0.000000,  3.117667,  6.235335,  9.353002, 12.470669, 15.588337, 
                                             18.706004, 21.823671, 24.941339, 28.059006, 31.176673, 34.294341,
                                             37.412008, 40.529675, 43.647343, 46.765010, 49.882677, 53.000345))
  
  
  }

testthat::test_that("we can read an input file.", test_readRaw())
testthat::test_that("createPtrSet works", test_createSet())
testthat::test_that("defineKnots works", test_knots())

