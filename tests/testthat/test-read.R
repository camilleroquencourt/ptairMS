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
  ret <- readRaw(file,calibTIS = FALSE)
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

test_checkSet <- function(){
  
  library(ptairData)
  directory <- system.file("extdata/mycobacteria", package = "ptairData")
  files <- list.files(directory,full.names = TRUE,recursive = TRUE,pattern = "\\.h5*")
  check<- checkSet(files,mzCalibRef =c(21.022,59.049), fracMaxTIC = 0.6 )
  testthat::expect_is(check,"list")
  
  }


testthat::test_that("we can read an input file.", test_readRaw())
testthat::test_that("createPtrSet works", test_createSet())
testthat::test_that("checkSet works", test_checkSet())
