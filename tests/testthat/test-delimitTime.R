testthat::context("Delimitation in time dimension")

test_delimitTime<-function(){
  
  # No matrix in input 
  testthat::expect_error(timeLimits(NULL))
  testthat::expect_error(timeLimits(rep(0,1)))

  # Type and Row names
  library(ptairData)
  filePath <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
  file<-readRaw(filePath, calibTIS = FALSE)
  res<-timeLimits(file)
  resExp<-res$exp
  resBg<-res$backGround
  
  testthat::expect_is(resExp,'matrix')
  testthat::expect_identical(c("start","end"),row.names(resExp))
  testthat::expect_equal(c(resExp),c(11,22,31,45))
  
  testthat::expect_is(resBg,'integer')
  testthat::expect_equal(resBg,c(1,2,3,4,5,6,7,25,26,27,28,47,48,49))
  
  # Bad ratio input 
  testthat::expect_error(timeLimits(file,fracMaxTIC = 2))
  
}

testthat::test_that("timeLimits", test_delimitTime())
