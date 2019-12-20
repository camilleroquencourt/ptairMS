testthat::context("Delimitation in time dimension")

test_delimitTime<-function(){
  
  # No matrix in input 
  testthat::expect_error(timeLimits(NULL))
  testthat::expect_error(timeLimits(rep(0,1)))

  # Type and Row names
  filePath <- system.file("extdata/exhaledAir/ind1", "ind1-1.h5", package = "ptairData")
  file<-readRaw(filePath, calibTIS = FALSE)
  res<-timeLimits(file,plot=FALSE)
  
  testthat::expect_is(res,'matrix')
  testthat::expect_identical(c("start","end"),row.names(res))
  testthat::expect_equal(c(res),c(8,19,55,69))
  
  # Bad ratio input 
  testthat::expect_error(timeLimits(file,intensityFrac = 2))
  
}

testthat::test_that("timeLimits", test_delimitTime())
