testthat::context("Test calibration")

test_calib<-function(){
  # Error
  testthat::expect_error(calibration(NULL))

  library(ptairData)
  filePath <-  system.file("extdata/exhaledAir/ind1", "ind1-1.h5", 
                           package = "ptairData")
  file <- readRaw(filePath,calib = FALSE)

  # Warning because mzRef Calib by default or out of the mass axis of example file
  fileCalibrate <- calibration(file,mzCalibRef = c(21.022,59.049),
                               calibrationPeriod=60)

  # check class
  testthat::expect_is(fileCalibrate, 'ptrRaw')
  
  #check format 
  testthat::expect_equal(dim(getCalibrationInfo(fileCalibrate)$calibCoef[[1]]), c(2,1))
  testthat::expect_equal(rownames(getCalibrationInfo(fileCalibrate)$calibCoef[[1]]), c("a","b"))
  
  #check values
  testthat::expect_equal(c(round(getCalibrationInfo(fileCalibrate)$calibCoef[[1]])), c(8839,-219))
  testthat::expect_true(all(getCalibrationInfo(fileCalibrate)$calibError<20))
  
}

testthat::test_that("we can make calibration", test_calib())
