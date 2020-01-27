# vi: fdm=marker ts=4 et cc=80 tw=80

testthat::context("Test alignment functions.")

test_findEqualGreater <- function() {
  ret <- findEqualGreaterM(c(1,2,2,4), c(3))
  testthat::expect_identical(ret,4)
}

test_align <- function() {
  testthat::expect_error(align(NULL))
  wrongMat1 <- matrix(0,ncol=2)
  testthat::expect_error(align(wrongMat1))
  wrongMat2 <- cbind(mz=c("21","22"),group=c(1,2))
  testthat::expect_error(align(wrongMat2))

  goodMat<- cbind(Mz=c(21.002,21.0021,21.0019),quanti=c(100,120,115),group=c(1,2,3))
  res <- align(goodMat)
  testthat::expect_equal(length(res),1)
  testthat::expect_is(res,'list')
  testthat::expect_equal(dim(res[[1]]),c(3,3))
}

test_alignExpirations <- function(){
  
  Ex<- cbind(Mz=c(21.002,21.0021,21.0019,21.002),quanti=c(100,120,115,40), 
                       delta_mz=rep(0.5,4), resolution = rep(5000,4),group=c(1,2,3,0))
  res <- alignExpirations(as.data.frame(Ex))
  testthat::expect_equal(dim(res)[1],1)
  testthat::expect_is(res,'data.table')  
}

test_alignSamples <- function(){
  directory <-  system.file("extdata/mycobacteria",  package = "ptairData")
  dirSet <- createPtrSet(directory, setName = "test", mzCalibRef =c(21.022,59.049))
  dirSet <- detectPeak(dirSet, mzNominal = c(21,59))
  eset <- alignSamples(dirSet, bgThreshold = 0)
  
  # output is expression set
  testthat::expect_is(eset,'ExpressionSet')
  
  #deux pics aligné
  testthat::expect_equal(nrow(Biobase::exprs(eset)),2)
  testthat::expect_equal(ncol(Biobase::exprs(eset)),6)
  
  #test filer
  eset <- alignSamples(dirSet,bgThreshold = 2)
  testthat::expect_equal(ncol(Biobase::fData(eset)),6)
  testthat::expect_equal(nrow(Biobase::exprs(eset)),1)
  }

testthat::test_that("findEqualGreater() works correctly.", test_findEqualGreater())
testthat::test_that("align() works correctly.", test_align())
testthat::test_that("alignExpirations() works correctly.", test_alignExpirations())
testthat::test_that("alignSamples() works correctly.", test_alignSamples())
