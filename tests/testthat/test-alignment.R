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


test_alignSamples <- function(){
  library(ptairData)
  directory <-  system.file("extdata/mycobacteria",  package = "ptairData")
  dirSet <- createPtrSet(directory, setName = "test", mzCalibRef =c(21.022,59.049))
  dirSet <- detectPeak(dirSet, mzNominal = c(21,59),smoothPenalty = 0)
  eset <- alignSamples(X = dirSet,quanti="ppb",pValGreaterThres = 1)
  
  # output is expression set
  testthat::expect_is(eset,'ExpressionSet')
  
  #deux pics aligné
  testthat::expect_equal(nrow(Biobase::exprs(eset)),2)
  testthat::expect_equal(ncol(Biobase::exprs(eset)),6)
  
  #test filer
  eset <- alignSamples(dirSet,pValGreaterThres = 0.05)
  testthat::expect_equal(ncol(Biobase::fData(eset)),7)
  testthat::expect_equal(nrow(Biobase::exprs(eset)),1)
}

test_impute<-function(){
  library(ptairData)
  directory <-  system.file("extdata/mycobacteria",  package = "ptairData")
  dirSet <- createPtrSet(directory, setName = "test", mzCalibRef =c(21.022,59.049))
  dirSet <- detectPeak(dirSet,mz=c(21,63,77),minIntensity =50,smoothPenalty = 0)
  eset <- alignSamples(dirSet, pValGreaterThres = 1,fracGroup = 0,bgCorrected =FALSE)
  testthat::expect_equal(sum(is.na(Biobase::exprs(eset))),3)
  X<-Biobase::exprs(eset)
  eset<-impute(eset,dirSet)
  Ximpute<-Biobase::exprs(eset)
  testthat::expect_equal(sum(is.na(Biobase::exprs(eset))),0)
  testthat::expect_true(round(Ximpute[2,1],2)<1)
}

testthat::test_that("findEqualGreater() works correctly.", test_findEqualGreater())
testthat::test_that("align() works correctly.", test_align())
testthat::test_that("alignSamples() works correctly.", test_alignSamples())
testthat::test_that("impute() works correctly.", test_impute())
