#' .vocRef
#'
#' reference VOC
#'
#' @param mzRgeVn range of mz values
#' @return data frame with the reference VOCs
.ptrms_vocRef <- function(mzRgeVn) {
  vocDF <- .ptrms_refTable("voc")
  if (length(mzRgeVn) == 1) {
    vocRefDF <- vocDF[vocDF[, "nominal"] == round(mzRgeVn), ]
  } else {
    vocRefDF <- vocDF[vocDF[, "mz"] >= mzRgeVn[1] & vocDF[, "mz"] <= mzRgeVn[2], ]
  }
  return(vocRefDF)
}


.ptrms_refTable <- function(refC = c("element",
                                     "voc",
                                     "deLacyCostello2014"),
                            dirC = "//10.0.238.33/Data/Phenostore/data/Exhalomics/reference/") {
  
  read.table(paste0(dirC, refC, ".tsv"),
             header = TRUE,
             quote = "\"",
             sep = "\t",
             stringsAsFactors = FALSE)
  
}
