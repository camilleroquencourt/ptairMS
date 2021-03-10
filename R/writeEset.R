


#### writing ####

#' Exporting an ExpressionSet instance into 3 tabulated files 'dataMatrix.tsv',
#' sampleMetadata.tsv', and 'variableMetadata.tsv'
#'
#' Note that the \code{dataMatrix} is transposed before
#' export (e.g., the samples are written column wise in the 'dataMatrix.tsv'
#' exported file).
#'
#' @param x An S4 object of class \code{ExpressionSet}
#' @param dirName Character: directory where the tables should be written
#' @param overwrite Logical: should existing files be overwritten?
#' @param verbose Logical: should messages be printed?
#' @return No object returned.
#' @rdname writeEset
#' @export
#' @examples
#' data(exhaledPtrset )
#' eset <- ptairMS::alignSamples(exhaledPtrset ) 
#'\dontrun{
#' writeEset(eset, dirName = file.path(getwd(), "processed_dataset"))
#'}
setMethod("writeEset", "ExpressionSet",
          function(x,
                   dirName,
                   overwrite = FALSE,
                   verbose = TRUE){
              
              if (!(file.exists(dirName) && file.info(dirName)[, "isdir"])) {
                if (verbose)
                  message("'", dirName, "' directory has been created")
                dir.create(dirName,
                           showWarnings = verbose)
              }
              
            tableFilesVc <- c(dataMatrix = file.path(dirName, 
                                                     "dataMatrix.tsv"),
                              sampleMetadata = file.path(dirName, 
                                                         "sampleMetadata.tsv"),
                              variableMetadata = file.path(dirName, 
                                                           "variableMetadata.tsv"))

            for (tableC in names(tableFilesVc)) {
              
              tableFileC <- tableFilesVc[tableC]
              
              if (file.exists(tableFileC) && !overwrite)
                stop("The following file already exists:\n", tableFileC,
                     "\nPlease choose another file name.")
              
            }
            
            ## Writing
            
            tdatMN <- Biobase::exprs(x)
            samDF <- Biobase::pData(x)
            varDF <- Biobase::fData(x)
            chkLs <- .checkTableFormat(t(tdatMN), samDF, varDF)
            
            if (!chkLs[["chkL"]]) {
              stop("Sample and/or variable names do not match 
                   between your tables.")
            } else if (chkLs[["ordL"]]) {
              tdatMN <- t(chkLs[["datMN"]])
            }
            
            datDF <- cbind.data.frame(dataMatrix = rownames(tdatMN),
                                      as.data.frame(tdatMN))
            
            utils::write.table(datDF,
                               file = tableFilesVc['dataMatrix'],
                               quote = FALSE,
                               row.names = FALSE,
                               sep = "\t")

              samDF <- cbind.data.frame(sampleMetadata = rownames(samDF),
                                        samDF)
              utils::write.table(samDF,
                                 file = tableFilesVc['sampleMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")

              varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                                        varDF)
              utils::write.table(varDF,
                                 file = tableFilesVc['variableMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")
 
            if (verbose)
              message("The following file(s) have been written:\n",
                  paste(tableFilesVc[!is.na(basename(tableFilesVc))], 
                        collapse = "\n"),
                  "\n")
            
          })

.checkTableFormat <- function(datMNw, samDFw, varDFw,
                            infCw = "interactive") {
  
  chkL <- TRUE
  ordL <- FALSE
  
  if (mode(datMNw) != "numeric") {
    cat("The dataMatrix is not of the 'numeric' type\n")
    chkL <- FALSE
  }
  
  if (!identical(rownames(datMNw), rownames(samDFw))) {
    ## checking sample names
    
    datSamDifVc <- setdiff(rownames(datMNw), rownames(samDFw))
    
    if (length(datSamDifVc)) {
      if (infCw != "none")
        cat("The following samples were found in the dataMatrix column names 
            but not in the sampleMetadata row names:\n", sep = "")
      print(cbind.data.frame(col = as.numeric(vapply(datSamDifVc,
                                                     function(samC) 
                                                       which(rownames(datMNw) == samC),
                                                     FUN.VALUE = 1)),
                             name = datSamDifVc))
      chkL <- FALSE
    }
    
    samDatDifVc <- setdiff(rownames(samDFw), rownames(datMNw))
    
    if (length(samDatDifVc)) {
      if (infCw != "none")
        cat("The following samples were found in the sampleMetadata row names 
            but not in the dataMatrix column names:\n",
            sep = "")
      print(cbind.data.frame(row = as.numeric(vapply(samDatDifVc, 
                                                     function(samC) which(rownames(samDFw) == samC),
                                                     FUN.VALUE =1 )),
                             name = samDatDifVc))
      chkL <- FALSE
    }
    
    if (nrow(datMNw) != nrow(samDFw)) {
      if (infCw != "none")
        cat("The dataMatrix has ", nrow(datMNw), " columns (ie samples) whereas 
            the sampleMetadata has ", nrow(samDFw), " rows\n",
            sep = "")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(datMNw)), rownames(samDFw))) {
      if (infCw != "none")
        cat("The dataMatrix column names start with an 'X' but not the 
            sampleMetadata row names\n", sep = "")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(samDFw)), rownames(datMNw))) {
      if (infCw != "none")
        cat("The sampleMetadata row names start with an 'X' but not the 
            dataMatrix column names\n", sep = "")
      chkL <- FALSE
    } else if (identical(sort(rownames(datMNw)), sort(rownames(samDFw)))) {
      if (infCw != "none")
        cat("Message: Re-ordering dataMatrix sample names to match 
            sampleMetadata\n")
      datMNw <- datMNw[rownames(samDFw), , drop = FALSE]
      stopifnot(identical(sort(rownames(datMNw)), sort(rownames(samDFw))))
      ordL <- TRUE
    } else {
      if (infCw != "none")
        cat("The dataMatrix column names and the sampleMetadata row names are 
            not identical:\n", sep = "")
      print(cbind.data.frame(indice = seq_len(nrow(datMNw)),
                             dataMatrix_columnnames = rownames(datMNw),
                             sampleMetadata_rownames = rownames(samDFw))[rownames(datMNw) != rownames(samDFw), , drop = FALSE])
      chkL <- FALSE
    }
    
  }
  
  if (!identical(colnames(datMNw), rownames(varDFw))) {
    ## checking variable names
    
    datVarDifVc <- setdiff(colnames(datMNw), rownames(varDFw))
    
    if (length(datVarDifVc)) {
      if (infCw != "none")
        cat("The following variables were found in the dataMatrix row names but 
            not in the variableMetadata row names:\n", sep = "")
      print(cbind.data.frame(row = as.numeric(vapply(datVarDifVc, 
                                                     function(varC) which(colnames(datMNw) == varC),
                                                     FUN.VALUE =1 )),
                             name = datVarDifVc))
      chkL <- FALSE
    }
    
    varDatDifVc <- setdiff(rownames(varDFw), colnames(datMNw))
    
    if (length(varDatDifVc)) {
      if (infCw != "none")
        cat("The following variables were found in the variableMetadata row names 
            but not in the dataMatrix row names:\n", sep = "")
      print(cbind.data.frame(row = as.numeric(vapply(varDatDifVc, 
                                                     function(varC) which(rownames(varDFw) == varC),
                                                     FUN.VALUE = )),
                             name = varDatDifVc))
      chkL <- FALSE
    }
    
    if (ncol(datMNw) != nrow(varDFw)) {
      if (infCw != "none")
        cat("The dataMatrix has ",
            nrow(datMNw),
            " rows (ie variables) whereas the variableMetadata has ",
            nrow(varDFw),
            " rows\n",
            sep = "")
      chkL <- FALSE
    } else if (identical(sort(colnames(datMNw)), sort(rownames(varDFw)))) {
      if (infCw != "none")
        cat("Message: Re-ordering dataMatrix variable names to match 
            variableMetadata\n")
      datMNw <- datMNw[, rownames(varDFw), drop = FALSE]
      stopifnot(identical(sort(colnames(datMNw)), sort(rownames(varDFw))))
      ordL <- TRUE
    } else {
      if (infCw != "none")
        cat("\n\nThe dataMatrix row names and the variableMetadata row names 
            are not identical:\n",
            sep = "")
      print(cbind.data.frame(row = seq_len(ncol(datMNw)),
                             dataMatrix_rownames = colnames(datMNw),
                             variableMetadata_rownames = rownames(varDFw))[colnames(datMNw) != rownames(varDFw), , drop = FALSE])
      chkL <- FALSE
    }
  }
  
  return(list(chkL = chkL,
              ordL = ordL,
              datMN = datMNw))
  
}


