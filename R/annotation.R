## annotateVOC (ExpressionSet) ----

#' Putatively annotate VOC mz by using the reference compilation from the literature
#'
#' Putatively annotate VOC mz by using the reference compilation from the literature
#'
#' @param x Expression set object (resp. data.frame) (resp. numeric vector) containing
#' the PTR-MS processed data (resp. containing a column with the mz values) (resp. containing the mz values)
#' @param mzColname Character: column name from the fData (resp. from the data.frame) containing
#' the mz values; if set to NA, featureNames from the Expression set (resp. row.names
#' from the data.frame) will be converted to numerics and used
#' @param ppm Numeric: tolerance
#' @param prefix Character: prefix for the new 'annotation' columns [default: 'vocDB_']
#' @param fields Characer vector: fields of the 'vocDB' database to be queried among:
#' 'mz_Hplus' [default], 'formula_Hplus' [default], 'name' [default], 
#' reference', 'cas', 'chebi', 'chebi.formula', 'chebi.monoisotopic.mass',
#' 'chebi.name', 'chebi.inchi', 'chebi.kegg.compound.id'
#' @param ... Additional parameters to be passed to the internal function
#' @return Returns the Expression set with additional columns in the fData (resp.
#' the data.frame with additional columns) (resp. a data.frame with columns)
#' containing the matched 'mz', 'formula', 'name', 'reference', and 'matrix' putative annotations
#' @examples
#' library(ptairMS)
#' library(ptairData)
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' bacteria.ptrset <- ptairMS::createPtrSet(directory, setName = "bacteria",
#' mzCalibRef = c(21.022,59.049))
#' bacteria.ptrset <- ptairMS::detectPeak(bacteria.ptrset)
#' bacteria.eset <- ptairMS::alignSamples(bacteria.ptrset)
#' bacteria.eset <- ptairMS::annotateVOC(bacteria.eset)
#' Biobase::fData(bacteria.eset)
#' @rdname annotation
#' @export
setMethod("annotateVOC", "ExpressionSet",
          function(x,
                   mzColname = NA,ppm = 20, prefix = "vocDB_",
                   fields = c("mz_Hplus",
                              "formula_Hplus",
                              "cas.name",
                              "reference",
                              "cas",
                              "chebi",
                              "chebi.formula",
                              "chebi.monoisotopic.mass",
                              "chebi.name",
                              "chebi.inchi",
                              "kegg.compound.id")[1:3],
                   ...) {
            
            fdataDF <- Biobase::fData(x)
            
            if (is.na(mzColname)) {
              if (is.null(rownames(fdataDF)))
                stop("The 'mz' column name is missing and rownames from fData are NULL.",
                     call. = FALSE)
              mzVn <- suppressWarnings(as.numeric(rownames(fdataDF)))
              if (any(is.na(mzVn)))
                stop("The 'mz' column name is missing and rownames from fData contain values which cannot be converted into numerics.",
                     call. = FALSE)
            } else {
              
              mzColnameI <- grep(mzColname, colnames(fdataDF))
              
              if (length(mzColnameI) != 1)
                stop("No or multiple columns found in the fdataDF with the '", mzColname, "' name.",
                     call. = FALSE)
              
              mzVn <- fdataDF[, mzColnameI]
              
            }
            
            annotateDF <- .annotate(mz_Hplus = mzVn, ppm, prefix,fields,...)
            
            for (annotateC in colnames(annotateDF))
              fdataDF[, annotateC] <- annotateDF[, annotateC]
            
            Biobase::fData(x) <- fdataDF
            
            x
            
          })

## annotateVOC (data.frame) ----

#' @rdname annotation
#' @export
setMethod("annotateVOC", "data.frame",
          function(x,
                   mzColname = NA,
                   ...) {
            
            if (is.na(mzColname)) {
              if (is.null(rownames(x)))
                stop("The 'mz' column name is missing and rownames are 'NULL'.",
                     call. = FALSE)
              mzVn <- suppressWarnings(as.numeric(rownames(x)))
              if (any(is.na(mzVn)))
                stop("The 'mz' column name is missing and rownames contain values which cannot be converted into numerics.",
                     call. = FALSE)
            } else {
              
              mzColnameI <- grep(mzColname, colnames(x))
              
              if (length(mzColnameI) != 1)
                stop("No or multiple columns found with the '", mzColname, "' name.",
                     call. = FALSE)
              
              mzVn <- x[, mzColnameI]
              
            }
            
            annotateDF <- .annotate(mz_Hplus = mzVn, ...)
            
            for (annotateC in colnames(annotateDF))
              x[, annotateC] <- annotateDF[, annotateC]
            
            x
            
          })

## annotateVOC (numeric) ----

#' @rdname annotation
#' @export
setMethod("annotateVOC", "numeric",
          function(x, ...) {
            
            if (!is.numeric(x))
              stop("'mz' values must be of 'numeric' mode.", call. = FALSE)
            
            .annotate(mz_Hplus = x, ...)
            
            
          })

.annotate <- function(mz_Hplus,
                      ppm = 20,
                      prefix = "vocDB_",
                      fields = c("mz_Hplus",
                                 "formula_Hplus",
                                 "cas.name",
                                 "reference",
                                 "cas",
                                 "chebi",
                                 "chebi.formula",
                                 "chebi.monoisotopic.mass",
                                 "chebi.name",
                                 "chebi.inchi",
                                 "kegg.compound.id")[1:3]) {
  
  if (!is.numeric(mz_Hplus))
    stop("'mz_Hplus' must be of 'numeric' mode.", call. = FALSE)
  
  vocdbDF <- .loadVocDB()
  
  
  matrixVc <- c("blood", "breath", "faeces", "milk", "saliva", "skin", "urine")
  
  fielddbVl <- fields %in% colnames(vocdbDF)
  names(fielddbVl) <- fields
  if ("matrix" %in% names(fielddbVl))
    fielddbVl["matrix"] <- TRUE
  if (sum(!fielddbVl) > 0) {
    warnings("The following fields were not found in the vocDB database and will be ignored:\n", paste(fields[!fielddbVl], collapse = ", "))
  }
  fields <- fields[fielddbVl]
  if ("matrix" %in% fields && length(matrixVc) == 0)
    fields <- setdiff(fields, "matrix")
  
  annotateDF <- data.frame(row.names = as.character(mz_Hplus),
                           stringsAsFactors = FALSE)
  
  for (fieldC in fields)
    annotateDF[, paste0(prefix, fieldC)] <- character(nrow(annotateDF))
  
  for (i in 1:nrow(annotateDF)) {
    
    mzN <- mz_Hplus[i]
    vocVi <- which(abs(mzN - vocdbDF[, "mz_Hplus"]) < ppm * 1e-6 * mzN)
    
    if (length(vocVi) > 0) {
      
      for (fieldC in fields) {
        
        if (fieldC == "matrix") {
          
          vocMatrixVc <- character()
          
          for (k in 1:length(vocVi)) {
            
            vocMatrixVl <- as.logical(vocdbDF[vocVi[k], matrixVc, drop = TRUE])
            vocMatrixVc <- c(vocMatrixVc, matrixVc[vocMatrixVl]) 
            
          }
          
          annotateDF[i, paste0(prefix, "matrix")] <- paste(sort(unique(vocMatrixVc)), collapse = ", ")
          
        } else {
          
          vocFieldVc <- character()
          
          for (k in 1:length(vocVi)) {
            
            vocdbVc <- as.character(vocdbDF[vocVi[k], fieldC, drop = TRUE])
            
            for (splitC in c("|", " / ", ", "))
              vocdbVc <- unlist(lapply(vocdbVc, function(vocdbC)
                unlist(strsplit(vocdbC,
                                split = splitC, fixed = TRUE))))
            
            vocFieldVc <- c(vocFieldVc, vocdbVc)
            
          }
          
          annotateDF[i, paste0(prefix, fieldC)] <- paste(sort(unique(vocFieldVc)), collapse = ", ")
          
          
        }
        
      }
      
    }
    
  }
  
  annotateDF
  
}


.loadVocDB <- function() {
  
  vocdbDF <- utils::read.table(file = system.file("extdata/reference_tables/vocDB.tsv",
                                           package = "ptairMS"),
                        comment.char = "",
                        header = TRUE,
                        quote = "",
                        sep = "\t",
                        stringsAsFactors = FALSE)
  
  vocdbDF[, "mz_Hplus"] <- floor(vocdbDF[, "mz_Hplus"] * 1e5) / 1e5
  
  vocdbDF
  
}



# .build_vocDB.tsv <- function(databases = "de_Lacy_Costello_2014") {
#   
#   if ("de_Lacy_Costello_2014" %in% databases) {
#     
#     # loading the de Lacy Costello database (2014)
#     # [CAS have been converted to CHEBIs by Pierrick Roger-Mele with the 'biodb' package]
#     cost.df <- utils::read.table(file = system.file("extdata/reference_tables/vocDB_tables/de_Lacy_Costello_2014.tsv",
#                                              package = "ptairMS"),
#                           comment.char = "",
#                           header = TRUE,
#                           quote = "",
#                           sep = "\t",
#                           stringsAsFactors = FALSE)
#     
#     message("Initial number annotations: ", nrow(cost.df))
#     
#     
#     # adding extra information from CTS conversion of CAS ids
#     # https://cts.fiehnlab.ucdavis.edu/batch
#     cts.df <- utils::read.table(file = system.file("extdata/reference_tables/vocDB_tables/de_Lacy_Costello_2014_cts-20191029142351.csv",
#                                             package = "ptairMS"),
#                          comment.char = "",
#                          header = TRUE,
#                          quote = "\"",
#                          sep = ",",
#                          stringsAsFactors = FALSE)
#     cts.df[cts.df[, "ChEBI"] == "No result", "ChEBI"] <- ""
#     
#     ## additional ChEBI to those already found by biodb
#     add_chebi.vi <- setdiff(which(cts.df[, "ChEBI"] != ""),
#                             which(cost.df[, "chebi.accession"] != ""))
#     add_chebi.vc <- cts.df[add_chebi.vi, "ChEBI"]
#     add_chebi.vc <- gsub("CHEBI:", "", gsub("\nCHEBI:", "|", add_chebi.vc))
#     add_chebi.vl <- !grepl("|", add_chebi.vc, fixed = TRUE)
#     add_chebi.vi <- add_chebi.vi[add_chebi.vl]
#     add_chebi.vc <- add_chebi.vc[add_chebi.vl]
#     
#     ## getting information for those new ChEBIs with biodb
#     mybiodb <- biodb::Biodb()
#     
#     chebi <- mybiodb$getFactory()$createConn('chebi')
#     
#     select.vl <- sapply(add_chebi.vc, function(chebi.c) {
#       entry.ls <- try(chebi$getEntry(chebi.c))
#       !inherits(entry.ls, "try-error")
#     }) # 4 ChEBI needs to be discarded because of 'UTF-8' errors
# 
#     add_chebi.vi <- add_chebi.vi[select.vl]
#     add_chebi.vc <- add_chebi.vc[select.vl]
#     
#     entries.ls <- chebi$getEntry(add_chebi.vc)
#     chebi.df <- mybiodb$entriesToDataframe(entries.ls,
#                                            fields = c("chebi.id", "formula", "monoisotopic.mass", "name", "inchi"))
#     
#     mybiodb$terminate()
#     
#     add_chebi.vi <- add_chebi.vi[-14] # weird formula (XXX)nH2O
#     chebi.df <- chebi.df[-14, ]
#     cost.df[add_chebi.vi, "chebi.accession"] <- chebi.df[, "chebi.id"]
#     cost.df[add_chebi.vi, "chebi.formula"] <- chebi.df[, "formula"]
#     cost.df[add_chebi.vi, "chebi.monoisotopic.mass"] <- chebi.df[, "monoisotopic.mass"]
#     cost.df[add_chebi.vi, "chebi.name"] <- chebi.df[, "name"]
#     cost.df[add_chebi.vi, "chebi.inchi"] <- chebi.df[, "inchi"]
#     
#     ## getting mz_Hplus and formula_Hplus for those new ChEBI
#     mz_formula.vn <- formula2mass(chebi.df[, "formula"],
#                                           protonizeL = TRUE)
#     cost.df[add_chebi.vi, "mz_Hplus"] <- mz_formula.vn
#     cost.df[add_chebi.vi, "formula_Hplus"] <- names(mz_formula.vn)
#     
#     
#     # removing all rows (molecules) without CAS converted ids (ie unkown mz and formula)
#     cost.df <- cost.df[!is.na(cost.df[, "mz_Hplus"]), ]
#     
#     message("Number of annotations with converted CAS ids: ", nrow(cost.df))
#     
#     # indexing isomers
#     
#     isomer.i <- 1
#     
#     isomer.vi <- numeric(nrow(cost.df))
#     
#     for (k in 1:nrow(cost.df)) {
#       
#       if (isomer.vi[k] == 0) {
#         
#         isomer.vi[which(cost.df[, "formula_Hplus"] == cost.df[k, "formula_Hplus"])] <- isomer.i
#         isomer.i <- isomer.i + 1
#         
#       }
#       
#     }
#     
#     cost.df[, "isomer.group"] <- isomer.vi
#     
#     cas.vi <- suppressWarnings(as.numeric(sapply(cost.df[, "cas.number"],
#                                                  function(cas.c) {
#                                                    gsub("-", "", cas.c)
#                                                  })))
#     
#     cost.df <- cost.df[order(cost.df[, "mz_Hplus"],
#                              cas.vi,
#                              cost.df[, "compound.name"]), ]
#     
#     # concatenating isomers
#     
#     cost.remaining.vi <- 1:nrow(cost.df)
#     cost.discarding.vi <- integer()
#     
#     while (length(cost.remaining.vi) > 0) {
#       
#       cost.selecting.vi <- which(cost.df[, "isomer.group"] == cost.df[cost.remaining.vi[1], "isomer.group"])
#       
#       if (length(cost.selecting.vi) > 1) {
#         
#         for (colname.c in setdiff(colnames(cost.df), c("mz_Hplus",
#                                                        "formula_Hplus",
#                                                        "reference",
#                                                        "chebi.formula",
#                                                        "chebi.monoisotopic.mass"))) {
#           
#           cost.df[cost.selecting.vi[1], colname.c] <- paste(cost.df[cost.selecting.vi, colname.c],
#                                                             collapse = "|")
#           
#         }
#         
#         cost.discarding.vi <- c(cost.discarding.vi, cost.selecting.vi[-1])
#         
#       }
#       
#       cost.remaining.vi <- setdiff(cost.remaining.vi, cost.selecting.vi)
#       
#     }
#     
#     cost.df <- cost.df[-cost.discarding.vi, ]
#     
#     message("Final number of annotations after concatenating isomers: ", nrow(cost.df))
#     
#     colnames(cost.df) <- gsub("cas.number", "cas", colnames(cost.df))
#     colnames(cost.df) <- gsub("compound.name", "cas.name", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.accession", "chebi", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.formula", "formula", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.monoisotopic.mass", "monoisotopic.mass", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.inchi", "inchi", colnames(cost.df))
#     cost.df[, "isomer.group"] <- NULL
#     
#     # building a 'summary' annotation
#     
#     summary.vc <- apply(cost.df, 1,
#                         function(voc.vc) {
#                           
#                           mz.n <- floor(as.numeric(voc.vc["mz_Hplus"]) * 1000)/1000
#                           formula.c <- voc.vc["formula_Hplus"]
#                           name.c <- voc.vc["cas.name"]
#                           if (grepl("|", name.c, fixed = TRUE)) {
#                             name.c <- unlist(strsplit(name.c, "|", fixed = TRUE))[1]
#                             suffix.c <- "..."
#                           } else
#                             suffix.c <- ""
#                           
#                           if (grepl(", ", name.c, fixed = TRUE)) {
#                             name.c <- unlist(strsplit(name.c, ", ", fixed = TRUE))[1]
#                             suffix.c <- "..."
#                           }
#                           
#                           paste0(mz.n, ", ",
#                                  formula.c, ", ",
#                                  name.c,
#                                  suffix.c)
#                           
#                         })
#     
#     cost.df[, "summary"] <- summary.vc
#     
#     cost.df[, "chebi"] <- paste0("CHEBI:", cost.df[, "chebi"])
#     
#     vocDB.df <- cost.df[, c("mz_Hplus",
#                             "formula_Hplus",
#                             "cas.name",
#                             "summary",
#                             "cas",
#                             "blood",
#                             "breath",
#                             "faeces",
#                             "milk",
#                             "saliva",
#                             "skin",
#                             "urine",
#                             "reference",
#                             "chebi",
#                             "monoisotopic.mass",
#                             "formula",
#                             "chebi.name",
#                             "inchi",
#                             "kegg.compound.id")]
#     
#     
#     
#   }
#   
#   utils::write.table(vocDB.df,
#               file = "vocDB.tsv",
#               quote = FALSE,
#               row.names = FALSE,
#               sep = "\t")
#   
#   # 'vocDB.tsv' file to be moved to the 'ptairMS/extdata/reference_tables/' directory
#   
#   invisible(vocDB.df)
#   
# }

#' Compute exact mass.
#'
#' Compute exact mass from an elemental formula
#'
#' @param formulaVc Vector of molecular formulas.
#' @param protonizeL Should a proton be added to the formula?
#' @return Vector of the corresponding (protonated) masses.
#' @export
#' @examples
#' formula2mass("CO2")
formula2mass <- function(formulaVc,
                         protonizeL = TRUE) {
  
  ## Gross, J. (2004). Mass spectrometry: a textbook (Springer). p70.
  ## http://www.ciaaw.org/atomic-masses.htm
  
  elemDF <- utils::read.table(system.file("extdata/reference_tables/elements.tsv",
                                   package = "ptairMS"),
                       header = TRUE,
                       quote = "\"",
                       sep = "\t",
                       stringsAsFactors = FALSE)
  
  elemVn <- elemDF[, "mass"]
  names(elemVn) <- elemDF[, "symbol"]
  
  sapply(formulaVc,
         function(formC) {
           .form2mass(formC, elemVn, protonizeL)
         },
         USE.NAMES = FALSE)
  
}

.form2mass <- function(frmC, eleVn, proL) {
  
  if (frmC == "[H6-16O2-18O+H]+") {
    masN <- 57.04
    names(masN) <- frmC
    return(masN)
  } else if (frmC == "[C2H6-34S+H]+") {
    masN <- 65.04362859832
    names(masN) <- frmC
    return(masN)
  }
  
  if (grepl("[", frmC, fixed = TRUE))
    frmC <- substr(frmC, 2, nchar(frmC))
  
  if (grepl("+H]+", frmC, fixed = TRUE))
    frmC <- gsub("+H]+", "H+", frmC, fixed = TRUE)
  
  if (substr(frmC, nchar(frmC), nchar(frmC)) == "+") {
    ## formula is already protonized (e.g. Herbig09)
    
    ## removing the final '+'
    frmC <- substr(frmC, 1, nchar(frmC) - 1)
    
    atoVn <- .findAtom(frmC, eleVn)
    
    ## removing one 'H'
    if (atoVn["H"] == 1) {
      atoVn <- atoVn[!(names(atoVn) == "H")]
    } else
      atoVn["H"] <- atoVn["H"] - 1
    
    proL <- TRUE
    
  } else {
    
    atoVn <- .findAtom(frmC, eleVn)
    
  }
  
  masN <- sum(eleVn[names(atoVn)] * atoVn)
  
  atoVc <- as.character(atoVn)
  atoVc[atoVc == "1"] <- ""
  frmC <- paste(paste0(names(atoVn),
                       atoVc),
                collapse = "")
  
  if (proL) {
    
    masN <- masN + eleVn["proton"]
    frmC <- paste0("[", frmC, "+H]+")
    
  }
  
  names(masN) <- frmC
  
  return(masN)
  
}

.findAtom <- function(frmC, eleVn) {
  
  frmSplVc <- unlist(strsplit(frmC, ""))
  atoVn <- numeric()
  atoI <- 0
  splI <- 1
  while (splI <= length(frmSplVc)) {
    atoI <- atoI + 1
    if (splI == length(frmSplVc)) {
      atoVn <- c(atoVn, 1)
      names(atoVn)[atoI] <- frmSplVc[splI]
      break
    } else if (splI == length(frmSplVc) - 1 &&
               frmSplVc[splI + 1] %in% letters) {
      atoVn <- c(atoVn, 1)
      names(atoVn)[atoI] <- paste0(frmSplVc[splI],
                                   frmSplVc[splI + 1])
      break
    } else {
      atoVn <- c(atoVn, 0)
      if (frmSplVc[splI + 1] %in% letters) {
        names(atoVn)[atoI] <- paste0(frmSplVc[splI],
                                     frmSplVc[splI + 1])
        splI <- splI + 1
      } else
        names(atoVn)[atoI] <- frmSplVc[splI]
      
      numI <- 1
      while ((splI + numI) <= length(frmSplVc) &&
             !(frmSplVc[splI + numI] %in% LETTERS))
        numI <- numI + 1
      numI <- numI - 1
      if (numI == 0) {
        atoVn[atoI] <- 1
      } else {
        atoVn[atoI] <- as.numeric(paste(frmSplVc[(splI + 1):(splI + numI)],
                                        collapse = ""))
      }
      if ((splI + numI) == length(frmSplVc)) {
        break
      } else
        splI <- splI + numI + 1
    }
  }
  atoVn <- table(rep(names(atoVn), times = atoVn))
  if (any(!(names(atoVn) %in% names(eleVn)))) {
    stop(paste(names(atoVn)[!(names(atoVn) %in% names(eleVn))],
               collapse = ", "), " mass(es) cannot be provided currently",
         call. = FALSE)
  }
  atoVn <- atoVn[names(eleVn)[names(eleVn) %in% names(atoVn)]]
  
  atoVn
  
}

# building the reference table of element masses
# source: http://www.ciaaw.org/atomic-masses.htm
.build_elements.tsv <- function() {
  
  ciaaw.df <- utils::read.table(system.file("extdata/reference_tables/vocDB_tables/ciaaw.tsv",
                                   package = "ptairMS"),
                       header = TRUE,
                       quote = "\"",
                       sep = "\t",
                       stringsAsFactors = FALSE)
  
  element.df <- ciaaw.df[!is.na(ciaaw.df[, "Z"]), ]
  
  element.df[, "mass"] <- sapply(element.df[, "mass"],
         function(mass.c) {
           
           if (mass.c == " 12(exact") {
             
             return(12)
             
           } else {
             
             mass.c <- substr(mass.c, 1, nchar(mass.c) - 3)
             mass_split.vc <- unlist(strsplit(mass.c, split = ""))
             mass_split.vl <- !is.na(as.numeric(mass_split.vc)) | mass_split.vc == "."
             mass_split.vc <- mass_split.vc[mass_split.vl]
             mass.n <- floor(as.numeric(paste(mass_split.vc, collapse = "")) * 1e6) / 1e6
             
             return(mass.n)
             
           }
           
         })
  
  # adding supplementary information about electron, proton, and neutron
  supp.df <- element.df[1:3, ]
  supp.df[, "Z"] <- rep(NA_integer_, 3)
  supp.df[, "symbol"] <- c("electron", "proton", "neutron")
  supp.df[, "element"] <- rep("", 3)
  supp.df[, "A"] <- rep("", 3)
  supp.df[, "mass"] <- c(0.000548, 1.0072765, 1.008665549)
  
  element.df <- rbind.data.frame(supp.df,
                                 element.df,
                                 stringsAsFactors = FALSE)
  
  utils::write.table(element.df,
              file = "elements.tsv",
              row.names = FALSE,
              sep = "\t")
  
}
