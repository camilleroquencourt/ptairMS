## annotateVOC (ExpressionSet) ----

#' Putatively annotate VOC mz by using the reference compilation from the literature
#'
#' Putatively annotate VOC mz by using the reference compilation from the literature, 
#' and detect isotope thanks to \code{findIsotope} function. 
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
#' library(ptairData)
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' bacteria.ptrset <- createPtrSet(directory, setName = "bacteria",
#' mzCalibRef = c(21.022,59.049))
#' bacteria.ptrset <- detectPeak(bacteria.ptrset)
#' bacteria.eset <- alignSamples(bacteria.ptrset)
#' bacteria.eset <- annotateVOC(bacteria.eset)
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
                              "kegg.compound.id")[c(1,2,3)],
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
            xIso<-try(findIsotope(x))
            if(!is.null(attr(xIso,"condition"))) 
              return(x) 
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
                                 "kegg.compound.id")[c(1,2,3)]) {
  
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
  
  for (i in seq_len(nrow(annotateDF))) {
    
    mzN <- mz_Hplus[i]
    vocVi <- which(abs(mzN - vocdbDF[, "mz_Hplus"]) < ppm * 1e-6 * mzN)
    
    if (length(vocVi) > 0) {
      
      for (fieldC in fields) {
        
        if (fieldC == "matrix") {
          
          vocMatrixVc <- character()
          
          for (k in seq_along(vocVi)) {
            
            vocMatrixVl <- as.logical(vocdbDF[vocVi[k], matrixVc, drop = TRUE])
            vocMatrixVc <- c(vocMatrixVc, matrixVc[vocMatrixVl]) 
            
          }
          
          annotateDF[i, paste0(prefix, "matrix")] <- paste(sort(unique(vocMatrixVc)), collapse = ", ")
          
        } else {
          
          vocFieldVc <- character()
          
          for (k in seq_along(vocVi)) {
            
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
  
  vapply(formulaVc,
         function(formC) {
           .form2mass(formC, elemVn, protonizeL)
         },FUN.VALUE = 1.01,
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

#'Isotope detection and validation of isotope clusters in PTR-TOF-MS peak table
#'
#'This function identify possible isotope cluster of C13,O17 and O18 atomes after alignment. It writes in the 
#'features data of the expressionSet (Biobase::fData()) the isotope m/z in the "isotope" column at the ligne 
#'of the more intenisf peak of the clusters.
#'@param eSet an expression set of PTR-TOF-MS data aligned 
#'@param ppm presision for the mz matching
#'@return an expresion with the column isotope added in teh features data
#' @examples
#' library(ptairData)
#' directory <- system.file("extdata/mycobacteria",  package = "ptairData")
#' bacteria.ptrset <- createPtrSet(directory, setName = "bacteria",
#' mzCalibRef = c(21.022,59.049))
#' bacteria.ptrset <- detectPeak(bacteria.ptrset)
#' bacteria.eset <- alignSamples(bacteria.ptrset,fracGroup=0.9)
#' bacteria.eset<- impute(bacteria.eset,bacteria.ptrset)
#' bacteria.eset <-findIsotope(bacteria.eset)
#' Biobase::fData(bacteria.eset)[,"isotope",drop=FALSE]
#' @export
findIsotope<-function(eSet,ppm=100){
  X<-Biobase::exprs(eSet)
  fDATA<-Biobase::fData(eSet)
  mz<-as.numeric(row.names(X))
  
  #FIND AND VALIDE ISOTOPE GROUP
  for (i in seq_along(mz)){
    iso <- ptairMS:::isotopeMzMatching(mz[i], mz[(i+1):length(mz)],ppm)
    if(length(iso)){
      if(validateGroup(c(mz[i],iso),X)){
        fDATA[i,"isotope"]<- paste(iso,collapse = "/")
      }
    }
  }
  Biobase::fData(eSet)<-fDATA
  return(eSet)
}


#find all possible isotope match in m/z dimension 
isotopeMzMatching<-function(m,mzSub,ppm,max=1){
  diffN <- 1.003355
  diffO<-c(1.0042,2.0042) #C13, O17, O18
  iso<-list("C13"=NULL,"O17"=NULL,"O18"=NULL)
  for(j in seq(max)){
    Iso<-lapply(c(diffN,diffO),function(diff) mzSub[which(abs(mzSub-(m+j*diff))*10^6/m < 100)])
    iso<-mapply(c, iso, Iso, SIMPLIFY=FALSE)
  }
  iso<-unique(Reduce(c,iso))
  return(iso)
}

#validate isotope group
validateGroup<-function(groupIso,X){
  
  #correlation inter sample
  testCorPval<-vapply(as.character(groupIso)[-1], function(y) cor.test(X[as.character(groupIso)[1],],
                                                                       X[y,],alternative = c("greater"))$p.value,1.1)
  testCor<- testCorPval < 0.01
  
  #ratio
  ratio<-X[as.character(groupIso)[-1],]/
    matrix(
      rep(X[as.character(groupIso)[1],,drop=F],2),
      nrow=length(as.character(groupIso)[-1]),byrow=TRUE)
  
  testRatio<- apply(ratio,1,function(x) median(x,na.rm = T)) < 0.5
  
  return(all(c(testCor,testRatio)))
}

