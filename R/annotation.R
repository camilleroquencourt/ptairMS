## annotateVOC (ExpressionSet) ----

#' @rdname annotation
#' @export
setMethod("annotateVOC", "ExpressionSet",
          function(x,
                   ionMassColname = "ion_mass", ppm = 20, prefix = "vocDB_",
                   fields = c("ion_mass",
                              "ion_formula",
                              "formula",
                              "mass_monoiso",
                              "name_iupac",
                              "pubchem_cid",
                              "inchi",
                              "inchikey",
                              "ref_year",
                              "ref_pmid",
                              "disease_name",
                              "disease_meshid")[c(1,2,5)]) {
            
            fdataDF <- Biobase::fData(x)
            
            ionMassColnameI <- which(colnames(fdataDF) == ionMassColname)
            
            if (length(ionMassColnameI) != 1)
              stop("No or multiple columns found in the fdataDF with the '", ionMassColname, "' name.",
                   call. = FALSE)
            
            ion_mass.vn <- fdataDF[, ionMassColnameI]
            
            annotateDF <- .annotate(ion_mass = ion_mass.vn,
                                    ppm = ppm,
                                    prefix = prefix,
                                    fields = fields)
            
            for (annotateC in colnames(annotateDF))
              fdataDF[, annotateC] <- annotateDF[, annotateC]
            
            Biobase::fData(x) <- fdataDF
            x <- findIsotope(x)
            
            x
            
          })

## annotateVOC (data.frame) ----

#' @rdname annotation
#' @export
setMethod("annotateVOC", "data.frame",
          function(x,
                   ionMassColname = "ion_mass", ppm = 20, prefix = "vocDB_",
                   fields = c("ion_mass",
                              "ion_formula",
                              "formula",
                              "mass_monoiso",
                              "name_iupac",
                              "pubchem_cid",
                              "inchi",
                              "inchikey",
                              "ref_year",
                              "ref_pmid",
                              "disease_name",
                              "disease_meshid")[c(1,2,5)]) {
            
            ionMassColnameI <- which(colnames(x) == ionMassColname)
            
            if (length(ionMassColnameI) != 1)
              stop("No or multiple columns found with the '", ionMassColname, "' name.",
                   call. = FALSE)
            
            ion_mass.vn <- x[, ionMassColnameI]
            
            annotateDF <- .annotate(ion_mass = ion_mass.vn)
            
            for (annotateC in colnames(annotateDF))
              x[, annotateC] <- annotateDF[, annotateC]
            
            x
            
          })

## annotateVOC (numeric) ----

#' @rdname annotation
#' @export
setMethod("annotateVOC", "numeric",
          function(x,
                   ionMassColname = "",
                   ppm = 20,
                   prefix = "vocDB_",
                   fields = c("ion_mass",
                              "ion_formula",
                              "formula",
                              "mass_monoiso",
                              "name_iupac",
                              "pubchem_cid",
                              "inchi",
                              "inchikey",
                              "ref_year",
                              "ref_pmid",
                              "disease_name",
                              "disease_meshid")[c(1,2,5)]) {
            
            if (!is.numeric(x))
              stop("Ion mass values must be of 'numeric' mode.", call. = FALSE)
            
            .annotate(ion_mass = x,
                      ppm = ppm,
                      prefix = prefix,
                      fields = fields)
            
            
          })

.annotate <- function(ion_mass,
                      ppm = 20,
                      prefix = "vocDB_",
                      fields = c("ion_mass",
                                 "ion_formula",
                                 "formula",
                                 "mass_monoiso",
                                 "name_iupac",
                                 "pubchem_cid",
                                 "inchi",
                                 "inchikey",
                                 "ref_year",
                                 "ref_pmid",
                                 "disease_name",
                                 "disease_meshid")[c(1,2,5)]) {
  
  if (!is.numeric(ion_mass))
    stop("'ion_mass' must be of 'numeric' mode.", call. = FALSE)
  
  
  vocdbDF <- .loadVocDB()
  
  
  fielddbVl <- fields %in% colnames(vocdbDF)
  names(fielddbVl) <- fields
  
  if (sum(!fielddbVl) > 0)
    warnings("The following fields were not found in the vocDB database and will be ignored:\n", paste(fields[!fielddbVl], collapse = ", "))
  
  fields <- fields[fielddbVl]
  
  annotateDF <- data.frame(row.names = as.character(ion_mass),
                           stringsAsFactors = FALSE)
  
  for (fieldC in fields)
    annotateDF[, paste0(prefix, fieldC)] <- character(nrow(annotateDF))
  
  for (i in seq_len(nrow(annotateDF))) {
    
    massN <- ion_mass[i]
    vocVi <- which(abs(massN - vocdbDF[, "ion_mass"]) < ppm * 1e-6 * massN)
    
    if (length(vocVi) > 0) {
      
      for (fieldC in fields) {
        
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
  
  annotateDF
  
}


.loadVocDB <- function() {
  
  vocdbDF <- utils::read.table(file = system.file("extdata/reference_tables/vocDB.tsv",
                                                  package = "ptairMS"),
                               check.names = FALSE,
                               comment.char = "",
                               header = TRUE,
                               quote = "\"",
                               sep = "\t",
                               stringsAsFactors = FALSE)
  
  vocdbDF[, "ion_mass"] <- floor(vocdbDF[, "ion_mass"] * 1e5) / 1e5
  
  vocdbDF
  
}





#' Compute exact mass.
#'
#' Compute exact mass from an elemental formula
#'
#' @param formula.vc Vector of molecular formulas.
#' @param protonate.l Should a proton be added to the formula?
#' @return Vector of the corresponding (protonated) masses.
#' @export
#' @examples
#' formula2mass("CO2")
formula2mass <- function(formula.vc,
                         protonate.l = TRUE) {
  
  ## Gross, J. (2004). Mass spectrometry: a textbook (Springer). p70.
  ## http://www.ciaaw.org/atomic-masses.htm
  
  atomic_weights.df <- utils::read.table(system.file("extdata/reference_tables/atomic_weights.tsv",
                                                     package = "ptairMS"),
                                         header = TRUE,
                                         quote = "\"",
                                         sep = "\t",
                                         stringsAsFactors = FALSE)
  
  atomic_weights.vn <- atomic_weights.df[, "mass"]
  names(atomic_weights.vn) <- atomic_weights.df[, "symbol"]
  
  mass_and_formula.mat <- vapply(formula.vc,
                                 function(formula.c) {
                                   .formula2mass(formula.c = formula.c,
                                                 atomic_weights.vn = atomic_weights.vn,
                                                 protonate.l = protonate.l)
                                 }, FUN.VALUE = list(mass = 12, formula = "C"),
                                 USE.NAMES = FALSE)
  
  mass_and_formula.vn <- unlist(mass_and_formula.mat[1, ])
  names(mass_and_formula.vn) <- unlist(mass_and_formula.mat[2, ])
  
  return(mass_and_formula.vn)
  
}

.formula2mass <- function(formula.c, atomic_weights.vn, protonate.l) {
  
  if (formula.c == "[H6-16O2-18O+H]+") {
    mass.n <- 57.04
    names(mass.n) <- formula.c
    return(mass.n)
  } else if (formula.c == "[C2H6-34S+H]+") {
    mass.n <- 65.04362859832
    names(mass.n) <- formula.c
    return(mass.n)
  }
  
  if (grepl("[", formula.c, fixed = TRUE))
    formula.c <- substr(formula.c, 2, nchar(formula.c))
  
  if (grepl("+H]+", formula.c, fixed = TRUE))
    formula.c <- gsub("+H]+", "H+", formula.c, fixed = TRUE)
  
  if (substr(formula.c, nchar(formula.c), nchar(formula.c)) == "+") {
    ## formula is already protonized (e.g. Herbig09)
    
    ## removing the final '+'
    formula.c <- substr(formula.c, 1, nchar(formula.c) - 1)
    
    atoms.vn <- .findAtom(formula.c = formula.c,
                          atomic_weights.vn = atomic_weights.vn)
    
    ## removing one 'H'
    if (atoms.vn["H"] == 1) {
      atoms.vn <- atoms.vn[!(names(atoms.vn) == "H")]
    } else
      atoms.vn["H"] <- atoms.vn["H"] - 1
    
    protonate.l <- TRUE
    
  } else {
    
    atoms.vn <- .findAtom(formula.c = formula.c,
                          atomic_weights.vn = atomic_weights.vn)
    
  }
  
  mass.n <- sum(atomic_weights.vn[names(atoms.vn)] * atoms.vn)
  
  atoms.vc <- as.character(atoms.vn)
  atoms.vc[atoms.vc == "1"] <- ""
  formula.c <- paste(paste0(names(atoms.vn),
                            atoms.vc),
                     collapse = "")
  
  if (protonate.l) {
    
    mass.n <- mass.n + atomic_weights.vn["proton"]
    formula.c <- paste0("[", formula.c, "+H]+")
    
  }
  
  # names(mass.n) <- formula.c
  
  return(list(mass = unname(mass.n),
              formula = formula.c))
  
}

.findAtom <- function(formula.c, atomic_weights.vn) {
  
  formula_split.vc <- unlist(strsplit(formula.c, ""))
  atoms.vn <- numeric()
  atom.i <- 0
  split.i <- 1
  while (split.i <= length(formula_split.vc)) {
    atom.i <- atom.i + 1
    if (split.i == length(formula_split.vc)) {
      atoms.vn <- c(atoms.vn, 1)
      names(atoms.vn)[atom.i] <- formula_split.vc[split.i]
      break
    } else if (split.i == length(formula_split.vc) - 1 &&
               formula_split.vc[split.i + 1] %in% letters) {
      atoms.vn <- c(atoms.vn, 1)
      names(atoms.vn)[atom.i] <- paste0(formula_split.vc[split.i],
                                        formula_split.vc[split.i + 1])
      break
    } else {
      atoms.vn <- c(atoms.vn, 0)
      if (formula_split.vc[split.i + 1] %in% letters) {
        names(atoms.vn)[atom.i] <- paste0(formula_split.vc[split.i],
                                          formula_split.vc[split.i + 1])
        split.i <- split.i + 1
      } else
        names(atoms.vn)[atom.i] <- formula_split.vc[split.i]
      
      num.i <- 1
      while ((split.i + num.i) <= length(formula_split.vc) &&
             !(formula_split.vc[split.i + num.i] %in% LETTERS))
        num.i <- num.i + 1
      num.i <- num.i - 1
      if (num.i == 0) {
        atoms.vn[atom.i] <- 1
      } else {
        atoms.vn[atom.i] <- as.numeric(paste(formula_split.vc[(split.i + 1):(split.i + num.i)],
                                             collapse = ""))
      }
      if ((split.i + num.i) == length(formula_split.vc)) {
        break
      } else
        split.i <- split.i + num.i + 1
    }
  }
  atoms.vn <- table(rep(names(atoms.vn), times = atoms.vn))
  if (any(!(names(atoms.vn) %in% names(atomic_weights.vn)))) {
    stop(paste(names(atoms.vn)[!(names(atoms.vn) %in% names(atomic_weights.vn))],
               collapse = ", "), " mass(es) cannot be provided currently",
         call. = FALSE)
  }
  atoms.vn <- atoms.vn[names(atomic_weights.vn)[names(atomic_weights.vn) %in% names(atoms.vn)]]
  
  atoms.vn
  
}

#' Isotope detection and validation of isotope clusters in PTR-TOF-MS peak table
#'
#' This function identify possible isotope cluster of C13,O17 and O18 atomes after alignment. It writes in the 
#' features data of the expressionSet (Biobase::fData()) the isotope m/z in the "isotope" column at the ligne 
#' of the more intenisf peak of the clusters.
#' @param eSet an expression set of PTR-TOF-MS data aligned 
#' @param ppm presision for the mass matching
#' @return an expresion with the column isotope added in teh features data
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
    iso <- isotopeMzMatching(mz[i], mz[(i+1):length(mz)],ppm)
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
  testCorPval<-vapply(as.character(groupIso)[-1], function(y) stats::cor.test(X[as.character(groupIso)[1],],
                                                                       X[y,],alternative = c("greater"))$p.value,1.1)
  testCor<- testCorPval < 0.01
  
  #ratio
  ratio<-X[as.character(groupIso)[-1],]/
    matrix(
      rep(X[as.character(groupIso)[1],,drop=F],2),
      nrow=length(as.character(groupIso)[-1]),byrow=TRUE)
  
  testRatio<- apply(ratio,1,median) < 0.5
  
  return(all(c(testCor,testRatio)))
}

