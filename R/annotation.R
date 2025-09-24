## annotateVOC (ExpressionSet) ----

#' @rdname annotation
#' @export
setMethod("annotateVOC", "ExpressionSet",
          function(x,
                   ionMassColname = "ion_mass", ppm = 50, 
                   ppm_formula=30,
                   top_formula=5,
                   extra_formula=TRUE,
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
            
            fdataDF <- Biobase::fData(x)
            
            ionMassColnameI <- which(colnames(fdataDF) == ionMassColname)
            
            
            if (length(ionMassColnameI) != 1)
              stop("No or multiple columns found in the fdataDF with the '", 
                   ionMassColname, "' name.",
                   call. = FALSE)
            
            ion_mass.vn <- fdataDF[, ionMassColnameI]
            
            annotateDF <- .annotate(ion_mass = ion_mass.vn,
                                    ppm = ppm,
                                    prefix = prefix,extra_formula=extra_formula,
                                    fields = fields,ppm_formula=ppm_formula,top_formula=top_formula)
            
            for (annotateC in colnames(annotateDF))
              fdataDF[, annotateC] <- annotateDF[, annotateC]
            
            Biobase::fData(x) <- fdataDF

            xIso<-try(findIsotope(eSet = x,ppm=ppm))
            if(!is.null(attr(xIso,"condition"))) return(x) else return(xIso)

            #x <- findIsotope(x)
            
            
          })

## annotateVOC (data.frame) ----

#' @rdname annotation
#' @export
setMethod("annotateVOC", "data.frame",
          function(x,
                   ionMassColname = "ion_mass", ppm = 50, 
                   ppm_formula=30,
                   top_formula=5,
                   extra_formula=TRUE,
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
            
            ionMassColnameI <- which(colnames(x) == ionMassColname)
            
            if (length(ionMassColnameI) != 1)
              stop("No or multiple columns found with the '", ionMassColname, 
                   "' name.",
                   call. = FALSE)
            
            ion_mass.vn <- x[, ionMassColnameI]
            
            annotateDF <- .annotate(ion_mass = ion_mass.vn,
                                    prefix = prefix,
                                    extra_formula=extra_formula,
                                    fields = fields,
                                    ppm_formula=ppm_formula,
                                    top_formula=top_formula)
            
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
                   ppm = 50,
                   ppm_formula=30,
                   top_formula=5,
                   extra_formula=TRUE,
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
                      extra_formula=extra_formula,
                      fields = fields,
                      ppm_formula=ppm_formula,
                      top_formula=top_formula)
            
            
          })

.annotate <- function(ion_mass,
                      ppm = 30,
                      ppm_formula=30,
                      top_formula=5,
                      extra_formula=TRUE,
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
  elements = c(Rdisop::initializeCHNOPS(),Rdisop::initializeElements(c("Cl","Br","F","I")))
  
  
  fielddbVl <- fields %in% colnames(vocdbDF)
  names(fielddbVl) <- fields
  
  if (sum(!fielddbVl) > 0)
    
    warnings("The following fields were not found in the vocDB database and 
             will be ignored:\n", paste(fields[!fielddbVl], collapse = ", "))
  
  fields <- fields[fielddbVl]
  
  annotateDF <- data.frame(row.names = as.character(ion_mass),
                           stringsAsFactors = FALSE)
  
  for (fieldC in fields)
    annotateDF[, paste0(prefix, fieldC)] <- character(nrow(annotateDF))
  
  for (i in seq_len(nrow(annotateDF))) {
    
    massN <- ion_mass[i]
    vocVi <- which(abs(massN - vocdbDF[, "ion_mass"]) < ppm * 1e-6 * massN)
    
    if (length(vocVi) > 0) {
      if(length(vocVi)>1)  vocVi <- which.min(abs(massN - vocdbDF[, "ion_mass"]))
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
        
        annotateDF[i, paste0(prefix, fieldC)] <- paste(sort(unique(vocFieldVc)), 
                                                       collapse = ", ")
        
      }
      
    } else {
        if(extra_formula){
            formula<-MassTools::calcMF(massN- 1.0072765, z=0,ppm=ppm_formula,top=top_formula,elements = elements) # M+H+
            if(!is.null(formula)){
                annotateDF[i, paste0(prefix, "ion_mass")]<- paste(round(formula$mz,4 )+ 1.0072765,collapse = "/")
                annotateDF[i, paste0(prefix, "ion_formula")]<- paste0("[",formula$MF,"+H]+",collapse = "/")
            }
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
  
  manualAdd <- utils::read.table(file = system.file("extdata/reference_tables/manualAdd.tsv",
                                                    package = "ptairMS"),
                                 check.names = FALSE,
                                 comment.char = "",
                                 header = TRUE,
                                 quote = "\"",
                                 sep = "\t",
                                 stringsAsFactors = FALSE)
  
  vocdbDF<-rbind(vocdbDF,manualAdd)
  
  vocdbDF<-vocdbDF[order(vocdbDF$ion_mass),]
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
  atoms.vn <- atoms.vn[names(atomic_weights.vn)[names(atomic_weights.vn) %in% 
                                                  names(atoms.vn)]]
  
  atoms.vn
  
}


## To label as isotope: 
## - difference un ion mass < ppm for atomic in "extdata/reference_tables/atomic_isotopes.tsv"
## - Correlation inter patient >0.8
## - ratio < 20 % or to the corresponding isotopic ratio if annotated formula
findIsotope<-function(eSet,ppm=50){
  X<-Biobase::exprs(eSet)
  fDATA<-Biobase::fData(eSet)
  mz<-as.numeric(row.names(X))
  
  #FIND AND VALIDE ISOTOPE GROUP
  for (i in seq_along(mz)){
      iso <- isotopeMzMatching(m = mz[i], mzSub = mz[(i+1):length(mz)],ppm = ppm)
      if(length(iso)){
        testIso<-validateGroup(groupIso = c(mz[i],iso),X = X,ppm)
        if(any(testIso)){
          fDATA[i,"isotope"]<- paste(iso[testIso],collapse = "/")
            #fDATA[as.character(iso[testIso]),"isotope"]<- mz[i]
        }
      }
    }
    Biobase::fData(eSet)<-fDATA
    return(eSet)
}


#find all possible isotope match in m/z dimension 
isotopeMzMatching<-function(m,mzSub,ppm,max=1){
  isotopes <- utils::read.table(system.file("extdata/reference_tables/atomic_isotopes.tsv",
                                            package = "ptairMS"),
                                header = TRUE,
                                quote = "\"",
                                sep = "\t",
                                stringsAsFactors = FALSE)
  
  
  isotopes<-isotopes[isotopes$abundance>0.001,]
  # anno<-annotateVOC(m,ppm=ppm,extra_formula = FALSE)
  # if(anno[,"vocDB_ion_formula"] != ""){
  #     element<-unique(isotopes$element)[vapply(unique(isotopes$element),function(e) grepl(pattern = e,x = anno[,"vocDB_ion_formula"]),FUN.VALUE = TRUE)]
  #    isotopes<-isotopes[isotopes$element %in% element,]
  # }
  diff<-lapply(split(isotopes, isotopes$element),function(x) {
    if(nrow(x)>1) x$mass[-1]-x$mass[1] else return(NULL) })
  diff<-Reduce(c,diff)
  Iso <- Reduce(c,lapply(diff,function(d) mzSub[which(abs(mzSub-(m+d))*10^6/m < ppm*1.3)]))
  iso<-unique(Iso)
  return(iso)
}

#validate isotope group
validateGroup<-function(groupIso,X,ppm){
  
  #correlation inter sample
  testCorPval<-vapply(as.character(groupIso)[-1], 
                      function(y) stats::cor.test(X[as.character(groupIso)[1],],
                                                  X[y,],alternative = c("greater"))$estimate,1.1)
  
  testCor<- testCorPval > 0.78
  
  #ratio
  ratio<-X[as.character(groupIso)[-1],]/
    matrix(
      rep(X[as.character(groupIso)[1],,drop=FALSE],
          length(as.character(groupIso)[-1])),
      nrow=length(as.character(groupIso)[-1]),byrow=TRUE)
  anno<-.annotate(groupIso[1],ppm=ppm,extra_formula = FALSE)
  # isotopes <- utils::read.table(system.file("extdata/reference_tables/atomic_isotopes.tsv",
  #                                          package = "ptairMS"),
  #                              header = TRUE,
  #                              quote = "\"",
  #                              sep = "\t",
  #                              stringsAsFactors = FALSE)
  
  isotopes <- utils::read.table(system.file("extdata/reference_tables/enviPat_isotopes.tsv",
                                package = "ptairMS"))
  #isotopes<-isotopes[!apply(isotopes,1,function(x) all(is.na(x[c(3,4,5)]))),]
  #data(isotopes)
  if( anno[,"vocDB_ion_formula"] != ""){
    formula <- anno[,"vocDB_ion_formula"]
    isoDistrib<-enviPat::isopattern(isotopes,chemforms = formula,threshold = 0.1,
                                    verbose = FALSE,charge=FALSE,emass=0.00054858)[[1]]
    testRatio<- vapply(seq_len(nrow(ratio)),function(x){
      index<-which(abs(isoDistrib[,"m/z"]-groupIso[x+1])*10^6/
                     round(groupIso[x+1]) < ppm*1.3)
      if(length(index>1)) index<- which.min(abs(isoDistrib[,"m/z"]-groupIso[x+1]))
      if(length(index)!=0){
         # testRatio<-abs(stats::median(ratio[x,],na.rm = TRUE)*100-isoDistrib[index,"abundance"])< 0.3*isoDistrib[index,"abundance"] # difference de ratio <30 %
        testRatio<- round(stats::median(ratio[x,],na.rm = TRUE)*100) <=
          max(ceiling(isoDistrib[index,"abundance"])+1,7)
    }else {
      testRatio<- stats::median(ratio[x,],na.rm = TRUE) < 0.2
    }
      },FUN.VALUE = TRUE) 
  }else {
    testRatio<- apply(ratio,1,function(x) stats::median(x,na.rm = TRUE)) < 0.2
  }

  return(apply(cbind(testCor,testRatio),1,all))
}


adduct<- function(eSet,ppm=50){
    X<- Biobase::exprs(eSet)
    Cormat<- cor(t(X))
    
    list<-list(NULL)
    for( i in seq(1,nrow(X))){
        list[[i]]<-rownames(X)[which(apply(X,1,function(x) cor(x,X[i,],method = "spearman")) >0.9)]
        
    }  
    
    names(list)<-rownames(X)
    groups<- list[which(unlist(lapply(list,length))>1)]
}
