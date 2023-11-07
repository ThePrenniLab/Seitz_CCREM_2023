#' do.findmain2
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param cmpd integer: vector defining compound numbers to annotated.  if NULL (default), all compounds
#' @param mode character: "positive" or "negative"
#' @param mzabs.error numeric: absolute mass deviation allowd, default = 0.01
#' @param ppm.error numeric: ppm mass error _added_ to mzabs.error, default = 10
#' @param ads character: vector of allowed adducts, i.e. c("[M+H]+"). if NULL, default positive mode values of H+, Na+, K+, and NH4+, as monomer, dimer, and trimer, are assigned. Negative mode include "[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-" as monomer, dimer, and trimer.
#' @param nls  character: vector of allowed neutral losses, i.e. c("[M+H-H2O]+").  if NULL, an extensive list derived from CAMERA's will be used. 
#' @param plot.findmain logical: should pdf polts be generated for evaluation? detfault = TRUE. PDF saved to working.directory/spectra
#' @param writeMat logical: should individual .mat files (for MSFinder) be generated in a 'mat' subdirectory in the 'spectra' folder? default = TRUE.
#' @details a partially annotated ramclustR object.  base structure is that of a standard R heirarchical clustering output, with additional slots described in ramclustR documentation (?ramclustR).  New slots added after using the interpretMSSpectrum functionality include those described below. 
#' @return    $M:  The inferred molecular weight of the compound giving rise to the each spectrum
#' @return    $M.ppm:  The ppm error of all the MS signals annotated, high error values should be considered 'red flags'. 
#' @return    $M.ann:  The annotated spectrum supporting the interpretation of M
#' @return    $use.findmain:  Logical vector indicating whether findmain scoring (TRUE) or ramclustR scoring (FALSE) was used to support inference of M.  By default, findmain scoring is used.  When ramclustR scoring differs from findmain scoring, the scoring metric which predicts higher M is selected. 
#' @return    $M.ramclustr:  M selected using ramclustR scoring
#' @return    $M.ppm.ramclustr:  ppm error of M selected using ramclustR scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto. 
#' @return    $M.ann.ramclustr:  annotated spectrum supporting M using ramclustR scoring
#' @return    $M.nann.ramclustr:  number of masses annotated using ramclustR scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto. 
#' @return    $M.space.ramclustr:  the 'space' of scores between the best and second best ramclustR scores. Calculated as a ratio. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto. 
#' @return    $M.findmain:  M selected using findmain scoring
#' @return    $M.ppm.findmain:  ppm error of M selected using findmain scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto. 
#' @return    $M.ann.findmain:  annotated spectrum supporting M using findmain scoring
#' @return    $M.nann.findmain:  number of masses annotated using findmain scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto. 
#' @return    $M.space.findmain:  the 'space' of scores between the best and second best findmain scores. Calculated as a ratio. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto. 
#' @references Jaeger C, ... Lisec J. Compound annotation in liquid chromatography/high-resolution mass spectrometry based metabolomics: robust adduct ion determination as a prerequisite to structure prediction in electrospray ionization mass spectra. Rapid Commun Mass Spectrom. 2017 Aug 15;31(15):1261-1266. doi: 10.1002/rcm.7905. PubMed PMID: 28499062.
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept findMain
#' @concept interpretMSSpectrum
#' @concept xcms
#' @author Corey Broeckling
#' @export 


do.findmain2 <- function (
    ramclustObj = NULL, 
    cmpd = NULL, 
    mode = "positive", 
    mzabs.error = 0.005, 
    ppm.error = 10,
    mainpkthr = 0.15,
    ads = NULL, 
    nls = NULL,
    plot.findmain = TRUE, 
    writeMat = TRUE,
    writeMS = FALSE,
    writeMGF = TRUE
) 
{
  
  if (!requireNamespace("InterpretMSSpectrum", quietly = TRUE)) {
    stop("The use of this function requires package 'InterpretMSSpectrum'.")
  }
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  
  if (is.null(ads)) {
    if (grepl("p", mode)) {
      ads <- c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", 
               "[2M+H]+", "[2M+Na]+", "[2M+K]+", "[2M+NH4]+",
               "[3M+H]+", "[3M+Na]+", "[3M+K]+", "[3M+NH4]+")
    }
    if (grepl("n", mode)) {
      ads <- c("[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-", 
               "[2M-H]-", "[2M+Na-2H]-", "[2M+K-2H]-", "[2M+CH2O2- H]-",
               "[3M-H]-", "[3M+Na-2H]-", "[3M+K-2H]-", "[3M+CH2O2- H]-")
    }
    if (is.null(ads)) {
      stop("please define adducts using 'ads' or set mode to either'positive' or 'negative'")
    }
  }
  if (is.null(nls)) {
    if (grepl("p", mode)) {
      nls <- c("[M+H-COCH2]+", "[M+H-C2H3NO]+", "[M+H-H2O]+", 
               "[M+H-NH3]+", "[M+H-HCOOH]-", "[M+H-C6H12O6]+", "[M+H-C5H10O5]+", 
               "[M+H-C12H22O11]+")
    }
    if (grepl("n", mode)) {
      nls <- c("[M-H-NH3]-", "[M-H-H2O]-", "[M-H-COCH2]-", 
               "[M-H-CO2]-", "[M-H-NH3-CO2]-", "[M-H-HCOOH]-", "[M-H-C6H12O6]-", 
               "[M-H-C5H10O5]-", "[M-H-C12H22O11]-")
    }
    if (is.null(nls)) {
      stop("please define neutral losses using 'nls' or set mode to either'positive' or 'negative'")
    }
  }
  
  
  findmain <- as.list(rep(NA, length(ramclustObj$cmpd)))
  names(findmain) <-  ramclustObj$cmpd
  
  if(is.null(cmpd)) {
    cmpd <- (1:max(ramclustObj$featclus))
  }
  
  ramclustObj$ms1.spectrum <- as.list(rep(NA, length(ramclustObj$ann)))
  ramclustObj$ms2.spectrum <- as.list(rep(NA, length(ramclustObj$ann)))
  ramclustObj$ms2.precursor <- rep(NA, length(ramclustObj$ann))
  ramclustObj$ms2.precursor.iontype <- rep(NA, length(ramclustObj$ann))
  
  for (cl in cmpd) {
    s <- data.frame(
      mz = ramclustObj$fm[which(ramclustObj$featclus == cl)], 
      int = ramclustObj$msint[which(ramclustObj$featclus == cl)])
    
    
    s2 <- data.frame(
      mz = ramclustObj$fm[which(ramclustObj$featclus == cl)], 
      int = ramclustObj$msmsint[which(ramclustObj$featclus == cl)])
    
    s <- s[order(s$mz),]
    s2 <- s2[order(s2$mz),]
    
    ramclustObj$ms1.spectrum[[cl]] <- s
    ramclustObj$ms2.spectrum[[cl]] <- s2
    
    out <- InterpretMSSpectrum::findMAIN(
      s, 
      rules = c(ads), 
      adducthyp = ads[grep("[M",ads, fixed = TRUE)], 
      ionmode = mode, 
      mzabs = mzabs.error, 
      ppm = ppm.error,
      mainpkthr = mainpkthr
    )
    summarytable <- summary(out)
    fm.out <- summarytable
    
    out <- InterpretMSSpectrum::findMAIN(
      s, 
      adductmz = NULL, 
      ionmode = mode, 
      rules = c(ads, nls), 
      adducthyp = ads[grep("[M", ads, fixed = TRUE)], 
      # ms2spec = s2, 
      mzabs = mzabs.error,
      ppm = ppm.error, 
      mainpkthr = mainpkthr, 
      collapseResults = FALSE)
    summarytable <- summary(out)
    rc.out <- data.frame('index' = 1:nrow(summarytable),  summarytable)
    
    sum.score <- merge(rc.out, fm.out, by = 'neutral_mass')
    sum.score <- data.frame(sum.score, 'total.score' = (sum.score$total_score.x + sum.score$total_score.y)/2)
    if(nrow(sum.score) == 0) {
      sum.score <- rc.out
    }
    sum.score <- sum.score[order(sum.score$total.score, decreasing = TRUE),]
    sum.score <- sum.score[which(sum.score$total.score >= (max(sum.score$total.score * 0.90))),]
    
    
    out.list <- as.list(rep(NA, 2))
    names(out.list) <- c("summary", "details")
    out.list[[1]] <- sum.score
    out.list[[2]] <- out[sum.score$index]
    
    if (100 * round(cl/100, digits = 0) == cl) {
      cat(cl, "of", max(ramclustObj$featclus), "\n")
    }
    
    findmain[[cl]] <- out.list
  }
  
  ramclustObj$findmain <- findmain
  
  # generate master findmain summary table
  sum.table <- data.frame(
    'cmpd' = vector(mode = 'character', length = 0),
    'hypothesis' = vector(mode = 'character', length = 0),
    'findmain.score' = vector(mode = 'numeric', length = 0)
  )
  for (cl in cmpd) {
    tmp <- findmain[[cl]]$summary
    sub.sum.table <- data.frame(
      'cmpd' = rep(ramclustObj$cmpd[cl], nrow(tmp)),
      'hypothesis' = paste0(ramclustObj$cmpd[cl], ".", formatC((1:nrow(tmp)), width = 2, flag = 0)),
      'findmain.score' = tmp$total.score
    )
    sum.table <- rbind(sum.table, sub.sum.table)
  }
  
  ramclustObj$findmain.summary <- sum.table
  
  ###############################################################################################################
  ## insert code here
  ## create list out output spectra first
  ## include new interpretMSspectrum findmain
  ## optionally convert multiple charged species to singly charged
  ## store spectra in ramclustObj
  ## writeMat, writeMS, writeMGF then gets much easier
  ###############################################################################################################
  
  if (writeMat) {
    if (!dir.exists("spectra")) {
      dir.create("spectra")
    }
    dir.create("spectra/mat")  # cmpd <- 1:50
    for (cl in cmpd) {
      for(i in 1:length(findmain[[cl]]$details)) {
        ms <- findmain[[cl]]$details[[i]]
        ms$int <- round(ms$int, digits = 3)
        prcr <- which(ms[, "adduct"] %in% ads)
        prcr <- prcr[which.max(ms[prcr, "int"])]
        prcmz <- ms[prcr, "mz"]
        prctype <- ms[prcr, "adduct"]
        out <- paste("NAME: ", ramclustObj$cmpd[cl], "\n", 
                     "RETENTIONTIME: ", round(ramclustObj$clrt[cl], 
                                              digits = 2), "\n", "PRECURSORMZ: ", prcmz, 
                     "\n", "PRECURSORTYPE: ", prctype, "\n", "IONTYPE: ", 
                     mode, "\n", "SPECTRUMTYPE: Centroid", "\n", 
                     if ((!is.null(ramclustObj$msmsint))) {
                       paste("COLLISIONENERGY: ", as.character(ramclustObj$ExpDes[[2]][which(row.names(ramclustObj$ExpDes[[2]]) == 
                                                                                               "CE2"), 1]), "\n", sep = "")
                     }, "MSTYPE: ", "MS1", "\n", "Num Peaks: ", nrow(ms), 
                     "\n", sep = "")
        for (j in 1:nrow(ms)) {
          
          out <- paste(out, ms[j, 1], " ", ms[j, 2], "\n", 
                       sep = "")
        }
        if (!is.null(ramclustObj$msmsint)) {
          feats <- which(ramclustObj$featclus == cl)
          if (length(feats) > 0) {
            msms <- data.frame(mz = ramclustObj$fm[feats], 
                               int = ramclustObj$msmsint[feats])
            msms <- msms[which(msms[, "mz"] <= (prcmz + 
                                                  3)), , drop = FALSE]
            msms <- msms[order(msms[, "int"], decreasing = TRUE), 
                         , drop = FALSE]
            msms$int <- round(1000*(msms$int/max(msms$int)), digits = 3)
            if (nrow(msms) > 0) {
              out <- paste(out, '\n', "MSTYPE:", " MS2", "\n", 
                           "Num Peaks: ", nrow(msms), "\n", sep = "")
              for (k in 1:nrow(msms)) {
                out <- paste(out, msms[k, 1], " ", msms[k, 
                                                        2], "\n", sep = "")
              }
            }
          }
        }
        else {
          do <- which(ramclustObj$featclus == cl)
          if (length(do) > 0) {
            msms <- data.frame(mz = ramclustObj$fm[do], 
                               int = ramclustObj$msint[do])
            msms <- msms[which(msms[, "mz"] <= (prcmz + 
                                                  3)), , drop = FALSE]
            msms <- msms[order(msms[, "int"], decreasing = TRUE), 
                         , drop = FALSE]
            msms$int <- round(1000*(msms$int/max(msms$int)), digits = 3)
            if (nrow(msms) > 0) {
              out <- paste(out, '\n', "MSTYPE:", " MS2", "\n", 
                           "Num Peaks: ", nrow(msms), "\n", sep = "")
              for (k in 1:nrow(msms)) {
                out <- paste(out, msms[k, 1], " ", msms[k, 
                                                        2], "\n", sep = "")
              }
            }
          }
        }
        write(out, file = paste0("spectra/mat/", ramclustObj$cmpd[cl], ".", formatC(i, width = 2, flag = 0),
                                 ".mat"))
      }
    }
  }
  
  
  ramclustObj$history$do.findmain <- paste(
    " Molecular weight was inferred from in-source spectra (Broeckling 2016) using the do.findmain function, which calls the ", 
    "interpretMSSpectrum package (Jaeger 2016). ", 
    "Parameters for do.findmain were set to: ", 
    "mode = ", mode, ", mzabs.error = ", mzabs.error, ", ppm.error = ", 
    ppm.error, ", ads = ", paste(ads, collapse = " "), ", nls = ", 
    paste(nls, collapse = " "), ".", sep = "")
  cat("finished", "\n")
  return(ramclustObj)
  
  
}


#' annotate.msfinder
#'
#' After running MSFinder on .mat or .msp files, import the formulas that were predicted and their scores 
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param mat.dir optional path to .mat directory
#' @param msp.dir optional path to .msp directory
#' @details this function imports the output from the MSFinder program to support annotation of the ramclustR object
#' @return new slot at $msfinder.formula.details
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export

annotate.msfinder <- function (ramclustObj = NULL, 
                               mat.dir = NULL,
                               priority.db = NULL,
                               priority.inchikey = NULL,
                               priority.db.factor = 0.9,
                               priority.inchikey.factor = 0.9
) 
{
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  home.dir <- getwd()
  
  if(!is.null(priority.db) & !is.null(priority.inchikey)) {
    warning("both inchikey and database priority set - ensure they are independent. ", '\n')
  }
  
  r <- grep("msfinder.formula", names(ramclustObj))
  if (length(r) > 0) {
    warning("removed previously assigned MSFinder formulas and structures", 
            "\n")
    ramclustObj <- ramclustObj[-r]
    r <- grep("msfinder.structure", names(ramclustObj))
    if(length(r)>0) {
      ramclustObj <- ramclustObj[-r]
    }
    rm(r)
  }
  if (is.null(mat.dir)) {
    mat.dir = paste0(getwd(), "/spectra/mat")
  }
  
  if(is.null(ramclustObj$history)) {
    ramclustObj$history <- ""
  }
  
  findmain.summary <- ramclustObj$findmain.summary
  
  ## determine most recent batchparam file
  msf.form.param.file <- list.files(mat.dir, pattern = "batchparam", full.names = TRUE)
  if(length(msf.form.param.file) > 0) {
    msf.form.param.file <- msf.form.param.file[which.max(file.info(msf.form.param.file)$ctime)]
    msf.params <- readLines(msf.form.param.file)
  }
  
  
  ## find all formula output files: 
  form.files <- list.files(mat.dir, pattern = '.fgt', recursive = FALSE, full.names = TRUE)
  
  ## find all structure output files: 
  struc.files <- list.files(mat.dir, pattern = '.sfd', recursive = TRUE, full.names = TRUE)
  
  ## separate out spectal search results
  is.spec.db <- grep("Spectral", struc.files)
  spec.files <- struc.files[is.spec.db]
  struc.files <- struc.files[-is.spec.db]

  
  ###############################
  ## read in formula results
  cat(" -- importing formula results", '\n')
  hypothesis <- gsub(".fgt", "", basename(form.files), ignore.case = TRUE)
  
  form.results <- data.frame(
    'hypothesis' = vector(mode = 'character', length = 0),
    'formula' = vector(mode = 'character', length = 0),
    'formula.score'   = vector(mode =  'numeric', length = 0)
  )
  
  for(i in 1:length(form.files)) {
    tmp <- readLines(form.files[i])
    tmp.names <- grep("NAME: ", tmp)
    if(length(tmp.names) > 0) {
      tmp.scores <- grep("TOTALSCORE:", tmp)
      if((length(tmp.names)>0) & (length(tmp.names) == length(tmp.scores))) {
        tmp.out <- data.frame(
          'hypothesis' = rep(hypothesis[i], length(tmp.names)),
          'formula' = gsub("NAME: ", "", tmp[tmp.names]),
          'formula.score'   = as.numeric(gsub("TOTALSCORE: ", "", tmp[tmp.scores]))
        )
      }
      tmp.rm <- grep("Spectral", tmp.out$formula)
      if(length(tmp.rm)>0) {
        tmp.out <- tmp.out[-tmp.rm,]
      }
      if(nrow(tmp.out)>0) {
        form.results <- rbind(form.results, tmp.out)
      }
    }
    
  }
  
  ## read in structure results
  cat(" -- importing structure results", '\n')
  struc.path <- lapply(1:length(struc.files), FUN = function(x) {
    unlist(strsplit(struc.files[x], "/", fixed = TRUE))
  })
  struc.path <- t(data.frame(struc.path))
  hypothesis <- as.character(struc.path[,(ncol(struc.path)-1)])
  hypothesis.n <- as.numeric(sapply(1:length(hypothesis), FUN = function(x) unlist(strsplit(hypothesis[x], '.', fixed = TRUE))[2]))
  formula    <- as.character(gsub('.sfd', '', struc.path[,ncol(struc.path)]))
  cmpd       <- sapply(1:length(hypothesis), FUN = function(x) {
    strsplit(hypothesis[x], ".", fixed = TRUE)[[1]][1]
  }
  )
  
  struc.results <- data.frame(
    'cmpd' = vector(mode = 'character', length = 0),
    'hypothesis' = vector(mode = 'character', length = 0),
    'retention.time' = vector(mode = 'character', length = 0),
    'assigned.M' = vector(mode = 'character', length = 0),
    'formula'    = vector(mode = 'character', length = 0),
    'compound.name' = vector(mode = 'character', length = 0),
    'inchikey'   = vector(mode = 'character', length = 0),
    'dbs'        = vector(mode = 'character', length = 0),
    'direct.parent.class'   = vector(mode = "character", length = 0),
    'chemontid'  = vector(mode = 'character', length = 0),
    'structure.score'  = vector(mode =  'numeric', length = 0)
  )
  cmpd.rts <- round(ramclustObj$clrt, digits = 1)
  names(cmpd.rts) <- ramclustObj$cmpd
  for(i in 1:length(struc.files)) {
    tmp <- readLines(struc.files[i])
    tmp.names <- grep("NAME: ", tmp)
    if(length(tmp.names) > 0) {
      tmp.scores <- grep("TotalScore:", tmp)
      tmp.inchikey <- grep("INCHIKEY: ", tmp)
      tmp.dbs     <- grep("RESOURCES: ", tmp)
      tmp.subclass<- grep("Ontology: ", tmp)
      tmp.chemontid<- grep("OntologyID: ", tmp)
      same.length <- all.equal(
        length(tmp.names),
        length(tmp.scores),
        length(tmp.inchikey),
        length(tmp.dbs),
        length(tmp.subclass),
        length(tmp.chemontid)
      )
      if((length(tmp.names)>0) & same.length) {
        tmp.out <- data.frame(
          'cmpd'                  = rep(cmpd[i], length(tmp.names)),
          'hypothesis'            = rep(hypothesis[i], length(tmp.names)),
          'retention.time'        = as.vector(rep(cmpd.rts[cmpd[i]], length(tmp.names))),
          'assigned.M'            = rep(ramclustObj$findmain[[cmpd[i]]]$summary$neutral_mass[hypothesis.n[i]], length(tmp.names)),
          'formula'               = rep(formula[i], length(tmp.names)),
          'compound.name'         = gsub("NAME: ", "", tmp[tmp.names]),
          'inchikey'              = gsub("INCHIKEY: ", "", tmp[tmp.inchikey]),
          'dbs'                   = gsub("RESOURCES: ", "", tmp[tmp.dbs]),
          'direct.parent.class'   = gsub("Ontology: ", "", tmp[tmp.subclass]),
          'chemontid'             = gsub("OntologyID: ", "", tmp[tmp.chemontid]),
          'structure.score'       = as.numeric(gsub("TotalScore: ", "", tmp[tmp.scores]))
        )
      }
      tmp.rm <- grep("Spectral", tmp.out$formula)
      if(length(tmp.rm)>0) {
        tmp.out <- tmp.out[-tmp.rm,]
      }
      if(nrow(tmp.out)>0) {
        struc.results <- rbind(struc.results, tmp.out)
      }
    }
  }
  
  ## bring formula score into structure table 
  cat(" -- associating formula scores to structures", '\n')
  struct.names <- paste(struc.results$hypothesis, struc.results$formula)
  form.names   <- paste(form.results$hypothesis, form.results$formula)
  struc.to.form <- match(
    struct.names,
    form.names
  )
  struc.results$formula.score <- form.results$formula.score[struc.to.form]
  
  ## bring findmain score into structure table 
  struc.to.findmain <- match(
    struc.results$hypothesis ,
    findmain.summary$hypothesis
  )
  struc.results$findmain.score <- findmain.summary$findmain.score[struc.to.findmain]
  
  
  ## calculate a total score for each annotation hypothesis
  cat(" -- assigning total score", '\n')
  struc.results$total.score <- round(
    2*struc.results$structure.score * struc.results$formula.score * struc.results$findmain.score, 
    digits = 2
    )
  
  
  ## read in spectral match results
  cat(" -- importing spectral search results", '\n')
  spec.path <- lapply(1:length(spec.files), FUN = function(x) {
    unlist(strsplit(spec.files[x], "/", fixed = TRUE))
  })
  spec.path <- t(data.frame(spec.path))
  hypothesis <- as.character(spec.path[,(ncol(spec.path)-1)])
  cmpd       <- sapply(1:length(hypothesis), FUN = function(x) {
    strsplit(hypothesis[x], ".", fixed = TRUE)[[1]][1]
  }
  )
  
  spec.results <- data.frame(
    'cmpd' = vector(mode = 'character', length = 0),
    'hypothesis' = vector(mode = 'character', length = 0),
    'compound.name' = vector(mode = 'character', length = 0),
    'inchikey'   = vector(mode = 'character', length = 0),
    'dbs'        = vector(mode = 'character', length = 0),
    'spectral.match.score'  = vector(mode =  'numeric', length = 0),
    'direct.parent.class'   = vector(mode = "character", length = 0),
    'chemontid'  = vector(mode = 'character', length = 0)
  )
  
  for(i in 1:length(spec.files)) {
    tmp <- readLines(spec.files[i])
    tmp.names <- grep("NAME: ", tmp)
    if(length(tmp.names) > 0) {
      tmp.scores <- grep("TotalScore:", tmp)
      tmp.inchikey <- grep("INCHIKEY: ", tmp)
      tmp.dbs     <- grep("RESOURCES: ", tmp)
      tmp.subclass<- grep("Ontology: ", tmp)
      tmp.chemontid<- grep("OntologyID: ", tmp)
      same.length <- all.equal(
        length(tmp.names),
        length(tmp.scores),
        length(tmp.inchikey),
        length(tmp.dbs),
        length(tmp.subclass),
        length(tmp.chemontid)
      )
      if((length(tmp.names)>0) & same.length) {
        tmp.out <- data.frame(
          'cmpd'       = rep(cmpd[i], length(tmp.names)),
          'hypothesis' = rep(hypothesis[i], length(tmp.names)),
          'compound.name'   = gsub("NAME: ", "", tmp[tmp.names]),
          'inchikey'   = gsub("INCHIKEY: ", "", tmp[tmp.inchikey]),
          'dbs'        = gsub("RESOURCES: ", "", tmp[tmp.dbs]),
          'spectral.match.score'      = as.numeric(gsub("TotalScore: ", "", tmp[tmp.scores])),
          'direct.parent.class'   = gsub("Ontology: ", "", tmp[tmp.subclass]),
          'chemontid'  = gsub("OntologyID: ", "", tmp[tmp.chemontid])
        )
      }
      if(nrow(tmp.out)>0) {
        spec.results <- rbind(spec.results, tmp.out)
      }
    }
  }
  
  struc.results <- data.frame(
    'assigned' = rep(FALSE, nrow(struc.results)),
    struc.results,
    'spectral.match.score' = rep(NA, nrow(struc.results)),
    'spectral.match.name'  = rep(NA, nrow(struc.results)),
    'spectral.match.inchikey' = rep(NA, nrow(struc.results))
  )
  
  ## update spectral match results table
  if(length(spec.results) > 0) {
    cmpds <- unique(spec.results$cmpd)
    if(length(cmpds) > 0) {
      for(i in 1:length(cmpds)) {
        tmp <- spec.results[which(spec.results$cmpd == cmpds[i]),]
        tmp <- tmp[which.max(tmp$spectral.match.score),]
        cmpd.match <- which(struc.results$cmpd == tmp$cmpd[1])
        if(length(cmpd.match)>0) {
          struc.results[cmpd.match, "spectral.match.score"] <- as.numeric(tmp$spectral.match.score[1])
          struc.results[cmpd.match, "spectral.match.name"] <- as.character(tmp$compound.name[1])
          struc.results[cmpd.match, "spectral.match.inchikey"] <- as.character(tmp$inchikey[1])
        }
      }
    }
  }
  
  ## assign database priority score, if applicable
  
  if(!is.null(priority.db)) {
    cat(" -- assinging database priority score", '\n')
    priority.factor.v <- rep(priority.db.factor, nrow(struc.results))
    for(i in 1:length(priority.db)) {
      db.match <- grep(priority.db[i], struc.results$dbs, ignore.case = TRUE)
      priority.factor.v[db.match] <- 1
    }
    struc.results$db.priority.factor <- priority.factor.v
    struc.results$total.score <- struc.results$total.score * priority.factor.v
  }
  
  
  ## assign inchikey priority score, if applicable
  if(!is.null(priority.inchikey)) {
    cat(" -- assinging inchikey priority score", '\n')
    priority.factor.v <- rep(priority.inchikey.factor, nrow(struc.results))
    short.inchikey.priority <- sapply(priority.inchikey, FUN = function(x) {
      unlist(strsplit(x, "-"))[1]
    })
    short.inchikey.struc <- sapply(struc.results$inchikey, FUN = function(x) {
      unlist(strsplit(x, "-"))[1]
    })
    inchi.match <- short.inchikey.struc %in% short.inchikey.priority
    priority.factor.v[inchi.match] <- 1
    struc.results$inchikey.priority.factor <- priority.factor.v
    struc.results$total.score <- struc.results$total.score * priority.factor.v
  }
  
  ## remove any results which have no findmain score (and therefore no total score)
  struc.results <- struc.results[which(!is.na(struc.results$findmain.score)),]
  
  ## reset annotations
  ramclustObj$ann  <- ramclustObj$cmpd
  
  ## update structure table
  ramclustObj$use.spectral.match <- rep(FALSE, length(ramclustObj$cmpd))
  
  cat(" -- annotating", '\n')
  for(i in 1:length(ramclustObj$ann)) {
    if(ramclustObj$ann[i] != ramclustObj$cmpd[i]) next
    cmpd.match <- which(struc.results$cmpd == ramclustObj$cmpd[i])
    if(length(cmpd.match)==0) next
    ann.sub <- struc.results[cmpd.match,]
    
    use.spectral.match <- any(!is.na(ann.sub$spectral.match.score))
    if(use.spectral.match) {
      ramclustObj$use.spectral.match[i] <- TRUE
      sp.m.index <- which.max(ann.sub$spectral.match.score)
      struc.results$assigned[cmpd.match[sp.m.index]] <- TRUE
      ramclustObj$ann[i] <-  ann.sub$spectral.match.name[sp.m.index]
      ramclustObj$inchikey <- ann.sub$spectral.match.inchikey[sp.m.index]
      ramclustObj$spectral.match.score[i] <- ann.sub$spectral.match.score[sp.m.index]
    } else {
      best.score <- which.max(ann.sub$total.score)
      sel.hyp <- as.numeric(gsub(paste0(ramclustObj$cmpd[i], "."), "", ann.sub$hypothesis[best.score]))
      fm <- ramclustObj$findmain[[i]]$details[[sel.hyp]]
      fm <- fm[order(fm[,1]),]
      fm.sum <- ramclustObj$findmain[[i]]$summary[sel.hyp,]
      struc.results$assigned[cmpd.match[best.score]] <- TRUE
      ramclustObj$ms1.spectrum[[i]] <- fm
      ramclustObj$ms2.spectrum[[i]] <- data.frame(
        ramclustObj$ms2.spectrum[[i]][order(ramclustObj$ms2.spectrum[[i]][,1]),],
        fm[3:ncol(fm)]
      )
      ramclustObj$ann[i] <- ann.sub$compound.name[best.score]
      ramclustObj$M[i] <- fm.sum$neutral_mass[1]
      ramclustObj$formula[i] <- ann.sub$formula[best.score]
      ramclustObj$inchikey[i] <- ann.sub$inchikey[best.score]
      ramclustObj$findmain.score[i] <- ann.sub$findmain.score[best.score]
      ramclustObj$formula.score[i] <- ann.sub$formula.score[best.score]
      ramclustObj$structure.score[i] <- ann.sub$structure.score[best.score]
      ramclustObj$total.score[i] <- ann.sub$total.score[best.score]
    }
  }
  
  ramclustObj$annotations.full <- struc.results
  ramclustObj$annotations.selected <- struc.results[which(struc.results$assigned),2:ncol(struc.results)]
  
  ramclustObj$history$msfinder <- paste(
    "MSFinder (Tsugawa 2016) was used for spectral matching,",
    "formula inference, and tentative structure assignment.",
    "Results were imported into the RAMClustR object.", 
    "A total score was calculated based on the product scores from the findmain function",
    "and the MSfinder formula and structure scores.", 
    "A total of", length(unique(ramclustObj$findmain.summary$hypothesis)), "annotation hypotheses were tested", 
    "for", length(ramclustObj$ann), "compounds.", 
    "A complete spreadsheet of all annotation hypothesis and scores can be found in the 'spectra/all.annotations.csv' file,",
    "and a subset of only those selected for annotation can be found in the 'spectra/assigned.annotations.csv' file.",
    "Spectra matches took precedence over computational inference based annotations.")
  
  
  
  if(!is.null(priority.db)) {
    ramclustObj$history$msfinder <-paste(
      as.character(ramclustObj$history$msfinder),
      "The following database(s) were assigned as 'priority': ", paste0(paste(priority.db, collapse = ', '), "."),
      "The database priority.factor was set to", priority.db.factor, "to decrease scores for compounds which tailed to match priority database(s)."
    )
  }
  
  if(!is.null(priority.inchikey)) {
    ramclustObj$history$msfinder <-paste(
      as.character(ramclustObj$history$msfinder),
      "The a list of",  length(priority.inchikey), "inchikeys was provided.", 
      "The inchikey priority.factor was set to", priority.inchikey.factor, "to decrease scores for compounds with non-matching inchikey(s)."
    )
  }
  ramclustObj$history$msfinder <-paste(
    as.character(ramclustObj$history$msfinder),
    "The highest total score was selected for each compound, considering all hypotheses."
  )
  
  ## write annotaion tables to 'spectra' directory
  out.dir <- gsub(basename(mat.dir), "",  mat.dir)
  write.csv(struc.results, file = paste0(out.dir, "all.annotations.csv"), row.names = FALSE)
  write.csv(ramclustObj$annotations.selected, file = paste0(out.dir, "assigned.annotations.csv"), row.names = FALSE)
  return(ramclustObj)
}

rc.run.msfinder <- function(
    msfinder.dir = "",
    input.dir  =  "",
    param.file = ""
) {
  
  if(substring(msfinder.dir, nchar(msfinder.dir)) != "/") {
    msfinder.dir <- paste0(msfinder.dir, "/")
  }
  
  if(substring(input.dir, nchar(input.dir)) != "/") {
    input.dir <- paste0(input.dir, "/")
  }
  
  cat(" ----- pre_mssearch", '\n')
  #  MsfinderConsoleApp.exe mssearch -i .\Data\ -o .\ -m .\MsfinderConsoleApp-Param.txt
  system(
    paste0(
      msfinder.dir, "MsfinderConsoleApp",
      " mssearch ",
      " -i ", input.dir,
      " -o ", input.dir,
      " -m ", param.file
    )
  )
  
  cat(" ----- post_mssearch", '\n')
  
  ## rename directories containing spectral search results to keep them from being overwritten
  spec.results <- list.files(
    paste0(input.dir),
    pattern = "MsSearch result",
    recursive = FALSE)
  versions <- as.numeric(gsub(".txt", "", gsub("MsSearch result-", "", spec.results)))
  spec.results <- spec.results[which.max(versions)]
  res <- read.delim(paste0(input.dir, spec.results))
  from.name <- unique(res$File.name)
  to.name <- paste0(from.name, ".mssearch")
  
  for(i in 1:length(from.name)) {
    dir.create(paste0(input.dir,  to.name[i]))
    files.to.move <- list.files(paste0(input.dir, from.name[i]), recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
    file.copy(files.to.move, to = paste0(input.dir,  to.name[i]), recursive = TRUE)
  }
  cat(" ----- post_modify.directory.names", '\n')
  
  system(
    paste0(
      msfinder.dir, "MsfinderConsoleApp",
      " predict ",
      " -i ", input.dir,
      " -o ", input.dir,
      " -m ", param.file
    )
  )
  
  cat (" ----- post_predict", '\n')
}


rc.run.msfinder.gc <- function(
    msfinder.dir = "",
    input.dir  =  "",
    param.file = ""
) {
  
  if(substring(msfinder.dir, nchar(msfinder.dir)) != "/") {
    msfinder.dir <- paste0(msfinder.dir, "/")
  }
  
  if(substring(input.dir, nchar(input.dir)) != "/") {
    input.dir <- paste0(input.dir, "/")
  }
  
  cat(" ----- pre_mssearch", '\n')
  #  MsfinderConsoleApp.exe mssearch -i .\Data\ -o .\ -m .\MsfinderConsoleApp-Param.txt
  system(
    paste0(
      msfinder.dir, "MsfinderConsoleApp",
      " mssearch ",
      " -i ", input.dir,
      " -o ", input.dir,
      " -m ", param.file
    )
  )
  
  cat(" ----- post_mssearch", '\n')
  
}


######################################################################
## get classyfire heirarchy for a given vector of classyfire chemonids

get.classyfire.local <- function(
    ramclustObj = NULL,
    chemontid = NULL,
    chemont.obo.location = NULL,
    circle.plot = TRUE
) {
  
  if(!is.null(ramclustObj) & is.null(chemontid)) {
    chemontid <- ramclustObj$annotations.full$chemontid
  }
  
  if(is.null(chemontid)) {stop("no chemontids submitted", '\n')}
  if(is.null(chemont.obo.location)) {stop("no chemont.obo file specified", '\n')}
  
  classyfire <- ontologyIndex::get_ontology(
    file = chemont.obo.location,
    propagate_relationships = "is_a",
    extract_tags = "minimal",
    merge_equivalent_terms = TRUE
  )
  
  cmpd.classyfire <- data.frame(
    'superclass'     = rep(NA, 0),
    'class'          = rep(NA, 0),
    'subclass'       = rep(NA, 0),
    'direct.parent'  = rep(NA, 0)
  )
  
  for(i in 1:length(chemontid)) {
    
    if(
      nchar(ramclustObj$annotations.full[i,'chemontid'])==0 | ramclustObj$annotations.full[i,'chemontid'] == "NA"
      ) {
      out <- data.frame(
        'superclass'     = NA,
        'class'          = NA,
        'subclass'       = NA,
        'direct.parent'  = NA
      )
      cmpd.classyfire <- rbind(
        cmpd.classyfire,
        out
      )
    } else {
      classyfication <- as.vector(classyfire$name[unlist(classyfire$ancestors[chemontid[i]])][3:5])
      out <- data.frame(
        'superclass'     = classyfication[1],
        'class'          = classyfication[2],
        'subclass'       = classyfication[3],
        'direct.parent'  = classyfire$name[chemontid[i]]
      )
      cmpd.classyfire <- rbind(
        cmpd.classyfire,
        out
      )
      
    }
    if(nrow(cmpd.classyfire)>i) {stop("too many rows. On: ", i, '\n')}
  }
  
  cmpd.classyfire[which(is.na(cmpd.classyfire), arr.ind = TRUE)] <- 'unassigned'
  
  if(is.null(ramclustObj)) {
    return(cmpd.classyfire)} else {
      ramclustObj$annotations.full <- 
        data.frame(
          ramclustObj$annotations.full,
          cmpd.classyfire
        )
      ramclustObj$annotations.selected <- ramclustObj$annotations.full[which(ramclustObj$annotations.full$assigned),2:ncol(ramclustObj$annotations.full)]
      return(ramclustObj)
    }
  
  if(circle.plot) {
    
    ## circle plot first
    require(dplyr)
    require(circlepackeR)
    require(data.tree)
    ct <- ramclustObj$annotations.selected[,c("superclass", "class", "subclass", "direct.parent", "compound.name")]
    ct <- replace(ct, is.na(ct), "unknown")
    ct <- ct %>% count(superclass, class, subclass, compound.name)
    # ct <- ct[sample(1:nrow(ct), 100),]
    ct <- data.frame(
      'root' = rep('classyfire', nrow(ct)),
      ct
    )
    ct$pathString <- paste("classyfire", ct$superclass, ct$class, ct$subclass, ct$compound.name, sep = "/")
    p <- circlepackeR::circlepackeR(data.tree::as.Node(ct), size = "n", color_min = "hsl(56,80%,80%)", color_max = "hsl(341,30%,40%)")
    p$sizingPolicy$browser$fill = TRUE
    p$sizingPolicy$knitr$figure = TRUE
    p$sizingPolicy$knitr$defaultHeight = 0.5
    p$sizingPolicy$knitr$defaultWidth = 1
    htmlwidgets::saveWidget(p, file=paste0( getwd(), "/spectra/annotation.classyfire.summary.html"), 
                            knitrOptions = list())
    
    ## sunburst plot
    ## may be able to add more info with 'hover' function to add things like stats results, or compound description.
    ct <- ramclustObj$annotations.selected[,c("superclass", "class", "subclass", "direct.parent", "compound.name")]
    ct <- replace(ct, is.na(ct), "unassigned")
    for(i in 1:ncol(ct)) {
      ct[,i] <- gsub("-", "_", ct[,i])
    }
    d <- data.frame(
      "ids" = rep(NA, 0),
      "labels" = rep(NA, 0),
      "parents" = rep(NA, 0),
      "values"= rep(0, 0)
    )
    superclasses <- unique(ct[,"superclass"])
    for(i in 1:length(superclasses)) {
      superclass <- superclasses[i]
      ## add root node element
      ctsub <- ct[which(ct$superclass == superclass),]
      out <- data.frame(
        "ids" = superclass,
        "labels" = superclass,
        "parents" = "",
        "values" = nrow(ctsub)
      )
      d <- rbind(d, out)
      
      ## for all classes within superclass i
      classes <- unique(ctsub[,"class"])
      for(j in 1:length(classes)) {
        class <- classes[j]
        ctsub2 <- ctsub[which(ctsub$class == class),]
        if(nrow(ctsub2) == 0) next
        out <- data.frame(
          "ids" = paste0(c(superclass, class), collapse = "-"),
          "labels" = class,
          "parents" = superclass,
          "values" = nrow(ctsub2)
        )
        d <- rbind(d, out)
        
        ## for all subclasses within class j
        subclasses <- unique(ctsub2[,"subclass"])
        for(k in 1:length(subclasses)) {
          subclass <- subclasses[k]
          ctsub3 <- ctsub2[which(ctsub2$subclass == subclass),]
          if(nrow(ctsub3) == 0) next
          if(class == "unassigned" & subclass == "unassigned") next
          out <- data.frame(
            "ids" = paste0(c(class, subclass), collapse = "-"),
            "labels" = subclass,
            "parents" = paste0(c(superclass, class), collapse = "-"),
            "values" = nrow(ctsub3)
          )
          d <- rbind(d, out)
          
          ## for all direct.parents within subclass k
          # direct.parents <- unique(ctsub3[,"direct.parent"])
          # for(l in 1:length(direct.parents)) {
          #   direct.parent <- direct.parents[l]
          #   ctsub4 <- ctsub3[which(ctsub3$direct.parent == direct.parent),]
          #   if(nrow(ctsub4) == 0) next
          #   if(subclass == "unassigned" & direct.parent == "unassigned") next
          #   out <- data.frame(
          #     "ids" = paste0(c(subclass, direct.parent), collapse = "-"),
          #     "labels" = direct.parent,
          #     "parents" = paste0(c(class, subclass), collapse = "-"),
          #     "values" = nrow(ctsub4)
          #   )
          #   d <- rbind(d, out)
            
            # # for all compounds within direct.parent l
            # compounds <- unique(ctsub4[,"compound.name"])
            # for(m in 1:length(compounds)) {
            #   compound <- compounds[m]
            #   ctsub5 <- ctsub4[which(ctsub4$compound.name == compound),]
            #   if(nrow(ctsub5) == 0) next
            #   if(direct.parent == "unassigned" & compound == "unassigned") next
            #   out <- data.frame(
            #     "ids" = paste0(c(direct.parent, compound), collapse = "-"),
            #     "labels" = compound,
            #     "parents" = paste0(c(subclass, direct.parent), collapse = "-"),
            #     "values" = nrow(ctsub5)
            #   )
            #   d <- rbind(d, out)
          
          # for all compounds within subclass k
          compounds <- unique(ctsub3[,"compound.name"])
          for(m in 1:length(compounds)) {
            compound <- compounds[m]
            ctsub5 <- ctsub3[which(ctsub3$compound.name == compound),]
            if(nrow(ctsub5) == 0) next
            if(subclass == "unassigned" & compound == "unassigned") next
            out <- data.frame(
              "ids" = paste0(c(subclass, compound), collapse = "-"),
              "labels" = compound,
              "parents" = paste0(c(class, subclass), collapse = "-"),
              "values" = nrow(ctsub5)
            )
            d <- rbind(d, out)
          
            # }
          }
        }
      }
    }
    
    
    # for(i in 1:nrow(d)) {
    #   if(d$parents[i] == "") next
    #   if(length(which(d$ids == d$parents[i])) != 1) stop()
    # }
    
    
    fig2 <- plotly::plot_ly(d, 
                            ids = ~ids, 
                            labels = ~labels, 
                            parents = ~parents, 
                            values = ~values,
                            branchvalues = "total",
                            insidetextorientation='radial',
                            maxdepth=2,
                            type = 'sunburst')
    # fig2
    
    
    htmlwidgets::saveWidget(fig2, file=paste0( getwd(), "/spectra/annotation.classyfire.sunburst.html"), 
                            knitrOptions = list())
    
  }
}

# get.puchem.from.inchikey.local <- function(
    #     inchikey = NULL,
#     pc.inchi.file = "R:/RSTOR-PMF/Software/db/pubchem/CID-InChI-Key"
# ) {
#   d <- read.delim(file = pc.inchi.file, header = FALSE)
#   d <- d[,c(1,3)]
#   names(d) <- c("cid", "inchikey")
#   
# }

rc.merge.split.clusters <- function(
    ramclustObj = NULL,
    merge.threshold = 0.7,
    cor.method = 'spearman'
    
) {
  
  if(is.null(ramclustObj)) stop('please provide a ramclustObj', '\n')
  if(!is.numeric(merge.threshold)) merge.threshold <- as.numeric(merge.threshold)
  if(merge.threshold < 0 | merge.threshold > 1) stop("'merge.threshold' must be between zero and one", '\n')
  
  
  ## for all clusters, see if there are clusters at nearby retentiontimes which highly correlate
  ## if there are, join them, renumber featclus, and regenerate SpecAbund file.  
  orig.cl.n <- max(ramclustObj$featclus)
  cls <- max(ramclustObj$featclus):2
  # ramclustObj <- RC
  for(i in cls) {
    if(max(diff(sort(unique(ramclustObj$featclus)))) > 1) stop("error 1")
    potential.merges <- which(
      (abs(ramclustObj$clrt[1:i] - ramclustObj$clrt[i])) < (3*ramclustObj$clrtsd[i])
    )
    potential.merges <- potential.merges[!potential.merges == i]
    if(length(potential.merges) == 0) next
    
    rval <- cor(ramclustObj$SpecAbund[,i], ramclustObj$SpecAbund[,potential.merges], method = cor.method)
    merges <- potential.merges[which(rval >= merge.threshold)]
    if(length(merges) == 0) next
    
    merges <- max(merges)
    
    ramclustObj$featclus[which(ramclustObj$featclus == i)] <- merges
    # if(max(diff(sort(unique(ramclustObj$featclus)))) > 1) cat(i, "max diff = ", max(diff(sort(unique(ramclustObj$featclus)))), '\n')
    old.featclus <- ramclustObj$featclus
    new.featclus <- old.featclus
    decend.by.one <- which(old.featclus > i)
    new.featclus[decend.by.one] <- (old.featclus[decend.by.one])-1
    if(max(diff(sort(unique(new.featclus)))) > 1) stop("error 2")
    ramclustObj$featclus <- new.featclus
    
  }
  # sort(unique(ramclustObj$featclus))
  # sort(unique(RC$featclus))
  
  ## store SpecAbund sample names
  sa.rn <- dimnames(ramclustObj$SpecAbund)[[1]]
  
  # 0.99888
  
  # collapse feature dataset into spectrum dataset
  wts<-colSums(ramclustObj$MSdata)
  ramclustObj$SpecAbund<-matrix(nrow=nrow(ramclustObj$SpecAbund), ncol=max(ramclustObj$featclus))
  for (ro in 1:nrow(ramclustObj$SpecAbund)) { 
    for (co in 1:ncol(ramclustObj$SpecAbund)) {
      ramclustObj$SpecAbund[ro,co]<- weighted.mean(ramclustObj$MSdata[ro,which(ramclustObj$featclus==co)], wts[which(ramclustObj$featclus==co)])
    }
  }
  dimnames(ramclustObj$SpecAbund)[[1]] <- sa.rn
  
  strl <- nchar(max(ramclustObj$featclus)) - 1
  ramclustObj$cmpd <- paste("C", formatC(1:max(ramclustObj$featclus), digits = strl, flag = 0 ) , sep="")
  ramclustObj$ann <- ramclustObj$cmpd
  
  clrt<-aggregate(ramclustObj$frt, by=list(ramclustObj$featclus), FUN="mean")
  ramclustObj$clrt<-clrt[which(clrt[,1]!=0),2]
  clrtsd<-aggregate(ramclustObj$frt, by=list(ramclustObj$featclus), FUN="sd")
  ramclustObj$clrtsd<-clrtsd[which(clrtsd[,1]!=0),2]
  ramclustObj$nfeat<-as.vector(table(ramclustObj$featclus)[2:max(ramclustObj$featclus)])
  ramclustObj$nsing<-length(which(ramclustObj$featclus==0))
  ramclustObj$annconf<-rep(4, length(ramclustObj$clrt))
  ramclustObj$annnotes<-rep("", length(ramclustObj$clrt))
  dimnames(ramclustObj$SpecAbund)[[1]] <- ramclustObj$sample_names
  
  new.cl.n <- max(ramclustObj$featclus)
  cat(paste("Original cluster number =", orig.cl.n, '\n', "New cluster number =", new.cl.n, '\n'))
  
  return(ramclustObj)
  
}



#' annotate.msfinder
#'
#' After running MSFinder on .mat or .msp files, import the formulas that were predicted and their scores 
#' @param ramclustObj R object - the ramclustR object which was used to write the .mat or .msp files
#' @param mat.dir optional path to .mat directory
#' @param msp.dir optional path to .msp directory
#' @details this function imports the output from the MSFinder program to support annotation of the ramclustR object
#' @return new slot at $msfinder.formula.details
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export

annotate.msfinder.gcei <- function (ramclustObj = NULL, 
                                    mat.dir = NULL,
                                    priority.db = NULL,
                                    priority.inchikey = NULL,
                                    priority.db.factor = 0.9,
                                    priority.inchikey.factor = 0.9,
                                    adj.score = TRUE
) 
{
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  home.dir <- getwd()
  
  if(!is.null(priority.db) & !is.null(priority.inchikey)) {
    warning("both inchikey and database priority set - ensure they are independent. ", '\n')
  }
  
  r <- grep("msfinder.formula", names(ramclustObj))
  if (length(r) > 0) {
    warning("removed previously assigned MSFinder formulas and structures", 
            "\n")
    ramclustObj <- ramclustObj[-r]
    r <- grep("msfinder.structure", names(ramclustObj))
    if(length(r)>0) {
      ramclustObj <- ramclustObj[-r]
    }
    rm(r)
  }
  if (is.null(mat.dir)) {
    mat.dir = paste0(getwd(), "/spectra/mat")
  }
  
  if(is.null(ramclustObj$history)) {
    ramclustObj$history <- ""
  }
  
  findmain.summary <- ramclustObj$findmain.summary
  
  ## find all structure output files: 
  struc.files <- list.files(mat.dir, pattern = '.sfd', recursive = TRUE, full.names = TRUE)
  
  ## separate out spectal search results
  is.spec.db <- grep("Spectral", struc.files)
  spec.files <- struc.files[is.spec.db]
  rm(struc.files)
  # struc.files <- struc.files[-is.spec.db]
  
  ## read in spectral match results
  spec.path <- lapply(1:length(spec.files), FUN = function(x) {
    unlist(strsplit(spec.files[x], "/", fixed = TRUE))
  })
  spec.path <- t(data.frame(spec.path))
  # hypothesis <- as.character(spec.path[,(ncol(spec.path)-1)])
  cmpd       <- as.vector(spec.path[,(ncol(spec.path)-1)])
  
  spec.results <- data.frame(
    'cmpd' = vector(mode = 'character', length = 0),
    'rt' = vector(mode =  'numeric', length = 0),
    'median.signal.intensity'= vector(mode =  'numeric', length = 0),
    'compound.name' = vector(mode = 'character', length = 0),
    'inchikey'   = vector(mode = 'character', length = 0),
    'dbs'        = vector(mode = 'character', length = 0),
    'spectral.match.score'  = vector(mode =  'numeric', length = 0),
    'direct.parent.class'   = vector(mode = "character", length = 0),
    'chemontid'  = vector(mode = 'character', length = 0)
  )
  
  for(i in 1:length(spec.files)) {
    tmp <- readLines(spec.files[i])
    tmp.names <- grep("NAME: ", tmp)
    if(length(tmp.names) > 0) {
      tmp.scores <- grep("TotalScore:", tmp)
      tmp.inchikey <- grep("INCHIKEY: ", tmp)
      tmp.dbs     <- grep("RESOURCES: ", tmp)
      tmp.subclass<- grep("Ontology: ", tmp)
      tmp.chemontid<- grep("OntologyID: ", tmp)
      same.length <- all.equal(
        length(tmp.names),
        length(tmp.scores),
        length(tmp.inchikey),
        length(tmp.dbs),
        length(tmp.subclass),
        length(tmp.chemontid)
      )
      if((length(tmp.names)>0) & same.length) {
        tmp.out <- data.frame(
          'cmpd'       = rep(cmpd[i], length(tmp.names)),
          'rt'         = rep(ramclustObj$clrt[as.integer(as.numeric(gsub("C", "", cmpd[i])))], length(tmp.names)),
          'median.signal.intensity' = rep(median(ramclustObj$SpecAbund[,as.integer(as.numeric(gsub("C", "", cmpd[i])))], na.rm = TRUE), length(tmp.names)),
          'compound.name'   = gsub("NAME: ", "", tmp[tmp.names]),
          'inchikey'   = gsub("INCHIKEY: ", "", tmp[tmp.inchikey]),
          'dbs'        = gsub("RESOURCES: ", "", tmp[tmp.dbs]),
          'spectral.match.score'      = as.numeric(gsub("TotalScore: ", "", tmp[tmp.scores])),
          'direct.parent.class'   = gsub("Ontology: ", "", tmp[tmp.subclass]),
          'chemontid'  = gsub("OntologyID: ", "", tmp[tmp.chemontid])
        )
      }
      if(nrow(tmp.out)>0) {
        spec.results <- rbind(spec.results, tmp.out)
      }
    }
  }
  
  spec.results <- data.frame(
    'assigned' = rep(FALSE, nrow(spec.results)),
    spec.results
  )
  
  if(adj.score) {
    spec.results$spectral.match.score <- spec.results$spectral.match.score^0.5
  }
  spec.results$spectral.match.score <- spec.results$spectral.match.score^0.5
  ## assign database priority score, if applicable
  if(!is.null(priority.db)) {
    priority.factor.v <- rep(priority.db.factor, nrow(spec.results))
    for(i in 1:length(priority.db)) {
      db.match <- grep(priority.db[i], spec.results$dbs, ignore.case = TRUE)
      priority.factor.v[db.match] <- 1
    }
    spec.results$db.priority.factor <- priority.factor.v
    spec.results$total.score <- spec.results$spectral.match.score * priority.factor.v
  }
  
  
  ## assign inchikey priority score, if applicable
  if(!is.null(priority.inchikey)) {
    priority.factor.v <- rep(priority.inchikey.factor, nrow(spec.results))
    short.inchikey.priority <- sapply(priority.inchikey, FUN = function(x) {
      unlist(strsplit(x, "-"))[1]
    })
    short.inchikey.struc <- sapply(spec.results$inchikey, FUN = function(x) {
      unlist(strsplit(x, "-"))[1]
    })
    inchi.match <- short.inchikey.struc %in% short.inchikey.priority
    priority.factor.v[inchi.match] <- 1
    spec.results$inchikey.priority.factor <- priority.factor.v
    spec.results$total.score <- spec.results$spectral.match.score * priority.factor.v
  }
  
  if(is.null(spec.results$total.score)) {
    spec.results$total.score <- spec.results$spectral.match.score
  }
  
  ## reset annotations
  ramclustObj$ann  <- ramclustObj$cmpd
  
  ## update structure table
  # ramclustObj$use.spectral.match <- rep(FALSE, length(ramclustObj$cmpd))
  
  for(i in 1:length(ramclustObj$ann)) {
    if(ramclustObj$ann[i] != ramclustObj$cmpd[i]) next
    cmpd.match <- which(spec.results$cmpd == ramclustObj$cmpd[i])
    if(length(cmpd.match)==0) next
    ann.sub <- spec.results[cmpd.match,]
    ramclustObj$use.spectral.match[i] <- TRUE
    sp.m.index <- which.max(ann.sub$total.score)
    spec.results$assigned[cmpd.match[sp.m.index]] <- TRUE
    ramclustObj$ann[i] <-  ann.sub$compound.name[sp.m.index]
    ramclustObj$inchikey <- ann.sub$spectral.match.inchikey[sp.m.index]
    ramclustObj$spectral.match.score[i] <- ann.sub$spectral.match.score[sp.m.index]
  }
  
  ramclustObj$annotations.full <- spec.results
  ramclustObj$annotations.selected <- spec.results[which(spec.results$assigned),2:ncol(spec.results)]
  
  ramclustObj$history$msfinder <- paste(
    "MSFinder (Tsugawa 2016) was used for spectral matching.",
    "Results were imported into the RAMClustR object.", 
    "A total of", nrow(spec.results), "spectral matches were were evaluated", 
    "for", length(unique(spec.results$cmpd)), "compounds.", 
    "A complete spreadsheet of all spectral matches can be found in the 'spectra/all.spectral.matches.csv' file,",
    "and a subset of only those selected for annotation can be found in the 'spectra/assigned.annotations.csv' file.")
  
  
  
  if(!is.null(priority.db)) {
    ramclustObj$history$msfinder <-paste(
      as.character(ramclustObj$history$msfinder),
      "The following database(s) were assigned as 'priority': ", paste0(paste(priority.db, collapse = ', '), "."),
      "The database priority.factor was set to", priority.db.factor, "to decrease spectral match scores for compounds which failed to match priority database(s)."
    )
  }
  
  if(!is.null(priority.inchikey)) {
    ramclustObj$history$msfinder <-paste(
      as.character(ramclustObj$history$msfinder),
      "The a list of",  length(priority.inchikey), "inchikeys was provided.", 
      "The inchikey priority.factor was set to", priority.inchikey.factor, "to decrease scores for compounds with non-matching inchikey(s)."
    )
  }
  ramclustObj$history$msfinder <-paste(
    as.character(ramclustObj$history$msfinder),
    "The highest spectral match score was selected for each compound, considering all matches and any relevent score adjustments for inchikey or database prioritization."
  )
  
  ## write annotaion tables to 'spectra' directory
  out.dir <- gsub(basename(mat.dir), "",  mat.dir)
  write.csv(ramclustObj$annotations.full, file = paste0(out.dir, "all.spectral.matches.csv"), row.names = FALSE)
  write.csv(ramclustObj$annotations.selected, file = paste0(out.dir, "assigned.annotations.csv"), row.names = FALSE)
  return(ramclustObj)
}


#' write.gcei.mat
#'
#' Export GC-MS EI spectra for spectral searching in MSFinder
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @details exports files to a directory called 'spectra'.  a new directory 'spectra/mat' is created to hold the individual mat files.  
#' @return nothing, just exports files to the working directory
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept interpretMSSpectrum
#' @concept xcms
#' @author Corey Broeckling
#' @importFrom methods is
#' @export 
#' 

write.gcei.mat <- function(
    ramclustObj = NULL
) {
  
  if(!is(ramclustObj, "hclus") & 
     ramclustObj$dist.method != "RAMClustR") {
    stop("this is not a RAMClustR object")
  }
  
  if(!dir.exists('spectra')) {
    dir.create('spectra')
  }
  
  
  if(!dir.exists('spectra/mat')) {
    dir.create('spectra/mat')
  }
  
  ion.mode <- "Positive"
  
  out.list <- as.list(rep(NA, length(ramclustObj$cmpd)))
  
  for(i in 1:length(ramclustObj$cmpd)) {
    ions <- which(ramclustObj$featclus == i)
    
    spectrum <- data.frame(
      'mz' = round(0.99888*ramclustObj$fmz[ions], digits = 0),
      'int' = round(ramclustObj$msint[ions], digits = 0)
    )
    
    spectrum <- spectrum[order(spectrum[,"mz"], decreasing = FALSE), ]
    
    out <- paste0(
      "NAME: ", ramclustObj$cmpd[i], '\n', 
      "PRECURSORTYPE: [M]+.", '\n',
      "IONMODE: ", ion.mode, '\n',
      "COLLISIONENERGY: 70", '\n',
      "SPECTRUMTYPE: Centroid", '\n',
      "RETENTIONTIME: ", round(ramclustObj$clrt[i], 2), '\n'
    )
    
    if(any(names(ramclustObj)=="clri")) {paste0(
      out,
      "RETENTIONINDEX: ", round(ramclustObj$clri[i], 2),  '\n'
    )
    }
    
    out <- paste0(out, 
                  "PRECURSORMZ: ", round(max(spectrum[,1]),1), '\n' 
    )
    
    out <- paste0(
      out,
      "INSTRUMENTTYPE: GC-EI-Q", '\n',  
      "INSTRUMENT: ", '\n',
      "Authors: ", '\n', 
      "License: ", '\n',
      "FORMULA: ", '\n',
      "ONTOLOGY: ", '\n',
      "SMILES: ", '\n',  
      "INCHIKEY: ", '\n',
      "INCHI: ", '\n', 
      "METABOLITENAME: ", '\n',
      "SCANNUMBER: -1 ", '\n',
      "RETENTIONTIME: 0", '\n',
      "RETENTIONINDEX: 0", '\n',
      "CCS: 0", '\n',
      "INTENSITY: 0", '\n',
      "#Specific field for labeled experiment", '\n',
      "IsMarked: False", '\n'
    )
    
    out <- paste0(out,
                  "MSTYPE: MS1", '\n', "Num Peaks: 0", '\n'
    )
    
    out <- paste0(out,
                  "MSTYPE: MS2", '\n',
                  "Num Peaks: ", nrow(spectrum), '\n'
    )
    for(j in 1:nrow(spectrum)) {
      out <- paste0(out, 
                    spectrum[j,"mz"], 
                    " ", 
                    round(spectrum[j,"int"]),
                    '\n'
      )
    }
    out.list[[i]] <- out
  }
  
  for(i in 1:length(out.list)) {
    sink(paste0("spectra/mat/", ramclustObj$cmpd[[i]], ".mat"))
    cat(out.list[[i]], '\n')
    sink()
  }
  
}


#' rc.expand.sample.names
#'
#' turn concatenated sample names into factors  
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param delim what delimiter should be used to separate names into factors?  '-' by default
#' @param factor.names logical or character vector.  if TRUE, user will enter names one by on in console.  If character vector (i.e. c("trt", "time")) names are assigned to table
#' @param quiet logical .  if TRUE, user will not be prompted to enter names one by on in console.
#' @details THis function only works on newer format ramclustObjects with a $phenoData slot.
#' @details This function will split sample names by a delimiter, and enable users to name factors
#' @return  ramclustR object with normalized data.   
#'  
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export

rc.expand.sample.names <- function(
    ramclustObj = NULL,
    delim = "-",
    factor.names = TRUE,
    quiet = FALSE
) {
  
  params <- c(
    "delim" = delim
  )
  
  
  if(!is.null(ramclustObj$phenoData$sample.names.sample_name)) {
    sn <- as.character(ramclustObj$phenoData$sample.names.sample_name)
  }
  if(!is.null(ramclustObj$phenoData$sample.names.sn)) {
    sn <- as.character(ramclustObj$phenoData$sample.names.sn)
  }
  
  if(!is.null(ramclustObj$phenoData$sample.names)) {
    sn <- as.character(ramclustObj$phenoData$sample.names)
  }
  
  cat(paste(sn, collapse = ", "), '\n')
  
  if(!any(ls()=="sn")) {
    stop('missing sample names in phenoData slot', '\n')
  }
  des <-  strsplit(sn, delim)
  l <- sapply(1:length(des), FUN = function(x) {
    length(des[[x]])
  })
  if(length(table(l)) != 1) {
    cat("delimited sample names have variable lengths ranging from: ", '\n',
        range(l)[1], " to ", range(l)[2], '\n')
    ch <- which(l != median(l))
    stop("please fix sample names:", '\n', paste(" ", sn[ch], sep = '\n'))
  }
  des <- data.frame(t(data.frame(des, check.names = FALSE)), 
                    stringsAsFactors = FALSE, check.names = FALSE)
  
  if(is.logical(factor.names) & !quiet) {
    if(factor.names){
      for(x in 1:ncol(des)) {
        cat(
          "column",x, "variables:",'\n',
          unique(des[,x]), '\n')
        fn <- readline(prompt=cat("Type name and press [enter] to continue:", '\n'))
        dimnames(des)[[2]][x] <- fn
      }
    }
  }
  
  if(is.character(factor.names)) {
    if(length(factor.names) != ncol(des)) {
      stop(length(factor.names), "factor names and", ncol(des), "factors - please correct", '\n')
    }
    dimnames(des)[[2]] <- factor.names
  }
  rn <- row.names(ramclustObj$phenoData)
  ramclustObj$phenoData <- cbind(ramclustObj$phenoData, des[,1:ncol(des)])
  
  if(is.null(ramclustObj$params)) {ramclustObj$params <- list()}
  ramclustObj$params$rc.expand.sample.names <- params
  
  row.names(ramclustObj$phenoData) <- rn
  return(ramclustObj)
  
}

define_samples <- function(ramclustObj, tag) {
  ## define samples in each set
  if(length(tag) == 0) {
    stop("no tag provided", "\n")
  }
  
  samples <- grep(tag[1], ramclustObj$phenoData$sample.names)
  samples <- samples[which(samples <= nrow(ramclustObj$MSdata))]
  
  if (length(samples) == 0) {
    stop("no ", tag, " samples found using the tag ", "'", tag, "'", "\n")
  }
  return(samples)
}


#' rc.qc
#'
#' summarize quality control for clustering and for quality control sample variation based on compound ($SpecAbund) and feature ($MSdata and $MSMSdata, if present)
#'
#' @param ramclustObj ramclustR object to analyze
#' @param qc.tag qc.tag character vector of length one or two.  If length is two, enter search string and factor name in $phenoData slot (i.e. c("QC", "sample.type"). If length one (i.e. "QC"), will search for this string in the 'sample.names' slot by default.  
#' @param remove.qc logical - if TRUE (default) QC injections will be removed from the returned ramclustObj (applies to $MSdata, $MSMSdata, $SpecAbund, $phenoData, as appropriate). If FALSE, QC samples remain.
#' @param npc number of Principle components to calcuate and plot
#' @param scale "pareto" by default: PCA scaling method used
#' @param outfile.basename base name of output files. Extensions added internally. default = "ramclustQC"
#' @param view.hist logical.  should histograms be plotted? 
#' @details plots a ramclustR summary plot.  first page represents the correlation of each cluster to all other clusters, sorted by retention time.  large blocks of yellow along the diaganol indicate either poor clustering or a group of coregulated metabolites with similar retention time.  It is an imperfect diagnostic, particularly with lipids on reverse phase LC or sugars on HILIC LC systems.  Page 2: histogram of r values from page 1 - only r values one position from the diagonal are used.  Pages 3:5 - PCA results, with QC samples colored red.  relative standard deviation calculated as sd(QC PC scores) / sd(all PC scores).  Page 6: histogram of CV values for each compound int he dataset, QC samples only.  
#' @return   new RC object. Saves output summary plots to pdf and .csv summary tables to new 'QC' directory. If remove.qc = TRUE, moves QC samples to new $QC slot from original position. 
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @importFrom pcaMethods pca
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export

rc.qc<-function(ramclustObj=NULL,
                qc.tag="QC",
                remove.qc = FALSE,
                npc=4,
                scale="pareto",
                outfile.basename ="ramclustQC",
                view.hist = TRUE
                
){
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(is.null(qc.tag)) {
    stop("qc.tag = NULL; qc.tag must be defined to enable QC variance examination.", '\n')
  }
  
  if(is.null(outfile.basename)) {
    outfile.basename <- "ramclustQC"
  }
  
  define_samples <- function(ramclustObj, tag) {
    ## define samples in each set
    if(length(tag) == 0) {
      stop("no tag provided", "\n")
    }
    
    samples <- grep(tag[1], ramclustObj$phenoData$sample.names)
    samples <- samples[which(samples <= nrow(ramclustObj$MSdata))]
    
    if (length(samples) == 0) {
      stop("no ", tag, " samples found using the tag ", "'", tag, "'", "\n")
    }
    return(samples)
  }
  
  
  do.sets <- c("MSdata", "SpecAbund")
  
  if(is.null(ramclustObj$SpecAbund)) {
    do.sets <- do.sets[!(do.sets %in% "SpecAbund")]
  } 
  
  if(is.null(ramclustObj$MSdata)) {
    do.sets <- do.sets[!(do.sets %in% "MSdata")]
  } 
  
  do.sets.rows <- sapply(
    c(do.sets, "phenoData"), 
    FUN = function(x) {
      nrow(ramclustObj[[x]])
    })
  
  if(!sd(do.sets.rows) == 0) {
    stop("number of rows in MSdata, SpecAbund, and phenoData sets are not identical.")
  }
  
  ## define QC samples in each set
  qc <- define_samples(ramclustObj = ramclustObj, tag = qc.tag)
  
  ## create directory
  dir.create("QC")
  
  ## if cv threshold for compounds has been applied, use it, else create 'all cmpds' vector
  if(!is.null(ramclustObj$SpecAbund)) {
    if(!is.null(ramclustObj$cmpd.use)) {
      cmpd.use <- ramclustObj$cmpd.use
    } else {
      cmpd.use <- rep(TRUE, length(ramclustObj$ann))
    }
  }
  
  #visualize clustering
  ## if clustering was perfect, we should see a normal distribution of 
  ## correlational r values 1 step from the diagonal
  ## imperfect clustering introduces right skew
  ## load("datasets/RCobject.Rdata")
  
  # if(!is.null(ramclustObj$clrt)) {
  #   
  #   ## create file to collect figures. 
  #   pdf(file=paste("QC/", "ramclust_clustering_diagnostic.pdf", sep=""), 
  #       useDingbats=FALSE, width=8, height=8)  
  #   o<-order(ramclustObj$clrt[cmpd.use])
  #   c<-cor(ramclustObj$SpecAbund[,cmpd.use][,o])
  #   d<-diag(as.matrix((c[2:(nrow(c)), 1:ncol(c)-1])))
  #   hist(d, breaks=50, main="")
  #   title(main="histogram of pearson's r for each cluster to its adjacent cluster (by time)", cex.main=0.8,
  #         sub=paste("skew =", round(e1071::skewness(d), digits=3), " :values near zero are better", '\n', 
  #                   'WARNING:metabolic relationships will confound interpretation of this plot'), cex.sub=0.6)
  #   
  #   # ideally heatmap will have a bright yellow diagonal with no yellow squares near the diagonal
  #   # this is slow for larger numbers of clusters
  #   gplots::heatmap.2(c^2, trace="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, main="pearsons r^2, clusters sorted by rt", cex.main=0.5,
  #                     cexRow=0.02 + 1/log10(length(o)), cexCol=0.02 + 1/log10(length(o)))
  #   dev.off()
  # }
  # 
  
  ## PCA of QC samples
  ## histogram of feature and/or compound CVs for QC samples
  
  # cols<-rep(8, length(ramclustObj$sample_names))
  cols<-rep(8, nrow(ramclustObj$phenoData))
  cols[qc]<-2
  
  for(x in do.sets) {
    
    ## PCA plot
    if(x == "SpecAbund") {
      td <- ramclustObj[[x]][,cmpd.use]
    } else {
      td <- ramclustObj[[x]]
    }
    
    # if(!is.null(ramclustObj$MSMSdata) & x == "MSdata") {
    #   td <- td + ramclustObj$MSMSdata
    # }
    
    if(min(dim(td)) < npc) {npc <- min(dim(td))}
    PCA<-pcaMethods::pca(td, scale=scale, nPcs=npc, center=TRUE)
    sc<-PCA@scores
    write.csv(sc, file = paste0("QC/", outfile.basename, "_", x, "_pcascores.csv"))
    pdf(file = paste0("QC/", outfile.basename, "_", x, "_qc_diagnostic.pdf"), useDingbats=FALSE, width=8, height=8)  
    
    ld<-PCA@loadings
    for(i in 1:(ncol(sc)-1)) {
      plot(sc[,i], sc[,i+1], col=cols, pch=19, main=paste(
        "PCA analysis:", x, if(x == "SpecAbund") {
          "(compounds)"
        } else {"(features)"}
      ),
      xlab=paste("PC", i, "::  r^2 =", round(PCA@R2[i], digits=2), "::  QC(rel sd) = ", 
                 round(sd(sc[qc,i])/sd(sc[,i]), digits=2) ),
      ylab=paste("PC", i+1, "::  r^2 =", round(PCA@R2[i+1], digits=2), "::  QC(rel sd) = ", 
                 round(sd(sc[qc,i+1])/sd(sc[,i+1]), digits=2) )
      )
      legend(qc.tag, text.col=2, x="topright", bty="n")
    }
    
    ## histogram of QC relative standard deviations for all compounds/clusters
    sds<-apply(td[qc,], 2, FUN="sd", na.rm=TRUE)
    #cat(sds, '\n')
    means<-apply(td[qc,], 2, FUN="mean", na.rm=TRUE)
    cvs<-sds/means
    
    if(x == "MSdata") {
      ramclustObj$qc.cv.feature <- cvs
      ramclustObj$qc.cv.feature.msdata <- cvs
      if(!is.null(ramclustObj$MSMSdata)) {
        sds<-apply(ramclustObj$MSMSdata[qc,], 2, FUN="sd", na.rm=TRUE)
        #cat(sds, '\n')
        means<-apply(ramclustObj$MSMSdata[qc,], 2, FUN="mean", na.rm=TRUE)
        msms.cvs<-sds/means
        ramclustObj$qc.cv.feature.msmsdata <- msms.cvs
        cvs <- pmin(ramclustObj$qc.cv.feature.msdata, msms.cvs)
        ramclustObj$qc.cv.feature <- cvs
      }
      
    } else {
      ramclustObj$qc.cv.cmpd <- cvs
    }
    qs<-quantile(cvs, probs=seq(0,1,0.2), na.rm=TRUE)
    hist(cvs, breaks=50, main="")
    title(paste("histogram of", x,  "CVs from QC samples"), line=2.7)
    title("20% quantiles in red on top axis", col.main =2, cex.main=0.7, line=2)
    axis(side=3, col=2, col.ticks=2, col.axis=2, round(qs, digits=3), labels=TRUE, las=2, cex.axis=0.4)
    dev.off()
    
    if(view.hist) {
      hist(cvs, breaks=50, main="")
      title(paste("histogram of", x,  "CVs from QC samples"), line=2.7)
      title("20% quantiles in red on top axis", col.main =2, cex.main=0.7, line=2)
      axis(side=3, col=2, col.ticks=2, col.axis=2, round(qs, digits=3), labels=TRUE, las=2, cex.axis=0.4)
      
    }
    
    
    if(x == "SpecAbund") {
      out <- data.frame(
        "cmpd" = ramclustObj$cmpd[cmpd.use],
        "annotation" = ramclustObj$ann[cmpd.use],
        "rt" = ramclustObj$clrt[cmpd.use], 
        "mean.int" = means,
        "cv" = cvs
      )
    } else {
      if(is.null(ramclustObj$fmz)) next
      out <- data.frame(
        "mz" = ramclustObj$fmz,
        "rt" = ramclustObj$frt,
        "mean.int" = means,
        "cv" = cvs
      )
      if(length(ramclustObj$labels) > 0) {
        out <- data.frame(
          out, 
          "feature" = ramclustObj$labels,
          "cluster" = ramclustObj$featclus
        )
      }
    }
    write.csv(out, file = paste0("QC/", outfile.basename, "_", x, "_cv_summary.csv"))
  }
  
  if(remove.qc) {
    ramclustObj$qc <- list()
    for(x in c("phenoData", do.sets)) {
      ramclustObj$qc[[x]] <- ramclustObj[[x]][qc,]
      ramclustObj[[x]] <- ramclustObj[[x]][-qc,]
    }
  }
  
  ramclustObj$history$qc.summary <- paste(
    
    "Variance in quality control samples was described using the",
    "rc.qc function within ramclustR. Summary statistics are provided",
    "including the relative standard deviation of QC samples to all",
    "samples in PCA space, as well as the relative standard deviation",
    "of each feature/compound in QC samples, plotted as a histogram.", 
    if(!is.null(ramclustObj$cmpd.use)) {" Only compounds which passed the CV filter are reported."}
  )
  
  return(ramclustObj)
}


