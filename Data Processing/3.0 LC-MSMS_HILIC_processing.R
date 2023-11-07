# metabolomics workflow to process raw mzML files from HILIC LC-MS/MS

###############################################
## VARIABLES: 
###############################################

## search terms for QC and Blank samples
qc.tag                    <- "qc" #find qc files (case sensitive)
blank.tag                 <- "prep.blank" #ControlWater used as blanks

###############################################
## RUN THIS LOCALLY BEFORE RUNNING ON LINUX SERVER: 
## after runing, lead line with # to comment out text. 
###############################################
#library(csu.pmf.tools)
#ExpDes <- csu.pmf.tools::pmfDefineExperiment()
#save(ExpDes, file = "ExpDes.Rdata")


###############################################
## load libraries
###############################################
library(xcms)
#library(RAMClustR)
#library(InterpretMSSpectrum)

###############################################
## session preparation
###############################################

## generate file processing list from seq.csv file.
dir.create("datasets0.20")
seq                   <- read.csv('seq.csv', header = TRUE, check.names = FALSE) # seq <- seq[1:2,]
mzml                  <- c(paste0(as.character(seq[,1]), "_01.mzML")) #MS1 files will have the _01 extension
factor.names          <- names(seq)[3:ncol(seq)]

## check to ensure all files expected are present
if(any(!file.exists(mzml))) {
  stop(cat("missing files:", '\n', paste(mzml[!file.exists(mzml)], collapse = ", "), '\n'))
}
pd                    <- data.frame("filename" = mzml, seq[,2])

## qc samples are used for processing group in xcms - impacts feature filter in groupChromPeaks
sample.groups         <-  rep(0, nrow(seq))
qc                    <-  which(grepl(qc.tag, seq[,2], ignore.case = TRUE))
blanks                <-  which(grepl(blank.tag, seq[,2], ignore.case = TRUE))

#defines sample groups. you can define these however you want. in theory you can set groups to your treatments. so we could do 
#hairyVetch.tag <- "HairyVetch" #find Hairy Vetch files (case sensitive)
#hairyVetch <-  which(grepl(hairyVetch.tag, seq[,2], ignore.case = TRUE))
#sample.groups[hairyVetch] <- 1

sample.groups[qc]     <-  1 #this sets anything labeled as QC to a 1
sample.groups[blanks] <-  2 #this sets anything labeled as a blank to a 2
sample.groups         <- c(sample.groups) #remove 3+sample.groups if you are running DDA (MSn) data, this will be a vector with the same length as the number of files


#this combines the real samples (both samples and qcs)
not.blanks            <- which(sample.groups %in% c(0,1)) # sets the samples groups where 0s are samples and 1s are QCs


###############################################
## load experimental design
## if not premade, need to make a new experimental design
###############################################
load("ExpDes.Rdata")


###############################################
#Section A XCMS - process raw data
#perform feature detection, alignment, grouping, filling
###############################################

###############################################
## define peak detection parameters
cw.param <- CentWaveParam(
  ppm = 30,
  peakwidth = c(2.2, 15),
  snthresh = 5,
  prefilter = c(3, 10),
  mzCenterFun = "wMean",
  integrate = 1L,
  mzdiff = 0.015,
  fitgauss = TRUE,
  noise = 2,
  verboseColumns = TRUE,
  roiList = list(),
  firstBaselineCheck = TRUE,
  roiScales = numeric(),
  extendLengthMSW = TRUE
)

###############################################
##read raw data and perform peak detection
  raw_data          <- readMSData(files = mzml, pdata = new("NAnnotatedDataFrame", pd),
                                  mode = "onDisk") # takes mzML file list we have, index this onto the disk to help with memory demands.
  h                 <- header(raw_data) #ignore warning
  orig.msLevel      <- raw_data@featureData@data$msLevel
  raw_data@featureData@data$msLevel           <- rep(1, length(orig.msLevel))
  x.data            <- findChromPeaks(raw_data, param = cw.param) #performs peak detection, takes a while (hours), computational intensive!
  x.data
  save(x.data, file = "datasets0.20/x.1.Rdata")
  cat('\n', '\n', " - finished peak detection", '\n', '\n')
  gc()
  rm(raw_data)

################################################
### XCMS feature grouping, pre-RT correction

  pdp <- PeakDensityParam(sampleGroups = sample.groups, binSize = 0.015,
                          minFraction = 0.20, bw = 4) #define processing parameters (minFrac .45 retains features if they are found in 45% of samples, bin size combines features with related peak width and accounts for a relatively wide RT shift)
  x.data <- groupChromPeaks(x.data, param = pdp) #group features from sequential injections if masses and RT are similar enough
  save(x.data, file = "datasets0.20/x.2.Rdata")
  cat('\n', '\n', " - finished feature grouping - pre-RT adjustment", '\n', '\n')
  gc()

################################################
### XCMS retention time adjustment, peak density

  #how often should a feature be found in the samples, .35 is 35% of samples
  pgp <- PeakGroupsParam(
    minFraction = 0.20
  )
  x.data <- adjustRtime(x.data, param = pgp) #apply pgp parameters to adjust for RT drift
  save(x.data, file = "datasets0.20/x.3.Rdata")
  cat('\n', '\n', " - finished RT adjustment", '\n', '\n')
  x.data
  gc()

################################################
### XCMS feature grouping, post-RT correction
  #adjust pdp bw to be much lower since we now about RT adjusted peaks so we can set a much tighter RT to define sameness across samples
  pdp <- PeakDensityParam(sampleGroups = sample.groups, binSize = 0.015,
                          minFraction = 0.20, bw = 1.75)
  
  x.data <- groupChromPeaks(x.data, param = pdp)
  save(x.data, file = "datasets0.20/x.4.Rdata")

################################################
##XCMS fillPeaks

  fpp               <- FillChromPeaksParam(expandMz = 0, expandRt = 0, ppm = 0) #fill missing values with local noise/signal that did not pass the xcms s/n threshold. no quality evaluation done with this data, just filling the 0s
  x.data            <- fillChromPeaks(x.data, param = fpp, msLevel = 1L)
  save(x.data, file="datasets0.20/x.5.Rdata")
  cat('\n', '\n', " - finished fill peaks", '\n', '\n')
  x.data
  
  
  cat('\n', '\n', " - HILIC workflow complete", '\n', '\n')
