###############################################
## VARIABLES: 
###############################################

## search terms for QC and Blank samples
qc.tag                    <- "qc" #find qc files (case sensitive)
blank.tag                 <- "prep.blank" #ControlWater used as blanks

seq                   <- read.csv('seq.csv', header = TRUE, check.names = FALSE) # seq <- seq[1:2,]
mzml                  <- c(paste0(as.character(seq[,1]), "_01.mzML")) #MS1 files will have the _01 extension
factor.names          <- names(seq)[3:ncol(seq)]

library(RAMClustR)

###############################################
## Section B Start - Build RC obect, remove features
## found at levels in blank samples similar to QC.
###############################################

load("ExpDes.Rdata")

#download all xcms objects from alpine and load them here:
load("datasets0.20/x.1.Rdata")
load("datasets0.20/x.2.Rdata")
load("datasets0.20/x.3.Rdata")
load("datasets0.20/x.4.Rdata")
load("datasets0.20/x.5.Rdata")


## build empty ramclustObj, pull out dataframe that pulls out all of the signal intensities for all of the features
RC <- RAMClustR::rc.get.xcms.data(
  xcmsObj = x.data,
  ExpDes = ExpDes,
  MStag = "_01.mzML", 
  mzdec = 4 #keep four decimal points
)

dimnames(RC$phenoData)[[2]][2] <- "sample.names" 

## turn sample names into factor names
RC <- RAMClustR::rc.expand.sample.names(
  ramclustObj = RC,
  delim = "-", 
  factor.names = c("sampleID", "species", "rep", "exoA", "timepoint", "run_order", "run_batch") #you have to manually add these and make sure they are the same as in ExpDes
)

## replace missing values with low level noise (optionally also zero values), sometimes fill peaks from before doesn't replace all the zeros
RC <- RAMClustR::rc.feature.replace.na(
  ramclustObj = RC,
  replace.int = 0.5, #replaces zeros with noise, 
  replace.noise = 0.5, #adds noise at half the intensity of the lowest intensity peak
  replace.zero = FALSE
)


## create 'qc' plots pre normalization for comparison, how consistent are the QCs before RT drift correction?
#this step won't work unless you have at least 3 reps of QC...
RC <- RAMClustR::rc.qc(
  ramclustObj = RC,
  outfile.basename = "preNorm",
  qc.tag = qc.tag,
  remove.qc = FALSE, 
  npc = 4, 
  view.hist = TRUE
)

save(RC, file = "datasets0.20/RCobject_0.Rdata")
load("datasets0.20/RCobject_0.Rdata")



##############################################################
## normalize feature values to QC samples by run order and batch effect, look for trends in run or batch order and correct for any patterns in drift that are noticed. apply fit to samples

RC$sample_names <- RC$phenoData[,2]

RC <- RAMClustR::rc.feature.normalize.qc(
  ramclustObj = RC,
  order = as.numeric(as.character(RC$phenoData$run_order)), #these must be specified in the ExpDes which variables in the sample names correspond to the run order and batch
  batch = as.numeric(as.character(RC$phenoData$run_batch)),
  qc.tag = qc.tag,
  rsq.cut = 0.1, #r squared cutoff= if value is below 0.1, AND pvalue for that fit is greater than .5, dont correct
  p.cut = 0.05,
  output.plot = TRUE
)

#Features were normalized by linearly regressing run order versus qc feature intensities to account for instrument signal intensity drift. Only features with a regression pvalue less than 0.05 and an r-squared greater than 0.1 were corrected.  Of 31787 features, 6801 were corrected for run order effects in at least one batch.  Batch effects were normalized to median intensity for each feature.

## remove 'system peak' features, you can remove all feature that are not at least 3 times as high in the qc as in the blank. 
## normally corey uses this when there are blanks included. if no blanks are used then comment out.
RC <- RAMClustR::rc.feature.filter.blanks(
  ramclustObj = RC,
  qc.tag = qc.tag,
  blank.tag = blank.tag,   ## tag to be recognized in sample name
  sn = 3  ## if qc samples are not this times as high in signal than blanks, feature is removed
)

#43.4% of features move forward
#ma MSdata 
#df phenoData 
#Features which failed to demonstrate signal intensity of at least 3 fold greater in QC samples than in blanks were removed from the feature dataset. 17987 of 31787 features were removed.
##############################################################


## filter out noisy features
RC <- RAMClustR::rc.feature.filter.cv(
  ramclustObj = RC,
  qc.tag = qc.tag,
  max.cv = 1 #CV greater than .5 for any feature, it will be thrown out.
)

#MSdata :  18759 passed the CV filter 
#Features were filtered based on their qc sample CV values. Only features with CV vaules less than or equal to 
#0.5 in MSdata set were retained. 13028 of 31787 features were removed.


RC <- rc.feature.normalize.tic(
  ramclustObj = RC
)




MSdata_TIC_normFeat<- t(RC$MSdata)
write.csv(MSdata_TIC_normFeat, file= "MSdata_TIC_normFeat.csv")
write.csv(RC$fmz, file= "mz__TIC_normFeat.csv")
write.csv(RC$frt, file= "RT__TIC_normFeat.csv")

load("datasets0.20/RCobject_1_Keep_unk_TICnorm.Rdata")

#CLUSTERING: CLUSTER features into compounds, this may take some time...

RC <- RAMClustR::rc.ramclustr(
  ramclustObj = RC,
  st = NULL,
  sr = 0.5,
  maxt = NULL,
  deepSplit = FALSE,
  blocksize = 2000,
  mult = 5,
  hmax = NULL,
  collapse = TRUE,
  minModuleSize = 3,
  linkage = "average",
  cor.method = "pearson",
  rt.only.low.n = TRUE)

#RAMClust has condensed 31787 features into 1172 spectra 

source("server.functions.R")
## make post-clustering joining function for peaks which appear to have been split into two artificially

RC <- rc.merge.split.clusters(ramclustObj = RC)

#Original cluster number = 1172 
#New cluster number = 609


save(RC, file = "datasets0.20/RCobject_2_QCpostclusteringTICnorm.Rdata")

RC <- RAMClustR::rc.ramclustr(
  ramclustObj = RC,
  st = NULL,
  sr = 0.5,
  maxt = NULL,
  deepSplit = FALSE,
  blocksize = 2000,
  mult = 5,
  hmax = NULL,
  collapse = TRUE,
  minModuleSize = 1,
  linkage = "average",
  cor.method = "pearson",
  rt.only.low.n = TRUE)

#RAMClust has condensed 31787 features into 21778 spectra 

save(RC, file = "datasets0.20/RCobject_2_QCpostclusteringTICnormminmodulesize1.Rdata")

#everything above jackie did for me


source("server.functions.R")
## make post-clustering joining function for peaks which appear to have been split into two artificially

RC <- rc.merge.split.clusters(ramclustObj = RC)

#Original cluster number = 21778 
#New cluster number = 20766

save(RC, file = "datasets0.20/RCobject_2_QCpostclusteringTICnormMergeClustersminmodsize1.Rdata")

MSdata_TIC_normFeat_clust_Merge_minMod1<- t(RC$MSdata)
write.csv(MSdata_TIC_normFeat_clust_Merge_minMod1, file= "MSdata_TIC_normFeat_clust_MergeminMod1.csv")
write.csv(RC$fmz, file= "mz__TIC_normFeatminMod1.csv")
write.csv(RC$frt, file= "RT__TIC_normFeatminMod1.csv") # in HILIC dir


RAMClustR::exportDataset(ramclustObj=RC, which.data="SpecAbund") #rename manually in folder, in datasets as SpecAbundTICNormwithQC

# create 'qc' plots post normalization for comparison, remove QC
RC <- RAMClustR::rc.qc(
  ramclustObj = RC,
  outfile.basename = "post_Clustering",
  qc.tag = qc.tag,
  remove.qc = TRUE
)

save(RC, file = "datasets0.20/RCobject_2_NoQC_postclustering.Rdata")
dim(RC$SpecAbund) #spec abund matrix of normalized feature intensity and compound.

## Export SpecAbund Dataset - contains quantitative signal intensity values
dir.create("datasets")
RAMClustR::exportDataset(ramclustObj=RC, which.data="SpecAbund") #this is pre-fixing the split clusters, so if you wantt this one then manually chagne the file name before doing the next step


## make post-clustering joining function for peaks which appear to have been split into two artificially WITH QC
source("server.functions.R")
load("datasets0.20/RCobject_2_QCpostclustering.Rdata")
RC <- rc.merge.split.clusters(ramclustObj = RC)


#export the new SpecAbund with collapsed clsuters, in datasets named: SpecAbundTICNormwithQCafterMerge
RAMClustR::exportDataset(ramclustObj=RC, which.data="SpecAbund")

#export the .MSP for manual MS1 feature matching to MSMS/MS1 feature matching

rc.export.msp.rc(ramclustObj = RC, one.file = TRUE, mzdec=4) #called it; CCREM_newJackie_keepUnknowns.rc
save(RC, file = "datasets0.20/RCobject_2_QCpostclusteringTICnormMergeClustersminMod1.Rdata")

###############################################
## Section C: FindMain feature annotation
###############################################


library(InterpretMSSpectrum) #make sure this is the most updated version
source("do.findmain.R")

#val starting here 08142023

load("datasets0.20/RCobject_2_QCpostclusteringTICnormMergeClustersminMod1.Rdata")

RC <- do.findmain(
  ramclustObj = RC,
  mode = "positive",
  mzabs.error = 0.005,
  ppm.error = 10,
  ads = NULL,
  nls = NULL,
  scoring = "auto",
  plot.findmain = FALSE,
  writeMat = TRUE,
  writeMS = FALSE
)

save(RC, file = "datasets/RCobject_4_QCpostclusteringTICnormMergeClusters.Rdata")

load("datasets/RCobject_4_QCpostclusteringTICnormMergeClusters.Rdata")

#this will output folders that we will use as input to MSFinder

###############################################
## Section E (part2): use MSFINDER 
###############################################

# complete outside of R


