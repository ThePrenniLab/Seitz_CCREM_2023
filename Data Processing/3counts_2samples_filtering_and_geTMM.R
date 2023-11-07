library(writexl)
library(readxl)
library(tidyr)
library(dplyr)
library(edgeR)

#RAN ON OLD FEATURE COUNT FILES

# read in featurecounts table and gene key
counts = read.delim("metaT_output_feature_counts_revstranded_07292023_95perc_NEW_NEW.out.txt",sep="\t",header = TRUE) #this file has one line of sample header
stats=read.delim("metaT_output_feature_counts_revstranded_07292023_95perc_NEW_NEW.out.summary.txt",header = TRUE,sep="\t")
row.names(counts)=counts[,1]
counts=counts[,-1]

gene_info=counts[,1:5] # get the first 5 column of gene information
gene_info$name=row.names(gene_info)
counts=counts[,6:45] # change these numbers based on your sample numbers

# moving on to filtering
colSums(counts)
counts$sum=rowSums(counts) # make sum column
range(counts$sum) # range 0 to 3954780

# filter out genes that are 0 across all samples
filtered_counts=counts %>%
  filter(sum>0)
# check to make sure range starts at 1 now
range(filtered_counts$sum) # range 1 to 3954780

#  making read counts less than 3 = 0
counts3 = filtered_counts
counts3[counts3==1|counts3==2] <-0

# convert to presence/absence so we can do the sample number filtering
pa_counts=ifelse(counts3>0,1,0)
pa_counts=as.data.frame(pa_counts)
pa_counts=pa_counts[,-41] # get rid of old sum column in filtered counts - change this number to be the column number of your sum column
pa_counts$n=rowSums(pa_counts)

# just has to be in 2 samples = change n>=X to be X samples
n2 = pa_counts %>%
  filter(n>=2)

x = row.names(n2)
x=as.data.frame(x)
x$n=n2$n
filtered_counts$name=row.names(filtered_counts)

filtered_1x = inner_join(filtered_counts,x,by=c("name"="x"))
colnames(filtered_1x)
# move name column to the front and then have the sample columns
filtered_1x = filtered_1x[,c(42,1:40)]

write_xlsx(filtered_1x,"filtered_3counts_2samples.xlsx")




# convert to geTMM
count=filtered_1x
count=left_join(count,gene_info,by="name")%>%select(-Strand,-Chr,-Start,-End)
rpk = ((count[,2:41]*10^3)/count$Length) # change count[,2:X] to be your column numbers
row.names(rpk)=count$name
group <- c(rep("A",ncol(rpk)))

rpk.norm <- DGEList(counts=rpk,group=group)

# get library sizes
colSums(stats[,-1])
# needs to be in same order as your columns in filtered_1x
lib_size=c(27036356,34584227,34265674,34281102,21578297,106533045, 29068247, 41321554, 50373307, 36788973, 66649744, 78392547,20379447, 37459801, 29070161, 35559341, 77804219, 28113329, 38098917, 35682822, 17246911, 75682920, 23220604, 30713441, 29141994, 18188316, 33917706, 113864775, 18753895, 30674000, 32073439, 41229002, 81953255, 25179485, 35780983, 30051844, 21509705, 88037893, 30717205, 36824969)
lib_size=as.data.frame(lib_size)
lib_size$sample=c(colnames(filtered_1x[,2:41])) # change these numbers to be your columns
row.names(lib_size)=lib_size$sample
lib_size=select(lib_size,-sample)

rpk.norm$samples$lib.size = lib_size[rownames(rpk.norm$samples),]
rpk.norm$samples$lib.size
rpk.norm$samples

rpk.norm <- calcNormFactors(rpk.norm)
getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df) #0 to 101710.1

getmms_df$gene=row.names(getmms)
colnames(getmms_df)
# rearrange columns so gene is first
getmms_df = getmms_df[,c(41,1:40)]
colnames(getmms_df)
write_xlsx(getmms_df, 'getmms_ge3counts_2samples.xlsx')



########################################
## adding kayla's code

########################################

#join in annotations
anno=read.delim("../annotations.tsv") #from DRAM
row.names(anno)=anno[,1]
anno=anno[,-1]
getmms_df_withanno=merge(getmms_df,anno, by = 'row.names')
write_xlsx(getmms_df_withanno,'geTMM_norm.counts.rpk_edger_genes_getmms_ge3counts_2samples_withannoNEW.xlsx')

#now need to collapse by genome
library(dplyr)

# make row names its own column so they can be extractable by tibble::rownames_to_column
#row_names <- getmms_df$Row.names
#getmms_df <- getmms_df %>% select(-Row.names)
#rownames(getmms_df) <- row_names

# make a new df with a duplicated column which we will make into a genome column so that we can sum across just the genome
norm.counts.rpk_edger.zerosremoved_for_bins <- as.data.frame(getmms_df) #make DF
norm.counts.rpk_edger.zerosremoved_for_bins <- tibble::rownames_to_column(norm.counts.rpk_edger.zerosremoved_for_bins, "rn") # make genes the row names, call it rn
norm.counts.rpk_edger.zerosremoved_for_bins$rn<-gsub('_k121(.*)','',norm.counts.rpk_edger.zerosremoved_for_bins$rn) #remove anything with _k121* and replace with nothing
norm.counts.rpk_edger.zerosremoved_for_bins$rn<-gsub('_scaffold(.*)','',norm.counts.rpk_edger.zerosremoved_for_bins$rn) #remove anything with _scaffold* and replace with nothing
norm.counts.rpk_edger.zerosremoved_for_bins_nogeneNames <- norm.counts.rpk_edger.zerosremoved_for_bins %>%
  select(-gene)

#collapse rows and get the mean gene count for that genome across samples

norm.counts.rpk_edger.bins_mean<-norm.counts.rpk_edger.zerosremoved_for_bins_nogeneNames%>%
  group_by(rn)%>%
  summarise_all(list(mean=mean, mean = ~ mean(.x, na.rm = TRUE)))

norm.counts.rpk_edger.bins_mean<-norm.counts.rpk_edger.zerosremoved_for_bins_nogeneNames%>%
  group_by(rn)%>%
  summarise_all(list(mean=mean, mean = ~ mean(.x, na.rm = TRUE)))

norm.counts.rpk_edger.bins_mean<-norm.counts.rpk_edger.zerosremoved_for_bins_nogeneNames%>%
  group_by(rn)%>%
  summarise_all(funs(count))#count num genes ON, filter later for genomes that have at least 20 ON to be considered doing something

norm.counts.rpk_edger.bins_sum<-norm.counts.rpk_edger.zerosremoved_for_bins_nogeneNames%>%
  group_by(rn)%>%
  summarise_all(list(sum=sum))

norm.counts.rpk_edger.bins_median<-norm.counts.rpk_edger.zerosremoved_for_bins_nogeneNames%>%
  group_by(rn)%>%
  summarise_all(list(median=median))


write.csv(norm.counts.rpk_edger.bins_mean,"norm.counts.rpk_edger.bins_mean.csv")



