library(readxl)
library(writexl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
library(vegan)
library(pheatmap)
library(Maaslin2)
library(RColorBrewer)

setwd("/Users/mcgivern.9/Desktop/Projects/agribiome/2023_Val_CoverCrops/analyses/")


# read in geTMM file
getmm=read_excel("data_from_val/geTMM_norm.counts.rpk_edger_genes_getmms_ge3counts_3samples_withannoNEW.xlsx")

# split up this file into genes and annotations
genes=getmm[,1:43]
annos=getmm[,c(2,43:69)]
#update annos with K00371 tree calls
annos2=annos%>%mutate(ko_id=ifelse(gene=="C_No_T0_bin.34.renamed_C_No_T0_bin.34_scaffold_176_1","nxrB",
                            ifelse(gene=="L_E1_T20_COA1_bin.17.renamed_L_E1_T20_COA1_bin.17_scaffold_18132_9","narY",
                                   ifelse(gene=="CerealRye_bin.12_CerealRye_bin.12_k121_2186179_1","narH",
                                          ifelse(gene=="C_No_T20_B_bin.39.renamed_C_No_T20_B_bin.39_scaffold_18357_5","nxrB",
                                                 ifelse(gene=="C_No_T20_C_bin.46.renamed_C_No_T20_C_bin.46_scaffold_1621_13","narY",
                                                        ifelse(gene=="C_No_T20_C_bin.140.renamed_C_No_T20_C_bin.140_scaffold_4078_12","narH",
                                                        ifelse(gene=="HairyVetch_bin.17_HairyVetch_bin.17_k121_3981816_20","narH",
                                                               ifelse(gene=="L_E1_T20_C_bin.131.renamed_L_E1_T20_C_bin.131_scaffold_15570_4","narY",
                                                                      ifelse(gene=="C_No_T20_A_bin.54.renamed_C_No_T20_A_bin.54_scaffold_137_23","narH",
                                                                             ko_id))))))))))
####################
## part 1. key mags
####################

# convert geTMM from gene table to MAG table (using sum geTMM per mag for mags with 10+ genes)
mag_sum_10genes=genes%>%select(-Row.names)%>%gather(-gene,-fasta,key="sample",value="geTMM")%>%
  filter(geTMM>0)%>%group_by(fasta,sample)%>%summarise(mag_sum=sum(geTMM),n=n())%>%filter(n>=10)

## using geTMM
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]
#write_xlsx(wide_mag_geTMM,"mag_metaT_sum.geTMM_10genes.xlsx")

## using abundance
mag_total=mag_sum_10genes%>%ungroup()%>%group_by(sample)%>%summarise(total=sum(mag_sum))
mag_sum_10genes_relabun=mag_sum_10genes%>%left_join(.,mag_total,by="sample")%>%ungroup()%>%group_by(sample,fasta)%>%summarise(abund=mag_sum/total)
wide_mag_abundance=mag_sum_10genes_relabun%>%spread(key=sample,value=abund)
wide_mag_abundance[is.na(wide_mag_abundance)]<-0
wide_mag_abund_mat=as.data.frame(wide_mag_abundance)
row.names(wide_mag_abund_mat)=wide_mag_abund_mat[,1]
wide_mag_abund_mat=wide_mag_abund_mat[,-1]
#write_xlsx(wide_mag_abundance,"mag_metaT_abund_10genes.xlsx")

## read in gtdb
gtdb=read_excel("data_from_val/MAG_RelAbund_taxonomy_allMAGsinDatabaseNEW.xlsx",sheet="taxonomy")

## plotting area plot of MAGs in metaT
mag_ave_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),sep="_")%>%group_by(fasta,treat,time)%>%summarise(ave=mean(mag_sum))
mag_ave_geTMM_tot=mag_ave_geTMM%>%ungroup()%>%group_by(treat,time)%>%summarise(sum=sum(ave))
mag_ave_geTMM%>%left_join(.,mag_ave_geTMM_tot,by=c("treat","time"))%>%ungroup()%>%group_by(fasta,treat,time)%>%summarise(ave_abun=ave/sum)%>%mutate(time=ifelse(time=="T0",0,ifelse(time=="T5",5,21)))%>%
  ggplot()+
  geom_area(aes(x=time,y=ave_abun,fill=fasta),show.legend = FALSE)+
  theme_classic()+
  facet_grid(~treat)

## plotting area plot of genera in metaT
genus_ave_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),sep="_")%>%left_join(.,gtdb,by=c("fasta"="MAG"))%>%separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%group_by(genus,treat,rep,time)%>%summarise(g_sum=sum(mag_sum))%>%ungroup()%>%group_by(genus,treat,time)%>%summarise(g_ave=mean(g_sum))
genus_ave_geTMM_tot=genus_ave_geTMM%>%ungroup()%>%group_by(treat,time)%>%summarise(sum=sum(g_ave))

genus_ave_abund=genus_ave_geTMM%>%left_join(.,genus_ave_geTMM_tot,by=c("treat","time"))%>%ungroup()%>%group_by(genus,treat,time)%>%summarise(ave_abun=g_ave/sum)
gen=genus_ave_abund%>%mutate(plot_abun=ifelse(ave_abun>=0.01,genus,"rare"))%>%ungroup()%>%group_by(plot_abun,treat,time)%>%summarise(abun=sum(ave_abun))%>%ungroup()%>%select(plot_abun)%>%distinct()

#only genera >1%
genus_ave_abund%>%mutate(plot_abun=ifelse(ave_abun>=0.01,genus,"rare"))%>%ungroup()%>%group_by(plot_abun,treat,time)%>%summarise(abun=sum(ave_abun))%>%mutate(time=ifelse(time=="T0",0,ifelse(time=="T5",5,21)))%>%
  ggplot()+
  geom_bar(aes(x=ordered(as.character(time),levels=c("0","5","21")),y=abun*100,fill=plot_abun),color="black",position="stack",stat="identity",show.legend = TRUE)+
  theme_classic()+
  scale_fill_manual(values=genus_colors)+
  facet_grid(~ordered(treat,levels=c("ControlSoil","CerealRye","HairyVetch","Rapeseed","Sorghum")))+ylab("metaT relative abundance (%)")+xlab("Day")+theme(legend.position="bottom",legend.text = element_text(size=2))+guides(fill=guide_legend(ncol=7))

g=genus_ave_geTMM%>%ungroup()%>%select(genus)%>%distinct()
g%>%separate(genus,into=c("d","p","c","o","f","g"),sep=";")%>%group_by(p)%>%summarise(n=n())

genus_colors=c("#006d2c","#2ca25f","#66c2a4","#b2e2e2",#p__Thermoproteota
"#e36c44","#787060","#dc3e25","#e9ccb5","#df375e","#eabc38","#d25e7e","#c09836",#p__Actinomycetota
"#80b6b5","#458ee8","#3fe7c0","#517ec2","#759ee5","#68e2d9","#5f86b7","#b4e3e0","#3a9ac9",#p__Bacteroidota
"#810f7c","#8856a7", #p__Bdellovibrionota
"#78c679","#c2e699",#p__Chloroflexota
"#fb9a99",#p__Desulfobacterota_B
"#980043",#p__Nitrospirota
"#fe9929",#p__Planctomycetota;
"#9c52a7","#e0c1d3","#734ce7","#6e595f","#a845e4","#7b656b","#ce33cf","#a192a9",#p__Pseudomonadota
"black")#rare

# all of them
genus_ave_geTMM%>%left_join(.,genus_ave_geTMM_tot,by=c("treat","time"))%>%ungroup()%>%group_by(genus,treat,time)%>%summarise(ave_abun=g_ave/sum)%>%mutate(time=ifelse(time=="T0",0,ifelse(time=="T5",5,21)))%>%
  ggplot()+
  geom_bar(aes(x=ordered(as.character(time),levels=c("0","5","21")),y=ave_abun*100,fill=genus),color="black",position="stack",stat="identity",show.legend = TRUE)+
  theme_classic()+
  scale_fill_manual(values=genus_colors)+
  facet_grid(~ordered(treat,levels=c("ControlSoil","CerealRye","HairyVetch","Rapeseed","Sorghum")))+ylab("metaT relative abundance (%)")+xlab("Day")+theme(legend.position="bottom",legend.text = element_text(size=2))+guides(fill=guide_legend(ncol=7))

#######
## maaslin2 between exudate and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%select(fasta,sample,mag_sum)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]

t_wide_mag_getmm_mat=t(wide_mag_getmm_mat)

#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_mag_getmm_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_mag_getmm_mat,
         input_metadata = metadata,
         output="maaslin_mags_groups_d5_11.7.23",
         fixed_effects = c("group"))

#plotting significant ones
mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%
  filter(fasta=="C_No_T20_C_bin.131.renamed"|fasta=="C_No_T20_C_bin.35.renamed"|fasta=="Rapeseed_bin.53"|fasta=="CerealRye_bin.5"|fasta=="HairyVetch_bin.34"|fasta=="coassembly2_bin.227"|fasta=="C_No_T0_bin.34.renamed"|fasta=="C_No_T20_C_bin.15.renamed"|fasta=="coassembly2_bin.208")%>%left_join(.,gtdb,by=c("fasta"="MAG"))%>%separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%mutate(genus=paste(f,g,sep=";"))%>%select(-d,-p,-c,-o,-f,-g,-s)%>%
  ggplot()+
  geom_jitter(aes(x=treat,y=mag_sum,color=genus))+
  geom_boxplot(aes(x=treat,y=mag_sum),fill="NA")+
  facet_wrap(~fasta,scales="free")+
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("treatment")+ylab("MAG geTMM")

#######
## maaslin2 between exudate and controls at day 21
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T21")%>%select(fasta,sample,mag_sum)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]

t_wide_mag_getmm_mat=t(wide_mag_getmm_mat)

#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_mag_getmm_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_mag_getmm_mat,
         input_metadata = metadata,
         output="maaslin_mags_groups_d21",
         fixed_effects = c("group"))

#######
## maaslin2 between exudate and controls at day 0
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T0")%>%select(fasta,sample,mag_sum)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]

t_wide_mag_getmm_mat=t(wide_mag_getmm_mat)

#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_mag_getmm_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_mag_getmm_mat,
         input_metadata = metadata,
         output="maaslin_mags_groups_d0",
         fixed_effects = c("group"))


#######
## maaslin2 between cereal rye and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="ControlSoil"|treat=="CerealRye")%>%select(fasta,sample,mag_sum)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]

t_wide_mag_getmm_mat=t(wide_mag_getmm_mat)

#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_mag_getmm_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_mag_getmm_mat,
         input_metadata = metadata,
         output="maaslin_mags_CvsCR_d5",
         fixed_effects = c("group"))

#######
## maaslin2 between hairy vetch and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="ControlSoil"|treat=="HairyVetch")%>%select(fasta,sample,mag_sum)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]

t_wide_mag_getmm_mat=t(wide_mag_getmm_mat)

#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_mag_getmm_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_mag_getmm_mat,
         input_metadata = metadata,
         output="maaslin_mags_CvsHV_d5",
         fixed_effects = c("group"))


#######
## maaslin2 between rapeseed and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="ControlSoil"|treat=="Rapeseed")%>%select(fasta,sample,mag_sum)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]

t_wide_mag_getmm_mat=t(wide_mag_getmm_mat)

#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_mag_getmm_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_mag_getmm_mat,
         input_metadata = metadata,
         output="maaslin_mags_CvsR_d5",
         fixed_effects = c("group"))

#######
## maaslin2 between sorghum and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_mag_geTMM=mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="ControlSoil"|treat=="Sorghum")%>%select(fasta,sample,mag_sum)%>%spread(key=sample,value=mag_sum)
wide_mag_geTMM[is.na(wide_mag_geTMM)]<-0
wide_mag_getmm_mat=as.data.frame(wide_mag_geTMM)
row.names(wide_mag_getmm_mat)=wide_mag_getmm_mat[,1]
wide_mag_getmm_mat=wide_mag_getmm_mat[,-1]

t_wide_mag_getmm_mat=t(wide_mag_getmm_mat)

#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_mag_getmm_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_mag_getmm_mat,
         input_metadata = metadata,
         output="maaslin_mags_CvsS_d5",
         fixed_effects = c("group"))

#plotting significant ones from sorghum and hairy vetch
mag_sum_10genes%>%select(-n)%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%
  filter(fasta=="C_No_T20_B_bin.90.renamed"|fasta=="CerealRye_bin.3")%>%
  ggplot()+
  geom_jitter(aes(x=treat,y=mag_sum))+
  geom_boxplot(aes(x=treat,y=mag_sum),fill="NA")+
  facet_wrap(~fasta,scales="free")+
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####################
## part 2. key functions
####################

## making annotations file for DRAM - using genes in mags with 10+ genes, that appear in at least one replicate
tidy_gene=genes%>%select(-Row.names)%>%gather(-gene,-fasta,key="sample",value="geTMM")%>%
  filter(geTMM>0)
genes_per_sample=mag_sum_10genes_relabun%>%left_join(.,tidy_gene,by=c("fasta","sample"))%>%separate(sample,into=c("treat","rep","time"),sep="_")%>%mutate(group=paste(treat,time,sep="_"))%>%filter(geTMM>0)%>%select(gene,group)%>%distinct()
new_annos=genes_per_sample%>%left_join(.,annos,by="gene")%>%select(-fasta)
#write.table(new_annos,"treat_time_annotations.txt",sep="\t",quote=FALSE,na="",row.names = FALSE)

## making function table for NMDS and heatmap
fx_getmm_tidy=mag_sum_10genes_relabun%>%left_join(.,tidy_gene,by=c("fasta","sample"))%>%left_join(.,annos2,by=c("gene","fasta"))%>%
  mutate(final_fx=ifelse(is.na(cazy_best_hit)==0,cazy_best_hit,
          ifelse(is.na(peptidase_family)==0,peptidase_family,ko_id)))%>%
  mutate(final_fx=str_replace_all(final_fx,".hmm",""))%>%separate(final_fx,into=c("fx","junk"),sep="_")%>%
  ungroup()%>%group_by(sample,fx)%>%summarise(fx_total=sum(geTMM))%>%filter(fx!="")
wide_fx=fx_getmm_tidy%>%spread(key=sample,value=fx_total)
wide_fx[is.na(wide_fx)]<-0
wide_fx_mat=as.data.frame(wide_fx)
row.names(wide_fx_mat)=wide_fx_mat[,1]
wide_fx_mat=wide_fx_mat[,-1]
#write_xlsx(wide_fx,"metaT_fx_table.xlsx")

# NMDS at function id level
t_data=t(wide_fx_mat)
set.seed(13)
NMDS_Bray_data <-metaMDS(t_data, distance = "bray", k=2,
          autotransform = FALSE, noshare = 0.1, trace = 1,trymax=100)
bray_dist = metaMDSdist(t_data, distance = "bray", k=2,
         autotransform = FALSE, noshare = 0.1, trace = 1)
# stress 0.06
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data, display="sites"))
# bare plot
ggplot(ord.data)+
  geom_point(aes(x=NMDS1,y=NMDS2)) +
  theme_classic()

#adding factors
ord.data$sample=row.names(ord.data)
ord.data=ord.data%>%separate(sample,into=c("treat","rep","time"))
ord.data$group=paste(ord.data$treat,ord.data$time,sep="_")
ggplot(ord.data)+
  geom_point(aes(x=NMDS1,y=NMDS2,color=treat),size=3,alpha=0.6) +
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=time))+
  scale_color_manual(values=c("#255068","#cfcfcf","#ebab51","#b572a6","#5e894f","#431a54","#fde61e","#25838f"))+
  theme_classic()
adonis2(t_data~ord.data$treat*ord.data$time, dist="bray",perm=999)


# doing on just fx ids in dram distillate
genome_summary=read.delim("genome_summary_form_wCAMPER.tsv",sep="\t")
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
#write_xlsx(wide_dram_fx,"fx_id_dram_metaT_sum.geTMM.xlsx")

# NMDS at function id level (only ids in dram distillate)
t_data=t(wide_dram_fx_mat)
set.seed(13)
NMDS_Bray_data <-metaMDS(t_data, distance = "bray", k=2,
          autotransform = FALSE, noshare = 0.1, trace = 1,trymax=100)
bray_dist = metaMDSdist(t_data, distance = "bray", k=2,
         autotransform = FALSE, noshare = 0.1, trace = 1)
# stress 0.06
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data, display="sites"))
# bare plot
ggplot(ord.data)+
  geom_point(aes(x=NMDS1,y=NMDS2)) +
  theme_classic()

#adding factors
ord.data$sample=row.names(ord.data)
ord.data=ord.data%>%separate(sample,into=c("treat","rep","time"))
ord.data$group=paste(ord.data$treat,ord.data$time,sep="_")
ggplot(ord.data)+
  geom_point(aes(x=NMDS1,y=NMDS2,color=treat),size=3,alpha=0.6) +
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=time))+
  scale_color_manual(values=c("#255068","#cfcfcf","#ebab51","#b572a6","#5e894f","#431a54","#fde61e","#25838f"))+
  theme_classic()
adonis2(t_data~ord.data$treat*ord.data$time, dist="bray",perm=999)


#######
## maaslin2 between exudate and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
t_wide_dram_fx_mat=t(wide_dram_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_dram_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_dram_fx_mat,
         input_metadata = metadata,
         output="maaslin_dram.ids_group_T5",
         fixed_effects = c("group"))
pheatmap(log2(t_wide_dram_fx_mat+1),scale="column",cluster_rows = TRUE)

#######
## maaslin2 between exudate and controls at day 21
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T21")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
t_wide_dram_fx_mat=t(wide_dram_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_dram_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_dram_fx_mat,
         input_metadata = metadata,
         output="maaslin_dram.ids_group_T21",
         fixed_effects = c("group"))

#######
## maaslin2 between cereal rye and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="CerealRye"|treat=="ControlSoil")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
t_wide_dram_fx_mat=t(wide_dram_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_dram_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_dram_fx_mat,
         input_metadata = metadata,
         output="maaslin_dram.ids_CRvsC_T5",
         fixed_effects = c("group"))

#######
## maaslin2 between cereal rye and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="HairyVetch"|treat=="ControlSoil")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
t_wide_dram_fx_mat=t(wide_dram_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_dram_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_dram_fx_mat,
         input_metadata = metadata,
         output="maaslin_dram.ids_HVvsC_T5",
         fixed_effects = c("group"))

#######
## maaslin2 between hairy vtech and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="HairyVetch"|treat=="ControlSoil")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
t_wide_dram_fx_mat=t(wide_dram_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_dram_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_dram_fx_mat,
         input_metadata = metadata,
         output="maaslin_dram.ids_HVvsC_T5",
         fixed_effects = c("group"))

#######
## maaslin2 between rapeseed and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="Rapeseed"|treat=="ControlSoil")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
t_wide_dram_fx_mat=t(wide_dram_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_dram_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_dram_fx_mat,
         input_metadata = metadata,
         output="maaslin_dram.ids_RvsC_T5",
         fixed_effects = c("group"))

#######
## maaslin2 between sorghum and controls at day 5
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_dram_fx=fx_getmm_tidy%>%left_join(.,genome_summary,by=c("fx"="gene_id"),relationship = "many-to-many")%>%filter(is.na(header)==0)%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%filter(treat=="Sorghum"|treat=="ControlSoil")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_dram_fx[is.na(wide_dram_fx)]<-0
wide_dram_fx_mat=as.data.frame(wide_dram_fx)
row.names(wide_dram_fx_mat)=wide_dram_fx_mat[,1]
wide_dram_fx_mat=wide_dram_fx_mat[,-1]
t_wide_dram_fx_mat=t(wide_dram_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_dram_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_dram_fx_mat,
         input_metadata = metadata,
         output="maaslin_dram.ids_SvsC_T5",
         fixed_effects = c("group"))


#######
## maaslin2 between exudate and controls at day 5 = all function ids
#######
#Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
wide_fx=fx_getmm_tidy%>%ungroup()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%
  select(sample,fx,fx_total)%>%distinct()%>%spread(key=sample,value=fx_total)
wide_fx[is.na(wide_fx)]<-0
wide_fx_mat=as.data.frame(wide_fx)
row.names(wide_fx_mat)=wide_fx_mat[,1]
wide_fx_mat=wide_fx_mat[,-1]
t_wide_fx_mat=t(wide_fx_mat)
#Metadata file
#Formatted with features as columns and samples as rows.
metadata=row.names(t_wide_fx_mat)
metadata=as.data.frame(metadata)
colnames(metadata)=c("sample")
metadata=metadata%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)
metadata=metadata%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))
row.names(metadata)=metadata$sample
metadata=metadata[,-1]
Maaslin2(input_data = t_wide_fx_mat,
         input_metadata = metadata,
         output="maaslin_ids_nar.nxr_group_T5",
         fixed_effects = c("group"))


# #######
# ## making a heatmap of significant functions between exudates and controlat day 5
# #######
fx_of_interest=read_excel("../day5_maaslin2_functions.xlsx",sheet="maaslin_sig_results")
wide_fx_of_interests=fx_of_interest%>%select(feature)%>%left_join(.,fx_getmm_tidy,by=c("feature"="fx"))%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%
  select(feature,sample,fx_total)%>%spread(key=sample,value=fx_total)
wide_fx_of_interests[is.na(wide_fx_of_interests)]<-0
wide_fx_of_interests=as.data.frame(wide_fx_of_interests)
row.names(wide_fx_of_interests)=wide_fx_of_interests[,1]
wide_fx_of_interests=wide_fx_of_interests[,-1]


wide_fx_of_interests=fx_of_interest%>%select(feature,category,group)%>%left_join(.,fx_getmm_tidy,by=c("feature"="fx"))%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%filter(time=="T5")%>%mutate(row=paste(group,category,feature,sep="_"))%>%
  select(row,sample,fx_total)%>%spread(key=sample,value=fx_total)

annotation=fx_of_interest$category
annotation=as.data.frame(annotation)
row.names(annotation)=paste(fx_of_interest$group,fx_of_interest$category,fx_of_interest$feature,sep="_")

cols <- brewer.pal(12, "Set3")
cols <- colorRampPalette(cols)(25)
mycolors <- cols[1:21]
names(mycolors) <- unique(annotation$annotation)
mycolors <- list(annotation = mycolors)
pheatmap(t(wide_fx_of_interests),scale="column",annotation_col = annotation,annotation_colors=mycolors,cluster_cols = FALSE,labels_col =fx_of_interest$feature)

#####
## getting all N genes from maaslin
#####
n_ids=annos2 %>% filter(ko_id=="K03320"|ko_id=="K15576"|ko_id=="K15577"|ko_id=="K15578"|ko_id=="K01430"|ko_id=="K01725"|ko_id=="K01915"|ko_id=="K04751"|ko_id=="K04752"|ko_id=="K00362"|ko_id=="K00363"|ko_id=="K00368"|ko_id=="K10535"|ko_id=="K10945"|ko_id=="K02448"|ko_id=="nxrB")%>%select(gene,ko_id)%>%left_join(.,tidy_gene,by="gene")%>%separate(sample,into=c("treat","rep","time"),sep="_",remove = FALSE)%>%filter(time=="T5")%>%left_join(.,gtdb,by=c("fasta"="MAG"))%>%separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))%>%select(gene,genus,ko_id,fasta,sample,group,geTMM)%>%group_by(genus,ko_id,sample,group)%>%summarise(total=sum(geTMM))%>%ungroup()%>%group_by(genus,ko_id,group)%>%summarise(ave=mean(total))%>%spread(key=ko_id,value=ave)
#write_xlsx(n_ids,"nitrogen_id_genera.xlsx")

annos2 %>% filter(ko_id=="K03320"|ko_id=="K15576"|ko_id=="K15577"|ko_id=="K15578"|ko_id=="K01430"|ko_id=="K01725"|ko_id=="K01915"|ko_id=="K04751"|ko_id=="K04752"|ko_id=="K00362"|ko_id=="K00363"|ko_id=="K00368"|ko_id=="K10535"|ko_id=="K10945"|ko_id=="K02448"|ko_id=="nxrB")%>%select(gene,ko_id)%>%left_join(.,tidy_gene,by="gene")%>%separate(sample,into=c("treat","rep","time"),sep="_",remove = FALSE)%>%filter(time=="T5")%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))%>%group_by(sample,treat,group,ko_id)%>%summarise(total=sum(geTMM))%>%
  ggplot()+
  geom_jitter(aes(x=ordered(treat,levels=c("ControlSoil","CerealRye","HairyVetch","Rapeseed","Sorghum")),y=total))+
  geom_boxplot(aes(x=ordered(treat,levels=c("ControlSoil","CerealRye","HairyVetch","Rapeseed","Sorghum")),y=total),fill="NA")+
  facet_wrap(~ko_id,scales="free")+
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("treatment")+ylab("geTMM")

annos2 %>% filter(ko_id=="K02298"|ko_id=="K15066")%>%select(gene,ko_id)%>%left_join(.,tidy_gene,by="gene")%>%separate(sample,into=c("treat","rep","time"),sep="_",remove = FALSE)%>%filter(time=="T21")%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))%>%group_by(sample,treat,time,group,ko_id)%>%summarise(total=sum(geTMM))%>%
  ggplot()+
  geom_jitter(aes(x=ordered(treat,levels=c("ControlSoil","CerealRye","HairyVetch","Rapeseed","Sorghum")),y=total))+
  geom_boxplot(aes(x=ordered(treat,levels=c("ControlSoil","CerealRye","HairyVetch","Rapeseed","Sorghum")),y=total),fill="NA")+
  facet_wrap(~ko_id,scales="free_y")+
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("treatment")+ylab("geTMM")

#######
## ACC deaminase
#######
# acdS = K01505
# amoA = K10944 = confirmed by MAG taxonomy these are amoA

acc=read_excel("data_from_val/ACC_clean.xlsx",sheet="ACC")

# good acdS = C_No_T20_C_bin.46.renamed_C_No_T20_C_bin.46_scaffold_57508_3
tidy_gene%>%filter(gene=="C_No_T20_C_bin.46.renamed_C_No_T20_C_bin.46_scaffold_57508_3"|gene=="C_No_T20_C_bin.45.renamed_C_No_T20_C_bin.45_scaffold_3103_8")%>%left_join(.,acc,by=c("sample"="Sample"))%>%mutate(group=ifelse(treat=="ControlSoil","control","exudate"))%>%
  ggplot()+
  geom_jitter(aes(x=time,y=geTMM,color=treat))+
  theme_classic()+
  facet_grid(group~gene)

# mags expressing K00166, K09699 (2-oxobutanoate -> propanoyl-coA)
# mags expressing K01652, K01653, K11258 (2-oxobutanoate -> isoleucine)
tidy_gene%>%filter(fasta=="C_No_T20_C_bin.46.renamed"|fasta=="C_No_T20_C_bin.45.renamed")%>%left_join(.,new_annos,by="gene",relationship = "many-to-many")%>%filter(ko_id=="K00166"|ko_id=="K09699"|ko_id=="K01652"|ko_id=="K01653"|ko_id=="K11258"|ko_id=="K01505"|ko_id=="K00052")%>%select(gene,fasta,sample,geTMM,ko_id)%>%distinct()%>%separate(sample,into=c("treat","rep","time"),remove=FALSE)%>%
  group_by(ko_id,fasta,sample,treat,time,rep)%>%summarise(ko_sum=sum(geTMM))%>%spread(key=ko_id,value=ko_sum)%>%mutate(K01652=ifelse(is.na(K01652)==1,0,K01652))%>%mutate(K01505=ifelse(is.na(K01505)==1,0,K01505))%>%mutate(K00166=ifelse(is.na(K00166)==1,0,K00166))%>%mutate(K09699=ifelse(is.na(K09699)==1,0,K09699))%>%mutate(K01653=ifelse(is.na(K01653)==1,0,K01653))%>%mutate(K00052=ifelse(is.na(K00052)==1,0,K00052))%>%
  ggplot()+
  geom_point(aes(x=K01505,y=K00052,color=treat))+
  theme_classic()+
  facet_wrap(~fasta,scales = "free_y")

tidy_gene%>%filter(fasta=="C_No_T20_C_bin.46.renamed"|fasta=="C_No_T20_C_bin.45.renamed")%>%left_join(.,new_annos,by="gene",relationship = "many-to-many")%>%filter()












