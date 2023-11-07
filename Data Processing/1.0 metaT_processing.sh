#metaT processing


################################# 
#### Trimming and filtering ####
################# ################

#!/bin/bash
for element in $(<$1)
do
#unzip
gunzip "$element".fastq.gz
# Trim adapters and quality trim
bbduk.sh in=../Raw_Data/"$element".fastq out="$element"_trimmed.fastq ktrim=r k=23 Mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 Maq=10 t=7 ref=/opt/bbtools/bbmap/resources/adapters.fa 

# Zip reads for rqcfilter
gzip -c "$element"_trimmed.fastq > "$element"_trimmed.fastq.gz 

# Filter reads using rqcfilter
rqcfilter2.sh jni=t in="$element"_trimmed.fastq.gz path="$element"_filter  rqcfilterdata=/home/opt/RQCFilterData rna=t trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t removemicrobes=t sketch kapa=t clumpify=t barcodefilter=f trimpolyg=5 usejni=f
done

################# 
#### Mapping ####
################# 

#!/bin/bash
for element in $(<$1)
do
# Unzip reads into mapping dir
gunzip -c /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/trimmed_reads/"$element"_filter/"$element"_trimmed.anqrpht.fastq.gz > /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/"$element"_trimmed.anqrpht.fastq 

# Deinterleave trimmed and filtered fastqs
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping
/ORG-Data/scripts/deinterleave_fastq.sh < "$element"_trimmed.anqrpht.fastq "$element"_trimmed.anqrpht_R1.fastq "$element"_trimmed.anqrpht_R2.fastq
done

#MAPPING AT 95% ID
#!/bin/bash
for element in $(<$1)
do
#currently using the 95%id bowtie
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB/99perMAG_DB -S "$element"_bowtie.sam -1 "$element"_trimmed.anqrpht_R1.fastq -2 "$element"_trimmed.anqrpht_R2.fastq 
# Sam to bam file
samtools view -@ 30 -bS "$element"_bowtie.sam > "$element"_bowtie.bam 
# Reformat bam file to retain mapped reads at 97% id, the error was here. the in file said in=<sampleid>_bowtie.sorted.bam when it should just be bowtie.bam
reformat.sh in="$element"_bowtie.bam out="$element"_bowtie.reformat95.bam minidfilter=0.95 primaryonly=t pairedonly=f 
# Name-sort bam files
samtools sort -n -o "$element"_bowtie.reformat95_NAMESORT.bam "$element"_bowtie.reformat95.bam -@ 30 
done

################# 
#### Counting ####
################# 

# Feature counts, gff file from dram
#95%
featureCounts -T 20 -t CDS -g ID -s 2 -p -a /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.gff -o metaT_output_feature_counts_revstranded_07132023_95perc.out *_bowtie.reformat95_NAMESORT.bam

# Filtered and transformed data tables
# Remove first line of file before putting into R
sed -i '1d' metaT_output_feature_counts_revstranded_07132023_95perc.out

# Run R script 01_metaT_counts_to_geTMM.R to obtain geTMMs






