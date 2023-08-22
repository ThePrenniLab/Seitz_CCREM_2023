# to do: 

1. delete the unnessary directories (keep: Sequencing_QC_Reports, Raw_Data, Metagenome_Report_Tables). combine metaG reports and seq QCs into one directory as JGI stats, done
2. rename the directories, done
3. create a JGI_QC_Reports and move the Sequencing_QC_Reports and Metagenome_Report_Tables into it. 
4. rename the raw reads as sample name
5. trimming


#Map metaT reads to metaG database

cd /home/projects/Agribiome/RootExudate_Prenni/metaT

# raw reads in: 
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads

################################# 
#### Directory Cleaning 
################# ################

# Delete the unnecessary directories from the metaT files; keep: Sequencing_QC_Reports, Raw_Data, Metagenome_Report_Tables

#make list file for looping: 
use ls >> list.txt to list all the fa directories needed

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"
rm -r IMG_Data
rm -r Filtered_Raw_Data
rm -r QC_and_Genome_Assembly
rm -r QC_Filtered_Raw_Data
rm -r RNASeq_Data
done

#save
clean_metaT_directories_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash clean_metaT_directories_loop.sh list.txt

#save
clean_metaT_directories.sh


#create a JGI_QC_Reports and move the Sequencing_QC_Reports and Metagenome_Report_Tables into it. 

#!/bin/bash
for element in $(<$1)
do
cp -r /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element" /home/projects/Agribiome/RootExudate_Prenni/metaT/JGI_QC_Reports

done
#save
copy_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash copy_loop.sh list.txt

#save
copy.sh

#then remove the Sequencing_QC_Reports and Metagenome_Report_Tables from the raw reads dir

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"
rm -r Sequencing_QC_Reports
rm -r Metagenome_Report_Tables
done

#save
clean_loop1.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash clean_loop1.sh list.txt

#save
clean1.sh

#then remove the raw reads dir from the JGI_QC_Reports dir 
#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/JGI_QC_Reports/"$element"
rm -r Raw_Data

done

#save
clean_loop2.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash clean_loop2.sh list.txt

#save
clean2.sh

#rename raw reads files as sample names

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/Raw_Data
mv *.fastq.gz "$element".fastq.gz

done

#save
rename_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash rename_loop.sh list.txt

#save
rename.sh



################################# 
#### Trimming and filtering ####
################# ################


#!/bin/bash
for element in $(<$1)
do

#unzip
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/Raw_Data
gunzip "$element".fastq.gz

cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"
mkdir trimmed_reads
cd trimmed_reads

# Trim adapters and quality trim
bbduk.sh in=../Raw_Data/"$element".fastq out="$element"_trimmed.fastq ktrim=r k=23 Mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 Maq=10 t=7 ref=/opt/bbtools/bbmap/resources/adapters.fa 

# Zip reads for rqcfilter
gzip -c "$element"_trimmed.fastq > "$element"_trimmed.fastq.gz 

# Filter reads using rqcfilter
rqcfilter2.sh jni=t in="$element"_trimmed.fastq.gz path="$element"_filter  rqcfilterdata=/home/opt/RQCFilterData rna=t trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t removemicrobes=t sketch kapa=t clumpify=t barcodefilter=f trimpolyg=5 usejni=f
done

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash trim_loop.sh list.txt

#save, submitted in metaT/raw_reads/slurm june 26 11:30am, 
trimming.sh


################# 
#### Mapping ####
################# 

cd /home/projects/Agribiome/RootExudate_Prenni/metaT/
mkdir mapping
cd mapping

#!/bin/bash
for element in $(<$1)
do
# Unzip reads into mapping dir
gunzip -c /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/trimmed_reads/"$element"_filter/"$element"_trimmed.anqrpht.fastq.gz > /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/"$element"_trimmed.anqrpht.fastq 

# Deinterleave trimmed and filtered fastqs
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping

/ORG-Data/scripts/deinterleave_fastq.sh < "$element"_trimmed.anqrpht.fastq "$element"_trimmed.anqrpht_R1.fastq "$element"_trimmed.anqrpht_R2.fastq
done

#save
unzip_deinterleave_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash unzip_deinterleave_loop.sh list.txt

#save, submitted in metaT/slurm july 5
unzip_deinterleave.sh


# Mapping all metaT with bowtie to the gff derived from DRAM annotations of MAG database
#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping

#currently using the 95%id bowtie
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB/99perMAG_DB -S "$element"_bowtie.sam -1 "$element"_trimmed.anqrpht_R1.fastq -2 "$element"_trimmed.anqrpht_R2.fastq 


# Sam to bam file
samtools view -@ 30 -bS "$element"_bowtie.sam > "$element"_bowtie.bam 


# Reformat bam file to retain mapped reads at 97% id, the error was here. the in file said in=<sampleid>_bowtie.sorted.bam when it should just be bowtie.bam
reformat.sh in="$element"_bowtie.bam out="$element"_bowtie.reformat97.bam minidfilter=0.97 primaryonly=t pairedonly=f 


# Name, sort bam files
samtools sort -n -o "$element"_bowtie.reformat97_NAMESORT.bam "$element"_bowtie.reformat97.bam -@ 30 

done

#save
mapping_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash mapping_loop.sh list.txt

#save, submitted in metaT/slurm july 6 11am, finished around mignight
mapping.sh

#MAPPING AT 95% ID

# Mapping all metaT with bowtie to the gff derived from DRAM annotations of MAG database
#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping

#currently using the 95%id bowtie
#bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB/99perMAG_DB -S "$element"_bowtie.sam -1 "$element"_trimmed.anqrpht_R1.fastq -2 "$element"_trimmed.anqrpht_R2.fastq 


# Sam to bam file
#samtools view -@ 30 -bS "$element"_bowtie.sam > "$element"_bowtie.bam 


# Reformat bam file to retain mapped reads at 97% id, the error was here. the in file said in=<sampleid>_bowtie.sorted.bam when it should just be bowtie.bam
reformat.sh in="$element"_bowtie.bam out="$element"_bowtie.reformat95.bam minidfilter=0.95 primaryonly=t pairedonly=f 


# Name, sort bam files
samtools sort -n -o "$element"_bowtie.reformat95_NAMESORT.bam "$element"_bowtie.reformat95.bam -@ 30 

done

#save
mapping_loop95%.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash mapping_loop95%.sh list.txt

#save, submitted in metaT/slurm july 11 9am, finished 3pm
mapping95%.sh



################# 
#### Counting ####
################# 

# Feature counts, gff file from dram
#95%
featureCounts -T 20 -t CDS -g ID -s 2 -p -a /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.gff -o metaT_output_feature_counts_revstranded_07132023_95perc.out *_bowtie.reformat95_NAMESORT.bam
#97%
featureCounts -T 20 -t CDS -g ID -s 2 -p -a /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.gff -o metaT_output_feature_counts_revstranded_070723_95perc.out *_bowtie.reformat97_NAMESORT.bam



#Server location of feature counts output metaT_output_feature_counts_revstranded_070722_97perc.out (3991706 lines)
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/metaT_output_feature_counts_revstranded_07132023_95perc.out
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/metaT_output_feature_counts_revstranded_070723_95perc.out


# Filtered and transformed data tables

1. metaT gene table transformed to geTMM (not filtered to genes in 3 or more samples)
# all genes

# Remove first line of file before putting into R
sed -i '1d' metaT_output_feature_counts_revstranded_070723_97perc.out

sed -i '1d' metaT_output_feature_counts_revstranded_07132023_95perc.out

# do both counting, use #2 unless you're looking for something in particular then look at 1

# Now metaT_output_feature_counts_revstranded_07192022_97perc.out is 3991705 lines
# Run R script 01_metaT_counts_to_geTMM.R to obtain geTMMs

# Output on W2 server
/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/USA/geTMM_tables/geTMM_norm.counts.rpk_edger_genes_nofilter_092122.csv  



2. metaT gene table transformed to geTMM filtered to genes in 3 or more samples
#genes that show up in 3+ samples

# Filter feature count output to genes that are in 3 or more samples

with open("./metaT_output_feature_counts_revstranded_07192023_95perc_NEW.out", 'r') as counts:
with open("./metaT_output_feature_counts_revstranded_07192023_95perc_sample_filterd.csv", 'w') as out_data:
    for i, l in enumerate(counts):
        if  i < 2:
            _ = out_data.write(l)
            continue
        if len([j for j in l.split('\t')[6:] if float(j)> 0]) > 5:
            ln = out_data.write(l) 

#Server location of filtered feature counts output metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd.csv (41829 lines)
/home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/metaT_output_feature_counts_revstranded_07192023_95perc_sample_filterd.csv 

#Remove first line of file before putting into R
sed -i '1d' metaT_output_feature_counts_revstranded_07192023_95perc_sample_filterd.csv

# Run R script 01_metaT_counts_to_geTMM.R to obtain geTMMs


#done





# REDO RAPESEED 2 T0 (MISSING FILE) AND RAPESEED 1 T5

# first clean and rename the rapeseed 2 file / directory

################################# 
#### Trimming and filtering ####
################# ################

#!/bin/bash
for element in $(<$1)
do

#unzip
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/Raw_Data
gunzip "$element".fastq.gz

# Trim adapters and quality trim in the trimmed reads dir
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/trimmed_reads

bbduk.sh in=../Raw_Data/"$element".fastq out="$element"_trimmed.fastq ktrim=r k=23 Mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 Maq=10 t=7 ref=/opt/bbtools/bbmap/resources/adapters.fa 

# Zip reads for rqcfilter
gzip -c "$element"_trimmed.fastq > "$element"_trimmed.fastq.gz 

# Filter reads using rqcfilter
rqcfilter2.sh jni=t in="$element"_trimmed.fastq.gz path="$element"_filter  rqcfilterdata=/home/opt/RQCFilterData rna=t trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t removemicrobes=t sketch kapa=t clumpify=t barcodefilter=f trimpolyg=5 usejni=f


################# 
#### Mapping ####
################# 
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping
# Unzip reads into mapping dir
gunzip -c /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/trimmed_reads/"$element"_filter/"$element"_trimmed.anqrpht.fastq.gz > /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/"$element"_trimmed.anqrpht.fastq 

# Deinterleave trimmed and filtered fastqs
/ORG-Data/scripts/deinterleave_fastq.sh < "$element"_trimmed.anqrpht.fastq "$element"_trimmed.anqrpht_R1.fastq "$element"_trimmed.anqrpht_R2.fastq

#currently using the 95%id bowtie
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB/99perMAG_DB -S "$element"_bowtie.sam -1 "$element"_trimmed.anqrpht_R1.fastq -2 "$element"_trimmed.anqrpht_R2.fastq 

# Sam to bam file
samtools view -@ 30 -bS "$element"_bowtie.sam > "$element"_bowtie.bam 

# Reformat bam file to retain mapped reads at 95% id
reformat.sh in="$element"_bowtie.bam out="$element"_bowtie.reformat95.bam minidfilter=0.95 primaryonly=t pairedonly=f 


# Name, sort bam files
samtools sort -n -o "$element"_bowtie.reformat95_NAMESORT.bam "$element"_bowtie.reformat95.bam -@ 30 

done

#save
redo_rapeseed_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash redo_rapeseed_loop.sh list_RS.txt

#save, submitted in metaT/raw_reads/slurm 
redo_rapeseed.sh

#redo the feature counts with rapeseed files
featureCounts -T 20 -t CDS -g ID -s 2 -p -a /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.gff -o metaT_output_feature_counts_revstranded_07192023_95perc_NEW.out *_bowtie.reformat95_NAMESORT.bam


#Server location of feature counts output metaT_output_feature_counts_revstranded_070722_97perc.out (3991706 lines)
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/metaT_output_feature_counts_revstranded_07192023_95perc_NEW.out

#remove first line
sed -i '1d' metaT_output_feature_counts_revstranded_07192023_95perc_NEW.out


2. metaT gene table transformed to geTMM filtered to genes in 2 or more samples

#put this into a python script called sample_count_filter.py in mapping. submit using "python sample_count_filter.py"
#no filter	
with open("./metaT_output_feature_counts_revstranded_07192023_95perc_NEW.out", 'r') as counts:
    with open("./metaT_output_feature_counts_revstranded_07192023_95perc_sample_NOfilter.csv", 'w') as out_data:
        for i, l in enumerate(counts):
            if  i < 0:
                _ = out_data.write(l)
                continue
            if len([j for j in l.split('\t')[6:] if float(j)> 0]) > 5:
                ln = out_data.write(l)
    
#filter to genes in 2+ samples            
with open("./metaT_output_feature_counts_revstranded_07192023_95perc_NEW.out", 'r') as counts:
    with open("./metaT_output_feature_counts_revstranded_07192023_95perc_sample_filtered.csv", 'w') as out_data:
        for i, l in enumerate(counts):
            if  i < 2:
                _ = out_data.write(l)
                continue
            if len([j for j in l.split('\t')[6:] if float(j)> 0]) > 5:
                ln = out_data.write(l)


#Server location of filtered feature counts output metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd.csv (41829 lines)
/home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/metaT_output_feature_counts_revstranded_07192023_95perc_sample_filterd.csv 

#Remove first line of file before putting into R
sed -i '1d' metaT_output_feature_counts_revstranded_07192023_95perc_sample_filterd.csv

# Run R script 01_metaT_counts_to_geTMM.R to obtain geTMMs




# CLEAN UP
1. zip raw reads
2. delete all SAM files and redundant BAM files.
3. delete the 97% mapping dir?


1. zip raw reads

#!/bin/bash
for element in $(<$1)
do

#zip raw reads 
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/Raw_Data
gzip "$element".fastq


done

#save
clean_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash clean_loop.sh list.txt

#!/bin/bash
for element in $(<$1)
do

#zip raw reads 
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/"$element"/Raw_Data
rm "$element".fastq

done

#save
rm_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
bash rm_loop.sh list.txt


2. delete sams
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#delete all sam files
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping
rm *.sam















#redo the sorting steps for RS, this failed last time because the overwrite = false was the default, so i deleted the files and rerunning.
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping

reformat.sh in=Rapeseed_2_T0_bowtie.bam out=Rapeseed_2_T0_bowtie.reformat95.bam minidfilter=0.95 primaryonly=t pairedonly=f     

# Name, sort bam files
samtools sort -n -o Rapeseed_2_T0_bowtie.reformat95_NAMESORT.bam Rapeseed_2_T0_bowtie.reformat95.bam -@ 30





#redo RS_1_t5

################################# 
#### Trimming and filtering ####
################# ################

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#unzip
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/Rapeseed_1_T5/Raw_Data
gunzip Rapeseed_1_T5.fastq.gz

# Trim adapters and quality trim in the trimmed reads dir
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/Rapeseed_1_T5/trimmed_reads

bbduk.sh in=../Raw_Data/Rapeseed_1_T5.fastq out=Rapeseed_1_T5_trimmed.fastq ktrim=r k=23 Mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 Maq=10 t=7 ref=/opt/bbtools/bbmap/resources/adapters.fa 

# Zip reads for rqcfilter
gzip -c Rapeseed_1_T5_trimmed.fastq > Rapeseed_1_T5_trimmed.fastq.gz 

# Filter reads using rqcfilter
rqcfilter2.sh jni=t in=Rapeseed_1_T5_trimmed.fastq.gz path=Rapeseed_1_T5_filter  rqcfilterdata=/home/opt/RQCFilterData rna=t trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t removemicrobes=t sketch kapa=t clumpify=t barcodefilter=f trimpolyg=5 usejni=f


################# 
#### Mapping ####
################# 
cd /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping
# Unzip reads into mapping dir
gunzip -c /home/projects/Agribiome/RootExudate_Prenni/metaT/raw_reads/Rapeseed_1_T5/trimmed_reads/Rapeseed_1_T5_filter/Rapeseed_1_T5_trimmed.anqrpht.fastq.gz > /home/projects/Agribiome/RootExudate_Prenni/metaT/mapping/Rapeseed_1_T5_trimmed.anqrpht.fastq 

# Deinterleave trimmed and filtered fastqs
/ORG-Data/scripts/deinterleave_fastq.sh < Rapeseed_1_T5_trimmed.anqrpht.fastq Rapeseed_1_T5_trimmed.anqrpht_R1.fastq Rapeseed_1_T5_trimmed.anqrpht_R2.fastq

#currently using the 95%id bowtie
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB/99perMAG_DB -S Rapeseed_1_T5_bowtie.sam -1 Rapeseed_1_T5_trimmed.anqrpht_R1.fastq -2 Rapeseed_1_T5_trimmed.anqrpht_R2.fastq 

# Sam to bam file
samtools view -@ 30 -bS Rapeseed_1_T5_bowtie.sam > Rapeseed_1_T5_bowtie.bam 

# Reformat bam file to retain mapped reads at 95% id
reformat.sh in=Rapeseed_1_T5_bowtie.bam out=Rapeseed_1_T5_bowtie.reformat95.bam minidfilter=0.95 primaryonly=t pairedonly=f 


# Name, sort bam files
samtools sort -n -o Rapeseed_1_T5_bowtie.reformat95_NAMESORT.bam Rapeseed_1_T5_bowtie.reformat95.bam -@ 30 

#save metaT/raw_reads/slurm 
redo_rapeseed_1_t5.sh






#run DRAM adjectives because we need the output for the getTMM script
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40 
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

/opt/Miniconda2/miniconda2/envs/scripts/bin/rule_adjectives /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/annotations.tsv adjectives.tsv


#save MetaG/slurm
adjectives.sh
sbatch --dependency=afterok:61068 adjectives.sh











