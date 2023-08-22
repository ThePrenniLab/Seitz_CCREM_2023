###########################################################################################
###########################################################################################
# CCREM MetaG Sample Processing
# Valerie Seitz
# Jan 16 2023 - March X
###########################################################################################
###########################################################################################


# Working Directory
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG

# raw reads deposited here: 
cd /ORG-Data-wrighton/Agribiome/RootExudate_Prenni/MetaG/raw_reads

###########################################################################################
######## Make working directories
mkdir MetaG

mkdir slurm #or call it scripts

#make directories for each sample...
mkdir CerealRye_T3
mkdir CerealRye_T5
# etc.. 

###########################################################################################
######## Make the txt file for looping: 

## IMPORTANT: when you make this file, you must make it in command line (do NOT make the file in word/excel and then transfer onto the server!)
cd slurm
nano

# TYPE (do not copy/paste) your sample names
CerealRye_T10
Control_T0
Control_T7
# etc...

###########################################################################################
######## Make subdirectories within each sample directory:
cd slurm
nano 

#!/bin/bash
#$1 list of names of the paired reads that you want to make dirs.
for element in $(<$1)
do
cd "$element"
mkdir raw_reads
mkdir assembly
mkdir trimmed_reads
cd ../
done

#save
makdir_loop.sh

nano
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash makdir_loop.sh list.txt

#save
mkdir.sh

sbatch mkdir.sh

###########################################################################################
######## Copy the files from ORG to your directory’s raw reads folder
cd slurm
cd nano

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#change to the directory where your reads are located
cd ../raw_reads 

# copy reads from ORG-DATA location 
cp /ORG-Data-wrighton/Agribiome/RootExudate_Prenni/*_S*_L002_R* .

# Rename reads to sample name
mv 1_S6_L002_R1_001.fastq.gz Control_T0_R1.fastq.gz
mv 1_S6_L002_R2_001.fastq.gz Control_T0_R2.fastq.gz
mv 2_S7_L002_R1_001.fastq.gz Control_T3_R1.fastq.gz
mv 2_S7_L002_R2_001.fastq.gz Control_T3_R2.fastq.gz
mv 3_S8_L002_R1_001.fastq.gz Rapeseed_T3_R1.fastq.gz
# and so on, see "MAG Data Processing.xlsx" file for all file name changes.

###########################################################################################
######## Unzip fastq files
cd slurm/
nano

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../raw_reads

gunzip *.fastq.gz 

# Submit the job
sbatch unzip.sh

###########################################################################################
######## FastQC on raw reads

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../raw_reads

fastqc *.fastq #perform QC 

#submited jan 16 at 12:45pm
sbatch fastqc.sh 

###########################################################################################
######## Copy reads from MetaG/raw_reads to individual sample directories

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/raw_reads/CerealRye_T10_R1.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/CerealRye_T10/raw_reads
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/raw_reads/Control_T0_R1.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/Control_T0/raw_reads
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/raw_reads/Control_T0_R2.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/Control_T0/raw_reads
# ...See "MAG Data Processing" for the rest of the "copy" code.

# This is how it should have been done:
# Remove the raw data files from the original transfer to ../MetaG/raw_reads

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
rm /home/projects/Agribiome/RootExudate_Prenni/MetaG/raw_reads/*.fastq

###########################################################################################
######## Sickle Trimming on raw reads

##  I did NOT do this step (but this is what the looping would look like)
### I individually did sickle, see "sickle_info" for that code and output

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/raw_reads

sickle pe -f "$element"_R1.fastq  -r "$element"_R2.fastq   -t sanger -o "$element"_R1_trimmed.fastq -p "$element"_R2_trimmed.fastq -s discared_R1R2.fastq
rm R1R2_singles.fastq
mv *trimmed.fastq ../trimmed_reads
../../
done

#save
sickle_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash sickle_loop.sh list.txt

#save
sickle.sh

###########################################################################################
######## FastQC on trimmed reads

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/trimmed_reads

#rerun fastqc
fastqc "$element"_R1_trimmed.fastq "$element"_R2_trimmed.fastq

cd ../../
done

#save
trimmed_fastqc_loop.sh


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash trimmed_fastqc_loop.sh list.txt

#save and submit
trimmed_fastqc.sh

sbatch trimmed_fastqc.sh


###########################################################################################
###### Megahit Assembly

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/assembly

mkdir megahit
cd megahit

megahit -1 ../../trimmed_reads/"$element"_R1_trimmed.fastq -2 ../../trimmed_reads/"$element"_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 -m 0.4  -t 50

cd ../../../
done

#save
megahit_loop.sh

nano

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=28-00:00:00 #changed this to 28 
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash megahit_loop.sh list.txt

#save and submit - job restarted on feb 21 2023 with increased time limit to 28d, increased tasks too. finished marc 3 5pm
megahit.sh

###########################################################################################
###### Megahit Assembly Stats

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/assembly/megahit/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o "$element"_final.contigs_STATS
cd ../../../
done

#save
assembly_stats_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash assembly_stats_loop.sh list.txt

#save
assembly_stats.sh

sbatch assembly_stats.sh #submitted march 3 8pm, finished 10pm

# add assembly stats to the file: assembly_stats.xlsx file, 
## found out that two samples did not assembly so reruning those individually (see end of page for that code)

###########################################################################################
######## BINNING

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/assembly/megahit/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o "$element"_final.contigs_2500.fa

# map to these scaffolds to get sam file

bbmap.sh -Xmx48G threads=20 overwrite=t ref="$element"_final.contigs_2500.fa in1=../../../trimmed_reads/"$element"_R1_trimmed.fastq in2=../../../trimmed_reads/"$element"_R2_trimmed.fastq out="$element"_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS "$element"_final.contigs_2500_mapped.sam > "$element"_final.contigs_2500_mapped.bam
samtools sort -T "$element"_final.contigs_2500_mapped.sorted -o "$element"_final.contigs_2500_mapped.sorted.bam "$element"_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in="$element"_final.contigs_2500_mapped.sorted.bam out="$element"_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

#bin
runMetaBat.sh "$element"_final.contigs_2500.fa "$element"_final.contigs_2500_mapped99per.sorted.bam

cd ../../../../
done

#save
map_bin_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

bash map_bin_loop.sh list.txt

sbatch map_bin.sh #submitted march 6 9am, finished march 7 5am

###########################################################################################
######## QC on bins 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/assembly/megahit/megahit_out/"$element"_final.contigs_2500.fa.metabat-bins
checkm lineage_wf -t 20 -x fa . checkm

cd checkm
checkm qa -o 2 -f results.txt --tab_table -t 20 lineage.ms .

cd ../../../../../
done

checkm_loop.sh 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

bash checkm_loop.sh list.txt

# save, submitted march 7 10 am, finished: 1pm
checkm.sh
sbatch checkm.sh


###########################################################################################
####### Evaluate binning results (manual step): 

# use this to go directly to extracting the M and HQ bins
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/CerealRye_T10/assembly/megahit/megahit_out/CerealRye_T10_final.contigs_2500.fa.metabat-bins/checkm/results.txt

# Record good bins in new excel file: "MAG processing manual code"

#EXAMPLE (SEE MAG processing manual code FOR ALL MANUAL CODE)
# make MAG directories in each sample folder
mkdir MAGs
# Copy the HQ/MQ MAGs to your MAGs folder (make the MAG directory in the same sample directory. potential path would be: ../Control_T3/MAGs)
cd MAGs
# Copy MQ + HQ MAGs to this directory (while in MAGs directory):
cp ../assembly/megahit/megahit_out/CerealRye_T10_final.contigs_2500.fa.metabat-bins/bin.4.fa .
cp ../assembly/megahit/megahit_out/CerealRye_T10_final.contigs_2500.fa.metabat-bins/bin.6.fa .
# Rename the bins
mv bin.4.fa CerealRye_T10_bin.4.fa 
mv bin.6.fa CerealRye_T10_bin.6.fa 

###########################################################################################
######### Evaluating MAG taxonomy

# here is where you can rename the gtdb with at least the date ran and version (if possible)

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs

#this will create SAM summary files for the taxonomy assigned for the bins
gtdbtk classify_wf -x fa --genome_dir . --out_dir gtdb_03102023 --cpus 20

cd ../../
done

#Exit, rename, and run: 
sbatch gtdb_loop.sh 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash gtdb_loop.sh list.txt

#save, submitted march 10 at noon. job ID: 58002, done around 5
gtdb.sh

###########################################################################################
######### Clean up files

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/assembly/megahit/megahit_out

# Cleaning up
rm "$element"_final.contigs_2500_mapped.sam
rm "$element"_final.contigs_2500_mapped.sorted.bam
rm "$element"_final.contigs_2500_mapped99per.sorted.bam
gzip "$element"_final.contigs_2500_mapped.bam

cd ../../../../
done

#SAVE
clean_finalcontigsfiles_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash clean_finalcontigsfiles_loop.sh list.txt

#save, submitted march 13 at 10:30
clean_finalcontigfiles.sh

###########################################################################################
####### Annotations using DRAM

#evidently i deleted the dram original run thinking i didn't need it on each assembly (i do) so rerunning may 50 2023 so we can run prodigal on all the genes.

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs

#add date behind DRAM_1.4.4
DRAM.py annotate -i '*fa' -o  ../DRAM_1.4.4_05052023 --min_contig_size 2500 --threads 20
cd ../../
done

#save
dram_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

bash dram_loop.sh list.txt

#save, submitted may 5 at 12 noon, finished around midnight
dram.sh

#rename dram directories: 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs

mv DRAM_1.2.4_output DRAM_1.4.4_05052023
cd ../../
done

clean_dram_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash clean_dram_loop.sh list.txt

#save, 
clean_dram.sh


### distill the annotations 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs

DRAM.py distill -i DRAM_1.4.4_05052023/annotations.tsv -o DRAM_1.4.4_05052023/distill

cd ../../
done

dram_distill_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

bash dram_distill_loop.sh list.txt

#save, submitted may 15 at 930am -10am
dram_distill.sh


###########################################################################################
####### build a genome database

### dRep

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs

dRep dereplicate dRep_v2.6.2 -p 20 -comp 50 -con 10 -g ./*fa
cd ../../
done

#save
sbatch dRep_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash dRep_loop.sh list_dRep.txt

#save, began march 20 11am
dRep.sh

#rename the dRep_v2.6.2 file.directory to indicate how many MAGs were used.

#CerealRye_T10
mv dRep_v2.6.2 dRep_v2.6.2_2MAGs
#CerealRye_T21
mv dRep_v2.6.2 dRep_v2.6.2_6MAGs
#Control_T7
mv dRep_v2.6.2 dRep_v2.6.2_2MAGs
#Control_T10
mv dRep_v2.6.2 dRep_v2.6.2_2MAGs
#Control_T21
mv dRep_v2.6.2 dRep_v2.6.2_6MAGs
#HairyVetch_T10
mv dRep_v2.6.2 dRep_v2.6.2_4MAGs
#HairyVetch_T21
mv dRep_v2.6.2 dRep_v2.6.2_7MAGs
#Rapeseed_T21
mv dRep_v2.6.2 dRep_v2.6.2_7MAGs
#Sorghum_T10
mv dRep_v2.6.2 dRep_v2.6.2_3MAGs
#Sorghum_T21
mv dRep_v2.6.2 dRep_v2.6.2_7MAGs


# The resulting dereplicated MAG set will be in dRep_v2.6.2_# OF MAGs/dereplicated_genomes. 
# To count the number of genomes in this directory use 
ls -dq *fa | wc -l

#CerealRye_T10
0 MAGs?
#CerealRye_T21
6MAGs
#Control_T7
2MAGs
#Control_T10
0 MAGs ?
#Control_T21
6MAGs
#HairyVetch_T10
4MAGs
#HairyVetch_T21
7MAGs
#Rapeseed_T21
7MAGs
#Sorghum_T10
3 MAGs
#Sorghum_T21
7 MAGs



###########################################################################################
#####Assessing the quality of your MAG database

### Calculate percent reads mapped to MAG database
# This will give us a sense of how much of our data we are using with our MAGs. 
#To calculate the percentage of metagenomic reads that map to the MAG database, users 
# should 1) rename MAG contig headers to match bin names (we need to be able to match the 
#contigs to the bins so we know where the contig came from) 2) concatenate MAG contig 
#fasta files 3) map to concatenated contig fasta file


# Renaming scaffolds (they have to have the same name in DRAM) 
# will rename the fasta headers within the file to the name of the bin file in the same 
# manner as DRAM. This doesn't need to be ran in slurm unless more than 500 genomes. 
# We need to rename the headers because we are compiling MAGs from different assemblies, 
# so there may be multiple scaffolds with the same names--the mappers cant deal with this.

# Whenever you rename files, always put the new renamed files in a NEW directory. 
# You should always check the files are the same except for the renaming (so you can check 
# file size or length)

rename_bins_like_dram.py -i '*fa' -o genome_renamed

# to be renamed:
# not sure why dRep didnt work on control T10 and CR t10, so i just copied the MAGs from the
# MAG directory into the dereplicated genomes directory. and then ran the rename command

#CRt10 and control t10 were manually moved
#CerealRye_T21
#Control_T7
#Control_T21
#HairyVetch_T10
#HairyVetch_T21
#Rapeseed_T21
#Sorghum_T10
#Sorghum_T21


### Concatenate MAG database and map a set of metagenomic reads

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs/dRep_v2.6.2_*MAGs/dereplicated_genomes/genome_renamed
cat *fa > cat_MAGs_dRep.fa
bbmap.sh -Xmx48G threads=20 overwrite=t ref=cat_MAGs_dRep.fa in1=../../../../trimmed_reads/"$element"_R1_trimmed.fastq in2=../../../../trimmed_reads/"$element"_R2_trimmed.fastq outm="$element"_final.contigs_2500_mapped_interleaved.fastq #reads that mapped to the genome database

cd ../../../../../
done

#wouldn't we want to add the MAG files that only had one MAG from that sample? going to stop here with regular assembly and just do coassembly for now....

#save
map_to_mags_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash map_to_mags_loop.sh list_map_to_mags.txt

#SAVE
sbatch map_to_mags.sh

#To count the number of reads that mapped use 
wc -l WCRC_304_final.contigs_2500_mapped_interleaved.fastq

#fastq files have 4 lines per read, so divide the number by 4 that gives the number of 
# reads that mapped. Divide that number by the total reads in your trimmed reads to get the 
# % of reads that mapped to your genome database.

#Divide this number by four to obtain the number of reads that mapped. The percent reads 
# that mapped to your database is then calculated as number of reads that mapped/total 
# trimmed reads that went into the assembly.

#Important: Delete these reads that mapped (WCRC_304_final.contigs_2500_mapped_interleaved.fastq) after recording the percentage information.

#Stop and think:
## What is the percent reads mapped? Is this percentage typical of what we see in the environment you are assembling from?
## What is the percent reads mapped distribution across your samples? Are any particularly lower than the rest?


###########################################################################################
##### Use singleM to profile microbial lineages represented in metagenomes

source /opt/Miniconda2/miniconda2/bin/activate singlem

singlem pipe --sequences <assembly> --otu_table <output OTU table> --output_extra --threads <#>



###########################################################################################
##### Mapping for determining abundance with metagenomic reads

#The first thing we need to do is build a bowtie2 database from our MAG database. Then we 
# will map our metaG reads.

# before you make bowtie database, make sure you count the number of scaffolds in every file to make sure they are all mapping and you don’t miss anything

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../

# make directory for bowtie database and build it
mkdir bowtie_DB
cd bowtie_DB
bowtie2-build ../../MAGs/DRAM_v1.4_111122/scaffolds.fna WCRC_99perMAG_DB --threads 15

# now map reads to the db and output a SAM file
cd ../
bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 15 -x bowtie_DB/WCRC_99perMAG_DB -S WCRC_304_mapped_99perMAGs.sam -1 ../WCRC_304/processed_reads/WCRC_304_R1_trimmed.fastq -2 ../WCRC_304/processed_reads/WCRC_304_R2_trimmed.fastq


sbatch bowtie_map.sh


#Step 2: Convert SAM to BAM, filter, sort
# Next we need to make the SAM into a BAM file. Then we will filter the BAM, and finally sort it.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../
# convert SAM to BAM
samtools view -@ 15 -bS WCRC_304_mapped_99perMAGs.sam > WCRC_304_mapped_99perMAGs.bam

# filter BAM to mappings 97% id, removing multi-mappings and keeping only proper pairs
reformat.sh in=WCRC_304_mapped_99perMAGs.bam out=WCRC_304_97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t

# position sort filtered BAM
samtools sort -@ 15 -o WCRC_304_97peridmapped_99perMAGs_POSSORT.bam WCRC_304_97peridmapped_99perMAGs.bam


sbatch sam2bam_filter_sort.sh

Step 2: Convert SAM to BAM, filter, sort
Next we need to make the SAM into a BAM file. Then we will filter the BAM, and finally sort it.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../
# convert SAM to BAM
samtools view -@ 15 -bS WCRC_304_mapped_99perMAGs.sam > WCRC_304_mapped_99perMAGs.bam

# filter BAM to mappings 97% id, removing multi-mappings and keeping only proper pairs
reformat.sh in=WCRC_304_mapped_99perMAGs.bam out=WCRC_304_97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t

# position sort filtered BAM
samtools sort -@ 15 -o WCRC_304_97peridmapped_99perMAGs_POSSORT.bam WCRC_304_97peridmapped_99perMAGs.bam


sbatch sam2bam_filter_sort.sh

#Stop and Think:
 ## is tossing multimappings best for your system? how might this affect your ultimate results?
 ## how many reads do you have mapped?

# Step 3: Determine MAG coverage.
# We will use a tool called coverM to calculate MAG coverage. We will run coverM three 
# times, outputting different metrics that we will then combine to determine MAG coverage.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../ # use MAG directory if your MAGs are not named fa

# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0 -m reads_per_base --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt

# output MAGs with min-covered_fraction >0.97
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0.75 --output-file coverm_min75.txt &> min75_stats.txt

# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 -m trimmed_mean --output-file coverm_trimmed_mean.txt &> trimmed_mean_stats.txt


sbatch coverm.sh

#Convert the coverm_reads_per_base.txt values to coverage by multiplying the number by read length (likely 151bp). Then, in Excel, combine the 3 coverm_* output files to have the trimmed mean of MAGs (>3X coverage across 75% MAG).
## you have to output the file yourself
## you can report the trimmed mean of MAG with at least 3X coverage across 75% of the MAGs.

# Finally, when you are satisfied with your data, delete all SAM files and redundant BAM 
# files.

# Next step would be mapping metaT reads


###########################################################################################
####### Other

## re-run megahit assembly on Control T10 and HairyVetch T3
# Control T10

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00 #changed this to 28 
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/Control_T10/assembly

mkdir megahit_rerun
cd megahit_rerun

megahit -1 ../../trimmed_reads/Control_T10_R1_trimmed.fastq -2 ../../trimmed_reads/Control_T10_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 -m 0.4  -t 50

#save & run: submitted march 4 930am
megahit_rerun_Control_T10.sh


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00 
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/Control_T10/assembly/megahit_rerun/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o Control_T10_final.contigs_STATS
#save & run: submitted march 4 930am
assembly_stats_C_T10_megahit_rerun.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/Control_T10/assembly/megahit_rerun/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o Control_T10_final.contigs_2500.fa

# map to these scaffolds to get sam file
bbmap.sh -Xmx48G threads=20 overwrite=t ref=Control_T10_final.contigs_2500.fa in1=../../../trimmed_reads/Control_T10_R1_trimmed.fastq in2=../../../trimmed_reads/Control_T10_R2_trimmed.fastq out=Control_T10_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS Control_T10_final.contigs_2500_mapped.sam > Control_T10_final.contigs_2500_mapped.bam
samtools sort -T Control_T10_final.contigs_2500_mapped.sorted -o Control_T10_final.contigs_2500_mapped.sorted.bam Control_T10_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=Control_T10_final.contigs_2500_mapped.sorted.bam out=Control_T10_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

#bin
runMetaBat.sh Control_T10_final.contigs_2500.fa Control_T10_final.contigs_2500_mapped99per.sorted.bam

#save and run, submitted march 9 930am, finished around 5
sbatch map_bin_ControlT10.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/Control_T10/assembly/megahit_rerun/megahit_out/Control_T10_final.contigs_2500.fa.metabat-bins
checkm lineage_wf -t 20 -x fa . checkm

cd checkm
checkm qa -o 2 -f results.txt --tab_table -t 20 lineage.ms .

#save, began match 10 at 9am
checkM_CT10.sh

awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/Control_T10/assembly/megahit_rerun/megahit_out/Control_T10_final.contigs_2500.fa.metabat-bins/checkm/results.txt

bin.1
bin.10

mkdir MAGs
# Copy the HQ/MQ MAGs to your MAGs folder (make the MAG directory in the same sample directory. potential path would be: ../Control_T3/MAGs)
cd MAGs
# Copy MQ + HQ MAGs to this directory (while in MAGs directory):
cp ../assembly/megahit_rerun/megahit_out/Control_T10_final.contigs_2500.fa.metabat-bins/bin.1.fa .
cp ../assembly/megahit_rerun/megahit_out/Control_T10_final.contigs_2500.fa.metabat-bins/bin.10.fa .
# Rename the bins
mv bin.1.fa Control_T10_bin.1.fa 
mv bin.10.fa Control_T10_bin.10.fa 


# Hairy Vetch T3

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00 #changed this to 28 
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T3/assembly

mkdir megahit_rerun
cd megahit_rerun

megahit -1 ../../trimmed_reads/HairyVetch_T3_R1_trimmed.fastq -2 ../../trimmed_reads/HairyVetch_T3_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 -m 0.4  -t 50

#save & run: submitted march 5 9am 
megahit_rerun_HairyVetch_T3.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00 
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T3/assembly/megahit_rerun/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o HairyVetch_T3_final.contigs_STATS
#save & run: submitted march 6 9 am, 
assembly_stats_HV_T3_megahit_rerun.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T3/assembly/megahit_rerun/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o HairyVetch_T3_final.contigs_2500.fa

# map to these scaffolds to get sam file
bbmap.sh -Xmx48G threads=20 overwrite=t ref=HairyVetch_T3_final.contigs_2500.fa in1=../../../trimmed_reads/HairyVetch_T3_R1_trimmed.fastq in2=../../../trimmed_reads/HairyVetch_T3_R2_trimmed.fastq out=HairyVetch_T3_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS HairyVetch_T3_final.contigs_2500_mapped.sam > HairyVetch_T3_final.contigs_2500_mapped.bam
samtools sort -T HairyVetch_T3_final.contigs_2500_mapped.sorted -o HairyVetch_T3_final.contigs_2500_mapped.sorted.bam HairyVetch_T3_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=HairyVetch_T3_final.contigs_2500_mapped.sorted.bam out=HairyVetch_T3_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

#bin
runMetaBat.sh HairyVetch_T3_final.contigs_2500.fa HairyVetch_T3_final.contigs_2500_mapped99per.sorted.bam

#save and run, submitted march 9 940, finished around 5
sbatch map_bin_HairyVetch_T3.sh


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T3/assembly/megahit_rerun/megahit_out/HairyVetch_T3_final.contigs_2500.fa.metabat-bins
checkm lineage_wf -t 20 -x fa . checkm

cd checkm
checkm qa -o 2 -f results.txt --tab_table -t 20 lineage.ms .

#save, began match 10 at 9am
checkM_HVT3.sh

awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T3/assembly/megahit_rerun/megahit_out/HairyVetch_T3_final.contigs_2500.fa.metabat-bins/checkm/results.txt
#no medium or HQ bins even after megahit rerun



#### Other
sickle pe -f CerealRye_T10_R1.fastq  -r CerealRye_T10_R2.fastq   -t sanger -o CerealRye_T10_R1_trimmed.fastq -p CerealRye_T10_R2_trimmed.fastq -s discared_R1R2.fastq

fastqc CerealRye_T10_R1_trimmed.fastq CerealRye_T10_R2_trimmed.fastq

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

###### FIXING MY MISTAKE OF NOT COMBINING ALL SAMPLE MQ AND HQ MAGS INTO ONE DIRECTORY AND RUNNING
# DRAM, DREP, GTDB ON THAT DIRECTORY. 

# REMOVE THE DRAM, DREP, GTDB DIRECTORIES FROM SAMPLE FOLDERS. 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs/

rm -r dRep_v2.6.2_*MAGs
rm -r DRAM_1.4.4_03142023
rm -r gtdb_03102023

cd ../../
done

#save
clean_mags_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash clean_mags_loop.sh list.txt

#SAVE
sbatch clean_mags.sh


## MOVE MQ AND HQ MAGS TO A NEW DIRECTORY

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG
mkdir MedHighQualityMAGs

#!/bin/bash
for element in $(<$1)
do
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/MAGs/*.fa /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityMAGs

cd ../../
done

copy_MAGs_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash copy_MAGs_loop.sh list.txt

sbatch copy_MAGs.sh


##########################################################################################
######### Evaluating MAG taxonomy (on MedHighQualityMAGs)


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityMAGs

#this will create SAM summary files for the taxonomy assigned for the two bins
gtdbtk classify_wf -x fa --genome_dir . --out_dir gtdb_03282023 --cpus 20

#save, submitted in MetaG/MedHighQualityBins/slurm march 
gtdb-tk.sh

###########################################################################################
####### Annotations using DRAM


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityMAGs
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
DRAM.py annotate -i '*fa' -o  ../DRAM_1.4.4_03282023 --min_contig_size 2500 --threads 20

#save, submitted march
dram.sh

#rename dram directories: 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityBins

mv DRAM_1.2.4_output DRAM_1.4.4_DATE
cd ../../
done

clean_dram_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash clean_dram_loop.sh list.txt

#save, 
clean_dram.sh


### distill the annotations 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityBins

DRAM.py distill -i DRAM_1.4.4_03142023/annotations.tsv -o DRAM_1.4.4_03142023/distill

cd ../../
done

dram_distill_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

bash dram_distill_loop.sh list.txt

#save
dram_distill_coA.sh




###########################################################################################
####### build a genome database

### dRep

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityBins

dRep dereplicate dRep_v2.6.2 -p 20 -comp 50 -con 10 -g ./*fa

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityMAGs

dRep dereplicate dRep_v2.6.2 -p 20 -comp 50 -con 10 -g ./*fa

#save, began april 3  
dRep.sh

#rename the dRep_v2.6.2 file.directory to indicate how many MAGs were used.


# The resulting dereplicated MAG set will be in dRep_v2.6.2_# OF MAGs/dereplicated_genomes. 
# To count the number of genomes in this directory use 
ls -dq *fa | wc -l

# 54 input
# 30 after dRep

###########################################################################################
#####Assessing the quality of your MAG database

rename_bins_like_dram.py -i '*fa' -o genome_renamed


### Concatenate MAG database and map a set of metagenomic reads

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MedHighQualityBins/dRep_v2.6.2_*MAGs/dereplicated_genomes/genome_renamed
cat *fa > cat_MAGs_dRep.fa
bbmap.sh -Xmx48G threads=20 overwrite=t ref=cat_MAGs_dRep.fa in1=../../../../trimmed_reads/"$element"_R1_trimmed.fastq in2=../../../../trimmed_reads/"$element"_R2_trimmed.fastq outm="$element"_final.contigs_2500_mapped_interleaved.fastq #reads that mapped to the genome database

cd ../../../../../
done

#wouldn't we want to add the MAG files that only had one MAG from that sample? going to stop here with regular assembly and just do coassembly for now....

#save
map_to_mags_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash map_to_mags_loop.sh list_map_to_mags.txt

#SAVE
sbatch map_to_mags.sh

#To count the number of reads that mapped use 
wc -l WCRC_304_final.contigs_2500_mapped_interleaved.fastq

#Important: Delete these reads that mapped (WCRC_304_final.contigs_2500_mapped_interleaved.fastq) after recording the percentage information.

###########################################################################################
##### Use singleM to profile microbial lineages represented in metagenomes

source /opt/Miniconda2/miniconda2/bin/activate singlem
singlem pipe --sequences <assembly> --otu_table <output OTU table> --output_extra --threads <#>

###########################################################################################
##### Mapping for determining abundance with metagenomic reads

#The first thing we need to do is build a bowtie2 database from our MAG database. Then we 
# will map our metaG reads.

# before you make bowtie database, make sure you count the number of scaffolds in every file to make sure they are all mapping and you don’t miss anything

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../

# make directory for bowtie database and build it
mkdir bowtie_DB
cd bowtie_DB
bowtie2-build ../../MAGs/DRAM_v1.4_111122/scaffolds.fna WCRC_99perMAG_DB --threads 15

# now map reads to the db and output a SAM file
cd ../
bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 15 -x bowtie_DB/WCRC_99perMAG_DB -S WCRC_304_mapped_99perMAGs.sam -1 ../WCRC_304/processed_reads/WCRC_304_R1_trimmed.fastq -2 ../WCRC_304/processed_reads/WCRC_304_R2_trimmed.fastq
sbatch bowtie_map.sh


#Step 2: Convert SAM to BAM, filter, sort
# Next we need to make the SAM into a BAM file. Then we will filter the BAM, and finally sort it.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../
# convert SAM to BAM
samtools view -@ 15 -bS WCRC_304_mapped_99perMAGs.sam > WCRC_304_mapped_99perMAGs.bam

# filter BAM to mappings 97% id, removing multi-mappings and keeping only proper pairs
reformat.sh in=WCRC_304_mapped_99perMAGs.bam out=WCRC_304_97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t

# position sort filtered BAM
samtools sort -@ 15 -o WCRC_304_97peridmapped_99perMAGs_POSSORT.bam WCRC_304_97peridmapped_99perMAGs.bam


sbatch sam2bam_filter_sort.sh

Step 2: Convert SAM to BAM, filter, sort
Next we need to make the SAM into a BAM file. Then we will filter the BAM, and finally sort it.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../
# convert SAM to BAM
samtools view -@ 15 -bS WCRC_304_mapped_99perMAGs.sam > WCRC_304_mapped_99perMAGs.bam

# filter BAM to mappings 97% id, removing multi-mappings and keeping only proper pairs
reformat.sh in=WCRC_304_mapped_99perMAGs.bam out=WCRC_304_97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t

# position sort filtered BAM
samtools sort -@ 15 -o WCRC_304_97peridmapped_99perMAGs_POSSORT.bam WCRC_304_97peridmapped_99perMAGs.bam

sbatch sam2bam_filter_sort.sh

# Step 3: Determine MAG coverage.
# We will use a tool called coverM to calculate MAG coverage. We will run coverM three 
# times, outputting different metrics that we will then combine to determine MAG coverage.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../ # use MAG directory if your MAGs are not named fa

# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0 -m reads_per_base --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt

# output MAGs with min-covered_fraction >0.97
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0.75 --output-file coverm_min75.txt &> min75_stats.txt

# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 -m trimmed_mean --output-file coverm_trimmed_mean.txt &> trimmed_mean_stats.txt


sbatch coverm.sh

# Next step would be mapping metaT reads




###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

###### MOVING ALL BINS TO NEW DIRECTORY TO RUN DREP/GTDB AND COMPARE RESULTS TO COASSEMBLY ALL BINS

# manually renamed all bins to have sample name, now move to all_bins

#!/bin/bash
for element in $(<$1)
do
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/assembly/megahit/megahit_out/"$element"_final.contigs_2500.fa.metabat-bins/*.fa /home/projects/Agribiome/RootExudate_Prenni/MetaG/all_bins

cd ../../../../../
done

copy_mags_to_allbins_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash copy_mags_to_allbins_loop.sh list.txt

sbatch copy_mags_to_allbins.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/all_bins

dRep dereplicate dRep_v2.6.2 -p 20 -comp 10 -con 10 -g ./*fa

#submitted MetaG/all_bins/slurm april 3 10:30am - 11:10am
sbatch dRep.sh


# The resulting dereplicated MAG set will be in dRep_v2.6.2_109 OF MAGs/dereplicated_genomes. 
# To count the number of genomes in this directory use 
ls -dq *fa | wc -l

#input MAGS was 213, output was 109

#rename the dRep_v2.6.2 to have the number of MAGs used = dRep_v2.6.2_109


##########################################################################################
######### Evaluating MAG taxonomy (on all_bins)


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/all_bins

#this will create SAM summary files for the taxonomy assigned for the two bins
gtdbtk classify_wf -x fa --genome_dir . --out_dir gtdb_04032023 --cpus 20

#save, submitted in MetaG/all_bins/slurm april 3  
gtdb.sh

sbatch --dependency=afterok:58534 gtdb.sh
















