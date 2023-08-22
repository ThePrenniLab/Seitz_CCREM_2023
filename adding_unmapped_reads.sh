### Concatenate MAG database and map a set of metagenomic reads, output the unmapped reads

## Idea here is that we will get a set of the unmapped reads, then we can assemble those
## and see if we can bin anymore scaffolds. we will then check and map those new bins to the 
## genome database and see if we get much of an increase in % reads mapped. 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed

bbmap.sh -Xmx48G threads=20 overwrite=t ref=MAGs_concatenated_dRep.fa in1=../../../trimmed_reads/concat_R1_coassembly_trimmed.fastq in2=../../../trimmed_reads/concat_R2_coassembly_trimmed.fastq outu1=R1_final.contigs_2500_unmapped.fastq outu2=R2_final.contigs_2500_unmapped.fastq #reads that were not mapped to the genome database

#SAVE, submitted april 5 10 am - april 9 5pm
sbatch extract_unmapped_reads.sh
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed
wc -l R1_final.contigs_2500_unmapped.fastq
# how many reads were unmapped? 2771689220/4 = 692,922,305 


###########################################################################################
###### Megahit Assembly on unmapped reads

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/unmapped_reads_assembly

mkdir megahit
cd megahit

megahit -1 ../../coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed/R1_final.contigs_2500_unmapped.fastq -2 ../../coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed/R2_final.contigs_2500_unmapped.fastq --k-min 31 --k-max 121 --k-step 10 -m 0.4  -t 50

#save in unmapped_reads_assembly/slurm began april 10 at 10am
megahit.sh


###########################################################################################
###### Megahit Assembly Stats

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/unmapped_reads_assembly/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o unmapped_final.contigs_STATS


#save, submitted 11am-1pm may 12
assembly_stats.sh

# add assembly stats to the file: assembly_stats.xlsx file 

###########################################################################################
######## BINNING

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/unmapped_reads_assembly/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o final.contigs_2500.fa

# map to these scaffolds to get sam file

bbmap.sh -Xmx48G threads=20 overwrite=t ref=final.contigs_2500.fa in1=../../coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed/R1_final.contigs_2500_unmapped.fastq in2=../../coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed/R2_final.contigs_2500_unmapped.fastq out=final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS final.contigs_2500_mapped.sam > final.contigs_2500_mapped.bam
samtools sort -T final.contigs_2500_mapped.sorted -o final.contigs_2500_mapped.sorted.bam final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=final.contigs_2500_mapped.sorted.bam out=final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

#bin
runMetaBat.sh final.contigs_2500.fa final.contigs_2500_mapped99per.sorted.bam

#submitted may 12 12:00, done may 19
sbatch map_bin.sh 

#roughly 250 bins

###########################################################################################
######## QC on bins 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/unmapped_reads_assembly/megahit_out/final.contigs_2500.fa.metabat-bins
checkm lineage_wf -t 20 -x fa . checkm

cd checkm
checkm qa -o 2 -f results.txt --tab_table -t 20 lineage.ms .

# save, submitted may 19 930am-1pm
checkm.sh




###########################################################################################
####### Evaluate binning results (manual step): 

# use this to go directly to extracting the M and HQ bins, anything with greater than 49% completion and less than 11 contamination
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/unmapped_reads_assembly/megahit_out/final.contigs_2500.fa.metabat-bins/checkm/results.txt

# Record good bins in new excel file: "MAG processing manual code"

#EXAMPLE (SEE MAG processing manual code FOR ALL MANUAL CODE)
# make MAG directories in each sample folder
mkdir MAGs

# Copy the HQ/MQ MAGs to your MAGs folder (make the MAG directory in the same sample directory. potential path would be: ../Control_T3/MAGs)
cd MAGs

# Copy MQ + HQ MAGs to this directory (while in MAGs directory):
# example:
cp ../megahit_out/final.contigs_2500.fa.metabat-bins/bin.126.fa .

# Rename the bins: this is an iterative coassembly, so names bins as coasssembly2
mv bin.108.fa coassembly2_bin.108.fa

# move these M/HQ mags to the full database, then drep that database again.

cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/unmapped_reads_assembly/MAGs/*fa /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022

###########################################################################################
######### Clean up files

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/unmapped_reads_assembly/megahit_out

# Cleaning up
rm final.contigs_2500_mapped.sam
rm final.contigs_2500_mapped.sorted.bam
rm final.contigs_2500_mapped99per.sorted.bam
gzip final.contigs_2500_mapped.bam


#save
clean.sh


# then rerun map_to_mags to see if the % of reads mapping increased. 



####### dereplicate full genome database

### dRep
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022

dRep dereplicate dRep_v2.6.2 -p 20 -comp 49 -con 10 -g ./*fa

#saved in MedHighQualityBins/slurm, began march 27 9pm - midnight
sbatch dRep.sh

#rename the dRep_v2.6.2 to have the number of MAGs used = dRep_v2.6.2_43MAGs

ls -dq *fa | wc -l
# there were 342 M and HQ bins before dRep
# there are 326 bins after dRep

mv dRep_v2.6.2 dRep_v2.6.2_326


use ls *fa >> list.txt to list all the fa files needed

wc -l use to check all files are there


###########################################################################################
####### DRAM on individual bins


#!/bin/bash
 
#$1 list of names of the paired reads that you want to trim.
 
for element in $(<$1)
do

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/
 
DRAM.py annotate -i "$element" -o  "$element"_DRAM_1.4.4_052023 --min_contig_size 2500 --threads 20
DRAM.py distill -i "$element"_DRAM_1.4.4_052023/annotations.tsv -o "$element"_DRAM_1.4.4_052023/distill

done

dram_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=128gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash dram_loop.sh list.txt



###########################################################################################
####### Annotations using DRAM + custom cover crop root exudate module on MedHighQualityBins

# 	WAITING TO DO THIS STEP UNTIL WE HAVE METAB DATA


#!/bin/bash
#$1 list of names of the paired reads that you want to trim.
 
for element in $(<$1)
do

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/
 
DRAM.py distill -i "$element"_DRAM_1.4.4_052023/annotations.tsv -o "$element"_DRAM_1.4.4_052023/distill_wCustomDistillate --custom_distillate custom_distillate_coverCrop_exudates.txt

done

#save, submitted april 3 4pm
dram_custom_distillate_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash dram_custom_distillate_loop.sh list.txt

###########################################################################################
####### merged dram output files


#!/bin/bash
#$1 list of names of the paired reads that you want to trim.
 
for element in $(<$1)
do
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes
DRAM.py merge_annotations -i "$element"_DRAM_1.4.4_052023/annotations.tsv -o merged_annotation


done

#save
dram_merge_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash dram_merge_loop.sh list.txt

#save, submitted in microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/slurm on may 22 245pm
dram_merge.sh


###########################################################################################
####### Annotations using DRAM on full database (no need to merge when done this way)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes

DRAM.py annotate -i '*fa' -o  DRAM_1.4.4_05222023 --min_contig_size 2500 --threads 20

# distill the annotations 
DRAM.py distill -i DRAM_1.4.4_05222023/annotations.tsv -o DRAM_1.4.4_05222023/distill

#saved in MedHighQualityBins/slurm, started may 22 545pm
dram_no_loop.sh

##########################################################################################
######### Rename Bins


### Renames bins to match DRAM, cd into MedHighQualityBins/dRep_v2.6.2_*MAGs/dereplicated_genomes/ first. 
### The scaffold names (the headers in the files) will be changed from  (e.g.): k121_2164525 to HairyVetch_bin.11_k121_2164525
rename_bins_like_dram.py -i '*fa' -o genomes_renamed


##########################################################################################
######### Evaluating MAG taxonomy

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/

#this will create SAM summary files for the taxonomy assigned for the two bins
# old gtdb code: gtdbtk classify_wf -x fa --genome_dir . --out_dir gtdb_03282023 --cpus 20

gtdbtk classify_wf --extension fa --genome_dir . --out_dir ./gtdb_classify_v2.3.0_release214 --skip_ani_screen

#save, submitted in microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/slurm
gtdb.sh




run dram on each bin from the dReped genome database from the added 250 and the old database
- dram in a loop over all files 
- run gtdb normally, there is a new verion
- then dram merge anotations
- bowtie? 
-  rerun map_to_mags to see if the % of reads mapping increased. 
	- 
MetaG: Bowtie needs to be run on EACH sample so that you can map MAG abundances to each timepoint/sample. Bowtie is mapping MAGs back to the original samples to tell you abundance and then running coverM after that gives you coverage. 

MetaT: use scaffold.fna from the merged annotations from dram to map metaT to that file. you will get more mapping this way. 
- we will use feature counts which takes genes.ff files which has genes and where the scaffold it came from and calculates where they map on the scaffold file and counts. 

running bowtie on each sample reads to dRep bins. 




Priority

1. 
--
checkM - done
rename and move MQ/HQ bins to the 2019_2022 database - done
dRep all those MAGs one more time. - done 
run DRAM on the dereplicated MAGs, run on each bin so each one will have its own annotation, done
run GTDB on the dReped bins, run on the full database, not individual bins. done 

#work on making the metadata file for metagenomes/bins
--

2. Mapping for metageome abundance (bowtie)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low


cd /home/projects/Agribiome/RootExudate_Prenni/MetaG

# make directory for bowtie database and build it using the merged dram file from the full database
mkdir bowtie_DB
cd bowtie_DB

#build a database of scaffolds from the dram scaffolds file, bowtie2-build outputs a set of 6 files 
#with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2. In the case of a 
#large index these suffixes will have a bt2l termination. These files together constitute the index: 
#they are all that is needed to align reads to that reference.

#the 99perMAGs means that we are using a 99% dereplicated MAG database
bowtie2-build ../microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/scaffolds.fna 99perMAG_DB --threads 15

#submitted metaG/slurm at 11am-11:15 june 1
build_bowtie.sh


# map metagenome reads to bowtie database

#!/bin/bash
for element in $(<$1)
do

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG

# now map individual metagenome reads to the db and output a SAM file
# map to trimmed reads for each sample to determine abundance of each MAG in each sample.
bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 15 -x bowtie_DB/99perMAG_DB -S bowtie_DB/"$element"_mapped_99perMAGs.sam -1 "$element"/trimmed_reads/"$element"_R1_trimmed.fastq -2 "$element"/trimmed_reads/"$element"_R2_trimmed.fastq

done

#submitted in MetaG/slurm
bowtie_map_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

#submit in MetaG/slurm
bash bowtie_map_loop.sh list.txt

#submit 1115 on june 1
sbatch --dependency=afterok:59723 bowtie_map.sh

#will have a sam file for each metagenome. so check you have 26 at the end



# convert SAM to BAM


#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB

# convert SAM to BAM
samtools view -@ 15 -bS "$element"_mapped_99perMAGs.sam > "$element"_mapped_99perMAGs.bam

# filter BAM to mappings 97% id, removing multi-mappings (multiple possible locations a read can map, so mapper choses BOTH, we want to remove this) 
# and keeping only proper pairs
# we are filtering out bad mappings by specifying that we want to map at 97% identity to our 99% deprelicated database
reformat.sh in="$element"_mapped_99perMAGs.bam out="$element"_97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t

# position sort filtered BAM, position sorting the BAM re-orders the lines of the BAM file according to the reference positions (ie. the contig positions).
samtools sort -@ 15 -o "$element"_97peridmapped_99perMAGs_POSSORT.bam "$element"_97peridmapped_99perMAGs.bam

done

#submit in metaG/slurm
sam2bam_filter_sort_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

bash sam2bam_filter_sort_loop.sh list.txt

#submit
sam2bam_filter_sort.sh

sbatch --dependency=afterok:59725 sam2bam_filter_sort.sh

# CoverM - genome abundance table

#then use coverM to get the coverage of each bin from each sample...abundnace data like OTU table

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB

# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0 -m reads_per_base coverm_reads_per_base.txt &> reads_per_base_stats.txt
# output MAGs with min-covered_fraction >0.75
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0.75 --output-file coverm_min75.txt &> min75_stats.txt
# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 -m trimmed_mean --output-file coverm_trimmed_mean.txt &> trimmed_mean_stats.txt




# REDO WITH 95% ID to see if we get more bins
# USE 95% ID, MORE BINS

# convert SAM to BAM

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB

# filter BAM to mappings 95% id, removing multi-mappings (multiple possible locations a read can map, so mapper choses BOTH, we want to remove this) 
# and keeping only proper pairs
# we are filtering out bad mappings by specifying that we want to map at 95% identity to our 99% deprelicated database
reformat.sh in="$element"_mapped_99perMAGs.bam out="$element"_95peridmapped_99perMAGs.bam minidfilter=0.95 primaryonly=t pairedonly=t

# position sort filtered BAM, position sorting the BAM re-orders the lines of the BAM file according to the reference positions (ie. the contig positions).
samtools sort -@ 15 -o "$element"_95peridmapped_99perMAGs_POSSORT.bam "$element"_95peridmapped_99perMAGs.bam

done

#submit in metaG/slurm
sam2bam_filter_sort_loop_95%ID.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

bash sam2bam_filter_sort_loop_95%ID.sh list.txt

#submited june 30 10am
sam2bam_filter_sort_95%ID.sh

# then move all those bam files into bowtie_95%ID dir and run coverM from there

# CoverM - genome abundance table

# then use coverM to get the coverage of each bin from each sample...abundnace data like OTU table
# to troubleshoot coverM errors (since they don't print in the slurm output, instead head or less into the min_75_stats or trimmed_mean_stats.txt files to see what went wrong.)
# this time coverM wasn't working because the concatenated MAG file was messing it up since it was counting the same bin twice (once ini genomes_renamed and once in the MAG_concatenated file.)
# moved the concatenated file into "not_fa_files" for now. 


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=400gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/bowtie_DB

# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAG_database/dRep_v2.6.2_326/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 -m reads_per_base --output-file coverm_reads_per_base_95%ID.txt &> reads_per_base_stats_95%ID.txt
# output MAGs with min-covered_fraction >0.75
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAG_database/dRep_v2.6.2_326/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0.75 --output-file coverm_min75_95%ID.txt &> min75_stats_95%ID.txt
# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAG_database/dRep_v2.6.2_326/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.95 -m trimmed_mean --output-file coverm_trimmed_mean_95%ID.txt &> trimmed_mean_stats_95%ID.txt

#save and run, took about 10 mins
coverM_95%ID.sh
















