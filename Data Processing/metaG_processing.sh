# metaG processing for 26 metagenomes

## Note: "$element" is a list of sample names.


###########################################################################################
######## FastQC on raw reads

fastqc *.fastq #perform QC 


###########################################################################################
######## Trim raw reads using Sickle

#!/bin/bash
for element in $(<$1)
do
sickle pe -f "$element"_R1.fastq  -r "$element"_R2.fastq   -t sanger -o "$element"_R1_trimmed.fastq -p "$element"_R2_trimmed.fastq -s discared_R1R2.fastq
done

###########################################################################################
######## FastQC on trimmed reads

#!/bin/bash
for element in $(<$1)
do
fastqc "$element"_R1_trimmed.fastq "$element"_R2_trimmed.fastq
done

###########################################################################################
###### Megahit (v1.2.9) Assembly, done for each metaG, coassembly (combined treatments), and iterative coassembly for unmapped reads

#!/bin/bash
for element in $(<$1)
do
megahit -1 ../../trimmed_reads/"$element"_R1_trimmed.fastq -2 ../../trimmed_reads/"$element"_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 -m 0.4  -t 50
done

###########################################################################################
###### Megahit Assembly Stats

#!/bin/bash
for element in $(<$1)
do
/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o "$element"_final.contigs_STATS
done


###########################################################################################
######## BINNING

#!/bin/bash
for element in $(<$1)
do
pullseq.py -i final.contigs.fa -m 2500 -o "$element"_final.contigs_2500.fa

# map to these scaffolds to get sam file
bbmap.sh -Xmx48G threads=20 overwrite=t ref="$element"_final.contigs_2500.fa in1=../../../trimmed_reads/"$element"_R1_trimmed.fastq in2=../../../trimmed_reads/"$element"_R2_trimmed.fastq out="$element"_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS "$element"_final.contigs_2500_mapped.sam > "$element"_final.contigs_2500_mapped.bam
samtools sort -T "$element"_final.contigs_2500_mapped.sorted -o "$element"_final.contigs_2500_mapped.sorted.bam "$element"_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in="$element"_final.contigs_2500_mapped.sorted.bam out="$element"_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

# bin
runMetaBat.sh "$element"_final.contigs_2500.fa "$element"_final.contigs_2500_mapped99per.sorted.bam
done

###########################################################################################
######## QC on bins 

#!/bin/bash
for element in $(<$1)
do
checkm lineage_wf -t 20 -x fa . checkm
cd checkm
checkm qa -o 2 -f results.txt --tab_table -t 20 lineage.ms .
done

# combine all M/HQ MAGs from each assembly method into one directory as your MAG database.

###########################################################################################
######## dereplicate full genome database

#!/bin/bash
dRep dereplicate dRep_v2.6.2 -p 20 -comp 49 -con 10 -g ./*fa

###########################################################################################
####### DRAM on individual bins


#!/bin/bash
for element in $(<$1)
do
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
DRAM.py annotate -i "$element" -o  "$element"_DRAM_1.4.4_052023 --min_contig_size 2500 --threads 20
DRAM.py distill -i "$element"_DRAM_1.4.4_052023/annotations.tsv -o "$element"_DRAM_1.4.4_052023/distill
done

##########################################################################################
######### Rename Bins
rename_bins_like_dram.py -i '*fa' -o genomes_renamed


##########################################################################################
######### Evaluating MAG taxonomy

#!/bin/bash
gtdbtk classify_wf --extension fa --genome_dir . --out_dir ./gtdb_classify_v2.3.0_release214 --skip_ani_screen

##########################################################################################
########## Concatenate MAG database and map a set of metagenomic reads

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/genomes_renamed
cat *fa > cat_MAGs_dRep.fa
bbmap.sh -Xmx48G threads=20 overwrite=t ref=cat_MAGs_dRep.fa in1=../../../../individual_assembly/"$element"/trimmed_reads/"$element"_R1_trimmed.fastq in2=../../../../individual_assembly/"$element"/trimmed_reads/"$element"_R2_trimmed.fastq outm="$element"_final.contigs_2500_mapped_interleaved.fastq #reads that mapped to the genome database = interleaved file
done


###########################################################################################
##### Mapping for determining abundance with metagenomic reads

#!/bin/bash
mkdir bowtie_DB
cd bowtie_DB

#build a database of scaffolds from the dram scaffolds file
bowtie2-build ../microcosm_MAG_database_2019_2022/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/scaffolds.fna 99perMAG_DB --threads 15

#submitted in /home/projects/Agribiome/RootExudate_Prenni/metaG/slurm at 11am-11:15 june 1
build_bowtie.sh

#!/bin/bash
for element in $(<$1)
do
# now map individual metagenome reads to the db and output a SAM file
# map to trimmed reads for each sample
bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 15 -x bowtie_DB/99perMAG_DB -S bowtie_DB/"$element"_mapped_99perMAGs.sam -1 "$element"/trimmed_reads/"$element"_R1_trimmed.fastq -2 "$element"/trimmed_reads/"$element"_R2_trimmed.fastq
done

# convert SAM to BAM

#!/bin/bash
for element in $(<$1)
do
# convert SAM to BAM
samtools view -@ 15 -bS "$element"_mapped_99perMAGs.sam > "$element"_mapped_99perMAGs.bam
# filter BAM to mappings 97% id, removing multi-mappings and keeping only proper pairs
reformat.sh in="$element"_mapped_99perMAGs.bam out="$element"_97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t
# position sort filtered BAM
samtools sort -@ 15 -o "$element"_97peridmapped_99perMAGs_POSSORT.bam "$element"_97peridmapped_99perMAGs.bam
done

###########################################################################################
##### coverM

#!/bin/bash
# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0 -m reads_per_base --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt
# output MAGs with min-covered_fraction >0.97
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0.75 --output-file coverm_min75.txt &> min75_stats.txt
# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ../MAGs/dRep_v2.6.2_2MAGs/dereplicated_genomes/genomes_renamed --bam-files *POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 -m trimmed_mean --output-file coverm_trimmed_mean.txt &> trimmed_mean_stats.txt
