# CCREM MetaG Sample Processing
# Valerie Seitz
# Jan 21 2023

# Working Directory
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly


###########################################################################################
####### Create Directory Infrastructure: 

# Create a coassembly directory
cd coassembly
# Create directory for each treatment
mkdir Control
# etc...
mkdir slurm
cd slurm
nano

#create text file for looping
Control
Rapeseed
CerealRye
Sorghum
HairyVetch
#save as list.txt

# Make subdirectories

#!/bin/bash
for element in $(<$1)
do
cd "$element"
mkdir assembly
mkdir trimmed_reads
cd ../
done

#save
mkdir_loop.sh

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

# Transfer the trimmed reads (put into each treatment folder)
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T3/trimmed_reads/HairyVetch_T3_R1_trimmed.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/trimmed_reads
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T3/trimmed_reads/HairyVetch_T3_R2_trimmed.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/trimmed_reads

cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T5/trimmed_reads/HairyVetch_T5_R1_trimmed.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/trimmed_reads
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/HairyVetch_T5/trimmed_reads/HairyVetch_T5_R2_trimmed.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/trimmed_reads
#...etc

###########################################################################################
###### Megahit Co-Assembly & Assembly Stats

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Control/assembly
#concatenate forward and reverse trimmed reads files 
cat ../trimmed_reads/Control_T0_R1_trimmed.fastq ../trimmed_reads/Control_T3_R1_trimmed.fastq ../trimmed_reads/Control_T5_R1_trimmed.fastq ../trimmed_reads/Control_T7_R1_trimmed.fastq ../trimmed_reads/Control_T10_R1_trimmed.fastq ../trimmed_reads/Control_T21_R1_trimmed.fastq > concat_R1_trimmed.fastq
cat ../trimmed_reads/Control_T0_R2_trimmed.fastq ../trimmed_reads/Control_T3_R2_trimmed.fastq ../trimmed_reads/Control_T5_R2_trimmed.fastq ../trimmed_reads/Control_T7_R2_trimmed.fastq ../trimmed_reads/Control_T10_R2_trimmed.fastq ../trimmed_reads/Control_T21_R2_trimmed.fastq > concat_R2_trimmed.fastq
#coassemble
megahit -1 concat_R1_trimmed.fastq -2 concat_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 --mem-flag 1 -m 429496729600 -t 20 -o megahit_out

#save and run
sbatch megahit_coA_control.sh

### Get assembly stats
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Control/assembly/megahit_out
/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o Control_final.contigs_STATS

#submit job
sbatch assembly_stats_control.sh


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --dependency:afterok=<jobid>

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Rapeseed/assembly
#concatenate forward and reverse trimmed reads files 
cat ../trimmed_reads/Rapeseed_T3_R1_trimmed.fastq ../trimmed_reads/Rapeseed_T5_R1_trimmed.fastq ../trimmed_reads/Rapeseed_T7_R1_trimmed.fastq ../trimmed_reads/Rapeseed_T10_R1_trimmed.fastq ../trimmed_reads/Rapeseed_T21_R1_trimmed.fastq > concat_R1_trimmed.fastq
cat ../trimmed_reads/Rapeseed_T3_R2_trimmed.fastq ../trimmed_reads/Rapeseed_T5_R2_trimmed.fastq ../trimmed_reads/Rapeseed_T7_R2_trimmed.fastq ../trimmed_reads/Rapeseed_T10_R2_trimmed.fastq ../trimmed_reads/Rapeseed_T21_R2_trimmed.fastq > concat_R2_trimmed.fastq
#coassemble
megahit -1 concat_R1_trimmed.fastq -2 concat_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 --mem-flag 1 -m 429496729600 -t 20 -o megahit_out

sbatch megahit_coA_rapeseed.sh

### Get assembly stats

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --dependency:afterok=<jobid>

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Rapeseed/assembly/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o rapeseed_final.contigs_STATS

#submit job
sbatch assembly_stats_rapeseed.sh


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --dependency:afterok=<jobid>

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/CerealRye/assembly
#concatenate forward and reverse trimmed reads files 
cat ../trimmed_reads/CerealRye_T3_R1_trimmed.fastq ../trimmed_reads/CerealRye_T5_R1_trimmed.fastq ../trimmed_reads/CerealRye_T7_R1_trimmed.fastq ../trimmed_reads/CerealRye_T10_R1_trimmed.fastq ../trimmed_reads/CerealRye_T21_R1_trimmed.fastq > concat_R1_trimmed.fastq
cat ../trimmed_reads/CerealRye_T3_R2_trimmed.fastq ../trimmed_reads/CerealRye_T5_R2_trimmed.fastq ../trimmed_reads/CerealRye_T7_R2_trimmed.fastq ../trimmed_reads/CerealRye_T10_R2_trimmed.fastq ../trimmed_reads/CerealRye_T21_R2_trimmed.fastq > concat_R2_trimmed.fastq
#coassemble
megahit -1 concat_R1_trimmed.fastq -2 concat_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 --mem-flag 1 -m 429496729600 -t 20 -o megahit_out

sbatch megahit_coA_CerealRye.sh


### Get assembly stats

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low


cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/CerealRye/assembly/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o cerealRye_final.contigs_STATS

#submit job
sbatch assembly_stats_cerealRye.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --dependency:afterok=<jobid>

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Sorghum/assembly
#concatenate forward and reverse trimmed reads files 
cat ../trimmed_reads/Sorghum_T3_R1_trimmed.fastq ../trimmed_reads/Sorghum_T5_R1_trimmed.fastq ../trimmed_reads/Sorghum_T7_R1_trimmed.fastq ../trimmed_reads/Sorghum_T10_R1_trimmed.fastq ../trimmed_reads/Sorghum_T21_R1_trimmed.fastq > concat_R1_trimmed.fastq
cat ../trimmed_reads/Sorghum_T3_R2_trimmed.fastq ../trimmed_reads/Sorghum_T5_R2_trimmed.fastq ../trimmed_reads/Sorghum_T7_R2_trimmed.fastq ../trimmed_reads/Sorghum_T10_R2_trimmed.fastq ../trimmed_reads/Sorghum_T21_R2_trimmed.fastq > concat_R2_trimmed.fastq
#coassemble
megahit -1 concat_R1_trimmed.fastq -2 concat_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 --mem-flag 1 -m 429496729600 -t 20 -o megahit_out

sbatch megahit_coA_Sorghum.sh
#submitted jan 24 2pm 

### Get assembly stats

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Sorghum/assembly/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o sorghum_final.contigs_STATS

#submit job
sbatch assembly_stats_sorghum.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --dependency:afterok=<jobid>

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/assembly
#concatenate forward and reverse trimmed reads files 
cat ../trimmed_reads/HairyVetch_T3_R1_trimmed.fastq ../trimmed_reads/HairyVetch_T5_R1_trimmed.fastq ../trimmed_reads/HairyVetch_T7_R1_trimmed.fastq ../trimmed_reads/HairyVetch_T10_R1_trimmed.fastq ../trimmed_reads/HairyVetch_T21_R1_trimmed.fastq > concat_R1_trimmed.fastq
cat ../trimmed_reads/HairyVetch_T3_R2_trimmed.fastq ../trimmed_reads/HairyVetch_T5_R2_trimmed.fastq ../trimmed_reads/HairyVetch_T7_R2_trimmed.fastq ../trimmed_reads/HairyVetch_T10_R2_trimmed.fastq ../trimmed_reads/HairyVetch_T21_R2_trimmed.fastq > concat_R2_trimmed.fastq
#coassemble
megahit -1 concat_R1_trimmed.fastq -2 concat_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 --mem-flag 1 -m 429496729600 -t 20 -o megahit_out

sbatch megahit_coA_HairyVetch.sh

### Get assembly stats

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/assembly/megahit_out

/ORG-Data/scripts/quicklooks/contig_stats.pl -i final.contigs.fa -o hairyVetch_final.contigs_STATS

#submit job
sbatch assembly_stats_hairyVetch.sh


###########################################################################################
######## BINNING

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Control/assembly/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o controlcoA_final.contigs_2500.fa

# map to your scaffolds over 2500bp to get sam file with all of them in one file
bbmap.sh -Xmx48G threads=20 overwrite=t ref=controlcoA_final.contigs_2500.fa in1=../concat_R1_trimmed.fastq in2=../concat_R2_trimmed.fastq out=controlcoA_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS controlcoA_final.contigs_2500_mapped.sam > controlcoA_final.contigs_2500_mapped.bam
samtools sort -T controlcoA_final.contigs_2500_mapped.sorted -o controlcoA_final.contigs_2500_mapped.sorted.bam controlcoA_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=controlcoA_final.contigs_2500_mapped.sorted.bam out=controlcoA_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

# bin our newly-assembled contigs/scaffolds into metagenome-assembled genomes (MAGs)
runMetaBat.sh controlcoA_final.contigs_2500.fa controlcoA_final.contigs_2500_mapped99per.sorted.bam

#Exit, name the script, and run: began feb 27 10 am, finished march 3

sbatch map_bin_control.sh 


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Rapeseed/assembly/megahit_out

# pull scaffolds >= 2500 bp
pullseq.py -i final.contigs.fa -m 2500 -o rapeseed_final.contigs_2500.fa

# map to these scaffolds to get sam file
bbmap.sh -Xmx48G threads=20 overwrite=t ref=rapeseed_final.contigs_2500.fa in1=../concat_R1_trimmed.fastq in2=../concat_R2_trimmed.fastq out=rapeseed_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS rapeseed_final.contigs_2500_mapped.sam > rapeseed_final.contigs_2500_mapped.bam
samtools sort -T rapeseed_final.contigs_2500_mapped.sorted -o rapeseed_final.contigs_2500_mapped.sorted.bam rapeseed_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=rapeseed_final.contigs_2500_mapped.sorted.bam out=rapeseed_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t overwrite=true

#bin
runMetaBat.sh rapeseed_final.contigs_2500.fa rapeseed_final.contigs_2500_mapped99per.sorted.bam

#Exit, name the script, and run: began march 4 9am - march 5 noon
sbatch map_bin_rapeseed.sh 


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/CerealRye/assembly/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o cerealrye_final.contigs_2500.fa

# map to these scaffolds to get sam file

bbmap.sh -Xmx48G threads=20 overwrite=t ref=cerealrye_final.contigs_2500.fa in1=../concat_R1_trimmed.fastq in2=../concat_R2_trimmed.fastq out=cerealrye_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS cerealrye_final.contigs_2500_mapped.sam > cerealrye_final.contigs_2500_mapped.bam
samtools sort -T cerealrye_final.contigs_2500_mapped.sorted -o cerealrye_final.contigs_2500_mapped.sorted.bam cerealrye_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=cerealrye_final.contigs_2500_mapped.sorted.bam out=cerealrye_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

#bin
runMetaBat.sh cerealrye_final.contigs_2500.fa cerealrye_final.contigs_2500_mapped99per.sorted.bam

#Exit, name the script, and run: 
sbatch map_bin_cerealrye.sh 


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Sorghum/assembly/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o sorghum_final.contigs_2500.fa

# map to these scaffolds to get sam file

bbmap.sh -Xmx48G threads=20 overwrite=t ref=sorghum_final.contigs_2500.fa in1=../concat_R1_trimmed.fastq in2=../concat_R2_trimmed.fastq out=sorghum_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS sorghum_final.contigs_2500_mapped.sam > sorghum_final.contigs_2500_mapped.bam
samtools sort -T sorghum_final.contigs_2500_mapped.sorted -o sorghum_final.contigs_2500_mapped.sorted.bam sorghum_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=sorghum_final.contigs_2500_mapped.sorted.bam out=sorghum_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

#bin
runMetaBat.sh sorghum_final.contigs_2500.fa sorghum_final.contigs_2500_mapped99per.sorted.bam

#Exit, name the script, and run: submitted march 7 at 2pm
sbatch map_bin_sorghum.sh 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/assembly/megahit_out

pullseq.py -i final.contigs.fa -m 2500 -o hairyVetch_final.contigs_2500.fa

# map to these scaffolds to get sam file

bbmap.sh -Xmx48G threads=20 overwrite=t ref=hairyVetch_final.contigs_2500.fa in1=../concat_R1_trimmed.fastq in2=../concat_R2_trimmed.fastq out=hairyVetch_final.contigs_2500_mapped.sam

# convert sam to bam, and sort
samtools view -@ 20 -bS hairyVetch_final.contigs_2500_mapped.sam > hairyVetch_final.contigs_2500_mapped.bam
samtools sort -T hairyVetch_final.contigs_2500_mapped.sorted -o hairyVetch_final.contigs_2500_mapped.sorted.bam hairyVetch_final.contigs_2500_mapped.bam -@ 20

# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g idfilter=0.99 in=hairyVetch_final.contigs_2500_mapped.sorted.bam out=hairyVetch_final.contigs_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t

#bin
runMetaBat.sh hairyVetch_final.contigs_2500.fa hairyVetch_final.contigs_2500_mapped99per.sorted.bam

#Exit, name the script, and run: submitted march 7 at 9pm
sbatch map_bin_hairyVetch.sh

###########################################################################################
######## QC on bins 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/"$element"/assembly/megahit_out/"$element"_final.contigs_2500.fa.metabat-bins
checkm lineage_wf -t 20 -x fa . checkm

cd checkm
checkm qa -o 2 -f results.txt --tab_table -t 20 lineage.ms .

cd ../../../../
done

checkM_loop.sh 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu 
#SBATCH --partition=wrighton-hi,wrighton-low

bash checkm_loop.sh list.txt

# save, finished march 10 11am
checkM.sh
sbatch checkM.sh

###########################################################################################
####### Evaluate binning results (manual step): 

# use this to go directly to extracting the M and HQ bins
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/CerealRye/assembly/megahit_out/CerealRye_final.contigs_2500.fa.metabat-bins/checkm/results.txt
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Control/assembly/megahit_out/Control_final.contigs_2500.fa.metabat-bins/checkm/results.txt
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Rapeseed/assembly/megahit_out/Rapeseed_final.contigs_2500.fa.metabat-bins/checkm/results.txt
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Sorghum/assembly/megahit_out/Sorghum_final.contigs_2500.fa.metabat-bins/checkm/results.txt
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/HairyVetch/assembly/megahit_out/HairyVetch_final.contigs_2500.fa.metabat-bins/checkm/results.txt


# Record good bins in new excel file: "MAG processing manual code"


##########################################################################################
######### Evaluating MAG taxonomy (on MedHighQualityBins)


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins

#this will create SAM summary files for the taxonomy assigned for the two bins
gtdbtk classify_wf -x fa --genome_dir . --out_dir gtdb_03282023 --cpus 20


#save, submitted in MedHighQualityBins/slurm march 28 at 11am 
gtdb.sh

###########################################################################################
######### Clean up files
#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/"$element"/assembly/megahit_out

# Cleaning up
rm "$element"_final.contigs_2500_mapped.sam
rm "$element"_final.contigs_2500_mapped.sorted.bam
rm "$element"_final.contigs_2500_mapped99per.sorted.bam
gzip "$element"_final.contigs_2500_mapped.bam

cd ../../../
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

bash clean_finalcontigfiles_loop.sh list.txt

#save, submitted march 13 at 11am
clean_finalcontigfiles.sh


###########################################################################################
####### Annotations using DRAM on combined MQ and HQ bins (CORRECT WAY)

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=128gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins

DRAM.py annotate -i '*fa' -o  DRAM_1.4.4_03272023 --min_contig_size 2500 --threads 20

# distill the annotations 
DRAM.py distill -i DRAM_1.4.4_03272023/annotations.tsv -o DRAM_1.4.4_03272023/distill

#saved in MedHighQualityBins/slurm, started march 27 at 3:30pm - 5am 
dram.sh

###########################################################################################
####### Annotations using DRAM + custom cover crop root exudate module on MedHighQualityBins


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

DRAM.py distill -i DRAM_1.4.4_03272023/annotations.tsv -o DRAM_1.4.4_03272023/distill_wCustomDistillate --custom_distillate custom_distillate_coverCrop_exudates.txt

#save, submitted april 3 4pm
dram_custom_distillate.sh


###########################################################################################
####### build a genome database (in MedHighQualityBins)

### dRep
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins

dRep dereplicate dRep_v2.6.2 -p 20 -comp 50 -con 10 -g ./*fa

#saved in MedHighQualityBins/slurm, began march 27 9pm - midnight
sbatch dRep.sh

#rename the dRep_v2.6.2 to have the number of MAGs used = dRep_v2.6.2_43MAGs

ls -dq *fa | wc -l
# there were 89 M and HQ bins before dRep
# there are 43 bins after dRep


###########################################################################################
#####Assessing the quality of your MAG database (in MedHighQualityBins)

### Renames bins to match DRAM, cd into MedHighQualityBins/dRep_v2.6.2_*MAGs/dereplicated_genomes/ first. 
### The scaffold names (the headers in the files) will be changed from  (e.g.): k121_2164525 to HairyVetch_bin.11_k121_2164525
rename_bins_like_dram.py -i '*fa' -o genome_renamed


# for the next step, all files need to be concatenated, so i just concatenated the coassmebly files
#copy the trimmed coassembly reads into a new trimmed reads folder in MedHighQualityBins
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/Control/assembly/concat_R*_trimmed.fastq /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/trimmed_reads

#concatenate the forward and reverse reads
cat Control_concat_R1_trimmed.fastq Rapeseed_concat_R1_trimmed.fastq CerealRye_concat_R1_trimmed.fastq Sorghum_concat_R1_trimmed.fastq HairyVetch_concat_R1_trimmed.fastq > concat_R1_coassembly_trimmed.fastq
cat Control_concat_R2_trimmed.fastq Rapeseed_concat_R2_trimmed.fastq CerealRye_concat_R2_trimmed.fastq Sorghum_concat_R2_trimmed.fastq HairyVetch_concat_R2_trimmed.fastq > concat_R2_coassembly_trimmed.fastq

###########################################################################################
### Concatenate MAG database and map a set of metagenomic reads
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed
cat *fa > MAGs_concatenated_dRep.fa # cat this file with 2022 file and drep again, 

bbmap.sh -Xmx48G threads=20 overwrite=t ref=MAGs_concatenated_dRep.fa in1=../../../trimmed_reads/concat_R1_coassembly_trimmed.fastq in2=../../../trimmed_reads/concat_R2_coassembly_trimmed.fastq outm=final.contigs_2500_mapped_interleaved.fastq #reads that mapped to the genome database

#SAVE, submitted march 28 at 2pm - march 30 5am
sbatch map_to_mags.sh

#To count the number of reads that mapped use the following command, 
# see map_toMags_%reads_Mapped excel file for calculations, percent mapped was 22%
wc -l final.contigs_2500_mapped_interleaved.fastq


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

# put outu1 and outu2 into megahit again and assemble
#SAVE, submitted april 5 10 am
sbatch extract_unmapped_reads.sh

# ... see the adding_unmapped_reads.sh file for the rest of the code used for this data analysis



###########################################################################################
##### Use singleM to profile microbial lineages represented in metagenomes (uses raw reads as specified by singleM github) on MedHighQualityBins
https://wwood.github.io/singlem/usage/pipe

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/

mkdir singleM
cd singleM

source /opt/Miniconda2/miniconda2/bin/activate singlem
singlem pipe --forward ../trimmed_reads/concat_R1_coassembly_trimmed.fastq --reverse ../trimmed_reads/concat_R2_coassembly_trimmed.fastq --otu_table singleM_output_table.tsv --output_extras --threads 20

#SAVE, submitted to coassembly/MedHighQualityBins/slurm march 30 6am
sbatch singleM.sh


###########################################################################################
##### Mapping to determine abundance with metagenomic reads (in MedHighQualityBins)

# you will use these steps to see where your MAGs are coming from, this needs to be done on every sample.
# i.e., this allows you to map bins back to samples to get timeseries data and see ALL samples where the bin in located. 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins

# make directory for bowtie database and build it
mkdir bowtie_DB
cd bowtie_DB
bowtie2-build ../DRAM_1.4.4_03272023/scaffolds.fna 99perMAG_DB --threads 15

# now map reads to the db and output a SAM file
cd ../

# map
bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 15 -x bowtie_DB/99perMAG_DB -S bowtie_DB/mapped_99perMAGs.sam -1 trimmed_reads/concat_R1_coassembly_trimmed.fastq -2 trimmed_reads/concat_R2_coassembly_trimmed.fastq

#submit, march 29 9pm - 5am
bowtie_map.sh

#Step 2: Convert SAM to BAM, filter, sort
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/bowtie_DB

# convert SAM to BAM
samtools view -@ 15 -bS mapped_99perMAGs.sam > mapped_99perMAGs.bam

# filter BAM to mappings 97% id, removing multi-mappings and keeping only proper pairs
reformat.sh in=mapped_99perMAGs.bam out=97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t

# position sort filtered BAM
samtools sort -@ 15 -o 97peridmapped_99perMAGs_POSSORT.bam 97peridmapped_99perMAGs.bam

#submit, march 30 1045am- 1pm 
sam2bam_filter_sort.sh


# Step 3: Determine MAG coverage.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins

# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed --bam-files bowtie_DB/*POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0 -m reads_per_base --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt

# output MAGs with min-covered_fraction >0.97
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed --bam-files bowtie_DB/*POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0.75 --output-file coverm_min75.txt &> min75_stats.txt

# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed --bam-files bowtie_DB/*POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 -m trimmed_mean --output-file coverm_trimmed_mean.txt &> trimmed_mean_stats.txt

# submited to MedHighQualityBins/slurm april 3 1:30pm 
coverM.sh

# move the coverM files to their own directory. Talk to Kayla about the coverage calculations.
# Finally, when you are satisfied with your data, delete all SAM files and redundant BAM files.


###########################################################################################
##### Mapping to determine abundance with metagenomic reads (in MedHighQualityBins)
##### REDO with bowtie mapping to EACH SAMPLE

# you will use these steps to see where your MAGs are coming from, this needs to be done on every sample.
# i.e., this allows you to map bins back to samples to get timeseries data and see ALL samples where the bin in located. 

#!/bin/bash
for element in $(<$1)
do
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins

# map
bowtie2 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 -p 15 -x bowtie_DB/99perMAG_DB -S bowtie_DB/"$element"_mapped_99perMAGs.sam -1 ../../"$element"/trimmed_reads/"$element"_R1_trimmed.fastq -2 ../../"$element"/trimmed_reads/"$element"_R2_trimmed.fastq

cd .
done

#save
bowtie_map_redo_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

bash bowtie_map_redo_loop.sh list.txt

#before submitting, move the list.txt file to the MedHighQualityBins/slurm
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/slurm/list.txt /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/slurm

#submit, april 5 11am 
bowtie_map_redo.sh

#Step 2: Convert SAM to BAM, filter, sort

# need to edit still

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/bowtie_DB

# convert SAM to BAM
samtools view -@ 15 -bS mapped_99perMAGs.sam > mapped_99perMAGs.bam

# filter BAM to mappings 97% id, removing multi-mappings and keeping only proper pairs
reformat.sh in=mapped_99perMAGs.bam out=97peridmapped_99perMAGs.bam minidfilter=0.97 primaryonly=t pairedonly=t

# position sort filtered BAM
samtools sort -@ 15 -o 97peridmapped_99perMAGs_POSSORT.bam 97peridmapped_99perMAGs.bam

#submit, march 30 1045am- 1pm 
sam2bam_filter_sort.sh


# Step 3: Determine MAG coverage.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins

# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed --bam-files bowtie_DB/*POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0 -m reads_per_base --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt

# output MAGs with min-covered_fraction >0.97
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed --bam-files bowtie_DB/*POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0.75 --output-file coverm_min75.txt &> min75_stats.txt

# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed --bam-files bowtie_DB/*POSSORT.bam --threads 15 --min-read-percent-identity-pair 0.97 -m trimmed_mean --output-file coverm_trimmed_mean.txt &> trimmed_mean_stats.txt

# submited to MedHighQualityBins/slurm april 3 1:30pm 
coverM.sh



###########################################################################################
##### Mapping metatranscriptomics with metagenomic 

# Next step would be mapping metaT reads




###########################################################################################
##### Extract amino acid fasta file for metaproteomics (combine the sorghum microcosms from 
##### 2019 with these new enrichments)

#copy all the coassembly genes.faa files into one location so i can concatenate into one file 
cp /home/projects/Agribiome/Lindstrom_McGivern_2022_RE_encrichments/MetaG/DRAM_3.18.21/genes.faa /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/metaP/gene_files_2019_2022_microcosms_combined
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/metaP/gene_files/genes_concat.faa /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/metaP/gene_files_2019_2022_microcosms_combined
mv genes.faa genes_2019_microcosms.faa
cat genes_2019_microcosms.faa genes_concat.faa > genes_concat_2019_2022_microcosms.faa

#once the below code was finished, i deleted the individual gene files and only left the genes_concat.faa file.


#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=128gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/metaP/gene_files_2019_2022_microcosms_combined
cd-hit -i genes_concat_2019_2022_microcosms.faa -o gene_database_2019_2022_microcosms.faa -c 1 -M 128000 -T 20 -d 0 -B 1 -bak 1

#save, took about 3 mins
gene_database_extract_2019_2022_microcosms.sh


###########################################################################################
# dREP on ALL bins for a sunburst plot like the one in the paper. 

#copy all bins into new dir

#!/bin/bash
for element in $(<$1)
do
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/"$element"/assembly/megahit_out/"$element"_final.contigs_2500.fa.metabat-bins/*.fa /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/all_bins

cd ../../../../
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

sbatch copy_MAGs.sh

### dRep, dereplicate ALL bins (low, medium, and high quality. so change -comp from 50 to 10% completion)

## needs to be rerun after the USDA CC meeting

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/all_bins
dRep dereplicate dRep_v2.6.2 -p 20 -comp 49 -con 10 -g ./*fa

sbatch dRep.sh

#rename the dRep_v2.6.2 to have the number of MAGs used = 
ls -dq *fa | wc -l
# 277 before dRep (when using the 10% completion parameter)
# 108 after dRep (when using the 10% completion parameter)

# assign taxonomy to all dereplicated bins

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --time=28-00:00:00
#SBATCH --mem=400gb 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/all_bins/dRep_v2.6.2/dereplicated_genomes

#this will create SAM summary files for the taxonomy assigned for the two bins
gtdbtk classify_wf -x fa --genome_dir . --out_dir gtdb_03292023 --cpus 20


#save, submitted in MedHighQualityBins/slurm march 28 at 11am 
gtdb.sh















