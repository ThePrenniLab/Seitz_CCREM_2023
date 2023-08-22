#########################################################################################
##### RUN PRODIGAL ON EACH ASSEMBLY TO EXTRACT/CALL GENES, EXPORT THOSE GENES 


# we need to run Prodigal on the list of genes from the genome database. this will be run on each gene, from EACH individual assembly (so yes each individual sample)
# you need to run on the contig.genes file for each sample.this can be run in a loop by Kayla didn't have an example to give. this will take some time, so start early

# Flags
# -a: Write protein translations to the selected file. 
# -d:  Write nucleotide sequences of genes to the selected file.



#!/bin/bash
#$1 list of individual sample assembly names

for element in $(<$1)
do

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/"$element"/assembly/megahit/megahit_out
prodigal -i "$element"_final.contigs_2500.fa -a "$element"_contig.genes.faa -d "$element"_contigs.genes.fna -p meta -m -o "$element"_prodigal_output
cd ../../../../

done

#save
prodigal_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=28-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#Call genes using prodigal for metagenomes (2 tasks, 50gb). this process will take a while since prodigal cannot be multithreaded so we only use 2 tasks.

bash prodigal_loop.sh list.txt

#save, submitted may 15 945am, done by 12? not sure that worked correctly
prodigal.sh

#submitted in 
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/slurm


#########################################################################################
##### RUN PRODIGAL ON COASSEMBLY TO EXTRACT/CALL GENES, EXPORT THOSE GENES 


## will need to be rerun on the iterative coaasmbly too...

#!/bin/bash
#$1 list of individual sample assembly names

for element in $(<$1)
do

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/"$element"/assembly/megahit_out
prodigal -i "$element"_final.contigs_2500.fa -a "$element"_contig.genes.faa -d "$element"_contigs.genes.fna -p meta -m -o "$element"_prodigal_output
cd ../../../

done

#save
prodigal_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=28-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

#Call genes using prodigal for metagenomes (2 tasks, 50gb). this process will take a while since prodigal cannot be multithreaded so we only use 2 tasks.

bash prodigal_loop.sh list.txt

#save, submitted may 16 230pm, 4pm
prodigal.sh

#submitted in 
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/slurm

#########################################################################################
##### RUN PRODIGAL ON ITERATIVE COASSEMBLY TO EXTRACT/CALL GENES, EXPORT THOSE GENES 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=28-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/iterative_coassembly/megahit_out

prodigal -i final.contigs_2500.fa -a iterative_coA_contig.genes.faa -d iterative_coA_contigs.genes.fna -p meta -m -o iterative_coA_prodigal_output


#save, submitted june 2 as 145
prodigal.sh


#submitted in 
cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/iterative_coassembly/slurm



############


#NEEDS TO BE RUN


############


#########################################################################################
##### ANNOTATE GENE DATABASE GENES USING DRAM

# we need to run DRAM on the list of bins from the genome database. this will be run on each bin, so it can be started now
# and rerun when we add more bins from the second assembly

# can wait on this (low priority)

# make a list of bin names
bin_list.txt #(need to make still)



#DRAM loop for annotation of BINS in MAG database

#!/bin/bash
#$1 list of bin names

for element in $(<$1)
do
DRAM.py annotate -i "$element" -o  "$element"_DRAM --min_contig_size 2500 --threads 20
done

#save
dram_annotate_genes_loop.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=128gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

bash dram_annotate_genes_loop.sh bin_list.txt

#########################################################################################
# Why are we doing this? 

# so we run prodigal on all possible genes from the individual sample assemblies (so every possible gene we assembled), but remember that some of those 
# genes didn't make it a) into a genome, or b) into a genome that passed the bin QCs to make it in the genome database. We therefore cannot say that a gene or a process "isn't"
# there just because it didn't make it into the genome database. Therefor we use prodigal to extract all possible genes, then we will annotate those genes separately (not part 
# of a M/HQ bin) to give an idea of all the possible genes we have so that when a reviewer asks for it, we already have it. Some reviewers are picky about presense/absense	
# of genes/processes so this is out way of proving these genes are present or not.




















