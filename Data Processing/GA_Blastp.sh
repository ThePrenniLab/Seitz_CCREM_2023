#BLAST GA operons to MAGs

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/blastp_gibberellin

#blastp
makeblastdb -in /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.faa -dbtype prot -title gibberellin
blastp -query GA_genes.faa -db /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.faa -out gibberellin_genes.txt -outfmt 6
 