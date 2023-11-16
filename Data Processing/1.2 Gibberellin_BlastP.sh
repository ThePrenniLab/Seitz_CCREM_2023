#BLAST GA operons to MAGs


#!/bin/bash

#used GA_genes.faa as query --> see 1.3 GA_genes.faa for that file

#blastp
makeblastdb -in /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.faa -dbtype prot -title gibberellin
blastp -query GA_genes.faa -db /home/projects/Agribiome/RootExudate_Prenni/MetaG/MAG_database/dRep_v2.6.2_326/dereplicated_genomes/DRAM_1.4.4_05222023/genes.faa -out gibberellin_genes.txt -outfmt 6


