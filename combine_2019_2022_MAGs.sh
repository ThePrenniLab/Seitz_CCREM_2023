#combine MAGs from 2019 and 2022 M/HQ bins to make a new database 

mkdir microcosm_MAG_database_2019_2022
cd microcosm_MAG_database_2019_2022/

#copy bins
cp /home/projects/Agribiome/RootExudate_Prenni/MetaG/coassembly/MedHighQualityBins/dRep_v2.6.2_43MAGs/dereplicated_genomes/genome_renamed/*fa /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022

cp /home/projects/Agribiome/Lindstrom_McGivern_2022_RE_encrichments/MetaG/all_MAGs_48_10/drep/dRep_out_SH/dereplicated_genomes/*fa /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022


#dRep all bins

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022
dRep dereplicate dRep_v2.6.2 -p 20 -comp 49 -con 10 -g ./*fa

sbatch dRep.sh #submitted april 11 5pm

#rename the dRep_v2.6.2 to have the number of MAGs used = 
ls -dq *fa | wc -l

#  before dRep (when using the 49% completion parameter)
#  after dRep (when using the 49% completion parameter)

###########################################################################################
### Concatenate MAG database for Gustavo
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/Agribiome/RootExudate_Prenni/MetaG/microcosm_MAG_database_2019_2022
cat *fa > 2019_2022_MAGs_concatenated_dRep.fa 




















