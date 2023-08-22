#!/bin/bash
#SBATCH --job-name=HILIC_0.20minFrac_CCREM_MS1
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu

export TMPDIR=/scratch/alpine/lindsval@colostate.edu/tmp
#export TMP=$TMPDIR

# add directory where files are located
cd /scratch/alpine/lindsval@colostate.edu/HILIC

# clear any existing modules
module purge
source /curc/sw/anaconda3/latest
conda activate lc-msms

# add R script
Rscript --vanilla xcms.ramclustr.mse_val_MS1_HILIC_0.20minFrac.R
