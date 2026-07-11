#!/bin/bash

#SBATCH --account=def-zjunxi_cpu

#SBATCH --job-name=conformal_sim

#SBATCH --time=3:02:00

#SBATCH --cpus-per-task=48

#SBATCH --mem=192G

#SBATCH --array=1-41:10

#SBATCH --output=sim_log_%A_%a.txt



module load StdEnv/2023

module load r/4.5.0



cd ~/scratch/CAMS/simulations


SETTING_LIST="cams_vs_vanilla_lower_tail_hd_mild,cams_vs_vanilla_lower_tail_hd_main,cams_vs_vanilla_lower_tail_hd_strong"


Rscript run_sim.R $SETTING_LIST $SLURM_ARRAY_TASK_ID 1 0