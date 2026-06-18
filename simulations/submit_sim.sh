#!/bin/bash

#SBATCH --account=def-zjunxi_cpu

#SBATCH --job-name=conformal_sim

#SBATCH --time=8:00:00

#SBATCH --cpus-per-task=192

#SBATCH --mem=748G

#SBATCH --array=1-49:2

#SBATCH --output=sim_log_%A_%a.txt



module load StdEnv/2023

module load r/4.5.0



cd ~/scratch/CAMS/simulations



Rscript run_sim.R "starve_hetero_high_dim" $SLURM_ARRAY_TASK_ID 0