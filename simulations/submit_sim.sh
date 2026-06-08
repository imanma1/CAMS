#!/bin/bash
#SBATCH --account=def-zjunxi_cpu
#SBATCH --job-name=conformal_sim
#SBATCH --time=10:30:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=512G
#SBATCH --array=1-46:5
#SBATCH --output=sim_log_%A_%a.txt

module load StdEnv/2023
module load r/4.5.0

cd ~/scratch/CAMS/simulations

Rscript run_sim.R $SLURM_ARRAY_TASK_ID