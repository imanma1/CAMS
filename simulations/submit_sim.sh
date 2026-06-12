#!/bin/bash

#SBATCH --job-name=conformal_sim

#SBATCH --time=00:30:00

#SBATCH --cpus-per-task=48

#SBATCH --mem=192G

#SBATCH --output=sim_log_%j.txt



module load StdEnv/2023

module load r/4.5.0



cd ~/scratch/CAMS/simulations



Rscript run_sim.R 1
