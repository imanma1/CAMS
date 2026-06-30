#!/bin/bash

#SBATCH --account=def-zjunxi_cpu

#SBATCH --job-name=conformal_sim

#SBATCH --time=6:15:00

#SBATCH --cpus-per-task=48

#SBATCH --mem=192G

#SBATCH --array=1-46:5

#SBATCH --output=sim_log_%A_%a.txt



module load StdEnv/2023

module load r/4.5.0



cd ~/scratch/CAMS/simulations


SETTING_LIST="homo_cens,cov_cens,prot_cens,heavy_prot_cens,heavy_inter_cens,surv_misspec,cens_misspec,simul_misspec,complex_surv,var_shift_heavy_cens,starve_hetero,starve_hetero_high_dim"


Rscript run_sim.R $SETTING_LIST $SLURM_ARRAY_TASK_ID 1 0