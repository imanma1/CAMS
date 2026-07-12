#!/bin/bash

#SBATCH --account=def-zjunxi_cpu

#SBATCH --job-name=conformal_sim

#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=48

#SBATCH --mem=192G

#SBATCH --array=1-46:5

#SBATCH --output=sim_log_%A_%a.txt



module load StdEnv/2023

module load r/4.5.0



cd ~/scratch/CAMS/simulations


SETTING_LIST="cams_vs_vanilla_lower_tail_hd_mild,cams_vs_vanilla_lower_tail_hd_main,cams_vs_vanilla_lower_tail_hd_strong"
SC_METHODS="oracle,km,aft_lognormal,km_x1,km"
AUGMENTATION_METHODS="correct,wrong"
HOMOSCEDASTIC_EVENT=1
SC_NTREE=1000
NUM_RUNS=5

Rscript run_sim.R \
  "$SETTING_LIST" \
  "$SLURM_ARRAY_TASK_ID" \
  "$SC_METHODS" \
  0 \
  "$AUGMENTATION_METHODS" \
  "$HOMOSCEDASTIC_EVENT" \
  "$SC_NTREE" \
  "$NUM_RUNS"
