#!/bin/bash

#SBATCH --account=def-zjunxi_cpu
#SBATCH --job-name=cams_intersection_screen
#SBATCH --time=3:02:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --array=1-100:10
#SBATCH --output=sim_log_%A_%a.txt

module load StdEnv/2023
module load r/4.5.0

cd ~/scratch/CAMS/simulations

SETTING_LIST="cov_cens_dr_tail_hetero"

SC_METHODS="power_oracle"
ONLY_CAMS=0
AUGMENTATION_METHODS="oracle_event"
HOMOSCEDASTIC_EVENT=0
SC_NTREE=1000
NUM_RUNS=10
USE_INTERSECTIONAL_R=1

Rscript run_sim.R \
  "$SETTING_LIST" \
  "$SLURM_ARRAY_TASK_ID" \
  "$SC_METHODS" \
  "$ONLY_CAMS" \
  "$AUGMENTATION_METHODS" \
  "$HOMOSCEDASTIC_EVENT" \
  "$SC_NTREE" \
  "$NUM_RUNS" \
  "$USE_INTERSECTIONAL_R"