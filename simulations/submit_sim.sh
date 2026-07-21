#!/bin/bash

#SBATCH --account=def-zjunxi_cpu
#SBATCH --job-name=cams_intersection_screen
#SBATCH --time=3:02:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --array=1-10:1
#SBATCH --output=sim_log_%A_%a.txt

module load StdEnv/2023
module load r/4.5.0

cd ~/scratch/CAMS/simulations

SETTING_LIST="cams_vs_vanilla_lower_tail_hd_mild,cams_vs_vanilla_lower_tail_hd_main,cams_vs_vanilla_lower_tail_hd_strong,anti_aligned_cens,weibull_aft_anti_cens,var_shift_heavy_cens,high_survival_heavy_cens,surv_misspec"

SC_METHODS="oracle"
ONLY_CAMS=0
AUGMENTATION_METHODS="same"
HOMOSCEDASTIC_EVENT=0
SC_NTREE=1000
NUM_RUNS=1
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