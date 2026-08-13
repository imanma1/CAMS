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


########################################
## Per-setting configuration
########################################

# Format:
#
# setting:n_train:n_calib:n_test:bernoulli_prob
#
# Separate multiple settings with ";"
#
# Here:
#   - first two settings use n_train = 2000
#   - rare setting uses n_train = 200
#
# Adjust n_calib / n_test / bernoulli_prob
# independently for any setting as needed.

SETTING_CONFIGS="\
intersection_scale_shift_hd:2000:8000:20000:0.10;\
intersection_early_event_mixture_hd:2000:8000:20000:0.10;\
rare_intersection_scale_shift_ld:200:800:20000:0.10"


########################################
## Other simulation configuration
########################################

SC_METHODS="aft_lognormal"

ONLY_CAMS=0

AUGMENTATION_METHODS="same"

HOMOSCEDASTIC_EVENT=0

SC_NTREE=1000

NUM_RUNS=10

USE_INTERSECTIONAL_R=1


########################################
## Run
########################################

Rscript run_sim.R \
  "$SETTING_CONFIGS" \
  "$SLURM_ARRAY_TASK_ID" \
  "$SC_METHODS" \
  "$ONLY_CAMS" \
  "$AUGMENTATION_METHODS" \
  "$HOMOSCEDASTIC_EVENT" \
  "$SC_NTREE" \
  "$NUM_RUNS" \
  "$USE_INTERSECTIONAL_R"