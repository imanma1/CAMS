#!/usr/bin/env bash
# Submit from the repository root: sbatch real_data/submit_support2.sh
# Native Experiment A only; one split seed per isolated array task.
#SBATCH --account=def-zjunxi_cpu
#SBATCH --job-name=cams_support2_A
#SBATCH --time=03:00:02
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --array=1-100%10
#SBATCH --output=support2_%A_%a.log

set -euo pipefail
task_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
task_repo_root="${CAMS_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$task_script_dir/.." && pwd)}}"
if [[ "${SUPPORT2_LOAD_MODULES:-1}" == "1" ]] && type module &>/dev/null; then
  module load StdEnv/2023
  module load "${SUPPORT2_R_MODULE:-r/4.5.0}"
fi
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/support2-library}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_seed="${TASK_SEED:-${SLURM_ARRAY_TASK_ID:-1}}"
if [[ ! "$task_seed" =~ ^[0-9]+$ ]]; then echo "Invalid task seed: $task_seed" >&2; exit 2; fi
task_results_root="${RESULTS_ROOT:-$task_repo_root/results/support2/experiment_A_${SLURM_ARRAY_JOB_ID:-local}}"
task_csv="${SUPPORT2_CSV:-$task_repo_root/real_data/data/support2.csv}"
if [[ ! -f "$task_csv" ]]; then
  echo "SUPPORT2 CSV missing: run bash real_data/setup_support2.sh on the login node first." >&2
  exit 2
fi
if [[ ! -f "$task_repo_root/real_data/run_support2.R" ]]; then
  echo "CAMS_ROOT/SLURM_SUBMIT_DIR must be the repository root." >&2; exit 2
fi
mkdir -p -- "$task_results_root/tasks"
"${RSCRIPT_BIN:-Rscript}" "$task_repo_root/real_data/run_support2.R" \
  --scenario=native --seeds="$task_seed" --csv="$task_csv" \
  --horizon="${SUPPORT2_HORIZON:-730}" --alpha="${SUPPORT2_ALPHA:-0.1}" \
  --feature-set="${SUPPORT2_FEATURE_SET:-entry}" --landmark-days="${SUPPORT2_LANDMARK_DAYS:-0}" \
  --probability-floor="${SUPPORT2_PROBABILITY_FLOOR:-0.05}" \
  --grid-step="${SUPPORT2_GRID_STEP:-0.001}" --integration-steps="${SUPPORT2_INTEGRATION_STEPS:-256}" \
  --dft-fixed-cutoff="${SUPPORT2_DFT_FIXED_CUTOFF:-auto}" \
  --aft-ridge-fallback="${SUPPORT2_AFT_RIDGE_FALLBACK:-true}" --aft-ridge-theta="${SUPPORT2_AFT_RIDGE_THETA:-0.1}" \
  --include-dft=true --include-cqr=true --include-subgroups=true --strict-methods=true \
  --output="$task_results_root/tasks/seed_$task_seed"
