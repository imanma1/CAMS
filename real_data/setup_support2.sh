#!/usr/bin/env bash
set -euo pipefail
task_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
task_repo_root="${CAMS_ROOT:-$(cd -- "$task_script_dir/.." && pwd)}"
if [[ "${SUPPORT2_LOAD_MODULES:-1}" == "1" ]] && type module &>/dev/null; then
  module load StdEnv/2023
  module load "${SUPPORT2_R_MODULE:-r/4.5.0}"
fi
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/support2-library}"
mkdir -p -- "$R_LIBS_USER"
"${RSCRIPT_BIN:-Rscript}" "$task_repo_root/real_data/setup_support2.R"
