# SUPPORT2 implementation

This pipeline reuses the existing CAMS calibration code. It replaces the
simulation data generator with a cohort loader and stratified patient splits,
and replaces latent-event evaluation with exact restricted-outcome or
censoring-aware evaluation. `simulations/simu.R` remains unchanged because its
DFT baselines and `T_test` evaluation cannot simply be reused for ordinary
right-censored records.

## Run

From the repository root, install the only required packages:

```r
install.packages(c("survival", "quantreg"))
```

Run Experiment A with original censoring and all comparisons:

```bash
Rscript real_data/run_support2.R --scenario=native --horizon=730 --seeds=1
```

Run five splits:

```bash
Rscript real_data/run_support2.R --scenario=native --horizon=730 --seeds=1,2,3,4,5
```

Data is downloaded from Vanderbilt into `real_data/data/support2.csv` if absent.
To use an existing **original** CSV:

```bash
Rscript real_data/run_support2.R --csv=/path/to/support2.csv --scenario=all --seeds=1
```

A new timestamped directory is created under `results/support2/` by default.
Use `--output=/path/to/results` to choose one. Existing seed/scenario directories
are never overwritten. The same output root can be resumed with **new seeds**;
previous aggregate summaries are retained.

```bash
Rscript real_data/run_support2.R --scenario=all --seeds=1 --output=results/support2/my_study
Rscript real_data/run_support2.R --scenario=all --seeds=2,3 --output=results/support2/my_study
```

For a quick smoke run, use `--grid-step=0.01 --integration-steps=128`.
The main defaults use the full 0.001 calibration grid and 256 integration
intervals. Keep these choices recorded when comparing runs.

## Defaults

The R CLI defaults to `native` at 180 days for compatibility with the earlier
restricted-horizon pilot. The HPC submission script explicitly selects native
Experiment A at **730 days**, using the original full observed follow-up for
fitting and calibration. A horizon caps candidates; it does not replace native
observations with administratively shortened records. Override it using
`SUPPORT2_HORIZON=180` or `--horizon=180`.

| Setting | Default |
|---|---|
| Cohort | Original 9,105-person SUPPORT2 CSV |
| Target event | All-cause death (`death`) |
| Observed time | `d.time`, in days |
| Prediction time | Study entry, using the entry predictor allowlist |
| Audit/decision groups | Female/male crossed with age <65 / >=65 |
| Split | 50% fitting, 25% calibration, 25% test, stratified by group |
| Target miscoverage | 0.10 |
| Lower-bound ceiling | 180 days |
| Candidate/event nuisance | Pooled Weibull AFT |
| Censoring nuisance | Lognormal AFT on `Surv(censored_T, 1-event)` |
| Censoring probability floor | 0.05, configurable |
| AIPCW integration | Optional continuous-hazard dense quadrature |

The entry allowlist is `age`, `sex`, `race`, `dzgroup`, `num.co`, `diabetes`,
and `dementia`. Encodings, numeric scales and imputation values are learned
only on the fitting split. Unseen categories map to the learned fitting
reference; missing categories are explicit if present during fitting.
Rank-deficient columns are removed using fitting data only. Original patient
IDs are preserved for reproducibility and are excluded from predictors.

To include day-3 physiology, request the explicit landmark:

```bash
Rscript real_data/run_support2.R --feature-set=day3 --landmark-days=3 --scenario=all --seeds=1
```

This retains records with follow-up beyond day 3 and subtracts three days from
follow-up. It changes the target cohort to people alive and observed at that
landmark. The pipeline interprets `d.time` as follow-up from study entry; verify
that origin and the measurement window against study documentation before
making a clinical interpretation. It excludes subsequent costs, disability,
prognostic scores and physician survival predictions.

## Three scenarios

1. `native`: original observed times and event indicators. Fitting uses the
   full original follow-up; candidates are capped at the chosen horizon during
   calibration and prediction.
2. `added_common`: restrict follow-up to the horizon and independently generate
   an exponential added censoring time with probability 0.8 of exceeding it.
3. `added_group`: same procedure, with that probability set to 0.8 for younger
   female patients, 0.65 for older female patients, 0.65 for younger male patients,
   and 0.5 for older male patients.

Added censoring is conditional on allowed group covariates, using fresh random
draws; it is not selected using survival outcomes. The saved original outcomes
are unavailable to fitting and calibration. Generated censoring times are saved
in the manifest for audit, but are not passed to CAMS or its nuisance models.
The probabilities above are **not** final observed censoring percentages:
administrative censoring at the horizon and earlier deaths affect those rates.

These are semi-synthetic early-censoring experiments on real patients, separate
from the native application. The loader rejects an added-censoring experiment
if original censored outcomes are not known through the requested horizon.
All original censored follow-up exceeds 180 days in this data version.

## Methods

- `Weibull-uncalibrated`: pooled alpha-quantile, capped at the horizon.
- `Marginal-IPCW`, `Marginal-AIPCW`: one pooled calibration level.
- `CAMS-IPCW`, `CAMS-AIPCW`: the **same pooled fitted candidate model**,
  with separate calibration levels in the four groups.
- `Separate-IPCW`, `Separate-AIPCW`: event and censoring models fitted separately
  in each group, followed by calibration in that group.
- `Adaptive DFT-T (joint)` and `Adaptive DFT-T (subgroup)`.
- `Adaptive DFT-CT (joint)` and `Adaptive DFT-CT (subgroup)`.
- `DFT-fixed (joint)` and `DFT-fixed (subgroup)`.
- `CQR (joint)` and `CQR (subgroup)`: the existing one-sided quantile-regression
  CQR fitted to observed time, with the opt-in empirical order-statistic correction.

Joint methods pool fitting and calibration records. Subgroup methods fit and
calibrate independently in each of the four sex/age cells. Constant and
rank-deficient predictors are removed using each group's fitting data only.
The event and censoring AFT fits use up to 1,000 iterations and explicit
convergence checks. An unusable unpenalized AFT fit is retried with the fixed
ridge penalty theta = 0.1, using the same fitting records and predictors.
This is a prespecified numerical fallback, decided without calibration/test
outcomes; its use and the original cause appear in method status and fitted
model metadata. It handles, for example, a rare missing-category indicator with
no observed censoring events in a subgroup. It is not an assumption that this
working model is correctly specified. Disable it with
`--aft-ridge-fallback=false`, or set `--aft-ridge-theta=0.1` explicitly for a
prespecified sensitivity study. A fit that remains unusable yields missing
predictions; it is never replaced by an uncalibrated or falsely successful zero
bound.

Separate methods remove constant/linearly dependent predictors within each
fitting group. Numerical failures are recorded and yield missing predictions,
not zero-bound successes. A zero bound from no feasible calibration candidate
is retained as a genuine fallback. CQR can emit a nonunique-solution warning
with tied outcomes; finite returned predictions remain recorded.

Use `--include-subgroups=false`, `--include-dft=false`, or `--include-cqr=false` to reduce the pilot.
The latter removes the `quantreg` requirement.

## DFT adaptation for original right censoring

The original DFT routines require potential censoring times C for all patients.
SUPPORT2 supplies only Y and the event indicator. This adapter uses a conditional
imputation step inspired by Algorithm 1 of Sesia and Svetnik (2025):

- If a record is censored, retain its observed C = Y.
- If a death is observed, draw C from the fitting-only lognormal censoring model,
  conditional on C > Y. In survival form, draw U uniformly on (0,1) and solve
  G(C | X) = U G(Y | X), using log probabilities for numerical stability.
- Never impute a test patient's censoring time or any patient's event time.

Joint DFT uses pooled nuisance models. Subgroup DFT uses each cell's nuisance
models. All three DFT comparisons within a scope/group share the same draws.
The default imputation seed derives from the split seed; `--imputation-seed=123`
allows an explicit sensitivity run. Original Y/event records remain unchanged
and are used for evaluation and all other methods.

These are **imputation-adapted DFT comparisons**, not unmodified known-C Type-I
DFT. Conditional independence and censoring-model accuracy matter. This code
does not claim the double robustness of DR-COSARC: it does not add that method's
extra uncalibrated-quantile cap. References:
[Gui et al.](https://arxiv.org/abs/2211.01227) and
[Sesia and Svetnik](https://proceedings.mlr.press/v267/sesia25a.html).

Adaptive DFT-T uses ordered Weibull event quantiles. DFT-CT additionally caps
each candidate at the fitted conditional (1 - probability-floor) censoring
quantile. Here the CT learner is **lognormal AFT**, rather than the quantile
forest in the simulation script. This choice keeps imputation, weights and
the CT cutoff coherent with one fitted censoring model, and requires only two
R packages. The candidate grid is explicit (default 0.001 through 0.999);
the historical weighted risk ratio and its monotone envelope are preserved.

DFT-fixed selects a constant cutoff on an internal split of the outer fitting
data, maximizing mean lower bounds over censoring-time quantile candidates.
Outer calibration/test outcomes never enter tuning. Use
`--dft-fixed-cutoff=180` for a prespecified constant. Calibration retains the
test-point weight at +Inf. Empty eligible sets yield a documented conservative
zero fallback. The fixed cutoff, effective sample sizes and clipping are saved.

The baseline adapter is separate from the original simulation utilities, so
their heavy automatic package installation and unrestricted parallelism are
not invoked in batch jobs.

## HPC: Experiment A

Copy the extracted repository to the cluster and run from its root:

```bash
# Login node, once: load R, install missing packages, and cache SUPPORT2.
bash real_data/setup_support2.sh

# 100 original-censoring splits, maximum 10 simultaneous jobs.
sbatch real_data/submit_support2.sh

# After the array finishes, replace JOBID with the returned array job ID.
module load StdEnv/2023 r/4.5.0
Rscript real_data/merge_support2.R \
  --input=results/support2/experiment_A_JOBID --expected-seeds=1:100
```

The submission defaults reuse the account and R module from `submit_sim.sh`:
`def-zjunxi_cpu`, `StdEnv/2023`, and `r/4.5.0`. Each seed requests one CPU,
8 GB RAM and one hour. Adjust the account/module if your allocation differs.
Submit from the repository root, or provide `CAMS_ROOT=/absolute/path/to/CAMS`.
Compute jobs require the cached CSV and installed packages; they never install
packages or download the dataset. Both scripts use `R_LIBS_USER`, defaulting
to `$HOME/R/support2-library`.

For three pilot array tasks:

```bash
sbatch --array=1-3%3 real_data/submit_support2.sh
```

Environment overrides include `RESULTS_ROOT`, `SUPPORT2_CSV`,
`SUPPORT2_HORIZON`, `SUPPORT2_ALPHA`, `SUPPORT2_FEATURE_SET`,
`SUPPORT2_LANDMARK_DAYS`, `SUPPORT2_PROBABILITY_FLOOR`,
`SUPPORT2_GRID_STEP`, `SUPPORT2_INTEGRATION_STEPS`,
`SUPPORT2_DFT_FIXED_CUTOFF`, `SUPPORT2_AFT_RIDGE_FALLBACK`,
`SUPPORT2_AFT_RIDGE_THETA` and `SUPPORT2_R_MODULE`.
Use `SUPPORT2_LOAD_MODULES=0` with an already configured R environment.

Run the same bash file locally, one seed at a time:

```bash
RESULTS_ROOT="$PWD/results/support2/local_A" TASK_SEED=1 \
  bash real_data/submit_support2.sh
RESULTS_ROOT="$PWD/results/support2/local_A" TASK_SEED=2 \
  bash real_data/submit_support2.sh
Rscript real_data/merge_support2.R \
  --input=results/support2/local_A --expected-seeds=1:2
```

Every task writes only its own `tasks/seed_N/` directory; no shared aggregate
is touched by concurrent workers. Seed/scenario directories are reserved before
fitting and never overwritten. A method failure saves available results and
returns a nonzero exit status under the default `--strict-methods=true`.
`--strict-methods=false` permits partial runs, while preserving failure records.

The merge step verifies matching code/data/settings and rejects duplicate seeds.
It writes completion counts, pooled and subgroup summaries, and method status.
With `--expected-seeds`, missing tasks are explicitly listed and produce a
nonzero exit status; available results are still saved. An interrupted task
without finalized files is marked incomplete. Re-run failed seeds into a fresh
result root, inspect their status, and keep duplicate seeds out of a combined
study. A result root from the earlier adapter version must not be resumed with
the expanded method list.

## Evaluation and interpretation

Every method is evaluated on the same test patients and audit groups. Reported
metrics include:

- `coverage_exact`: available when each original test outcome resolves whether
  T >= its predicted bound. It is **NA** if any original coverage indicator in
  that audit group is unknown. At the default 180-day ceiling it is directly
  observable for every method with finite predictions.
- `coverage_ipcw`, `coverage_aipcw`: common nuisance-based test evaluations,
  using pooled models fitted outside the test set. They are useful diagnostics,
  not observed truth for hidden outcomes. Estimation errors can produce values
  outside [0,1]; these values are not silently clipped.
- `coverage_observed_lower` and `coverage_observed_upper`: logical bounds on
  realized test coverage from the records given to methods.
- `unknown_indicator_fraction`: censored before the bound in those records.
- Mean/median bounds, zero/ceiling frequencies and relevant censoring clipping.

For added scenarios, exact evaluation uses original restricted outcomes while
the observable brackets describe the artificially censored records. For native
scenarios, both use the original records.

At 180 days the original cohort provides a clean calibration check, but has no
original early censoring. It cannot demonstrate censoring correction by itself;
that is the purpose of the added-censoring scenarios. To inspect the original
censoring application at a longer ceiling, use, for example:

```bash
Rscript real_data/run_support2.R --scenario=native --horizon=730 --seeds=1
```

Longer ceilings require censoring-support checks and may produce unknown exact
coverage indicators. No longer-horizon positivity guarantee is asserted here.

The four constraints are disjoint. Sex/age unions inherit their coverage; this
first adapter does **not** implement joint optimization for arbitrary overlapping
audit groups. Neither real-data event nor censoring working models are known
to be correctly specified. Repeated splits share patients and are not independent
clinical cohorts.

## What changed in existing utilities

`utils/cams.R` gains optional explicit group labels, a reusable fitted event
model, candidate cap, candidate grid, probability floor, and integration settings.
Historical defaults, hard-coded simulation groups and legacy integration remain.
The new runner opts into the new arguments.

The new dense integrator handles the observed censoring jump directly, and
integrates the compensator from zero to `min(Y,L)`. It uses the exact identity
`integral (1/G) dLambda = 1/G - 1` for continuous G and trapezoidal integration
for the other term. This covers hazard before the first censoring time, which
the legacy censoring-time-only grid can omit. It currently requires the smooth
lognormal AFT censoring model. It is not a general discrete-hazard/RSF extension.
Rounded observation times and the artificial administrative endpoint need a
consistent tie convention; the paper's continuous-time theorem does not itself
prove finite-sample validity of this discretized pilot.

`utils/cqr.R` gains `exact_quantile=TRUE`, selecting an empirical order statistic
from calibration scores plus Inf. Its default remains the previous interpolated
quantile, so old simulation calls retain their behavior.

## Saved outputs

Each scenario/seed directory contains:

- `metrics.csv`: overall and four-group results for all methods.
- `predictions.csv`: individual bounds and evaluation pseudo-losses.
- `calibration_diagnostics.csv`: selected levels, estimated risk and truncation.
- `method_status.csv`: success or failure for each comparison/group.
- `baseline_diagnostics.csv`: DFT selected levels/cutoffs, effective sample sizes,
  clipping, method runtimes, recorded warnings and conservative fallbacks.
- `dft_imputation_manifest.csv`: fitting/calibration censoring draws and whether
  each is observed or imputed; no test draws.
- `method_registry.csv`: method names, scopes, learners and censoring treatment.
- `split_manifest.csv`: IDs, split memberships, original and observed outcomes.
- `config.R`, `models.rds`, `sessionInfo.txt`: settings, fitted models and versions.

The output root contains `run_status.csv`, `metrics_all.csv` and `summary.csv`.
The summary reports completion counts alongside averages; methods with failures
are not silently treated as having the same comparison sample. A completed run
does not imply every separate-fitting or CQR comparison succeeded; inspect both
  status tables. Failed complete experiments and, by default, partial method
  comparisons return a nonzero exit code. The HPC merge writes `merged/summary.csv`,
  `summary_by_group.csv`, `metrics_all.csv`, `method_status_all.csv`, and `run_status.csv`.

## Interactive use

```r
source("real_data/support2.R")
result <- run_support2(
  repo_root = getwd(),
  scenario = "added_group",
  seed = 1,
  horizon = 180,
  output_dir = "results/support2/interactive_added_group"
)
result$metrics
result$diagnostics
```

For finer numerical checking, increase `integration_steps` to 512 or 1024.
Do not choose integration settings or audit groups by favorable test results.

## Validation

```bash
Rscript real_data/tests/test_support2.R
Rscript real_data/tests/test_support2_baselines.R
```

Checks cover analytic continuous AIPCW integration, zero bounds, compensation
before the first censoring jump, original outcome identification, split exhaustiveness
and reproducibility, preprocessing leakage, coherent added censoring, rejection of
unknown restricted outcomes, and explicit versus historical group labels.
Baseline checks compare DFT calibration with the historical formulas, validate
conditional censoring draws and random-stream isolation, and verify that
changing test outcomes cannot alter any of the eight baseline predictions.
They also reproduce a sparse-cell unpenalized failure, verify that the fixed
ridge fallback produces finite censoring predictions, and check that its use
is explicitly recorded.
Set `CAMS_ORIGINAL_CAMS_R` to an unmodified copy of `utils/cams.R` to additionally
compare existing simulation predictions and calibration diagnostics directly.

Tested on R 4.3.3, survival 3.5-8 and quantreg 5.97. Each run records its actual
package versions. The attached pilot is development evidence, not a final
statistical comparison.

## Sources

- https://hbiostat.org/data/
- https://hbiostat.org/data/repo/supportdesc
- https://hbiostat.org/data/repo/support2csv.zip
- https://archive.ics.uci.edu/dataset/880/support2

The original downloaded CSV used for the pilot has MD5
`9db505b0824c500f6feb746089da7259`.
For publication, cite the original SUPPORT study and acknowledge:
"Data obtained from http://hbiostat.org/data courtesy of the Vanderbilt University
Department of Biostatistics."
