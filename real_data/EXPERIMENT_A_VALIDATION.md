# Experiment A: validation of the expanded SUPPORT2 runner

All eight requested DFT/CQR comparisons are implemented in addition to the existing CAMS, marginal and separate-fitting methods. This report concerns original censoring only; no added event or censoring times replace native observations.

## Executed protocol

- Original public SUPPORT2: 9,105 patients; death event and d.time follow-up.
- Entry features; four sex by age-65 cells; stratified 50/25/25 splits.
- Seeds 1, 2 and 3: 4,551 fitting, 2,275 calibration, 2,279 test patients each.
- Alpha 0.10; bound ceiling 730 days; candidate grid 0.001 through 0.999.
- Lognormal censoring AFT; Weibull event AFT; probability floor 0.05.
- Dense AIPCW quadrature with 256 intervals.
- DFT conditional imputation of missing C; lognormal AFT CT cutoff learner.
- Training-only ridge fallback theta 0.1 for unusable unpenalized AFT fits.

Each seed was executed locally through the actual real_data/submit_support2.sh bash file with its module-loading option disabled and a configured R executable. All three returned exit code zero. All 15 methods completed: 75 metric rows and 33 method/group status rows per seed. Native Y/event fields remained unchanged; no test censoring times were imputed.

## Descriptive pilot results

Coverage is exact for these saved predictions because every original coverage indicator is resolved. This does not mean all survival times are observed or all 730-day restricted outcomes are known. At other splits/ceilings, the code reports unresolved exact coverage as NA.

The subgroup column averages the minimum realized subgroup coverage within each split. It is not a population coverage guarantee or an independent-cohort confidence interval.

| Method | Complete splits | Exact overall coverage | Mean split minimum subgroup coverage | Mean bound (days) |
|---|---:|---:|---:|---:|
| Weibull-uncalibrated | 3/3 | 0.9291 | 0.9121 | 12.19 |
| Marginal-IPCW | 3/3 | 0.9005 | 0.8825 | 17.45 |
| Marginal-AIPCW | 3/3 | 0.9005 | 0.8825 | 17.45 |
| CAMS-IPCW | 3/3 | 0.8995 | 0.8831 | 16.58 |
| CAMS-AIPCW | 3/3 | 0.8995 | 0.8831 | 16.58 |
| Separate-IPCW | 3/3 | 0.8983 | 0.8824 | 18.23 |
| Separate-AIPCW | 3/3 | 0.8983 | 0.8824 | 18.23 |
| Adaptive DFT-T (joint) | 3/3 | 0.9005 | 0.8825 | 17.45 |
| Adaptive DFT-CT (joint) | 3/3 | 0.9005 | 0.8825 | 17.45 |
| DFT-fixed (joint) | 3/3 | 0.9002 | 0.8779 | 14.75 |
| CQR (joint) | 3/3 | 0.8983 | 0.8806 | 14.51 |
| Adaptive DFT-T (subgroup) | 3/3 | 0.8983 | 0.8824 | 18.20 |
| Adaptive DFT-CT (subgroup) | 3/3 | 0.8983 | 0.8824 | 18.20 |
| DFT-fixed (subgroup) | 3/3 | 0.9010 | 0.8852 | 15.98 |
| CQR (subgroup) | 3/3 | 0.8988 | 0.8724 | 18.25 |

Overall coverages are near 0.90, but realized subgroup minima are lower. These three splits do not establish reliable subgroup control or a clear CAMS advantage. Native SUPPORT2 has no observed censoring before day 344, so typical short lower bounds do not strongly stress early-censoring correction. The joint adaptive T/CT predictions coincide in this pilot; this is expected when the CT cap and censoring selection do not change the selected candidates.

## Numerical failure found and resolved

Before adding the explicit ridge policy, seed 3 failed the older-female censoring AFT fit: a one-patient missing-race indicator had no observed censoring events and an undefined coefficient. The final helper retries an unusable unpenalized fit with a fixed ridge penalty on the same fitting data. It records the reason and penalty in fitted-model metadata and method status. Five seed-3 method/group status entries report this one shared subgroup fitting issue: the two separate calibrators and the three subgroup DFT methods. CQR is fitted independently. No test outcomes were used to choose the fallback or its penalty.

## Checks completed

- 19 existing identification, preprocessing and CAMS regression checks passed, including agreement with the original unmodified cams.R simulation behavior.
- 19 baseline checks passed: historical DFT T/CT risk formulas, weighted fixed-DFT calibration including test-point mass, conservative empty-selection fallback, conditional censoring draws, RNG isolation, untouched native records, exclusion of test IDs from imputation, no test-outcome influence on predictions, and recorded stabilization of the sparse-cell fit.
- Setup bash script completed using cached data and installed survival/quantreg packages.
- Both bash scripts passed syntax checks.
- A simulated Slurm-copied script resolved SLURM_SUBMIT_DIR, array seed, paths containing spaces, all method flags and one-thread settings correctly.
- Merger checks passed for complete, missing, interrupted, duplicate-seed and incompatible-code/settings cases.
- The final merger completed all three verified native tasks, preserving 15 method summaries and 75 group summaries.
- Saved outputs were checked for complete method status, coherent finite bounded predictions, preserved native records, observed C retained exactly, conditional C above Y, and absence of test imputation.

Validation used R 4.3.3, survival 3.5-8 and quantreg 5.97. The Slurm resource/module directives follow the supplied submit_sim.sh; no job was submitted to the actual HPC scheduler during local validation. Each run saves sessionInfo.txt for cluster-side reproducibility.

## Run on HPC

From the extracted repository root:

```bash
bash real_data/setup_support2.sh
sbatch real_data/submit_support2.sh
# After completion, replace JOBID with the returned array job ID:
module load StdEnv/2023 r/4.5.0
Rscript real_data/merge_support2.R --input=results/support2/experiment_A_JOBID --expected-seeds=1:100
```

The default array is 100 seeds with up to 10 concurrent tasks. For a cluster pilot use sbatch --array=1-3%3 real_data/submit_support2.sh and merge with --expected-seeds=1:3.

Begin with real_data/README.md for the protocol, adaptation details and overrides. Saved verification results are under results/support2/experiment_A_verified/merged/; detailed logs are under its validation/ directory. DFT imputation adaptations and AFT CT learner choice must remain explicit in any research report.
