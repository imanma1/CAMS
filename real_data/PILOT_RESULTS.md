# SUPPORT2 pilot results

This is a development pilot of the adapted pipeline, not a final statistical comparison.

Nine main experiments were run: seeds 1, 2 and 3 for each of native, common added censoring, and group-dependent added censoring. All pooled CAMS and marginal procedures completed. Each split used 4,551 fitting, 2,275 calibration and 2,279 test patients. The candidate grid was 0.001 through 0.999; alpha was 0.10; bounds were capped at 180 days; and the AIPCW integration used 256 intervals.

Exact test coverage is available here because all original censoring occurs after 180 days. The added scenarios retain original restricted outcomes only for evaluation.

## Selected results

Values are averages across the three splits. The worst-group column is the average of the minimum of the four realized subgroup coverages in each split, not the minimum of three-split average subgroup coverages. It is affected by test-sample variation and group selection.

| Scenario | Method | Exact overall coverage | Mean of split minimum group coverage | Mean bound (days) |
|---|---|---:|---:|---:|
| native | Weibull-uncalibrated | 0.9291 | 0.9121 | 12.19 |
| native | Marginal-IPCW | 0.9005 | 0.8825 | 17.44 |
| native | Marginal-AIPCW | 0.9005 | 0.8825 | 17.44 |
| native | CAMS-IPCW | 0.8995 | 0.8831 | 16.57 |
| native | CAMS-AIPCW | 0.8995 | 0.8831 | 16.57 |
| native | CQR-observed | 0.8983 | 0.8806 | 14.51 |
| added_common | Weibull-uncalibrated | 0.9160 | 0.9049 | 15.22 |
| added_common | Marginal-IPCW | 0.8979 | 0.8852 | 18.14 |
| added_common | Marginal-AIPCW | 0.8983 | 0.8852 | 18.06 |
| added_common | CAMS-IPCW | 0.8991 | 0.8807 | 17.72 |
| added_common | CAMS-AIPCW | 0.9004 | 0.8807 | 17.54 |
| added_common | CQR-observed | 0.9096 | 0.8927 | 11.82 |
| added_group | Weibull-uncalibrated | 0.9130 | 0.9005 | 15.74 |
| added_group | Marginal-IPCW | 0.8979 | 0.8857 | 18.32 |
| added_group | Marginal-AIPCW | 0.8981 | 0.8862 | 18.22 |
| added_group | CAMS-IPCW | 0.9008 | 0.8846 | 17.64 |
| added_group | CAMS-AIPCW | 0.9011 | 0.8846 | 17.57 |
| added_group | CQR-observed | 0.9191 | 0.9006 | 10.21 |


## What the pilot establishes

- The original dataset loads with the expected 9,105 patients and 6,201 observed deaths.

- The core pipeline runs using ordinary observed survival records, without latent event or censoring times supplied to CAMS.

- Pooled fitting is shared between marginal and CAMS calibration comparisons.

- Exact evaluation, reusable splits, model status and calibration diagnostics are saved.

- The existing simulation defaults match the unmodified CAMS implementation on the regression fixture.

- Continuous AIPCW quadrature matches an analytic exponential example, including below the first observed censoring time.

- On the saved added-group seed-1 predictions, increasing integration intervals from 256 to 512 changes average estimated test risk by at most 3.35e-7 across the eight methods. This is a numerical check on those predictions, not a universal error bound.

- Additional coarse-grid smoke runs completed for a 730-day native ceiling and for the day-3 landmark with physiology.


## Findings to carry into development

Overall CAMS coverage is close to 0.90 averaged over these three splits. The realized minimum subgroup coverage is lower, and the pilot does not establish reliable control in every subgroup or a clear CAMS advantage. Several individual subgroup results fall below 0.90. More splits, uncertainty reporting, model sensitivity and supported audit choices are needed before drawing conclusions.

CQR on observed time becomes more conservative in the added-censoring scenarios, with shorter bounds. These data do not reproduce all features of the high-dimensional synthetic examples.

One native seed-3 separate-fitting group, older female patients, produced nonfinite subgroup model coefficients. Both separate methods retain missing predictions for that group. Their native whole-cohort summaries use only two complete splits and explicitly report that count. Successful groups keep their own coverage estimates. Do not compare such averages to three-split pooled averages as if the samples were identical.

Quantile regression emitted nonunique-solution warnings in some runs, but returned finite predictions. Warnings and complete-case counts should remain visible in later research reporting.

Repeated splits share patients. No population confidence intervals or independent-cohort claims are made from this small pilot.


## Files

The complete raw outputs are under `results/support2/pilot/`. Begin with `summary.csv`, then inspect `metrics_all.csv`, per-run calibration diagnostics, and `method_status.csv`. `integration_check.csv` records the numerical refinement check.

Run new experiments with:

```bash
Rscript real_data/run_support2.R --scenario=all --seeds=4,5,6,7,8 --output=results/support2/development
```

Use the settings and caveats in `real_data/README.md`. The four audit groups are disjoint; this adapter does not yet implement general overlapping-group optimization.
