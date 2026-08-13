# Overview 
This repository contains the code to reproduce the numerical results 
in Section 4 and Section 5 of the paper: 
[Conformalized Survival Analysis with Adaptive Cutoffs](https://arxiv.org/abs/2211.01227).

# Folders
- `simulations/` contains the code for the simulation results in Section 4. 
- `real_data/` contains the code and data for the real data example in Section 5.
- `utils/` contains the helper functions.
- `results/` stores the results.

# Usage 
## Simulations
The script `run_sim.R` implements one run of the 
experiment. The users can specify the random seed when running the script. For example, to implement one run of the 
experiment in Section 4 with random seed 1, run the following command in your terminal:
```{r}
cd simulations 
Rscript run_sim.R homo_cens 1 oracle 1 correct 0 1000 1
```

### Censoring-model and augmentation robustness grid

`simulations/run_sim.R` also supports the controlled AIPCW experiment:

```text
Rscript run_sim.R <settings> <first-seed> <sc-methods> <only-cams> \
  <augmentation-methods> <homoscedastic-event> <rsf-ntree> <number-of-runs>
```

The censoring methods are `oracle`, `rsf`, `aft_lognormal`, `km_x1`, and
`km`. The augmentation methods are `correct`, `wrong`, and `same`. For the
definitive 2-by-2 double-robustness check, run:

```bash
cd simulations
Rscript run_sim.R \
  cams_vs_vanilla_lower_tail_hd_mild,cams_vs_vanilla_lower_tail_hd_main,cams_vs_vanilla_lower_tail_hd_strong \
  1 oracle,km 1 correct,wrong 1 1000 10
```

Setting `homoscedastic-event` to `1` uses a constant Weibull AFT scale, so the
full augmentation model is correctly specified. Results and plots are stored
under separate `sc_<method>/aug_<method>/event_<scale>` directories. Each
summary CSV reports means and standard deviations across seeds.

## Real data 
The script `pings_day_14.R` implements the real data example in Section 5.
To reproduce the results, run the following code in your terminal:
```{r}
cd real_data
Rscript pings_day_14.R
```
