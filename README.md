# Current DIERM workflow: MATLAB implementation

This package is the MATLAB conversion of the DIERM workflow used in the current manuscript. 

See `VALIDATION.md` for the completed checks, corrections, and remaining validation boundary.

## Model equation

For each site-year, the local response is

```text
log(ER) = beta0 + beta1 * (Ta - 12) + beta2 * (Ta - 12)^2
beta2 < 0
```

The centering temperature of 12 degrees C improves numerical stability; it is not an ecological threshold. The curve vertex is

```text
Topt = 12 - beta1 / (2 * beta2)
```

A HOT day is a day on which daily mean air temperature is greater than the same grid-cell/year Topt. Daily ER is reconstructed from the exponential response and integrated over 365 days.

## Site-year fitting rules

1. Retain paired daily observations with `-55 <= Ta <= 60 degrees C`, `ER > 0`, and finite values.
2. Require at least 300 valid daily pairs per site-year.
3. Aggregate data into 1-degree-C temperature bins, require at least three days per retained bin, and require at least six bins per site-year.
4. Fit the centered quadratic response by bounded nonlinear least squares from multiple starting values, with `beta2` constrained to `[-2, -1e-10]`.
5. Do not apply the legacy `R2 > 0.40` or `Topt <= 1.2*Tmax` filters. Vertex identifiability is recorded in `optimum_class` but is not an exclusion rule.

## Extra-Trees implementation

`fitExtraTreesRegressor.m` is a native MATLAB implementation of the key Extra-Trees rules used in the current analysis:

- 600 trees for the final model;
- no bootstrap resampling (each tree receives the complete training set);
- one uniformly random split threshold per candidate predictor at each node;
- 80% of predictors considered at each node;
- minimum leaf size of four; and
- a joint multi-output ensemble for the site-year OOF diagnostic, matching
  `ExtraTreesRegressor` in that step; and
- separate ensembles for `beta0`, `beta1`, and `beta2` in the final global
  model, matching the original `MultiOutputRegressor` configuration.

Predicted `beta2` is clipped to at most `-1e-10` before Topt and ER are calculated. Five-fold OOF evaluation is grouped by site, so all years from the same site stay in the same fold.

Because MATLAB and scikit-learn use different random-number generators and internal tree representations, a model retrained with this package will not be bit-for-bit identical to the previously serialized Python model. The equations, predictors, grouping rule, tree count, feature-sampling fraction, minimum leaf size, curvature constraint, and random seeds are preserved. Manuscript results should be regenerated as a complete set with one implementation rather than mixing maps or summaries from the two software implementations.

## Main files

| MATLAB file | Purpose |
|---|---|
| `refitPredictBetaSimulateER.m` | Fits site-year beta curves and performs site-grouped OOF evaluation |
| `trainGlobalBetaModel.m` | Trains and saves the final global beta model |
| `evaluateUnfilteredParameterPrediction.m` | Optional grouped-site parameter sensitivity diagnostics |
| `fitExtraTreesRegressor.m` | Native MATLAB Extra-Trees training |
| `predictExtraTreesRegressor.m` | Prediction from the native Extra-Trees model |
| `projectGlobalBetaCRUJRA.m` | Historical CRU-JRA/CEDAR-GPP reconstruction |
| `preaggregateCedarGPP.m` | Aggregates extracted CEDAR-GPP files from 0.05 to 0.5 degrees |
| `projectGlobalBetaSSP.m` | CESM2-WACCM SSP projections |
| `summarizeGlobalBetaSSP.m` | Baseline/future changes and Newey-West HAC trends |
| `makeFigure2CD.m` | Historical DIERM panels C and D of Figure 2 |
| `makeFutureBetaFigure.m` | Future HOT-day/ER figure |
| `selfTestDIERM.m` | Small synthetic-data test of equations, fitting, and Extra Trees |
| `selfTestWorkflow.m` | End-to-end synthetic test of fitting, grouped OOF prediction, ER reconstruction, and final model training |

The remaining `.m` files are documented helper functions for preprocessing, gridded I/O, area weighting, robust regression, and DIERM calculations.

## MATLAB requirements

- MATLAB R2022b or later is recommended.
- Statistics and Machine Learning Toolbox.
- Optimization Toolbox.
- NetCDF and HDF5 support included with MATLAB.
- Parallel Computing Toolbox is optional; the supplied code does not require it.

## Recommended execution order

Add the scripts directory to the MATLAB path:

```matlab
addpath('scripts');
```

### 1. Fit site-year curves and run grouped-site OOF evaluation

```matlab
results = refitPredictBetaSimulateER( ...
    "D:/data/fluxnet_daily.csv.gz", ...
    "D:/output/site_year_beta", ...
    Tref=12, NumTrees=500, Seed=42);
```

### 2. Train the final projection model

```matlab
bundle = trainGlobalBetaModel( ...
    "D:/data/fluxnet_daily.csv.gz", ...
    "D:/output/site_year_beta/all_site_year_beta_fits.csv", ...
    "D:/output/global_model", ...
    NumTrees=600, Seed=20260826);
```

The saved model is `global_beta_model.mat`.

### 3. Pre-aggregate CEDAR-GPP

Extract the CEDAR-GPP archive first, then run:

```matlab
preaggregateCedarGPP( ...
    "D:/data/CEDAR_GPP_extracted", ...
    "D:/output/cedar_gpp_0p5deg.mat", 1991, 2020);
```

### 4. Historical projection

```matlab
cfg = struct;
cfg.ModelFile = "D:/output/global_model/global_beta_model.mat";
cfg.CruDirectory = "D:/data/CRU_JRA";
cfg.CedarCache = "D:/output/cedar_gpp_0p5deg.mat";
cfg.OutputH5 = "D:/output/global_beta_crujra_cedar_1991_2020.h5";
cfg.OutputCSV = "D:/output/global_beta_crujra_cedar_1991_2020.csv";
cfg.StartYear = 1991;
cfg.EndYear = 2020;
historicalSummary = projectGlobalBetaCRUJRA(cfg);
```

### 5. SSP projection

```matlab
cfg = struct;
cfg.Scenario = "ssp585";
cfg.ModelFile = "D:/output/global_model/global_beta_model.mat";
cfg.DailyDirectory = "D:/data/CMIP6/daily";
cfg.MonthlyDirectory = "D:/data/CMIP6/monthly";
cfg.OutputH5 = "D:/output/global_beta_ssp585_2015_2100.h5";
cfg.OutputCSV = "D:/output/global_beta_ssp585_annual_summary.csv";
cfg.StartYear = 2015;
cfg.EndYear = 2100;
sspSummary = projectGlobalBetaSSP(cfg);
```

Run the function separately for `ssp126`, `ssp245`, `ssp370`, and `ssp585`.

### 6. Summaries and figures

```matlab
summarizeGlobalBetaSSP("D:/output/ssp_summaries", ...
    "D:/output/ssp_change_and_trends.csv");

makeFigure2CD( ...
    "D:/output/global_beta_crujra_cedar_1991_2020.csv", ...
    "D:/output/global_beta_crujra_cedar_1991_2020.h5", ...
    "D:/output/Figure2_CD");

makeFutureBetaFigure("D:/output/ssp_summaries", ...
    "D:/output/future_HOT_ER");
```

## Interpretation boundary

DIERM Topt is the vertex of a beta curve predicted from environmental covariates. It is different from the threshold used in the 94-site temporal holdout analysis, where a training-period MAT-Topt regression predicts validation-year Topt only when that training relationship is significant; otherwise, the training-period mean Topt is used. These two analyses serve different purposes and should not be combined in the Methods description.

Historical and future DIERM outputs are conditional simulations under the fitted unimodal response structure. They should not be presented as independent validation of the HOT-day hypothesis.
