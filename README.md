# Current manuscript version of DIERM: annotated code

This directory contains the production code for the Data-Informed Ecosystem Respiration Model (DIERM) used in the current manuscript. The added documentation does not change the model equations, filtering criteria, random seeds, or output definitions. It explains the statistical meaning and recommended execution order of each step.

## 1. Model definition

A centered unimodal temperature-response curve is first fitted for each site-year:

```text
log(ER) = beta0 + beta1 * (Ta - 12) + beta2 * (Ta - 12)^2
beta2 < 0
```

- `Ta`: daily mean air temperature (degrees C).
- `ER`: daily ecosystem respiration.
- `Tref = 12 degrees C`: a centering constant used to reduce parameter compensation and numerical correlation. It is not an ecological threshold.
- `beta0`: respiration level near 12 degrees C.
- `beta1`: first-order temperature response at 12 degrees C.
- `beta2`: curvature, constrained to be negative to impose a unimodal response.
- When the parabolic vertex can be calculated, `Topt = 12 - beta1 / (2 * beta2)`.
- `HOT days` is the annual number of days on which `Ta > Topt`.
- Daily ER is reconstructed as `exp(beta0 + beta1*x + beta2*x^2)`, and annual ER is the sum of the 365 daily estimates.

Before exponentiation, the linear predictor is clipped to `[-50, 50]` to prevent floating-point overflow. This is only a numerical safeguard and is not an additional ecological filter.

## 2. Site-year fitting rules

`scripts/refit_predict_beta_simulate_er.py` performs the following operations for every site-year:

1. Retain paired daily observations for which `-55 <= Ta <= 60 degrees C`, `ER > 0`, and both variables are finite.
2. Require at least 300 valid paired daily observations per site-year.
3. Aggregate observations into 1 degree C temperature bins, require at least three days per retained bin, and require at least six retained bins per site-year.
4. Fit `beta0`, `beta1`, and `beta2` with bounded nonlinear least squares using multiple initial values while constraining `beta2 < 0`.
5. Do not apply the older `R2 > 0.40` or `Topt <= 1.2*Tmax` filters. Responses with a fitted vertex above the observed temperature range, as well as nearly linear responses, remain in the ER parameter-model training set. `optimum_class` records vertex identifiability but is not used as an exclusion criterion.

## 3. Parameter model and independent validation

`scripts/train_global_beta_model.py` predicts the three beta coefficients from site-year air temperature, vapor-pressure deficit, precipitation, GPP, coordinates, and monthly statistics. The primary model is a multi-output Extra Trees regressor with:

- 600 trees;
- `min_samples_leaf = 4`;
- `max_features = 0.8`;
- median imputation of missing continuous predictors;
- joint prediction of `beta0`, `beta1`, and `beta2`; and
- clipping of predicted `beta2` to a maximum of `-1e-10` to preserve the unimodal structure.

Five-fold cross-validation is grouped by site using `GroupKFold`. All years from a given site occur on only one side of a training-test split. The out-of-fold results therefore constitute cross-site validation rather than a random split of site-year records. After out-of-fold evaluation, the final model is refitted using all eligible site-years and saved as `global_beta_model.joblib`.

## 4. Historical and future simulations

- `project_global_beta_crujra_cedar.py` uses daily or 6-hourly CRU-JRA meteorology and CEDAR-GPP to reconstruct global `beta`, `Topt`, HOT days, and ER from 1991 to 2020.
- `project_global_beta_ssp.py` uses CESM2-WACCM daily air temperature and humidity together with monthly GPP and precipitation for SSP1-2.6, SSP2-4.5, SSP3-7.0, and SSP5-8.5.
- Global summaries use grid-cell area weights. Annual ER is converted from `micromol CO2 m-2 s-1` to `g C m-2 yr-1`.
- HDF5 files store annual spatial fields, and CSV files store area-weighted annual summaries.

The historical and future results are conditional simulations under the fitted unimodal response structure. They should not be described as an independent validation of the HOT-day hypothesis.

## 5. File guide

| Script | Purpose |
|---|---|
| `dierm_core_annotated.py` | Annotated core equations and functions for Topt, HOT days, daily ER, and annual ER |
| `evaluate_unfiltered_parameter_prediction.py` | Constructs site-year environmental predictors and preprocessing pipelines |
| `refit_predict_beta_simulate_er.py` | Fits site-year beta curves and generates site-grouped out-of-fold results |
| `train_global_beta_model.py` | Trains and saves the final Extra Trees parameter model used for global projections |
| `preaggregate_cedar_gpp.py` | Aggregates CEDAR-GPP from 0.05-degree to the 0.5-degree CRU grid |
| `project_global_beta_crujra_cedar.py` | Performs the historical global reconstruction |
| `project_global_beta_ssp.py` | Performs future SSP projections |
| `summarize_global_beta_ssp.py` | Summarizes annual results across SSPs |
| `make_fig2_cd_global_beta_crujra_cedar.py` | Generates DIERM panels C and D of the current Figure 2 |
| `make_future_beta_figure.py` | Generates the future-scenario summary figure |

## 6. Recommended execution order

The following commands are path templates. Replace the values in angle brackets with the corresponding local paths.

```powershell
# 1) Fit local site-year responses and perform site-grouped OOF evaluation
python scripts/refit_predict_beta_simulate_er.py `
  --daily <fluxnet_daily.csv.gz> `
  --output-dir <beta_fit_output> `
  --tref 12

# 2) Train the final parameter model used for spatial projections
python scripts/train_global_beta_model.py `
  --daily <fluxnet_daily.csv.gz> `
  --fits <beta_fit_output/all_site_year_beta_fits.csv> `
  --output-dir <trained_model_output>

# 3) Optional: pre-aggregate CEDAR-GPP
python scripts/preaggregate_cedar_gpp.py `
  --zip <CEDAR_GPP.zip> `
  --output <cedar_gpp_0p5deg.h5> `
  --temp-dir <temporary_directory> `
  --start-year 1991 --end-year 2020

# 4) Historical reconstruction
python scripts/project_global_beta_crujra_cedar.py `
  --model <trained_model_output/global_beta_model.joblib> `
  --cru-dir <CRU_JRA_directory> `
  --cedar-cache <cedar_gpp_0p5deg.h5> `
  --output-h5 <global_beta_crujra_cedar_1991_2020.h5> `
  --output-csv <global_beta_crujra_cedar_1991_2020.csv> `
  --start-year 1991 --end-year 2020

# 5) Run each SSP separately; use ssp126, ssp245, ssp370, or ssp585
python scripts/project_global_beta_ssp.py `
  --scenario ssp585 `
  --model <trained_model_output/global_beta_model.joblib> `
  --daily-dir <CMIP6_daily_directory> `
  --monthly-dir <CMIP6_monthly_directory> `
  --large-output-dir <future_h5_output> `
  --summary-dir <future_csv_output> `
  --start-year 2015 --end-year 2100
```

## 7. Important distinction

DIERM derives `Topt` from the vertex of its predicted temperature-response curve. This differs from the threshold used in the time-holdout analysis of the 94 flux sites, where validation-year `Topt` is predicted from a significant training-period MAT relationship or otherwise represented by the mean training-period `Topt`. The two analyses serve different purposes and should not be conflated in the Methods.
