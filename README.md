# Shrimp SD-DRM Project

This repository contains the R code and data used to estimate the probability
of illness from pathogenic *Vibrio parahaemolyticus* in raw shrimp and to
evaluate how residual doxycycline exposure and antimicrobial resistance may
affect illness risk and the likely treatment outcome.

The analysis combines:

- a probabilistic exposure assessment for contamination-positive raw shrimp;
- estimation of doxycycline pharmacodynamic parameters from time-kill data;
- a simple-death dose-response model (SD-DRM);
- scenario analyses for residual doxycycline concentration and resistant
  fraction; and
- PAWN global sensitivity analysis.

## Repository structure

```text
Updated-shrimp-SD-DRM/
├── Identification of distribution of M/
│   ├── Distribution_fit.R
│   └── shrimp_consumption_days.csv
├── Final emax model-doxycycline/
│   ├── Antibiotic data with control.csv
│   ├── Emax code.R
│   ├── Fitting time kill data into PD model.pdf
│   └── Paper with time kill data.pdf
├── Exposure assessment code/
│   ├── Exposure assessment.R
│   └── raw_shrimp_mc2d_simulation_data.csv
└── SD-DRM and sensitivity analysis/
    ├── model_functions.R
    ├── raw_shrimp_mc2d_simulation_data.csv
    ├── Shrimp_SD-DRM.Rmd
    └── Sensitivity analysis.Rmd
```

## Model overview

### Exposure assessment

The pathogenic dose per shrimp-consumption day is calculated as:

$$
D = C_{retail} \times M \times F_{pathogenic}
$$

where:

- $C_{retail}$ is the concentration of *V. parahaemolyticus* in
  contamination-positive retail shrimp;
- $M$ is the amount of shrimp consumed per day; and
- $F_{pathogenic}$ is the fraction of the total concentration assumed to be
  pathogenic.

The exposure inputs are modelled as variability among consumption days:

| Input | Distribution | Parameters |
|---|---|---|
| Retail concentration, $C_{retail}$ | Log-uniform | 75 to 1,100 MPN/g |
| Consumption amount, $M$ | Gamma | Shape = 0.8760928; scale = 43.276368 g |
| Pathogenic fraction, $F_{pathogenic}$ | Uniform | 0.0018 to 0.0233 |

The Gamma distribution for consumption amount was fitted to weighted person-day
shrimp-consumption data. The exposure analysis uses 1,000 variability
simulations and seed 123.

### Pharmacodynamic model

Doxycycline pharmacodynamic parameters are estimated using an inhibitory
$E_{max}$ model fitted to log10-transformed bacterial concentrations from
time-kill data:

$$
\log_{10}(N) = a + \left(m_0 - \frac{E_{max}C}{EC_{50}+C}\right)t
$$

The fitted values used in the SD-DRM are:

| Parameter | Estimate | Standard error |
|---|---:|---:|
| $E_{max}$ | 1.476380 log10/hour | 0.10336972 |
| $EC_{50}$ | 0.166877 mg/L | 0.04584022 |

In the main SD-DRM analysis, $E_{max}$ and $EC_{50}$ are sampled jointly from
their fitted bivariate normal distribution using their variance-covariance
matrix. This preserves their estimated relationship. $E_{max}$ is converted
from log10/hour to natural-log/day before use in the SD-DRM.

### SD-DRM

The baseline dose-response relationship uses a beta-Poisson model with
$\alpha=0.60$ and $\beta=1.31\times10^6$. The total pathogenic dose is divided
into susceptible and resistant components according to the resistant fraction
$f_r$.

The analysis evaluates:

- residual doxycycline concentrations from 0% to 10% of the MIC;
- resistant fractions from 0 to 0.50;
- a doxycycline MIC of 0.13 micrograms/mL; and
- a food-safety time interval of one day.

The model calculates total illness risk and $\Delta P_{treat}$, defined as the
difference between the probabilities of more-likely-treatable and
less-likely-treatable illness outcomes. Positive values indicate that a
more-likely-treatable outcome is more probable.

### PAWN sensitivity analysis

PAWN sensitivity analysis evaluates the influence of the following inputs on
illness risk:

- retail concentration;
- consumption amount;
- pathogenic fraction;
- resistant fraction;
- residual doxycycline concentration;
- $E_{max}$; and
- $EC_{50}$.

The PAWN index is based on the Kolmogorov-Smirnov distance between the
unconditional output distribution and conditional output distributions.
Larger values indicate greater influence on predicted illness risk. $E_{max}$
and $EC_{50}$ are sampled independently in this analysis so that their
individual effects can be evaluated.

## Software requirements

The analyses were developed in R. Install the required packages with:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "patchwork",
  "mc2d",
  "MASS",
  "fitdistrplus",
  "scales",
  "knitr",
  "rmarkdown"
))
```

## Running the analyses

The scripts use relative file paths. Run each script from the folder in which
it is stored.

### 1. Fit the consumption distribution

Open `Identification of distribution of M/Distribution_fit.R`, set its folder
as the working directory, and run the script. It fits weighted Gamma and
lognormal distributions and reports their parameters and weighted AIC values.

### 2. Estimate the pharmacodynamic parameters

Open `Final emax model-doxycycline/Emax code.R`, set its folder as the working
directory, and run the script. It fits the inhibitory $E_{max}$ model and
reports the $E_{max}$ and $EC_{50}$ estimates, standard errors, and
variance-covariance matrix.

### 3. Run the exposure assessment

Open `Exposure assessment code/Exposure assessment.R`, set its folder as the
working directory, and run the script. It creates:

- `raw_shrimp_mc2d_simulation_data.csv`; and
- `dose_variability_mc2d.png`.

To rerun the SD-DRM using new exposure simulations, copy the updated
`raw_shrimp_mc2d_simulation_data.csv` into the
`SD-DRM and sensitivity analysis` folder.

### 4. Run the SD-DRM analysis

Open `SD-DRM and sensitivity analysis/Shrimp_SD-DRM.Rmd`, set its folder as the
working directory, and knit the document. The file imports the exposure
simulations and sources `model_functions.R` from the same folder.

The analysis produces conditional dose-response curves, expected-risk plots,
illness-risk cumulative distributions, and $\Delta P_{treat}$ plots. Generated
figures are saved in an `images` subfolder.

### 5. Run the PAWN sensitivity analysis

Open `SD-DRM and sensitivity analysis/Sensitivity analysis.Rmd`, set its folder
as the working directory, and knit the document. This produces the PAWN
sensitivity indices, parameter ranking, and sensitivity figure.

## Reproducibility notes

- The main simulation uses seed 123.
- The PAWN analysis uses seed 0.
- The main SD-DRM preserves the covariance between $E_{max}$ and $EC_{50}$.
- The PAWN analysis samples $E_{max}$ and $EC_{50}$ independently to evaluate
  their individual influence.
- Exposure variables represent variability only; pharmacodynamic parameter
  draws represent uncertainty.
- Generated figures are saved at 600 dpi.

## Author

Sumit Jyoti  
Atlantic Veterinary College, University of Prince Edward Island, Canada

## Citation

If you use this repository, please cite the associated publication when it
becomes available.
