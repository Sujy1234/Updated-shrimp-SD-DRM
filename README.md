# SD-DRM shrimp project

This repository contains the code and data used to develop and apply a simple-death dose-response model (SD-DRM) for antibiotic-resistant *Vibrio parahaemolyticus* associated with raw shrimp consumption.

## Repository structure

- `Identification of distribution of M`  
  Contains the shrimp-consumption data and R code used to fit and compare Gamma and lognormal distributions. The selected Gamma distribution is used to represent variability in the amount of shrimp consumed per day in the exposure assessment.

- `Final emax model-doxycycline`  
  Contains the doxycycline time-kill data, supporting publications, and R code used to estimate the pharmacodynamic parameters E<sub>max</sub> and EC<sub>50</sub>. The parameter estimates and their variance-covariance matrix are used to represent pharmacodynamic uncertainty in the SD-DRM.

- `Exposure assessment code`  
  Contains the R code used to estimate the dose of pathogenic *Vibrio parahaemolyticus* per raw-shrimp-consumption day. The simulated dose data generated in this folder are used as the exposure input in the SD-DRM analysis.

- `SD-DRM and sensitivity analysis`  
  Contains the main shrimp SD-DRM analysis, supporting model functions, exposure-input data, and PAWN sensitivity analysis. The model evaluates how residual doxycycline concentration and the fraction of resistant bacteria influence illness risk and the predicted treatability of illness. The sensitivity analysis evaluates the relative influence of the exposure and pharmacodynamic inputs on predicted illness risk.

The Gamma distribution fitted in the `Identification of distribution of M` folder is used in the exposure assessment. The pathogenic doses generated in the `Exposure assessment code` folder and the doxycycline pharmacodynamic estimates obtained in the `Final emax model-doxycycline` folder are then used as inputs in the main SD-DRM analysis.
