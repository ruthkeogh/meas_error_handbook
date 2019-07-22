# meas_error_handbook

Provides R code for performing analyses to correct for the impact of measurement error. The analyses are applied to example data from the
Third National Health and Nutrition Examination Survey (NHANES III), which is provided in a csv file.

The analysis of interest is a Cox regression with outcome of death due to cardiovascular disease and explanatory variables age, sex, diabetes status, smoking status, and systolic blood pressure (SBP). SBP is subject to measurement error and a repeated measurement is available in a subset of the cohort. Smoking status is subject to missing data,.  

The analyses performed are:
- a naive analysis ignoring measurement error in SBP and missing data in smoking status
- regression calibration for measurement error correction, ignoring missing data in smoking status
- Bayesian analysis, incorporating correction for measurement error and missing data
- Multiple imputation analysis using SMCFCS, incorporating correction for measurement error and missing data

# Reference

Ruth H Keogh & Jonathan W Bartlett. Measurement error as a missing data problem. In: Handbook of Measurement Error and Variable Selection. 2019. To appear. 

