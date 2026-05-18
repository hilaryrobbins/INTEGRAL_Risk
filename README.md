# INTEGRAL_Risk
Validation of the INTEGRAL-Risk protein-based lung cancer risk prediction model

These R scripts reproduce tables and figures in Zahed H, Feng X, et al. Biomarker-based eligibility for lung cancer screening - Validation of the protein-based INTEGRAL-Risk model in the Lung Cancer Cohort Consortium (LC3). JAMA 2026

Before using these scripts, we remind you of the following points:
1) Study design: The study design for both the model training and testing was a stratified case-cohort design that requires the application of appropriate sampling weights.
2) Population and time horizon: The study was limited to individuals with a smoking history, and the analysis is valid over a 3-year time horizon, because lung cancer cases diagnosed up to 3 years after blood draw were included.
3) Biomarker measurement: Blood measurements of 21 circulating proteins were measured for all study participants using the INTEGRAL panel. 
4) Access to the individual level data from this study is governed by the Lung Cancer Cohort Consortium (LC3) Access Policy: https://lc3.iarc.who.int/documents/lc3-access-policy.pdf.

The following scripts are included:
INTEGRAL_function_validation.R - This script specifies all functions. It must be loaded prior to executing any of the other scripts.
INTEGRAL_RISK_AIC.rds - This object file contains the fitted INTEGRAL-Risk model object.
INTEGRAL_validationscript.R - This reproduces main tables and figures related to the INTEGRAL-Risk model in Zahed et al.

We do not present scripts for the other considered models (best AUC, best AIC in Asia, best AUC in Aisa) since they were not selected as the main INTEGRAL-Risk model. We do not present scripts for eTables 1-5, since they are largely descriptive.

Key variables:
LC_event: lung cancer event (0,1)
t0: enrollment time
tevent: follow-up time until lung cancer diagnosis, death, or end of follow-up (3 years) 
