# JIA-polygenic-risk-and-cardiovascular-phenotypes
This R script relates to the manuscript "Juvenile idiopathic arthritis polygenic risk scores are associated with cardiovascular phenotypes in early adulthood". It details the ALSPAC variables used in this analysis, details of data cleaning and variable recoding, but excludes the variable extraction.

## Analysis pipeline
The analysis requires four main stages:
1. Extraction of ALSPAC data
2. Cleaning of ALSPAC data
3. Data exploration and summarisation
4. Regression analysis

### Extraction of ALSPAC data
A dataframe of ALSPAC offspring is generated from the ALSAPC cohort profile file. Analysis variables were also extracted from core ALSPAC files. Variables are coded as covariates ("covs"), cardiovascular ("cardio"), immunological ("immuno"), anthropometric ("anthro") or negative control ("negcon") variables. Polygenic risk scores and ALSPAC genetic principal components (PCs) were imported which were generated PRSice2.

### Cleaning of ALSPAC data
Data is coded according to participant age -  7yrs ("stage3"), 9yrs ("stage4"), 11yrs ("stage5), 13yrs ("stage6"), 15yrs ("stage7"), 17yrs ("stage8") and 24yrs ("stage9").
Extreme outlier values (greater than three times the interquartile range) were removed

### Data exploration and summarisation
Descriptive statistics for the cohort were generated and histograms for continuous variables were inspected

### Regression analysis
Linear and logistic regression analysis with and without covariates was undertaken and comprises:
1. Association between JIA PRS and positive/negative control variables using an unadjusted model, or models adjusted for genetic PCs, genetic PCs/age and genetic PC/age/sex.
2. Association between JIA PRS and cardiovascular phenotypes at age 24yrs
3. Sensitivity analyses:
    1. Association between JIA PRS and cardiovascular phenotypes at age 24yrs with MHC removal
    2. Association between JIA PRS and cardiovascular phenotypes at age 24yrs with participants with autoimmune disorders excluded
    3. Association between JIA PRS and cardiovascular phenotypes at age 24yrs using PRS scores with different P value thresholds
4. Longitudinal analysis of the association between JIA PRS and cardiovascular variables from age 7yrs to age 24yr
