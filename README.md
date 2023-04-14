# MS-phenotype-prediction

The code in this repository was used to generate the results presented in the article *Disease phenotype prediction in multiple sclerosis*.

## Background
Progressive multiple sclerosis (PMS) is currently diagnosed retrospectively. We investigated an alternative approach to evaluate the transition to PMS, using metabolomics and conformal prediction.

## Components
The analysis workflow consists of three main components:

1. The **data processing**, where data is cleaned and normalised.
2. The **variable selection**, where an elastic-net logistic regression model is used to select discriminatory metabolites.
3. The **conformal prediction**, where a support vector machine (SVM) model is built on the selected metabolites and supplemented with conformal prediction to enable estimation of prediction confidence in each new single patient prediction.
