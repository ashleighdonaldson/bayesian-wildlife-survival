# Bayesian Wildlife Survival Modelling in R

## Overview

This project investigates Bayesian approaches to estimating wildlife survival using capture–recapture data from European dippers (*Cinclus cinclus*).

The analyses examine annual apparent survival, differences between males and females, and detection probability. Maximum-likelihood and Bayesian Cormack–Jolly–Seber (CJS) models are fitted to account for imperfect detection and evaluate uncertainty in survival estimates.

All analyses are contained in a single R script.

## Data

The analyses use the `dipper` dataset included in the `marked` R package.

The dataset contains encounter histories for 294 individually marked European dippers monitored between 1981 and 1987, together with information on sex.

Encounter histories record whether each bird was detected during successive monitoring occasions. Individuals that were not encountered may have survived but remained undetected.

## Analyses

### Maximum-likelihood survival modelling

Cormack–Jolly–Seber models are fitted to estimate apparent survival and detection probability.

Five candidate models are considered:

- Constant survival and detection probability
- Sex-specific survival with constant detection
- Time-varying survival with constant detection
- Survival varying with sex and time
- Sex-specific survival with time-varying detection

Models are compared using Akaike's Information Criterion (AIC), with differences in AIC and model weights used to evaluate relative support.

### Bayesian survival modelling

A Bayesian CJS model is fitted using the `probitCJS` implementation in the `marked` package.

The model estimates sex-specific apparent survival and a common detection probability.

Posterior distributions are used to calculate mean survival probabilities, 95% credible intervals and the probability of a difference in survival between males and females.

### Model diagnostics

Four independent MCMC runs are used to evaluate sampling performance.

Convergence is assessed using Gelman–Rubin diagnostics, effective sample sizes and trace plots.

Maximum-likelihood estimates provide an additional reference for evaluating the Bayesian survival estimates.

## Methods

The analyses use capture–recapture modelling, maximum-likelihood estimation, Bayesian inference, Markov chain Monte Carlo sampling and information-theoretic model comparison.

The CJS framework estimates apparent survival while accounting for imperfect detection.

## R packages

- `marked`
- `ggplot2`
- `coda`

## Running the analysis

Install the required packages:

```r
install.packages(c("marked", "ggplot2", "coda"))
```

Run `bayesian_wildlife_survival.R` in RStudio.

The dataset is included in the `marked` package, so no additional data download is required.

The script creates `results` and `figures` directories containing model summaries, survival estimates, posterior distributions and diagnostic information.

## Outputs

The analysis produces:

- Annual encounter summaries
- AIC model comparison table
- Maximum-likelihood survival and detection estimates
- Bayesian sex-specific survival estimates and credible intervals
- Posterior survival distributions
- MCMC convergence diagnostics

## Limitations

The models estimate apparent survival rather than true survival, as mortality cannot be distinguished from permanent emigration.

The Bayesian model assumes constant detection probability across sexes and years, whereas the maximum-likelihood candidate set considers several alternative structures.

The Bayesian and maximum-likelihood models use different link functions, so regression coefficients are not directly comparable. Survival probabilities can be compared on the response scale.

The analysis uses established CJS modelling methods, but satisfactory convergence, model fit and parameter identifiability must be confirmed before drawing ecological conclusions. A formal goodness-of-fit assessment and posterior predictive checks are not included in the current script.

## Data source

European dipper capture–recapture dataset, available through the `marked` R package.

Package documentation: https://cran.r-project.org/package=marked

## Reproducibility

The analysis uses a fixed random seed and saves fitted models, numerical results and R session information to support reproducibility.
