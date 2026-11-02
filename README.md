# Bayesian Wildlife Survival Modelling

## Overview

This project explores Bayesian approaches to estimating wildlife survival using capture–recapture data from European dippers (*Cinclus cinclus*).

The analyses investigate annual apparent survival, differences between males and females, and variation in survival between monitoring years. Detection probability is estimated alongside survival to account for individuals that remain alive but are not encountered.

All analyses are contained in a single R script.

## Data

The analyses use the `dipper` dataset included in the `marked` R package.

The dataset contains encounter histories for 294 individually marked European dippers monitored between 1981 and 1987, together with information on sex.

Encounter histories indicate whether each individual was observed during successive monitoring occasions.

## Analyses

### Bayesian survival estimation

A Cormack–Jolly–Seber (CJS) model is used to estimate annual apparent survival and detection probability.

Posterior distributions are examined to calculate mean survival probabilities and 95% credible intervals.

### Sex differences in survival

A Bayesian CJS model is fitted to investigate whether apparent survival differs between males and females.

Posterior distributions are used to estimate sex-specific survival probabilities and the uncertainty surrounding their differences.

### Hierarchical Bayesian survival modelling

A hierarchical model is fitted to investigate variation in apparent survival between monitoring years.

Year-specific survival effects are modelled using a shared distribution, allowing survival estimates to vary between years while accounting for uncertainty.

### Model diagnostics

Markov chain Monte Carlo (MCMC) trace plots and convergence diagnostics are used to assess the sampling procedure.

## Methods

The analyses use Bayesian inference, capture–recapture modelling, hierarchical models and Metropolis–Hastings MCMC sampling.

## R packages

- `marked`
- `ggplot2`

## Running the analysis

Install the required packages:

```r
install.packages(c("marked", "ggplot2"))
```

Run `bayesian_wildlife_survival.R` in RStudio.

The dataset is included in the `marked` package, so no separate download is required.

The script generates posterior summaries, survival probability plots and MCMC diagnostics. Results are saved in the `results` and `figures` folders.

## Limitations

The models estimate apparent survival rather than true survival, as permanent emigration cannot be distinguished from mortality.

Detection probability is assumed to remain constant across sexes and monitoring years.

The hierarchical model is exploratory, and annual survival differences are not directly attributed to environmental conditions.

The MCMC implementation is intended for teaching purposes. Convergence and mixing should be evaluated before interpreting model estimates.

## Data source

European dipper capture–recapture dataset, available through the `marked` R package.

Package documentation: https://cran.r-project.org/package=marked
