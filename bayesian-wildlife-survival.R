
# Bayesian wildlife survival modelling
# European dipper capture-recapture data

# install.packages(c("marked", "ggplot2", "coda"))

library(marked)
library(ggplot2)
library(coda)

set.seed(123)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

data("dipper", package = "marked")

head(dipper)
str(dipper)
table(dipper$sex)


# Encounter histories

histories <- strsplit(as.character(dipper$ch), "")

stopifnot(
  all(lengths(histories) == 7),
  all(unlist(histories) %in% c("0", "1"))
)

encounters <- do.call(
  rbind,
  lapply(histories, as.integer)
)

encounter_summary <- data.frame(
  Year = 1981:1987,
  Encountered = colSums(encounters)
)

encounter_plot <- ggplot(
  encounter_summary,
  aes(Year, Encountered)
) +
  geom_line(
    colour = "steelblue",
    linewidth = 0.9
  ) +
  geom_point(
    colour = "steelblue",
    size = 2
  ) +
  theme_classic() +
  labs(
    title = "European dipper encounters",
    x = "Year",
    y = "Number of birds encountered"
  )

print(encounter_plot)


# Prepare CJS data

processed <- process.data(
  dipper,
  model = "CJS",
  groups = "sex",
  begin.time = 1981
)

design <- make.design.data(processed)


# Maximum-likelihood models

fit_cjs <- function(phi_formula, p_formula) {
  
  crm(
    processed,
    design,
    model.parameters = list(
      Phi = list(formula = phi_formula),
      p = list(formula = p_formula)
    ),
    hessian = TRUE
  )
  
}

model_constant <- fit_cjs(~1, ~1)
model_sex <- fit_cjs(~sex, ~1)
model_year <- fit_cjs(~time, ~1)
model_sex_year <- fit_cjs(~sex + time, ~1)
model_detection_year <- fit_cjs(~sex, ~time)

models <- list(
  Constant = model_constant,
  Sex = model_sex,
  Year = model_year,
  SexYear = model_sex_year,
  SexDetectionYear = model_detection_year
)


# Compare maximum-likelihood models

# The marked results object stores the information
# criterion in the AIC field.

aic_values <- vapply(
  models,
  function(x) {
    
    value <- x$results$AIC
    
    if (
      length(value) != 1 ||
      !is.numeric(value) ||
      !is.finite(value)
    ) {
      stop(
        "A valid AIC value was not found ",
        "for one of the fitted models."
      )
    }
    
    as.numeric(value)
  },
  numeric(1)
)

model_comparison <- data.frame(
  Model = names(models),
  AIC = unname(aic_values)
)

model_comparison <- model_comparison[
  order(model_comparison$AIC),
]

model_comparison$DeltaAIC <-
  model_comparison$AIC -
  min(model_comparison$AIC)

model_comparison$Weight <- exp(
  -0.5 * model_comparison$DeltaAIC
)

model_comparison$Weight <-
  model_comparison$Weight /
  sum(model_comparison$Weight)

rownames(model_comparison) <- NULL

print(model_comparison)


# Examine the best-supported model

best_model <- models[[
  model_comparison$Model[1]
]]

print(best_model)

print(best_model$results$beta)
print(best_model$results$reals)


# Maximum-likelihood survival estimates

# predict.crm is an S3 method called through
# the generic predict() function.

mle_survival <- predict(
  model_sex,
  ddl = design,
  parameter = "Phi",
  se = TRUE,
  merge = TRUE
)

print(mle_survival)


# Maximum-likelihood detection estimates

mle_detection <- predict(
  model_sex,
  ddl = design,
  parameter = "p",
  se = TRUE,
  merge = TRUE
)

print(mle_detection)


# Bayesian CJS model

# The probitCJS model uses a probit link and
# the established sampler supplied by marked.

bayes_parameters <- list(
  Phi = list(formula = ~sex),
  p = list(formula = ~1)
)

fit_bayesian <- function(seed) {
  
  set.seed(seed)
  
  crm(
    dipper,
    model = "probitCJS",
    groups = "sex",
    begin.time = 1981,
    model.parameters = bayes_parameters,
    burnin = 5000,
    iter = 15000
  )
  
}


# Fit four independent MCMC chains

bayes_fits <- lapply(
  1:4,
  function(i) fit_bayesian(123 + i)
)

bayes_model <- bayes_fits[[1]]

print(bayes_model)
print(bayes_model$results$beta)
print(bayes_model$results$reals)


# Extract posterior samples

get_samples <- function(fit, parameter) {
  
  samples <- fit$results$beta.mcmc[[parameter]]
  
  if (is.null(samples)) {
    stop(
      "Posterior samples were not found for ",
      parameter,
      ". Inspect names(fit$results)."
    )
  }
  
  as.matrix(samples)
}

phi_chains <- lapply(
  bayes_fits,
  get_samples,
  parameter = "Phi"
)

p_chains <- lapply(
  bayes_fits,
  get_samples,
  parameter = "p"
)

print(colnames(phi_chains[[1]]))
print(colnames(p_chains[[1]]))

stopifnot(
  all(vapply(phi_chains, ncol, integer(1)) == 2),
  all(vapply(p_chains, ncol, integer(1)) == 1)
)


# Check chain dimensions and parameter names

check_chains <- function(chains) {
  
  stopifnot(
    length(unique(
      vapply(chains, nrow, integer(1))
    )) == 1
  )
  
  reference <- colnames(chains[[1]])
  
  if (!is.null(reference)) {
    stopifnot(
      all(vapply(
        chains,
        function(x) identical(colnames(x), reference),
        logical(1)
      ))
    )
  }
  
  invisible(TRUE)
}

check_chains(phi_chains)
check_chains(p_chains)


# MCMC diagnostics

phi_mcmc <- mcmc.list(
  lapply(phi_chains, mcmc)
)

p_mcmc <- mcmc.list(
  lapply(p_chains, mcmc)
)

cat("Survival model convergence:\n")
print(gelman.diag(phi_mcmc, autoburnin = FALSE))

cat("Detection model convergence:\n")
print(gelman.diag(p_mcmc, autoburnin = FALSE))

cat("Survival effective sample sizes:\n")
print(effectiveSize(phi_mcmc))

cat("Detection effective sample sizes:\n")
print(effectiveSize(p_mcmc))

plot(phi_mcmc)
plot(p_mcmc)


# Combine posterior samples

phi_samples <- do.call(
  rbind,
  phi_chains
)

p_samples <- do.call(
  rbind,
  p_chains
)


# Posterior coefficient summaries

posterior_summary <- function(x) {
  
  data.frame(
    Parameter = colnames(x),
    Mean = colMeans(x),
    SD = apply(x, 2, sd),
    Lower = apply(
      x, 2, quantile, 0.025
    ),
    Upper = apply(
      x, 2, quantile, 0.975
    ),
    row.names = NULL
  )
}

phi_summary <- posterior_summary(phi_samples)
p_summary <- posterior_summary(p_samples)

print(phi_summary)
print(p_summary)


# Sex-specific survival probabilities

# The model uses an intercept and a coefficient
# for the second level of the sex factor.

sex_levels <- levels(factor(dipper$sex))

stopifnot(length(sex_levels) == 2)

phi_reference <- pnorm(
  phi_samples[, 1]
)

phi_other <- pnorm(
  phi_samples[, 1] +
    phi_samples[, 2]
)

summarise_probability <- function(x) {
  
  c(
    Mean = mean(x),
    Lower = unname(quantile(x, 0.025)),
    Upper = unname(quantile(x, 0.975))
  )
}

survival_summary <- data.frame(
  Sex = sex_levels,
  rbind(
    summarise_probability(phi_reference),
    summarise_probability(phi_other)
  ),
  row.names = NULL
)

print(survival_summary)


# Difference between sexes

survival_difference <- phi_other -
  phi_reference

cat(
  "Posterior survival difference:\n"
)

print(
  summarise_probability(survival_difference)
)

cat(
  "Probability of higher survival in",
  sex_levels[2], ":",
  mean(survival_difference > 0),
  "\n"
)


# Detection probability

detection_probability <- pnorm(
  p_samples[, 1]
)

cat("Posterior detection probability:\n")

print(
  summarise_probability(detection_probability)
)


# Bayesian survival plot

survival_plot <- ggplot(
  survival_summary,
  aes(Sex, Mean)
) +
  geom_point(
    size = 3,
    colour = "steelblue"
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.15,
    colour = "steelblue"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  labs(
    title = "Bayesian apparent survival by sex",
    x = "Sex",
    y = "Annual apparent survival probability"
  )

print(survival_plot)


# Posterior distributions

posterior_data <- rbind(
  data.frame(
    Survival = phi_reference,
    Sex = sex_levels[1]
  ),
  data.frame(
    Survival = phi_other,
    Sex = sex_levels[2]
  )
)

posterior_plot <- ggplot(
  posterior_data,
  aes(
    Survival,
    colour = Sex,
    fill = Sex
  )
) +
  geom_density(alpha = 0.15) +
  theme_classic() +
  labs(
    title = "Posterior survival distributions",
    x = "Annual apparent survival probability",
    y = "Density"
  )

print(posterior_plot)


# Bayesian and maximum-likelihood comparison

# Inspect the maximum-likelihood prediction table
# to confirm the sex labels and estimated values.

print(names(mle_survival))
print(head(mle_survival))

# Compare probability-scale estimates only after
# verifying the mapping of sex groups in both models.


# Save results

write.csv(
  encounter_summary,
  "results/encounter_summary.csv",
  row.names = FALSE
)

write.csv(
  model_comparison,
  "results/model_comparison.csv",
  row.names = FALSE
)

write.csv(
  mle_survival,
  "results/mle_survival.csv",
  row.names = FALSE
)

write.csv(
  mle_detection,
  "results/mle_detection.csv",
  row.names = FALSE
)

write.csv(
  survival_summary,
  "results/bayesian_survival.csv",
  row.names = FALSE
)

write.csv(
  phi_summary,
  "results/posterior_phi_coefficients.csv",
  row.names = FALSE
)

write.csv(
  p_summary,
  "results/posterior_detection_coefficients.csv",
  row.names = FALSE
)

saveRDS(
  best_model,
  "results/best_mle_model.rds"
)

saveRDS(
  bayes_fits,
  "results/bayesian_models.rds"
)

ggsave(
  "figures/encounter_histories.png",
  encounter_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/bayesian_survival.png",
  survival_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/posterior_distributions.png",
  posterior_plot,
  width = 8,
  height = 5,
  dpi = 300
)

writeLines(
  capture.output(sessionInfo()),
  "results/session_info.txt"
)
