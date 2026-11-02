
# Bayesian survival modelling
# European dipper capture-recapture data

install.packages(c("marked", "ggplot2"))

library(marked)
library(ggplot2)

set.seed(123)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

data("dipper", package = "marked")

head(dipper)
str(dipper)
table(dipper$sex, useNA = "ifany")


# Prepare encounter histories

data <- dipper[
  !is.na(dipper$ch) & !is.na(dipper$sex),
]

data$sex <- factor(data$sex)

stopifnot(nlevels(data$sex) == 2)

histories <- strsplit(as.character(data$ch), "")

stopifnot(
  all(lengths(histories) == 7),
  all(unlist(histories) %in% c("0", "1"))
)

y <- do.call(
  rbind,
  lapply(histories, as.integer)
)

sex <- as.integer(data$sex == levels(data$sex)[2])

n <- nrow(y)
n_years <- ncol(y)
n_intervals <- n_years - 1

first_capture <- apply(
  y,
  1,
  function(x) which(x == 1)[1]
)

cat("Number of birds:", n, "\n")
cat("Monitoring occasions:", n_years, "\n")

print(table(data$sex))
print(colSums(y))


# Encounter history summary

encounters <- data.frame(
  Year = 1981:1987,
  BirdsEncountered = colSums(y)
)

encounter_plot <- ggplot(
  encounters,
  aes(Year, BirdsEncountered)
) +
  geom_line(linewidth = 0.9, colour = "steelblue") +
  geom_point(size = 2, colour = "steelblue") +
  theme_classic() +
  labs(
    title = "Dipper encounter records",
    x = "Year",
    y = "Number of birds encountered"
  )

print(encounter_plot)


# Cormack-Jolly-Seber likelihood

cjs_loglik <- function(phi, p) {
  
  if (length(p) == 1) {
    p <- rep(p, n_intervals)
  }
  
  if (is.vector(phi)) {
    phi <- matrix(
      rep(phi, each = n),
      nrow = n
    )
  }
  
  alive <- rep(1, n)
  dead <- rep(0, n)
  loglik <- 0
  
  for (t in 2:n_years) {
    
    idx <- which(first_capture < t)
    
    if (length(idx) == 0) next
    
    survival <- phi[idx, t - 1]
    
    alive_new <- alive[idx] * survival
    
    dead_new <- dead[idx] +
      alive[idx] * (1 - survival)
    
    observed <- y[idx, t] == 1
    
    alive_new[observed] <- alive_new[observed] *
      p[t - 1]
    
    dead_new[observed] <- 0
    
    alive_new[!observed] <- alive_new[!observed] *
      (1 - p[t - 1])
    
    total <- alive_new + dead_new
    
    if (any(total <= 0)) {
      return(-Inf)
    }
    
    loglik <- loglik + sum(log(total))
    
    alive[idx] <- alive_new / total
    dead[idx] <- dead_new / total
  }
  
  loglik
}


# Bayesian model definitions

log_posterior <- function(theta, model) {
  
  alpha <- theta[1]
  detect <- theta[2]
  
  p <- plogis(detect)
  
  if (model == "constant") {
    
    phi <- rep(
      plogis(alpha),
      n_intervals
    )
    
    log_prior <- sum(
      dnorm(
        theta,
        mean = 0,
        sd = 2,
        log = TRUE
      )
    )
    
  } else if (model == "sex") {
    
    beta_sex <- theta[3]
    
    phi <- matrix(
      plogis(alpha + beta_sex * sex),
      nrow = n,
      ncol = n_intervals
    )
    
    log_prior <- sum(
      dnorm(
        theta,
        mean = 0,
        sd = 2,
        log = TRUE
      )
    )
    
  } else if (model == "hierarchical") {
    
    beta_sex <- theta[3]
    
    year_effect <- theta[4:9]
    
    log_sigma <- theta[10]
    sigma <- exp(log_sigma)
    
    phi <- sapply(
      seq_len(n_intervals),
      function(t) {
        plogis(
          alpha +
            beta_sex * sex +
            year_effect[t]
        )
      }
    )
    
    log_prior <- (
      dnorm(alpha, 0, 2, log = TRUE) +
        dnorm(detect, 0, 2, log = TRUE) +
        dnorm(beta_sex, 0, 2, log = TRUE) +
        sum(
          dnorm(
            year_effect,
            0,
            sigma,
            log = TRUE
          )
        ) +
        dnorm(sigma, 0, 1, log = TRUE) +
        log(2) +
        log_sigma
    )
    
  } else {
    
    stop("Unknown model")
    
  }
  
  log_prior + cjs_loglik(phi, p)
}


# Metropolis-Hastings sampler

run_mcmc <- function(
    model,
    iterations = 6000,
    burnin = 2000,
    thin = 2,
    chains = 4
) {
  
  parameter_names <- switch(
    model,
    constant = c("Intercept", "Detection"),
    sex = c(
      "Intercept",
      "Detection",
      "SexEffect"
    ),
    hierarchical = c(
      "Intercept",
      "Detection",
      "SexEffect",
      paste0("YearEffect", 1:6),
      "LogSigma"
    )
  )
  
  n_par <- length(parameter_names)
  
  keep <- seq(
    burnin + 1,
    iterations,
    by = thin
  )
  
  samples <- vector("list", chains)
  acceptance <- numeric(chains)
  
  for (chain in seq_len(chains)) {
    
    theta <- rep(0, n_par)
    
    theta[1] <- qlogis(0.6)
    theta[2] <- qlogis(0.6)
    
    theta <- theta +
      rnorm(n_par, 0, 0.1)
    
    current_lp <- log_posterior(
      theta,
      model
    )
    
    chain_samples <- matrix(
      NA_real_,
      nrow = length(keep),
      ncol = n_par
    )
    
    colnames(chain_samples) <- parameter_names
    
    accepted <- 0
    saved <- 1
    
    proposal_sd <- rep(0.12, n_par)
    
    for (i in seq_len(iterations)) {
      
      proposal <- theta +
        rnorm(n_par, 0, proposal_sd)
      
      proposal_lp <- log_posterior(
        proposal,
        model
      )
      
      if (
        is.finite(proposal_lp) &&
        log(runif(1)) <
        proposal_lp - current_lp
      ) {
        
        theta <- proposal
        current_lp <- proposal_lp
        accepted <- accepted + 1
        
      }
      
      if (i %in% keep) {
        
        chain_samples[saved, ] <- theta
        saved <- saved + 1
        
      }
      
    }
    
    samples[[chain]] <- chain_samples
    
    acceptance[chain] <- accepted / iterations
    
    cat(
      model,
      "- chain", chain,
      "- acceptance:",
      round(acceptance[chain], 3),
      "\n"
    )
    
  }
  
  list(
    samples = samples,
    acceptance = acceptance,
    model = model
  )
}


# Overall apparent survival

model_constant <- run_mcmc("constant")

posterior_constant <- do.call(
  rbind,
  model_constant$samples
)

survival_constant <- plogis(
  posterior_constant[, "Intercept"]
)

detection_constant <- plogis(
  posterior_constant[, "Detection"]
)

summarise_posterior <- function(x) {
  
  c(
    Mean = mean(x),
    Lower = unname(quantile(x, 0.025)),
    Upper = unname(quantile(x, 0.975))
  )
  
}

print(
  summarise_posterior(survival_constant)
)

print(
  summarise_posterior(detection_constant)
)


# Prior and posterior survival distributions

prior_survival <- plogis(
  rnorm(10000, 0, 2)
)

density_data <- rbind(
  data.frame(
    Probability = prior_survival,
    Distribution = "Prior"
  ),
  data.frame(
    Probability = survival_constant,
    Distribution = "Posterior"
  )
)

prior_plot <- ggplot(
  density_data,
  aes(
    Probability,
    colour = Distribution,
    fill = Distribution
  )
) +
  geom_density(alpha = 0.15) +
  theme_classic() +
  labs(
    title = "Bayesian apparent survival probability",
    x = "Annual apparent survival probability",
    y = "Density"
  )

print(prior_plot)


# Sex differences in apparent survival

model_sex <- run_mcmc("sex")

posterior_sex <- do.call(
  rbind,
  model_sex$samples
)

survival_group1 <- plogis(
  posterior_sex[, "Intercept"]
)

survival_group2 <- plogis(
  posterior_sex[, "Intercept"] +
    posterior_sex[, "SexEffect"]
)

sex_difference <- survival_group2 -
  survival_group1

sex_summary <- data.frame(
  Sex = levels(data$sex),
  rbind(
    summarise_posterior(survival_group1),
    summarise_posterior(survival_group2)
  ),
  row.names = NULL
)

print(sex_summary)

cat(
  "Posterior probability of a positive sex effect:",
  mean(posterior_sex[, "SexEffect"] > 0),
  "\n"
)

print(
  summarise_posterior(sex_difference)
)


# Plot sex-specific survival

sex_plot <- ggplot(
  sex_summary,
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
    title = "Apparent survival by sex",
    x = "Sex",
    y = "Annual apparent survival probability"
  )

print(sex_plot)


# Hierarchical Bayesian survival model

model_hier <- run_mcmc("hierarchical")

posterior_hier <- do.call(
  rbind,
  model_hier$samples
)

year_effects <- posterior_hier[
  ,
  paste0("YearEffect", 1:6),
  drop = FALSE
]

sigma_year <- exp(
  posterior_hier[, "LogSigma"]
)

print(
  summarise_posterior(sigma_year)
)


# Annual survival estimates

year_labels <- paste(
  1981:1986,
  1982:1987,
  sep = "-"
)

annual_results <- list()

for (s in 0:1) {
  
  for (t in 1:6) {
    
    survival_draws <- plogis(
      posterior_hier[, "Intercept"] +
        posterior_hier[, "SexEffect"] * s +
        year_effects[, t]
    )
    
    annual_results[[length(annual_results) + 1]] <-
      data.frame(
        Interval = year_labels[t],
        Sex = levels(data$sex)[s + 1],
        Mean = mean(survival_draws),
        Lower = quantile(
          survival_draws,
          0.025
        ),
        Upper = quantile(
          survival_draws,
          0.975
        )
      )
    
  }
  
}

annual_survival <- do.call(
  rbind,
  annual_results
)

annual_survival$Interval <- factor(
  annual_survival$Interval,
  levels = year_labels
)

print(annual_survival)


# Plot annual survival

annual_plot <- ggplot(
  annual_survival,
  aes(
    Interval,
    Mean,
    colour = Sex,
    group = Sex
  )
) +
  geom_ribbon(
    aes(
      ymin = Lower,
      ymax = Upper,
      fill = Sex
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    title = "Annual apparent survival",
    x = "Monitoring interval",
    y = "Apparent survival probability"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(annual_plot)


# Posterior detection probability

detection_hier <- plogis(
  posterior_hier[, "Detection"]
)

print(
  summarise_posterior(detection_hier)
)


# MCMC convergence diagnostics

gelman_rhat <- function(samples) {
  
  chains <- length(samples)
  iterations <- nrow(samples[[1]])
  
  chain_means <- sapply(
    samples,
    colMeans
  )
  
  chain_vars <- sapply(
    samples,
    function(x) apply(x, 2, var)
  )
  
  W <- rowMeans(chain_vars)
  
  B <- iterations *
    apply(chain_means, 1, var)
  
  var_hat <- (
    (iterations - 1) / iterations * W +
      B / iterations
  )
  
  sqrt(var_hat / W)
}

cat("Constant model R-hat:\n")
print(gelman_rhat(model_constant$samples))

cat("Sex model R-hat:\n")
print(gelman_rhat(model_sex$samples))

cat("Hierarchical model R-hat:\n")
print(gelman_rhat(model_hier$samples))


# Trace plots

trace_data <- do.call(
  rbind,
  lapply(
    seq_along(model_hier$samples),
    function(chain) {
      
      x <- model_hier$samples[[chain]]
      
      data.frame(
        Iteration = seq_len(nrow(x)),
        Chain = factor(chain),
        SurvivalIntercept = x[, "Intercept"]
      )
      
    }
  )
)

trace_plot <- ggplot(
  trace_data,
  aes(
    Iteration,
    SurvivalIntercept,
    colour = Chain
  )
) +
  geom_line(alpha = 0.7) +
  theme_classic() +
  labs(
    title = "MCMC trace plot",
    x = "Saved iteration",
    y = "Survival intercept"
  )

print(trace_plot)


# Save results

write.csv(
  encounters,
  "results/encounter_summary.csv",
  row.names = FALSE
)

write.csv(
  sex_summary,
  "results/sex_survival.csv",
  row.names = FALSE
)

write.csv(
  annual_survival,
  "results/annual_survival.csv",
  row.names = FALSE
)

saveRDS(
  model_constant,
  "results/model_constant.rds"
)

saveRDS(
  model_sex,
  "results/model_sex.rds"
)

saveRDS(
  model_hier,
  "results/model_hierarchical.rds"
)

ggsave(
  "figures/prior_posterior.png",
  prior_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/sex_survival.png",
  sex_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/annual_survival.png",
  annual_plot,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/mcmc_trace.png",
  trace_plot,
  width = 8,
  height = 5,
  dpi = 300
)
