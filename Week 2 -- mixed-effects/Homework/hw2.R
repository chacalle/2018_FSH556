
library(TMB)
library(data.table)
library(ggplot2)

rm(list = ls())
set.seed(30)

simulate_data <- function(n_sites = 10, n_counts_per_site = 10,
                          mu = 2, sigma_s = sqrt(1), sigma_y = sqrt(0.5)) {

  sims <- lapply(1:n_sites, function(s) {
    e_s <- rnorm(n = 1, mean = 0, sd = sigma_s)
    d_i <- rnorm(n = n_counts_per_site, mean = 0, sd = sigma_y)
    mean_y_i <- mu + d_i + e_s
    y_i <- rpois(n = n_counts_per_site, lambda = exp(mean_y_i))

    site_sims <- data.table(site = s, count = 1:n_counts_per_site,
                            epsilon_s = e_s, delta_i = d_i,
                            mean_y_i = mean_y_i, y_i = y_i)
    return(site_sims)
  })
  sims <- rbindlist(sims)
  return(sims)
}

sims <- simulate_data()

plot <- ggplot(data = sims, aes(x = y_i, fill = as.factor(site))) +
  geom_density(color = "black") +
  facet_wrap(~ site) +
  labs(title = "Site-specific simulations") +
  theme_bw()
print(plot)


# Compile model
compile("hw2.cpp")

fit_mod <- function(sim, among_site_variability = T, overdispersion = T) {

  n_sites <- max(sim$site)
  n_counts_per_site <- max(sim$count)

  # Build inputs
  sample <- dcast(sim, formula = site ~ count, value.var = "y_i")
  sample[, site := NULL]
  sample <- as.matrix(sample)

  Data = list(y_sc = sample)
  Parameters = list(mu = 5,
                    epsilon_s = rep(0, n_sites), log_sigma_s = 0,
                    delta_c = matrix(0, nrow = n_sites, ncol = n_counts_per_site), log_sigma_y = 0)
  Random = c("epsilon_s", "delta_c")

  # Define what variables to fix
  Map <- list()
  if (!among_site_variability) {
    Map[["epsilon_s"]] <- factor(rep(NA, n_sites))
    Map[["log_sigma_s"]] <- factor(NA)
    Parameters[["log_sigma_s"]] <- log(1e-10)
  }

  if (!overdispersion) {
    Map[["delta_c"]] <- factor(rep(NA, n_sites * n_counts_per_site))
    Map[["log_sigma_y"]] <- factor(NA)
    Parameters[["log_sigma_y"]] <- log(1e-10)
  }

  # Build object
  dyn.load(dynlib("hw2"))
  Obj = MakeADFun(data = Data, parameters = Parameters, random = Random, map = Map)

  # Optimize
  Opt = TMBhelper::Optimize(obj = Obj)

  # Get SEs and CIs

  means = as.list( Opt$SD, "Estimate" )
  std_error = as.list( Opt$SD, "Std. Error" )

  fit_results <- data.table(mean = means[["mu"]], se = std_error[["mu"]])
  fit_results[, lower := mean - (1.96 * se)]
  fit_results[, upper := mean + (1.96 * se)]

  return(fit_results)
}

result <- fit_mod(sims)

# create simulated data
n_simulated_datasets <- 100
full_simulated_data <- lapply(1:n_simulated_datasets, function(i) {
  sim <- simulate_data()
  sim[, sim := i]
})
full_simulated_data <- rbindlist(full_simulated_data)

fit_results <- lapply(1:4, function(mod) {
  mod_results <- lapply(1:n_simulated_datasets, function(i) {
    result <- fit_mod(full_simulated_data[sim == i],
                      among_site_variability = mod %in% c(2, 4),
                      overdispersion = mod %in% c(3, 4))
    result[, model_version := mod]
    result[, sim := i]
    return(result)
  })
  mod_results <- rbindlist(mod_results)
})
fit_results <- rbindlist(fit_results)

fit_results[, contains_truth := between(2, lower, upper)]
summary <- fit_results[, list(mu = mean(mean), coverage = sum(contains_truth) / .N), by = c("model_version")]
print(summary)
