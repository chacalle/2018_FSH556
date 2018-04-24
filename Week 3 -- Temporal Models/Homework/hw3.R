
library(TMB)
library(data.table)
library(ggplot2)

rm(list=ls())
set.seed(3)

# download data
cpue_data = FishData::download_catch_rates( survey="Eastern_Bering_Sea", species_set="Gadus chalcogrammus", error_tol=0.01, localdir=paste0(getwd(),"/") )
cpue_data <- data.table(cpue_data)
cpue_data <- cpue_data[, list(year = Year, cpue = Wt)]

# get rid of missing data points
cpue_data[cpue < 0, cpue := NA]
cpue_data <- cpue_data[!is.na(cpue)]

# calculate average cpue by year
cpue_data <- cpue_data[, list(cpue = mean(cpue)), by = "year"]

# Compile model
compile("hw3.cpp")
dyn.load(dynlib("hw3"))

fit_model <- function(fit_data, holdout_years = c(), ricker = F, quantile = 0.95) {

  original_years <- unique(fit_data$year)

  fit_data <- fit_data[!year %in% holdout_years]

  # Build inputs
  Data = list(log_y_t = log(fit_data$cpue), ricker = as.integer(ricker))
  Parameters = list(log_x_0 = 0, log_sigma = 1, theta = c(0, 0), log_x_t = rep(0, length(original_years)))
  Random = c("log_x_t")

  # Build object
  Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, DLL="hw3")

  # Optimize
  Opt = TMBhelper::Optimize(obj=Obj)

  # compile together estimates
  means = as.list( Opt$SD, "Estimate" )
  std_error = as.list( Opt$SD, "Std. Error" )

  fixed_variables <- c("log_x_0", "log_sigma", if (!ricker) c("alpha", "rho") else c("log_r", "log_K"))
  fixed_coefficients <- data.table(variable = fixed_variables, mean = unlist(means[c("log_x_0", "log_sigma", "theta")]), se = unlist(std_error[c("log_x_0", "log_sigma", "theta")]))
  random_coefficients <- data.table(variable = paste0("log_x_", original_years), mean = means[["log_x_t"]], se = std_error[["log_x_t"]])

  coefficients <- rbind(fixed_coefficients, random_coefficients)
  coefficients[, lower := mean - (qnorm((quantile + (1 - quantile) / 2)) * se)]
  coefficients[, upper := mean + (qnorm((quantile + (1 - quantile) / 2)) * se)]

  # compile together parameters
  parameters <- data.table(variable = c("log_x_0", "log_sigma", "theta1", "theta2"),
                           value = Opt$par)
  parameters[variable == "theta1", variable := "alpha"]
  parameters[variable == "theta2", variable := "rho"]

  return(list(Opt = Opt,
              coefficients = coefficients))
}

state_space_results <- fit_model(fit_data = cpue_data, ricker = F)
ricker_results <- fit_model(fit_data = cpue_data, ricker = T)

simulate_data <- function(parameters, years = 1982:2017) {
  log_x_t <- rep(NA, length(years))
  log_x_t[1] <- parameters[variable == "log_x_0", mean]
  for (y in 2:length(years)) {
    mean <- parameters[variable == "alpha", mean] + (parameters[variable == "rho", mean] * log_x_t[y - 1])
    log_x_t[y] <- rnorm(n = 1, mean = mean, sd = exp(parameters[variable == "log_sigma", mean]))
  }

  y_t <- rep(NA, length(years))
  for(y in 1:length(years)) {
    y_t[y] <- exp(rnorm(n = 1, mean = log_x_t[y], sd = exp(parameters[variable == "log_sigma", mean])))
  }

  sim_data <- data.table(year = years, cpue = y_t)
  return(sim_data)
}

simulation_experiment <- function(parameters, replicates = 100, ricker = F) {
  results <- lapply(1:replicates, function(i) {
    sims <- simulate_data(parameters)

    result <- fit_model(fit_data = sims, holdout_years = 2013:2017, ricker = ricker, quantile = 0.5)

    fixed_vars <- c("log_x_0", "log_sigma", "alpha", "rho", "log_r", "log_K")
    fixed_parameters <- state_space_results$coefficients[variable %in% fixed_vars, list(variable, fixed = mean)]

    sims <- sims[, list(variable = paste0("log_x_", year), data = log(cpue))]
    coef <- merge(result$coefficients, sims, by = "variable", all.x = T)
    coef <- merge(coef, fixed_parameters, by = "variable", all.x = T)
    coef[variable %in% fixed_vars, data := fixed]
    coef[, fixed := NULL]
    coef[, contains_truth := between(data, lower, upper)]
    coef[, replicate := i]

    return(coef)

  })
  results <- rbindlist(results)
  return(results)
}

sim_results <- simulation_experiment(state_space_results$coefficients)
sim_results_ricker <- simulation_experiment(state_space_results$coefficients, ricker = T)

plot_variables <- c("log_x_0", "log_sigma", "alpha", "rho")
ggplot(sim_results[variable %in% plot_variables], aes(x = replicate, y = mean, colour = contains_truth)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = data)) +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "replicate", y = "value") +
  theme_bw()
plot_variables <- c(paste0("log_x_", 2010:2017))
ggplot(sim_results[variable %in% plot_variables], aes(x = replicate, y = mean, colour = contains_truth)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "replicate", y = "value") +
  theme_bw()

coverage <- sim_results[, list(coverage = sum(contains_truth) / .N), by = c("variable")]
coverage <- coverage[variable %in% paste0("log_x_", 1982:2017)]
coverage[, year := as.integer(substr(variable, 7, 10))]
ggplot(coverage, aes(x = year, y = coverage)) +
  geom_bar(stat = "identity")

