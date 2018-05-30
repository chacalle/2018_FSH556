
library(data.table)
library(INLA)
library(RandomFields)
library(mvtnorm)
library(ggplot2)

sim_gompertz_fn <- function(ages = seq(15, 45, 5), n_per_year = 100,
                            years = 2016,
                            sd_tfr = 0.05, sd_alpha = 0.1, sd_beta = 0.1,
                            spatial_scale = 0.1, log_mu_tfr = log(4.0)) {

  # standard from Booth 1984
  y_s <- c(-1.07889, -0.31188, 0.35380, 1.05695,
           1.95343, 3.41302, 6.05569)

  # create simulation locations
  loc_xy <- data.table(s_i = 1:n_per_year, x = runif(n_per_year), y = runif(n_per_year))

  # simulate spatial random effects (tfr)
  model_tfr <- RandomFields::RMgauss(var=sd_tfr^2, scale=spatial_scale)
  model_alpha <- RandomFields::RMgauss(var=sd_alpha^2, scale=spatial_scale)
  model_beta <- RandomFields::RMgauss(var=sd_beta^2, scale=spatial_scale)

  sim_data <- data.table(s_i = 1:n_per_year, year = years)
  sim_data <- sim_data[, list(s_i,
                              log_tfr = RFsimulate(model=model_tfr, x=loc_xy$x, y=loc_xy$y)@data[,1],
                              alpha = RFsimulate(model=model_alpha, x=loc_xy$x, y=loc_xy$y)@data[,1],
                              log_beta = RFsimulate(model=model_beta, x=loc_xy$x, y=loc_xy$y)@data[,1]), by = "year"]
  sim_data[, beta := exp(log_beta)]
  sim_data[, log_mu_tfr := log_mu_tfr]
  sim_data[, tfr := exp(log_tfr + log_mu_tfr)]
  sim_data <- sim_data[, list(tfr, log_tfr, alpha, beta, age = ages, y_s = y_s), by = c("year", "s_i")]
  sim_data[, y_m := alpha + beta * y_s]
  sim_data[, F_m_hat := exp(-exp(-y_m))]
  sim_data[, F_m := tfr * F_m_hat]
  sim_data[, f := c(F_m[1], diff(F_m)) / 5, by = c("year", "s_i")]
  sim_data[, N := runif(.N, 5, 25)]
  sim_data[, B := rpois(.N, f * N)]

  result <- list(sim_data = sim_data,
                 loc_xy = loc_xy)
  return(result)
}

sim_fn <- function(ages = seq(15, 45, 5), n_per_year = 100,
                   years = 2000:2015, survey_years = c(2000, 2007, 2015),
                   sd_omega_l = 0.1, sd_epsilon_lt = 0.05, spatial_scale = 0.1,
                   sd_a = 0.4, rho_a = 0.9, sd_t = 0.7, rho_t = 0.1,
                   sd_observation = 0.02) {
  # log mean intercepts by age
  alpha_a <- data.table(age = ages, alpha_a = log(c(0.1, 0.2, 0.2, 0.15, 0.1, 0.05, 0.01)))

  # create simulation locations
  loc_xy <- data.table(s_i = 1:n_per_year, x = runif(n_per_year), y = runif(n_per_year))

  # simulate spatial random effects (omega and epsiolon)
  model_O <- RandomFields::RMgauss(var=sd_omega_l^2, scale=spatial_scale)
  model_E <- RandomFields::RMgauss(var=sd_epsilon_lt^2, scale=spatial_scale)

  # Simulate Omega
  omega <- data.table(s_i = 1:n_per_year)
  omega[, omega_l := RFsimulate(model = model_O, x=loc_xy$x, y=loc_xy$y)@data[,1]]

  # Simulate Epsilon
  epsilon <- data.table(s_i = 1:n_per_year)
  for(t in years){
    epsilon[, paste0(t) := list(RFsimulate(model=model_E, x=loc_xy$x, y=loc_xy$y)@data[,1])]
  }
  epsilon <- melt(epsilon, id.vars = c("s_i"), measure.vars = paste0(years), variable.name = "year", value.name = "epsilon_lt")
  epsilon[, year := as.integer(as.character(year))]

  # construct precision and covariance matrix for age
  Q_a <- construct_AR1_precision(n = length(ages), sigma = sd_a, rho = rho_a)
  C_a <- solve(Q_a)

  # construct precision and covariance matrix for time
  Q_t <- construct_AR1_precision(n = length(years), sigma = sd_t, rho = rho_t)
  C_t <- solve(Q_t)

  # combine covariance matrices together
  C_at <- kronecker(C_a, C_t)

  # simulate age-time random effect
  delta_at <- CJ(age = ages, year = years) # I think I have the ordering on this correct
  delta_at <- delta_at[, delta_at := as.vector(rmvnorm(1, mean = rep(0, length(ages) * length(years)), sigma = C_at))]

  # combine all the random effects into one data frame
  sim_data <- CJ(s_i = 1:n_per_year, year = years, age = ages)
  sim_data <- merge(sim_data, alpha_a, by = "age")
  sim_data <- merge(sim_data, omega, by = c("s_i"))
  sim_data <- merge(sim_data, epsilon, by = c("s_i", "year"))
  sim_data <- merge(sim_data, delta_at, by = c("age", "year"))


  sim_data[, mu_lat := alpha_a + omega_l + epsilon_lt + delta_at]
  sim_data[, log_f_lat := rnorm(n = 1, mean = mu_lat, sd = sd_observation), by = c("s_i", "age", "year")]
  sim_data[, f_lat := exp(log_f_lat)]

  sim_data <- merge(sim_data, loc_xy, by = "s_i")
  sim_data[year %in% survey_years, survey_year := T]
  return(sim_data)
}
