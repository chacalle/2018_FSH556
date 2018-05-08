library(TMB)
library(data.table)
library(ggplot2)

rm(list=ls())
set.seed(5)
setwd("/Users/chacalle/Documents/classes/2018_FSH556/Week 5 -- 1D spatial models/Homework/")
source("Homework_Week_5_simulator.R")
compile("hw5.cpp")
dyn.load(dynlib("hw5"))

run_simulation_experiment <- function(include_spatial_trend = T) {

  # generate simulated data
  sim_data <- data.table(Sim_Fn())
  setnames(sim_data, c("s_i", "a_i", "l_i", "true_linf_i")) # to be consistent with hw description
  sim_data <- sim_data[order(s_i)] # order data by latitude

  # Build inputs
  Data = list(s_i = sim_data$s_i, a_i = sim_data$a_i, l_i = sim_data$l_i)
  Parameters = list(log_l_0 = log(1), log_k = log(0.1), log_sigma2_measure = 1,
                    log_rho = log(0.1), log_sigma2_spatial = 1,
                    log_beta_0 = log(20), beta_s = 0,
                    epsilon_i = rep(0, nrow(sim_data)))
  Random = c("epsilon_i")

  # optionally include spatial trend
  Map <- list()
  if (!include_spatial_trend) {
    Map[["beta_s"]] <- factor(NA)
  }

  # Build object
  Obj = MakeADFun(data = Data, parameters = Parameters, random = Random, map = Map, DLL="hw5")

  # Obj$fn( Obj$par)
  # Obj$gr( Obj$par )
  # Obj$fn( Obj$env$last.par)
  # Obj$gr( Obj$par )

  # Optimize
  Opt = TMBhelper::Optimize(obj=Obj)

  # calculate rmse for asymptotic maximum size
  Report = Obj$report()
  sim_data[, pred_linf_i := Report$l_inf_i]
  Report$l_inf_i <- NULL
  Report$rmse <- sim_data[, sqrt(mean((true_linf_i - pred_linf_i) ^ 2))]
  setDT(Report)

  Report[, include_spatial_trend := include_spatial_trend]

  return(Report)
}

simulation_results <- lapply(1:10, function(i) {
  results <- lapply(c(F,T), function(v) {
    run_simulation_experiment(v)
  })
  results <- rbindlist(results)
  results[, experiment := i]
})
simulation_results <- rbindlist(simulation_results)
simulation_results <- melt(simulation_results, id.vars = c("include_spatial_trend", "experiment"))
mean_results <- simulation_results[, list(mean = mean(value)), by = c("include_spatial_trend", "variable")]

p <- ggplot(simulation_results, aes(x = value, fill = include_spatial_trend)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = mean_results, aes(xintercept = mean, colour = include_spatial_trend)) +
  facet_wrap( ~ variable, scales = "free") +
  theme_bw()
print(p)

