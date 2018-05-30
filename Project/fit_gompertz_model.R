rm(list = ls())

library(data.table)
library(INLA)
library(RandomFields)
library(mvtnorm)
library(ggplot2)
library(TMB)
library(TMBhelper)

source("simulator.R")

sim <- sim_gompertz_fn(n_per_year = 50)
sim_data <- sim$sim_data
loc_xy <- sim$loc_xy[, s_i := NULL]

# Make triangulated mesh
mesh = inla.mesh.create(loc_xy)
spde = inla.spde2.matern(mesh)

# Make inputs
Data = list("n_i" = nrow(sim_data), "n_x" = length(unique(sim_data$s_i)),
            "n_t" = length(unique(sim_data$year)), "n_a" = length(unique(sim_data$age)),
            "y_s" = c(-1.07889, -0.31188, 0.35380, 1.05695, 1.95343, 3.41302, 6.05569),
            "s_i" = sim_data[, s_i - min(s_i)], "t_i" = sim_data[, year - min(year)],
            "a_i" = sim_data[, (age - min(age)) %/% 5], "N_i" = sim_data[, N], "B_i" = sim_data[, B],
            "M0"=spde$param.inla$M0, "M1"=spde$param.inla$M1, "M2"=spde$param.inla$M2)
Params = list("log_tau_log_T" = log(1 / 0.05), "log_tau_alpha" = log(1 / 0.1), "log_tau_log_beta" = log(1 / 0.1),
              "log_kappa" = 0.0, "log_mu_T" = log(4.0),
              "log_T" = matrix(0, nrow=mesh$n, ncol=Data$n_t),
              "alpha" = matrix(0, nrow=mesh$n, ncol=Data$n_t),
              "log_beta" = matrix(0, nrow=mesh$n, ncol=Data$n_t))
Random = c("log_T", "alpha", "log_beta")

# Compile
Version = "gompertz_model"
compile(paste0(Version,".cpp"))
dyn.load(dynlib(Version))

# Build and run
Obj = MakeADFun( data=Data, parameters=Params, random=Random )

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1)
#Report = Obj$report()
#unlist( Report[c('Range','SigmaO','SigmaE')] )


