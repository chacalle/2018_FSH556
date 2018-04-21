#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {

  DATA_MATRIX(y_sc);

  PARAMETER(mu);

  PARAMETER_VECTOR(epsilon_s);
  PARAMETER(log_sigma_s);

  PARAMETER_MATRIX(delta_c);
  PARAMETER(log_sigma_y);

  int n_sites = y_sc.rows();
  int n_counts = y_sc.cols();

  // Objective funcction
  Type jnll = 0;

  // Probability of data conditional on fixed and random effect values
  matrix<Type> mean_y_sc;
  mean_y_sc.setZero(n_sites, n_counts);

  for (int s = 0; s < n_sites; s++) {
    for (int c = 0; c < n_counts; c++) {
      mean_y_sc(s, c) = exp(mu + epsilon_s(s)  + delta_c(s, c));
      jnll -= dpois(y_sc(s, c), mean_y_sc(s, c), true);
    }
  }

  // Probability of site specific random coefficients
  for (int s = 0; s < n_sites; s++) {
    jnll -= dnorm(epsilon_s(s), Type(0), exp(log_sigma_s), true);
  }

  // Probability of observation specific random coefficients (overdispersion)
  for (int s = 0; s < n_sites; s++) {
    for (int c = 0; c < n_counts; c++) {
      jnll -= dnorm(delta_c(s, c), Type(0), exp(log_sigma_y), true);
    }
  }

  return jnll;
}
