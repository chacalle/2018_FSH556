#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {

  // Data
  DATA_VECTOR(log_y_t);
  DATA_INTEGER(ricker);

  // Parameters
  PARAMETER(log_x_0);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(theta);

  PARAMETER_VECTOR(log_x_t);

  int data_years = log_y_t.size();
  int state_years = log_x_t.size();

  // Objective function
  Type jnll = 0;

  if (ricker == 1) {
    theta = exp(theta); // both r and K need to be positive
  }

  // Probability of random coefficients
  jnll -= dnorm(log_x_0, log_x_t(0), exp(log_sigma), true);
  for (int y = 1; y < state_years; y++) {
    Type current_log_x_t_value;
    if (ricker == 0) {
      current_log_x_t_value = theta(0) + theta(1) * log_x_t(y - 1);
    } else {
      current_log_x_t_value = exp(log_x_t(y - 1)) * exp(theta(0) * (1 - (exp(log_x_t(y - 1)) / theta(1))));
      current_log_x_t_value = log(current_log_x_t_value);
    }
    jnll -= dnorm(log_x_t(y), current_log_x_t_value, exp(log_sigma), true);
  }

  // Probability of data conditional on fixed and random coefficients
  for (int y = 0; y < data_years; y++) {
    jnll -= dnorm(log_y_t(y), log_x_t(y), exp(log_sigma), true);
  }

  REPORT(log_x_t);

  return(jnll);
}
