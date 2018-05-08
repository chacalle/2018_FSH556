
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() () {

  // Data
  DATA_VECTOR(s_i);
  DATA_VECTOR(a_i);
  DATA_VECTOR(l_i);

  // Parameters
  PARAMETER(log_l_0);
  PARAMETER(log_k);
  PARAMETER(log_sigma2_measure);
  PARAMETER(log_rho);
  PARAMETER(log_sigma2_spatial);
  PARAMETER(log_beta_0);
  PARAMETER(beta_s);

  // Random effects
  PARAMETER_VECTOR(epsilon_i);

  // Objective funcction
  int n_i = s_i.size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();

  Type l_0 = exp(log_l_0);
  Type k = exp(log_k);
  Type sigma2_measure = exp(log_sigma2_measure);
  Type rho = exp(log_rho);
  Type sigma2_spatial = exp(log_sigma2_spatial);
  Type beta_0 = exp(log_beta_0);

  // Probability of random effects
  vector<Type> dist_i(n_i);
  jnll_comp(1) -= dnorm(epsilon_i(0), Type(0), pow(sigma2_spatial, 0.5), true);
  for (int i=1; i < n_i; i++) {
    dist_i(i) = s_i(i) - s_i(i - 1);
    jnll_comp(1) -= dnorm(epsilon_i(i), pow(rho, dist_i(i)) * epsilon_i(i - 1), pow(sigma2_spatial*(1 - pow(rho, 2 * dist_i(i))), 0.5), true);
  }

  // Probability of data conditional on random effects
  vector<Type> l_inf_i(n_i);
  vector<Type> l_i_hat(n_i);
  for (int i = 0; i < n_i; i++) {
    l_inf_i(i) = beta_0 * exp(beta_s * s_i(i)) * exp(epsilon_i(i)); // this is exactly like the simulator function
    l_i_hat(i) = l_inf_i(i) - ((l_inf_i(i) - l_0) * exp(-k * a_i(i)));
    jnll_comp(0) -= dlognorm(l_i(i), log(l_i_hat(i)), pow(sigma2_measure, 0.5), true);
  }

  // Reporting
  Type jnll = jnll_comp.sum();

  REPORT(l_inf_i);
  REPORT(l_0);
  REPORT(k);
  REPORT(sigma2_measure);
  REPORT(rho);
  REPORT(sigma2_spatial);
  REPORT(beta_0);
  REPORT(beta_s);

  return jnll;
}
