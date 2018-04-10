
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() () {
  // Data
  DATA_INTEGER(model_version);
  DATA_VECTOR(y_i);
  DATA_MATRIX(X_ij);
  DATA_VECTOR(holdout_i);

  int n_data = y_i.size();
  int n_j = X_ij.row(0).size();

  // Parameters
  PARAMETER_VECTOR(b_j);

  PARAMETER_VECTOR(theta_z);
  Type zero_prob = invlogit(theta_z(0));

  // Linear predictor
  // Linear predictor
  vector<Type> log_linpred_i( n_data );
  for(int i = 0; i < n_data; i++){
    log_linpred_i(i) = 0;
    for(int j = 0; j < n_j; j++){
      log_linpred_i(i) += X_ij(i,j) * b_j(j);
    }
  }

  vector<Type> linpred_i(n_data);
  linpred_i= exp(log_linpred_i);



  // used in lognormal model
  Type lognorm_sd = exp(theta_z(1));

  // used in gamma model
  Type CV = exp(theta_z(1)); // log_CV
  Type shape = pow(CV, -2);
  vector<Type> scale = linpred_i * pow(CV, 2);

  // Objective function
  vector<Type> jnll_i(n_data);
  Type jnll = 0;
  Type pred_jnll = 0;


  // Probability of data conditional on fixed effect values
  for(int i = 0; i < n_data; i++){

    // when the catch is zero
    if (y_i(i) == 0) {
      jnll_i(i) -= log(zero_prob);
    } else {
      // when the catch is greater than zero
      if (model_version == 1 | model_version == 3) { // lognormal model

        jnll_i(i) -= log(1 - zero_prob) + dlognorm(y_i(i), log_linpred_i(i), lognorm_sd, true);

      } else if (model_version == 2) { // gamma model

        jnll_i(i) -= log(1 - zero_prob) + dgamma(y_i(i), shape, scale(i), true);

      }
    }

    // Running counter
    if(holdout_i(i) == 0) jnll += jnll_i(i);
    if(holdout_i(i) == 1) pred_jnll += jnll_i(i);
  }

  // Reporting
  REPORT( linpred_i );
  REPORT( pred_jnll );
  REPORT( jnll_i );


  return jnll;
}
