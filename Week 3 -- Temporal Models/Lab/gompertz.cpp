#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( nt );
  DATA_VECTOR( log_b_t );
  DATA_VECTOR( sigma2_lognormal );

  // Parameters
  PARAMETER( log_d0 );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaM );
  PARAMETER( alpha );
  PARAMETER( rho );
  PARAMETER_VECTOR( log_d_t );


  // Objective funcction
  Type jnll = 0;

  // Probability of random coefficients
  jnll -= dnorm( log_d_t(0), log_d0, exp(log_sigmaP), true );
  for( int t=1; t<nt; t++){
    jnll -= dnorm( log_d_t(t), alpha + rho*log_d_t(t-1), exp(log_sigmaP), true );
  }

  // Probability of data conditional on fixed and random effect values
  Type sigma2_M = exp(log_sigmaM * 2);
  for( int t=0; t<nt; t++){
    jnll -= dnorm( log_b_t(t), log_d_t(t), sigma2_M + sigma2_lognormal(t), true );
  }

  // Reporting
  Type sigmaP = exp(log_sigmaP);
  Type sigmaM = exp(log_sigmaM);

  REPORT( sigmaP );
  REPORT( sigmaM );
  REPORT( log_d_t );

  ADREPORT( sigmaP );
  ADREPORT( sigmaM );

  return jnll;
}
