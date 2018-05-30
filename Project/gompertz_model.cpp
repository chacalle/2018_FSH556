// Space time
#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Indices
  DATA_INTEGER(n_i); // Total number of observations
  DATA_INTEGER(n_x); // Number of vertices in SPDE mesh
  DATA_INTEGER(n_t); // Number of years
  DATA_INTEGER(n_a); // number of age groups

  // Data
  DATA_VECTOR(y_s); // Standard fertility age pattern (n_a)

  DATA_IVECTOR(s_i);
  DATA_IVECTOR(t_i);
  DATA_IVECTOR(a_i);
  DATA_VECTOR(N_i);
  DATA_VECTOR(B_i);

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Fixed effects
  PARAMETER(log_tau_log_T); // log inverse SD of TFR
  PARAMETER(log_tau_alpha); // log inverse SD of alpha
  PARAMETER(log_tau_log_beta); // log inverse SD of beta
  PARAMETER(log_kappa);
  PARAMETER(log_mu_T); // log mean of TFR

  // Random effects
  PARAMETER_ARRAY(log_T); // n_x x n_t
  PARAMETER_ARRAY(alpha); // n_x x n_t
  PARAMETER_ARRAY(log_beta); // n_x x n_t

  // objective function -- joint negative log-likelihood
  using namespace density;
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();

  // Derived quantities
  Type Range = sqrt(8) / exp(log_kappa);
  Type Sigma_log_T = 1 / sqrt(4 * M_PI * exp(2 * log_tau_log_T) * exp(2 * log_kappa));
  Type Sigma_alpha = 1 / sqrt(4 * M_PI * exp(2 * log_tau_alpha) * exp(2 * log_kappa));
  Type Sigma_log_beta = 1 / sqrt(4 * M_PI * exp(2 * log_tau_log_beta) * exp(2 * log_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4 * log_kappa) * M0 + Type(2.0)* exp(2*log_kappa) * M1 + M2;
  for ( int t=0; t<n_t; t++){
    jnll_comp(1) += SCALE(GMRF(Q), 1 / exp(log_tau_log_T))(log_T.col(t));
    jnll_comp(2) += SCALE(GMRF(Q), 1 / exp(log_tau_alpha))(alpha.col(t));
    jnll_comp(3) += SCALE(GMRF(Q), 1 / exp(log_tau_log_beta))(log_beta.col(t));
  }

  // Objects for derived values
  array<Type> T(n_x, n_t);
  array<Type> beta(n_x, n_t);
  array<Type> y_m(n_a, n_x, n_t);
  array<Type> F_m_hat(n_a, n_x, n_t);
  array<Type> F_m(n_a, n_x, n_t);
  array<Type> f(n_a, n_x, n_t);

  // transform random effects
  for (int x = 0; x < n_x; x++) {
    for (int t = 0; t < n_t; t++) {
      T(x, t) = exp(log_T(x, t) + log_mu_T);
      beta(x, t) = exp(log_beta(x, t));
    }
  }

  // gompertz model
  for (int a = 0; a < n_a; a++) {
    for (int x = 0; x < n_x; x++) {
      for (int t = 0; t < n_t; t++) {
        // std::cout << "a: " << a << std::endl;
        // std::cout << "x: " << x << std::endl;
        // std::cout << "t: " << t << std::endl;
        // std::cout << "alpha(x, t): " << alpha(x, t) << std::endl;
        // std::cout << "beta(x, t): " << beta(x, t) << std::endl;
        // std::cout << "y_s(a): " << y_s(a) << std::endl;

        y_m(a, x, t) = alpha(x, t) + (beta(x, t) * y_s(a));
        // std::cout << "y_m(a, x, t): " << y_m(a, x, t) << std::endl;
        F_m_hat(a, x, t) = exp(-exp(-y_m(a, x, t)));
        // std::cout << "F_m_hat(a, x, t): " << F_m_hat(a, x, t) << std::endl;
        F_m(a, x, t) = T(x, t) * F_m_hat(a, x, t);
        // std::cout << "T(x, t): " << T(x, t) << std::endl;
        // std::cout << "F_m(a, x, t): " << F_m(a, x, t) << std::endl;
        // transform from cumulative fertility to ASFR
        if (a == 0) {
          f(a, x, t) = F_m(a, x, t);
        } else {
          f(a, x, t) = F_m(a, x, t) - F_m(a - 1, x, t);
        }
        f(a, x, t) = f(a, x, t) / 5;
      }
    }
  }
  // std::cout << "Here is the matrix f:\n" << f << std::endl;


  // Probability of data conditional on random effects
  for (int i = 0; i < n_i; i++) {
    if (!isNA(N_i(i)) & !isNA(B_i(i))) {
      jnll_comp(0) -= dpois(B_i(i), N_i(i) * f(a_i(i), s_i(i), t_i(i)), true);
    }
  }

  // Objective function
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( jnll_comp );
  REPORT( jnll );

  return jnll;
}
