data {
  int<lower=0> N;         // Number of observations
  vector[N] y;            // Observations
  int<lower=0> N_rep;
}

parameters {
  real<lower=0> alpha;    // Shape parameter (must be positive)
  real<lower=0> beta;     // Rate parameter (must be positive)
}

model {
  // Likelihood
  y ~ gamma(alpha, beta);
}

generated quantities {
  vector[N] log_lik;      // Log-likelihood for each observation
  vector[N_rep] y_rep;        // Posterior predictive draws

  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(y[n] | alpha, beta);  // Log-likelihood
  }
  
  for(i in 1:N_rep) {
    y_rep[i] = gamma_rng(alpha, beta);            // Posterior predictive draw
  }
}