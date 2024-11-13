data {
  int<lower=0> N;         // Number of observations
  vector[N] y;            // Observations
  int<lower=0> N_rep;
}

parameters {
  real mu;                // Location parameter
  real<lower=0> sigma;    // Scale parameter (must be positive)
}

model {
  // Likelihood
  y ~ lognormal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;      // Log-likelihood for each observation
  vector[N_rep] y_rep;        // Posterior predictive draws

  for (n in 1:N) {
    log_lik[n] = lognormal_lpdf(y[n] | mu, sigma);  // Log-likelihood
  }
  
  for(i in 1:N_rep) {
    y_rep[i] = lognormal_rng(mu, sigma);            // Posterior predictive draw
  }
}