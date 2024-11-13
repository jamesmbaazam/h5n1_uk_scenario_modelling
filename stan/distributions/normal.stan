data {
  int<lower=0> N;         // Number of observations
  vector[N] y;            // Observations
}

parameters {
  real mu;                // Mean parameter
  real<lower=0> sigma;    // Standard deviation (must be positive)
}

model {
  // Likelihood
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;      // Log-likelihood for each observation
  vector[N] y_rep;        // Posterior predictive draws

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu, sigma);  // Log-likelihood
    y_rep[n] = normal_rng(mu, sigma);            // Posterior predictive draw
  }
}