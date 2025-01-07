data {
  int<lower=0> N;         // Number of observations
  vector[N] y;            // Observations
}

parameters {
  real<lower=0> lambda;   // Scale parameter (must be positive)
  real<lower=0> k;        // Shape parameter (must be positive)
}

model {
  // Likelihood
  y ~ weibull(k, lambda);
}

generated quantities {
  vector[N] log_lik;      // Log-likelihood for each observation
  vector[N] y_rep;        // Posterior predictive draws

  for (n in 1:N) {
    log_lik[n] = weibull_lpdf(y[n] | k, lambda);  // Log-likelihood
    y_rep[n] = weibull_rng(k, lambda);            // Posterior predictive draw
  }
}