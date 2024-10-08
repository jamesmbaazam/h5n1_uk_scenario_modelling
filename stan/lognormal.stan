data {
  int<lower=0> N;           // number of observations
  int N_rep;
  vector<lower=0>[N] y;     // observed intervals
}
parameters {
  real mu;                  // location parameter (mean of the log)
  real<lower=0> sigma;      // scale parameter (std dev of the log)
}
model {
  y ~ lognormal(mu, sigma); // likelihood
}
generated quantities {
  vector[N] y_rep; // Simulated data
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n] = lognormal_rng(mu, sigma); // Draws from the posterior predictive distribution
  }
  
  // Log-likelihood for LOO
  for (n in 1:N) {
    log_lik[n] = lognormal_lpdf(y[n] | mu, sigma);
  }
}