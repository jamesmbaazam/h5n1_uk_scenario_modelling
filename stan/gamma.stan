data {
  int<lower=0> N;           // number of observations
  int N_rep;
  vector<lower=0>[N] y;     // observed intervals
}
parameters {
  real<lower=0> alpha;      // shape parameter
  real<lower=0> beta;       // rate parameter
}
model {
  y ~ gamma(alpha, beta);   // likelihood
}
generated quantities {
  vector[N] y_rep; // Simulated data
  for (n in 1:N) {
    y_rep[n] = gamma_rng(alpha, beta); // Draws from the posterior predictive distribution
  }
}