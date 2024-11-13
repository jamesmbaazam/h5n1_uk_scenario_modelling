data {
  int<lower=1> n_vl;            // Number of viral load data points
  vector[n_vl] t_vl;            // Time points of viral load measurements since symptom onset
  vector[n_vl] y_vl;            // Viral load measurements (log_10 scale)

  int<lower=1> n_serial;        // Number of observed serial intervals
  vector[n_serial] y_serial;    // Observed serial intervals

  real<lower=0> detection_threshold;  // Detection threshold for PCR (log_10 scale)
  real<lower=0> max_time;             // Maximum time for the time steps
  int<lower=1> n_time_steps;          // Number of time steps for generated quantities
}

parameters {
  // Parameters for the incubation period (Lognormal distribution)
  real mu_log_inc;
  real<lower=0> sigma_log_inc;

  // Parameters for the serial interval (Gamma distribution)
  real<lower=0> alpha_serial;
  real<lower=0> beta_serial;

  // Parameters for time to peak viral load (Weibull distribution)
  real<lower=0> alpha_peak;
  real<lower=0> beta_peak;

  // Viral load parameters
  real<lower=0> vl_peak;        // Peak viral load (log_10 scale)
  real<lower=0> mu_vl_m2;       // Fall rate (positive)

  real<lower=0> sigma;          // Measurement error (log_10 scale)
}

transformed parameters {
  vector[n_vl] expected_vl;

  // Mean incubation period from Lognormal distribution
  real E_incubation_period = exp(mu_log_inc + 0.5 * square(sigma_log_inc));

  // Mean time to peak viral load from Weibull distribution
  real E_time_to_peak = beta_peak * tgamma(1 + 1 / alpha_peak);

  // Rise rate calculated from peak and time to peak
  real mu_vl_m1 = vl_peak / E_time_to_peak;

  for (i in 1:n_vl) {
    // Time since infection
    real t_since_inf = t_vl[i] + E_incubation_period;

    if (t_since_inf <= E_time_to_peak) {
      // Before the peak (rising phase)
      expected_vl[i] = mu_vl_m1 * t_since_inf;
    } else {
      // After the peak (declining phase)
      expected_vl[i] = vl_peak - mu_vl_m2 * (t_since_inf - E_time_to_peak);
    }
  }
}

model {
  // Priors for the incubation period parameters (Lognormal)
  mu_log_inc ~ normal(log(7), 0.1);
  sigma_log_inc ~ normal(log(0.5), 0.05);

  // Priors for the serial interval parameters
  alpha_serial ~ normal(9, 2);
  beta_serial ~ normal(1.2, 0.5);

  // Priors for time to peak viral load parameters
  alpha_peak ~ normal(2, 0.5);
  beta_peak ~ normal(5, 1);

  // Prior for peak viral load
  vl_peak ~ normal(5, 1);

  // Prior for fall rate
  mu_vl_m2 ~ normal(0.2, 0.1);

  // Measurement error prior
  sigma ~ normal(0.3, 0.1) T[0,];

  // Likelihoods
  // Observed serial intervals
  y_serial ~ gamma(alpha_serial, beta_serial);
  y_vl ~ normal(expected_vl, sigma);
}

generated quantities {
  real time_interval = max_time / n_time_steps;
  vector[n_time_steps] pcr_detection_prob;
  vector[n_time_steps] vl_sim;
  vector[n_time_steps] expected_vl_sim;

  for (t in 1:n_time_steps) {
    real t_current = (t - 1) * time_interval;  // Start at zero
    real t_since_inf = t_current;

    if (t_since_inf <= E_time_to_peak) {
      // Before the peak (rising phase)
      expected_vl_sim[t] = mu_vl_m1 * t_since_inf;
      vl_sim[t] = normal_rng(expected_vl_sim[t], sigma);
    } else {
      // After the peak (declining phase)
      expected_vl_sim[t] = vl_peak - mu_vl_m2 * (t_since_inf - E_time_to_peak);
      vl_sim[t] = normal_rng(expected_vl_sim[t], sigma);
    }

    // Calculate the PCR detection probability
     pcr_detection_prob[t] = inv_logit(expected_vl_sim[t] - detection_threshold);
  }
}