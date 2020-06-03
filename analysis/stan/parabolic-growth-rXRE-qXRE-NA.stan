// Parabolic growth model with random effects and covariates
// and partially unobserved initial length
// log(Lt) ~ N(log((L0 + r*dt)^q), sigma)

data {
  int<lower=1> N;           // number of observations
  int<lower=1> group[N];    // group IDs for random effects
  int<lower=1> K;           // total number of covariates (incl. intercept)
  matrix[N,K] X;            // covariates (first column is 1 for intercept)
  vector<lower=0>[N] L0;    // initial length
  vector<lower=0>[N] Lt;    // final length
  vector<lower=0>[N] dt;    // elapsed time
}

transformed data {
  int<lower=0> J;           // number of groups

  J = max(group);
}

parameters {
  vector[K] beta_log_r;         // regression coefs of log growth rate
  real<lower=0> sigma_log_r;    // hyper-SD of log growth rate
  vector[J] log_r_z;            // group-specific random effects on log r (z-scores)
  vector[K] beta_log_q;         // regression coefs of log growth curve shape
  real<lower=0> sigma_log_q;    // hyper-SD of log growth curve shape
  vector[J] log_q_z;            // group-specific random effects on log q (z-scores)
  real<lower=0> sigma;          // residual SD of log final length
}

model {
  vector[N] r;      // growth rate
  vector[N] q;      // growth curve shape

  // Priors
  beta_log_r ~ normal(0,5);
  sigma_log_r ~ normal(0,5);
  log_r_z ~ std_normal();  // implies log(r) ~ normal(mu_log_r, sigma_log_r)
  beta_log_q ~ normal(0,5);
  sigma_log_q ~ normal(0,5);
  log_q_z ~ std_normal();  // implies log(q) ~ normal(mu_log_q, sigma_log_q)
  sigma ~ normal(0,2);
  
  // Likelihood
  r = exp(X * beta_log_r + sigma_log_r * log_r_z[group]);
  q = exp(X * beta_log_q + sigma_log_q * log_q_z[group]);
  Lt ~ lognormal(q .* log(L0^(1/q) + r .* dt), sigma);
}
