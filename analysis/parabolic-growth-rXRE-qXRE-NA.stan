# Parabolic growth model with random effects and covariates
# and partially unobserved initial length
# Lt = (L0 + r*dt)^q

data {
  int<lower=1> N;           # number of observations
  int<lower=1> group[N];    # group IDs for random effects
  int<lower=1> K;           # total number of covariates (incl. intercept)
  matrix[N,K] X;            # covariates (first column is 1 for intercept)
  vector<lower=0>[N] L0;    # initial length
  vector<lower=0>[N] Lt;    # final length
  vector<lower=0>[N] dt;    # elapsed time
}

transformed data {
  int<lower=0> J;           # number of groups
  # vector[N] log_L0;         # log(L0)
  
  J = max(group);
  # log_L0 = log(L0);
}

parameters {
  vector[K] beta_log_r;         # regression coefs of log growth rate
  real<lower=0> sigma_log_r;    # hyper-SD of log growth rate
  vector[J] log_r_z;            # group-specific random effects on log r (z-scores)
  vector[K] beta_log_q;         # regression coefs of log growth curve shape
  real<lower=0> sigma_log_q;    # hyper-SD of log growth curve shape
  vector[J] log_q_z;            # group-specific random effects on log q (z-scores)
  real<lower=0> sigma;          # residual SD of log final length
}

model {
  vector[N] r;      # growth rate
  vector[N] log_q;  # log growth curve shape
  vector[N] q;      # growth curve shape
  vector[N] qinv;   # 1/q
  vector[N] L0qinv; # L0^(1/q)
  
  # Priors
  beta_log_r ~ normal(0,5);
  sigma_log_r ~ normal(0,5);
  log_r_z ~ normal(0,1);  # implies log(r) ~ normal(mu_log_r, sigma_log_r)
  beta_log_q ~ normal(0,5);
  sigma_log_q ~ normal(0,5);
  log_q_z ~ normal(0,1);  # implies log(q) ~ normal(mu_log_q, sigma_log_q)
  sigma ~ normal(0,2);
  
  # Likelihood
  r = exp(X * beta_log_r + sigma_log_r * log_r_z[group]);
  log_q = X * beta_log_q + sigma_log_q * log_q_z[group];
  q = exp(log_q);
  qinv = exp(-log_q);
  # L0qinv = exp(qinv .* log_L0);
  for(i in 1:N)
    L0qinv[i] = L0[i]^qinv[i];
  
  Lt ~ lognormal(q .* log(L0qinv + r .* dt), sigma);
}
