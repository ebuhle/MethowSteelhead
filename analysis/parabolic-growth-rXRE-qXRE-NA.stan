# Parabolic growth model with random effects and covariates
# and partially unobserved initial length

data {
  int<lower=1> N;           # number of observations
  int<lower=1> group[N];    # group IDs for random effects
  int<lower=1> K;           # total number of covariates (incl. intercept)
  matrix[N,K] X;            # covariates (first column is 1 for intercept)
  vector<lower=0>[N] L0     # initial length
  vector<lower=0>[N] Lt     # final length
  vector<lower=0>[N] dt     # elapsed time
}

transformed data {
  int<lower=0> J;           # number of groups
  
  J = max(group);
}

parameters {
  real beta_r;              # hyper-mean growth rate
  real<lower=0> sigma_r;    # hyper-SD growth rate
  vector[J] log_r_z;        # group-specific random effects on r (z-scores)
  real<lower=0> beta_q;     # growth curve shape
  
}

transformed parameters {
  
}

model {
  
}