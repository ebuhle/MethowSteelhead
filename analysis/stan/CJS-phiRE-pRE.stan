// Cormack-Jolly-Seber Model (aggregated array data format) with random effects on time-specific 
// survival and capture probabilities (possibly different grouping variables for each phi and p)

functions {
  int first_capture(int[] y_i) {
    for (t in 1:size(y_i))
      if (y_i[t])
        return t;
    return 0;
  }
  
  int last_capture(int[] y_i) {
    for (t_rev in 0:(size(y_i) - 1)) 
    {
      int t;
      t = size(y_i) - t_rev;
      if (y_i[t])
        return t;
    }
    return 0;
  }
  
  row_vector prob_uncaptured(int T, row_vector p, row_vector phi) {
    row_vector[T] chi;
    
    chi[T] = 1.0;
    for (t in 1:(T - 1)) 
    {
      int t_curr;
      int t_next;
      t_curr = T - t;
      t_next = t_curr + 1;
      chi[t_curr] = (1 - phi[t_curr]) + phi[t_curr] * (1 - p[t_next]) * chi[t_next];
    }
    return chi;
  }
}

data {
  int<lower=2> T;                  // number of capture events (includes marking)
  int<lower=0> M;                  // number of unique capture histories
  int<lower=1> group_phi[M,T-1];   // phi group IDs for each unique capture history
  int<lower=1> group_p[M,T];       // p group IDs for each unique capture history
  int<lower=0,upper=1> y[M,T];     // y[m,t]: history m captured at t
  int<lower=1> n[M];               // n[m]: number of individuals with capture history y[m,]
}

transformed data {
  int<lower=1> J_phi;              // number of groups for phi
  int<lower=1> J_p;                // number of groups for p
  int<lower=0,upper=T> first[M];   // first capture occasion
  int<lower=0,upper=T> last[M];    // last capture occasion
  int<lower=0,upper=T-1> last_minus_first[M];  // duh
  
  J_phi = max(to_array_1d(group_phi));
  J_p = max(to_array_1d(group_p));
  
  for (m in 1:M)
  {
    first[m] = first_capture(y[m,]);
    last[m] = last_capture(y[m,]);
    last_minus_first[m] = last[m] - first[m];
  }
}

parameters {
  vector<lower=0,upper=1>[T-1] mu_phi;  // inverse logit of mean(logit(phi[,t]))
  vector<lower=0>[T-1] sigma_phi;       // among-group SDs of logit(phi[,t])
  matrix[J_phi,T-1] logit_phi_z;        // group-specific random effects on phi (z-scores)
  vector<lower=0,upper=1>[T] mu_p;      // inverse logit of mean(logit(p[,t]))
  vector<lower=0>[T] sigma_p;           // among-group SDs of logit(p[,t])
  matrix[J_p,T] logit_p_z;              // group-specific random effects on p (z-scores)
}

transformed parameters {
  matrix<lower=0,upper=1>[M,T-1] phi_m;  // phi[,t]: Pr[alive at t + 1 | alive at t]
  matrix<lower=0,upper=1>[M,T] p_m;      // p[,t]: Pr[captured at t] (note p[,1] not used in model)
  matrix<lower=0,upper=1>[M,T] chi;      // chi[,t]: Pr[not captured >  t | alive at t]
  
  for(t in 1:(T-1))
    phi_m[,t] = inv_logit(logit(mu_phi[t]) + sigma_phi[t] * logit_phi_z[group_phi[,t],t]);
  
  for(t in 1:T)
    p_m[,t] = inv_logit(logit(mu_p[t]) + sigma_p[t] * logit_p_z[group_p[,t],t]);

  for(m in 1:M)
    chi[m,] = prob_uncaptured(T, p_m[m,], phi_m[m,]);
}

model {
  // implied uniform priors:
  // mu_phi ~ uniform(0,1)
  // mu_p ~ uniform(0,1)
  
  sigma_phi ~ normal(0,5);    // weakly informative
  to_vector(logit_phi_z) ~ normal(0,1);  // implies logit(phi[j,t]) ~ N(logit(mu_phi[t]), sigma_phi);
  sigma_p ~ normal(0,5);      // weakly informative
  to_vector(logit_p_z) ~ normal(0,1);    // implies logit(p[j,t]) ~ N(logit(mu_p[t]), sigma_p);
  
  // Likelihood of capture history
  // marginalized over discrete latent states
  for (m in 1:M) 
  {
    if (last_minus_first[m] > 0)  // if history m was recaptured
    {
      for(t in (first[m]+1):last[m])
      {
        target += n[m] * log(phi_m[m,t-1]);                 // survival from t - 1 to t
        target += n[m] * bernoulli_lpmf(y[m,t] | p_m[m,t]); // observation (captured or not)
      }
    }
    target += n[m] * log(chi[m,last[m]]);   // Pr[not detected after last[m]]
  }
}

generated quantities {
  matrix<lower=0,upper=1>[J_phi,T-1] phi; 
  matrix<lower=0,upper=1>[J_p,T] p; 
  
  for(t in 1:(T-1))
    for(j in 1:J_phi)
    {
      if(j <= max(group_phi[,t]))
        phi[j,t] = inv_logit(logit(mu_phi[t]) + sigma_phi[t] * logit_phi_z[j,t]);
      else
        phi[j,t] = 0;
    }
  
  p[,1] = rep_vector(1,J_p);  // set initial capture prob to 1
  for(t in 2:T)
    for(j in 1:J_p)
    {
      if(j <= max(group_p[,t]))
        p[j,t] = inv_logit(logit(mu_p[t]) + sigma_p[t] * logit_p_z[j,t]);
      else
        p[j,t] = 0;
    }
}
