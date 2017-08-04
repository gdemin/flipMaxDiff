data {
  int<lower=2> C; // Number of alternatives (choices) in each scenario
  int<lower=1> K; // Number of covariates of alternatives
  int<lower=1> R; // Number of respondents
  int<lower=1> S; // Number of scenarios per respondent
  int<lower=0> G; // Number of respondent covariates 
  int<lower=1,upper=C> YB[R, S]; // best choices
  int<lower=1,upper=C> YW[R, S]; // worst choices
  matrix[C, K] X[R, S]; // matrix of attributes for each obs
  matrix[G, R] Z; // vector of covariates for each respondent
}

parameters {
  matrix[K, R] Beta;
  matrix[K, G] Theta;
  corr_matrix[K] Omega;
  vector<lower=0>[K] tau;
}
transformed parameters {
  cov_matrix[K] Sigma = quad_form_diag(Omega, tau);
}
model {
  //priors
  to_vector(Theta) ~ normal(0, 10);
  tau ~ cauchy(0, 2.5); 
  Omega ~ lkj_corr(2);
  //likelihood
  for (r in 1:R) {
    Beta[,r] ~ multi_normal(Theta*Z[,r], Sigma);
    for (s in 1:S) {
      YB[r,s] ~ categorical_logit(X[r,s] * Beta[,r]);
      YW[r,s] ~ categorical_logit(-X[r,s] * Beta[,r]);
    }
  }
}

