//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/bym2.stan
//


data {
  int<lower=0> N;
  int<lower=0> y[N];                     // count outcomes
  vector<lower=0>[N] E;                  // exposure
  int<lower=1> K;                        // num covariates
  matrix[N, K] x;                        // design matrix

  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;            // intercept
  vector[K] betas;       // covariates
  real<lower=0> sigma;        // overall standard deviation
  vector[N] phi;         // spatial effects
}
transformed parameters {
// none
}
model {
  y ~ poisson_log(log_E + beta0 + x * betas + phi * sigma);  // co-variates
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  beta0 ~ normal(0.0, 1.0); // og = 1
  betas ~ normal(0.0, 1.0); // og = 1
  sigma ~ normal(0.0, 1.0); // og = 1
  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)
}
generated quantities {
  vector[N] eta = log_E + beta0 + x * betas + phi * sigma; // co-variates
  vector[N] mu = exp(eta);
  
  // http://mc-stan.org/loo/articles/loo2-with-rstan.html
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = poisson_log_lpmf(y[n] | eta[n]);
  }
}

