data {
  int<lower=0> n;
  int<lower=0> np;
  vector[n] wealth;
  vector[np] wealth_o;
  vector[np] wealth_p;
  matrix[n,3] D_mu_1;
  matrix[np,3] D_mu_o;
  matrix[np,3] D_mu_f;
  matrix[np,3] D_mu_m;
  vector[np] male;
}
parameters {
  vector[3] P_mu;
  real<lower=0> P_sigma;
}
model {
  P_mu ~ normal(0, 10);
  P_sigma ~ exponential(0.01);
  wealth ~ normal(D_mu_1 * P_mu, P_sigma);
}
generated quantities {
  vector[np] Dev_o;
  vector[np] Dev_f;
  vector[np] Dev_m;
  Dev_o = (wealth_o - (D_mu_o * P_mu)) / P_sigma;
  Dev_f = (wealth_p - (D_mu_f * P_mu)) / P_sigma;
  Dev_m = (wealth_p - (D_mu_m * P_mu)) / P_sigma;
}
