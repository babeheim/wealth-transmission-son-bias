data {
  int<lower=0> n_o;
  int<lower=0> n_p;
  int<lower=0> n_1;
  vector[n_o] RS_o;
  vector[n_p] RS_p;
  vector[n_1] RS_o1;
  vector[n_1] RS_f1;
  vector[n_1] RS_m1;
  matrix[n_o,3] D_mu_o;
  matrix[n_p,3] D_mu_p;
  matrix[n_1,3] D_mu_o1;
  matrix[n_1,3] D_mu_f1;
  matrix[n_1,3] D_mu_m1;
}
parameters {
  vector[3] P_mu_o;
  vector[3] P_mu_p;
  real<lower=0> sigma_o;
  real<lower=0> sigma_p;
}
model {
  RS_o ~ normal(D_mu_o * P_mu_o, sigma_o);
  RS_p ~ normal(D_mu_p * P_mu_p, sigma_p);
}
generated quantities {
  vector[n_1] Dev_o;
  vector[n_1] Dev_f;
  vector[n_1] Dev_m;
  Dev_f = (RS_f1 - (D_mu_f1 * P_mu_p)) / sigma_p;
  Dev_m = (RS_m1 - (D_mu_m1 * P_mu_p)) / sigma_p;
  Dev_o = (RS_o1 - (D_mu_o1 * P_mu_o)) / sigma_o;
}
