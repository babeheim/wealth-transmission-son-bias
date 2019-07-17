data {
  int<lower=0> np;
  int<lower=0> n1;
  int<lower=0> no;
  vector[np] rs_p;
  vector[n1] rs_1;
  vector[no] rs_f;
  vector[no] rs_m;
  vector[no] rs_o;
  matrix[np,2] D_mu_p;
  matrix[n1,3] D_1;
  matrix[no,2] D_mu_f;
  matrix[no,2] D_mu_m;
  matrix[no,3] D_o;
  vector[no] male;
}
parameters {
  vector[2] P_mu_p;
  real<lower=0> P_sigma_p;
  vector[3] P_mu_1;
  vector[3] P_sigma_1;
}
model {
  vector[n1] E_sigma_1;
  vector[no] Dev_o;
  vector[no] Dev_f;
  vector[no] Dev_m;
  P_mu_p ~ normal(0, 10);
  P_sigma_p ~ exponential(0.01);
  P_mu_1 ~ normal(0, 10);
  P_sigma_1 ~ normal(0, 10);
  rs_p ~ normal(D_mu_p * P_mu_p, P_sigma_p);
  E_sigma_1 = exp(D_1 * P_sigma_1);
  rs_1 ~ normal(D_1 * P_mu_1, E_sigma_1);
}
generated quantities {
  vector[no] Dev_o;
  vector[no] Dev_f;
  vector[no] Dev_m;
  Dev_o = (rs_o - (D_o * P_mu_1)) ./ exp(D_o * P_sigma_1);
  Dev_f = (rs_f - (D_mu_f * P_mu_p)) / P_sigma_p;
  Dev_m = (rs_m - (D_mu_m * P_mu_p)) / P_sigma_p;
}
