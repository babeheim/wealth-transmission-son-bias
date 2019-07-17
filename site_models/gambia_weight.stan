data {
  int<lower=0> n1;
  int<lower=0> n2;
  vector[n1] weight_d1;
  matrix[n1,8] D_mu_all;
  matrix[n1,4] D_sig_all;
  vector[n2] weight;
  vector[n2] weight_f;
  vector[n2] weight_m;
  matrix[n2,8] D_mu;
  matrix[n2,8] D_mu_f;
  matrix[n2,8] D_mu_m;
  matrix[n2,4] D_sig;
  matrix[n2,4] D_sig_f;
  matrix[n2,4] D_sig_m;
  vector[n2] male;
}
parameters {
  vector[8] P_mu;
  vector[4] P_sigma;
}
model {
  vector[n1] E_sig;
  P_mu ~ normal(0, 100);
  P_sigma ~ normal(0, 100);
  E_sig = exp(D_sig_all * P_sigma);  //exponentiate to make it positive
  weight_d1 ~ normal(D_mu_all * P_mu, E_sig);
}
generated quantities {
  vector[n2] Dev_o;
  vector[n2] Dev_f;
  vector[n2] Dev_m;
  Dev_o = (weight - (D_mu * P_mu)) ./ exp(D_sig * P_sigma);
  Dev_f = (weight_f - (D_mu_f * P_mu)) ./ exp(D_sig_f * P_sigma);
  Dev_m = (weight_m - (D_mu_m * P_mu)) ./ exp(D_sig_m * P_sigma);
}
