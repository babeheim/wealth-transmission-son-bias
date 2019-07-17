data {
  int<lower=0> n;
  int<lower=0> n2;
  vector[n] shares;
  vector[n2] shares_o;
  vector[n2] shares_p;
  matrix[n,3] D_mu_1;
  matrix[n,2] D_sig_1;
  matrix[n2,3] D_mu_o;
  matrix[n2,2] D_sig_o;
  matrix[n2,3] D_mu_f;
  matrix[n2,2] D_sig_f;
  matrix[n2,3] D_mu_m;
  matrix[n2,2] D_sig_m;
  vector[n2] male;
}
parameters {
  vector[3] P_mu;
  vector[2] P_sigma;
}
model {
  vector[n] E_sig;
  P_mu ~ normal(0, 10);
  P_sigma ~ normal(0, 10);
  E_sig = exp(D_sig_1 * P_sigma);
  shares ~ normal(D_mu_1 * P_mu, E_sig);
}
generated quantities {
  vector[n2] Dev_o;
  vector[n2] Dev_f;
  vector[n2] Dev_m;
  Dev_o = (shares_o - (D_mu_o * P_mu)) ./ exp(D_sig_o * P_sigma);
  Dev_f = (shares_p - (D_mu_f * P_mu)) ./ exp(D_sig_f * P_sigma);
  Dev_m = (shares_p - (D_mu_m * P_mu)) ./ exp(D_sig_m * P_sigma);
}
