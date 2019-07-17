data {
  int<lower=0> np;
  int<lower=0> n;
  int<lower=0> n1;
  vector[np] land_f;
  vector[n] land;
  vector[n1] land_o1;
  vector[n1] land_f1;
  matrix[np,2] D_mu_f;
  matrix[np,3] D_sig_f;
  matrix[n1,2] D_mu_f1;
  matrix[n1,3] D_sig_f1;
  matrix[n,3] D;
  matrix[n1,3] D1;
  vector[n1] male;
}
parameters {
  vector[2] P_mu_f;
  vector[3] P_sigma_f;
  vector[3] P_mu;
  vector[3] P_sigma;
}
model {
  vector[n] E_sig;
  vector[np] E_sig_f;
  P_mu_f ~ normal(0, 10);
  P_sigma_f ~ normal(0, 10);
  P_mu ~ normal(0, 10);
  P_sigma ~ normal(0, 10);
  E_sig_f = exp(D_sig_f * P_sigma_f);
  land_f ~ normal(D_mu_f * P_mu_f, E_sig_f);
  E_sig = exp(D * P_sigma);
  land ~ normal(D * P_mu, E_sig);
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  Dev_o = (land_o1 - (D1 * P_mu)) ./ exp(D1 * P_sigma);
  Dev_f = (land_f1 - (D_mu_f1 * P_mu_f)) ./ exp(D_sig_f1 * P_sigma_f);
}
