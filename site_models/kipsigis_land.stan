data {
  int<lower=0> n1;
  int<lower=0> n2;
  vector[n1] land;
  vector[n2] land_o;
  vector[n2] land_f;
  matrix[n1,5] D_mu;
  matrix[n1,2] D_sig;
  matrix[n2,5] D_mu_f;
  matrix[n2,2] D_sig_f;
  matrix[n2,5] D_mu_o;
  matrix[n2,2] D_sig_o;
  vector[n2] male;
}
parameters {
  vector[5] P_mu;
  vector[2] P_sigma;
}
model {
  vector[n1] E_sig;
  P_mu ~ normal(0, 10);
  P_sigma ~ normal(0, 10);
  E_sig = exp(D_sig * P_sigma);  //exponentiate to make it positive
  land ~ normal(D_mu * P_mu, E_sig);
}
generated quantities {
  vector[n2] Dev_o;
  vector[n2] Dev_f;
  Dev_o = (land_o - (D_mu_o * P_mu)) ./ exp(D_sig_o * P_sigma);
  Dev_f = (land_f - (D_mu_f * P_mu)) ./ exp(D_sig_f * P_sigma);
}
