data {
  int n1;
  int n2;
  matrix[n1,5] D_mu_1;
  matrix[n1,5] D_sig_1;
  matrix[n2,5] D_mu_o;
  matrix[n2,5] D_sig_o;
  matrix[n2,5] D_mu_f;
  matrix[n2,5] D_sig_f;
  matrix[n2,5] D_mu_m;
  matrix[n2,5] D_sig_m;
  vector[n1] ed;
  vector[n2] ed_o;
  vector[n2] ed_f;
  vector[n2] ed_m;
  vector[n2] male;
}
parameters {
  vector[5] P_mu;
  vector[5] P_sigma;
}
model {
  P_mu ~ normal(0, 10);
  P_sigma ~ normal(0, 10);
  ed ~ normal(exp(D_mu_1 * P_mu), exp(D_sig_1 * P_sigma));
}
generated quantities {
  vector[n2] Dev_o;
  vector[n2] Dev_f;
  vector[n2] Dev_m;
  Dev_o = (ed_o - exp(D_mu_o * P_mu)) ./ exp(D_sig_o * P_sigma);
  Dev_f = (ed_f - exp(D_mu_f * P_mu)) ./ exp(D_sig_f * P_sigma);
  Dev_m = (ed_m - exp(D_mu_m * P_mu)) ./ exp(D_sig_m * P_sigma);
}
