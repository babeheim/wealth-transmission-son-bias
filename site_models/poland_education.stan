data {
  int<lower=0> n;
  vector[n] ed;
  matrix[n,6] D_mu;
  matrix[n,2] D_sig;
}
parameters {
  vector[6] P_mu;
  vector[2] P_sigma;
}
model {
  P_mu ~ normal(0, 5);
  P_sigma ~ normal(0, 5);
  ed ~ normal(exp(D_mu * P_mu), exp(D_sig * P_sigma));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (weight_f - D_mu_f * P_mu) ./ exp(D_sig_f * P_sigma);
  Dev_m = (weight_m - D_mu_m * P_mu) ./ exp(D_sig_m * P_sigma);
  Dev_o = (weight_o - D_mu_o * P_mu) ./ exp(D_sig_o * P_sigma);
}
