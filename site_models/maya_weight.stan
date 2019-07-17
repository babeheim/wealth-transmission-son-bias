data {
  int<lower=0> n;
  int<lower=0> no;
  vector[n] weight;
  vector[no] weight_o;
  vector[no] weight_f;
  vector[no] weight_m;
  matrix[n,6] D_mu;
  matrix[n,2] D_sig;
  matrix[no,6] D_mu_o;
  matrix[no,2] D_sig_o;
  matrix[no,6] D_mu_f;
  matrix[no,2] D_sig_f;
  matrix[no,6] D_mu_m;
  matrix[no,2] D_sig_m;
}
parameters {
  vector[6] P_mu;
  vector[2] P_sigma;
}
model {
  P_mu ~ normal(0, 20);
  P_sigma ~ normal(0, 20);
  weight ~ normal(D_mu * P_mu, exp(D_sig * P_sigma));
}
generated quantities {
  vector[no] Dev_o;
  vector[no] Dev_f;
  vector[no] Dev_m;
  Dev_f = (weight_f - D_mu_f * P_mu) ./ exp(D_sig_f * P_sigma);
  Dev_m = (weight_m - D_mu_m * P_mu) ./ exp(D_sig_m * P_sigma);
  Dev_o = (weight_o - D_mu_o * P_mu) ./ exp(D_sig_o * P_sigma);
}
