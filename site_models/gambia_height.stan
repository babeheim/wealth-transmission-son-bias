data {
  int<lower=0> n;
  int<lower=0> no;
  vector[n] height;
  vector[no] height_o;
  vector[no] height_f;
  vector[no] height_m;
  matrix[n,4] D_mu;
  matrix[n,2] D_sig;
  matrix[no,4] D_mu_o;
  matrix[no,2] D_sig_o;
  matrix[no,4] D_mu_f;
  matrix[no,2] D_sig_f;
  matrix[no,4] D_mu_m;
  matrix[no,2] D_sig_m;
}
parameters {
  vector[4] P_mu;
  vector[2] P_sigma;
}
model {
  height ~ normal(D_mu * P_mu, exp(D_sig * P_sigma));
}
generated quantities {
  vector[no] Dev_o;
  vector[no] Dev_f;
  vector[no] Dev_m;
  Dev_f = (height_f - D_mu_f * P_mu) ./  exp(D_sig_f * P_sigma);
  Dev_m = (height_m - D_mu_m * P_mu) ./  exp(D_sig_m * P_sigma);
  Dev_o = (height_o - D_mu_o * P_mu) ./  exp(D_sig_o * P_sigma);
}
