data {
  int<lower=0> n;
  int<lower=0> n1;
  vector[n] hhgardsize;
  vector[n1] hhgardsize_o;
  vector[n1] hhgardsize_m;
  matrix[n,3] D_mu;
  matrix[n,2] D_sig;
  matrix[n1,3] D_mu_o;
  matrix[n1,2] D_sig_o;
  matrix[n1,3] D_mu_m;
  matrix[n1,2] D_sig_m;
}
parameters {
  vector[3] P_mu;
  vector[2] P_sigma;
}
model {
  P_mu ~ normal(0, 50);
  P_sigma ~ normal(0, 50);
  hhgardsize ~ normal(D_mu*P_mu, exp(D_sig*P_sigma));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_m;
  Dev_m = (hhgardsize_m - D_mu_m*P_mu) ./ exp(D_sig_m*P_sigma);
  Dev_o = (hhgardsize_o - D_mu_o*P_mu) ./ exp(D_sig_o*P_sigma);
}
