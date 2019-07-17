data {
  int<lower=0> n;
  int<lower=0> n1;
  vector[n] rs;
  vector[n1] rs1;
  vector[n1] rs1f;
  matrix[n,4] D_mu;
  matrix[n,3] D_sig;
  matrix[n1,4] D_1_mu;
  matrix[n1,3] D_1_sig;
  matrix[n1,4] D_1_mu_f;
  matrix[n1,3] D_1_sig_f;
}
parameters {
  vector[4] P_mu;
  vector[3] P_sigma;
}
model {
  P_mu ~ normal(0, 20);
  P_sigma ~ normal(0, 20);
  rs ~ normal(exp(D_mu * P_mu), exp(D_sig * P_sigma));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  Dev_f = (rs1f - exp(D_1_mu_f * P_mu)) ./ exp(D_1_sig_f * P_sigma);
  Dev_o = (rs1 - exp(D_1_mu * P_mu)) ./ exp(D_1_sig * P_sigma);
}