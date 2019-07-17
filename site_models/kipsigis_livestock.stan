data {
  int<lower=0> n1;
  int<lower=0> n2;
  vector[n1] livestock;
  vector[n2] livestock_o;
  vector[n2] livestock_f;
  matrix[n1,5] D_mu;
  matrix[n1,3] D_sig;
  matrix[n2,5] D_mu_f;
  matrix[n2,3] D_sig_f;
  matrix[n2,5] D_mu_o;
  matrix[n2,3] D_sig_o;
  vector[n2] male;
}
parameters {
  vector[5] P_mu;
  vector[3] P_sigma;
}
model {
  vector[n1] E_sig;
  P_mu ~ normal(0, 100);
  P_sigma ~ normal(0, 100);
  E_sig = exp(D_sig * P_sigma);  //exponentiate to make it positive
  livestock ~ normal(D_mu * P_mu, E_sig);
}
generated quantities {
  vector[n2] Dev_o;
  vector[n2] Dev_f;
  Dev_o = (livestock_o - (D_mu_o * P_mu)) ./ exp(D_sig_o * P_sigma);
  Dev_f = (livestock_f - (D_mu_f * P_mu)) ./ exp(D_sig_f * P_sigma);
}
