data {
  int<lower=0> n;
  int<lower=0> np;
  vector[n] ed;
  vector[np] ed_o;
  vector[np] ed_f;
  vector[np] ed_m;
  matrix[n,5] D;
  matrix[np,5] D_o;
  matrix[np,5] D_f;
  matrix[np,5] D_m;
  vector[np] male;
}
parameters {
  vector[5] P_mu;
  vector[5] P_sigma;
}
model {
  vector[n] E_sig;
  P_mu ~ normal(0, 10);
  P_sigma ~ normal(0, 10);
  E_sig = exp(D * P_sigma);
  ed ~ normal(exp(D * P_mu), E_sig);
}
generated quantities {
  vector[np] Dev_o;
  vector[np] Dev_f;
  vector[np] Dev_m;
  Dev_o = (ed_o - exp(D_o * P_mu)) ./ exp(D_o * P_sigma);
  Dev_f = (ed_f - exp(D_f * P_mu)) ./ exp(D_f * P_sigma);
  Dev_m = (ed_m - exp(D_m * P_mu)) ./ exp(D_m * P_sigma);
}
