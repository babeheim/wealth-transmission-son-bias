data {
  int<lower=0> n;
  int<lower=0> n1;
  vector[n] ed;
  vector[n1] ed_o;
  vector[n1] ed_f;
  vector[n1] ed_m;
  matrix[n,6] D1;
  matrix[n1,6] D_o;
  matrix[n1,6] D_f;
  matrix[n1,6] D_m;
}
parameters {
  vector[6] P_mu;
  vector[6] P_sigma;
}
model {
  P_mu ~ normal(0, 20);
  P_sigma ~ normal(0, 20);
  ed ~ normal(D1 * P_mu, exp(D1 * P_sigma));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (ed_f - D_f * P_mu) ./ exp(D_f * P_sigma);
  Dev_m = (ed_m - D_m * P_mu) ./ exp(D_m * P_sigma);
  Dev_o = (ed_o - D_o * P_mu) ./ exp(D_o * P_sigma);
}
