data {
  int<lower=0> n;
  int<lower=0> no;
  vector[n] height;
  vector[no] height_o;
  vector[no] height_f;
  vector[no] height_m;
  matrix[n,3] D_mu;
  matrix[no,3] D_mu_o;
  matrix[no,3] D_mu_f;
  matrix[no,3] D_mu_m;
}
parameters {
  vector[3] P_mu;
  real<lower=0> P_sigma;
}
model {
  P_mu ~ normal(0, 20);
  height ~ normal(D_mu * P_mu, P_sigma);
}
generated quantities {
  vector[no] Dev_o;
  vector[no] Dev_f;
  vector[no] Dev_m;
  Dev_f = (height_f - D_mu_f * P_mu) / P_sigma;
  Dev_m = (height_m - D_mu_m * P_mu) / P_sigma;
  Dev_o = (height_o - D_mu_o * P_mu) / P_sigma;
}
