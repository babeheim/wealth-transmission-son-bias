
data {
  int<lower=0> n;
  int<lower=0> n1;
  vector[n] land;
  vector[n1] land_o;
  vector[n1] land_f;
  matrix[n,3] D_mu;
  matrix[n1,3] D_mu_o;
  matrix[n1,3] D_mu_f;
}
parameters {
  vector[3] P_mu;
  real<lower=0> sigma;
}
model {
  land ~ normal(D_mu * P_mu, sigma);
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  Dev_f = (land_f - (D_mu_f * P_mu)) / sigma;
  Dev_o = (land_o - (D_mu_o * P_mu)) / sigma;
}
