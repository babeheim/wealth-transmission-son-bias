
data {
  int<lower=0> n;
  int<lower=0> no;
  vector[n] inc;
  vector[no] inc_o;
  vector[no] inc_f;
  matrix[n, 4] D_mu;
  matrix[no, 4] D_mu_o;
  matrix[no, 4] D_mu_f;
}
parameters {
  vector[4] P_mu;
  real<lower=0> sigma;
}
model {
  inc ~ normal(D_mu * P_mu, sigma);
}
generated quantities {
  vector[no] Dev_o;
  vector[no] Dev_f;
  Dev_f = (inc_f - (D_mu_f * P_mu)) / sigma;
  Dev_o = (inc_o - (D_mu_o * P_mu)) / sigma;
}
