data {
  int<lower=0> nf;
  int<lower=0> no;
  int<lower=0> n1;
  vector[nf] partners_f;
  vector[no] partners_o;
  vector[n1] partners_f1;
  vector[n1] partners_o1;
  matrix[no,2] D_o;
  matrix[n1,2] D_o1;
}
parameters {
  real<lower=0> P_mu_f;
  real<lower=0> P_sigma_f;
  vector[2] P_mu_o;
  real<lower=0> P_sigma_o;
}
model {
  P_mu_f ~ exponential(.01);
  P_sigma_f ~ exponential(.01);
  P_mu_o ~ normal(0, 50);
  P_sigma_o ~ exponential(.01);
  partners_f ~ normal(P_mu_f, P_sigma_f);
  partners_o ~ normal(exp(D_o*P_mu_o), P_sigma_o);
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  Dev_f = (partners_f1 - P_mu_f) / P_sigma_f;
  Dev_o = (partners_o1 - exp(D_o1*P_mu_o)) / P_sigma_o;
}
