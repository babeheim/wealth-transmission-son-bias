data {
  int<lower=0> np;
  int<lower=0> no;
  int<lower=0> n1;
  vector[np] wealth_p;
  vector[no] wealth_o;
  vector[n1] wealth_o1;
  vector[n1] wealth_f1;
  vector[n1] wealth_m1;
  matrix[np,2] D_p;
  matrix[n1,2] D_f1;
  matrix[n1,2] D_m1;
}
parameters {
  vector[2] P_mu_p;
  real<lower=0> P_sigma_p;
  real P_mu_o;
  real<lower=0> P_sigma_o;
}
model {
  P_mu_p ~ normal(0, 50);
  P_sigma_p ~ exponential(0.01);
  P_mu_o ~ normal(0, 50);
  P_sigma_o ~ exponential(0.01);
  wealth_p ~ normal(D_p * P_mu_p, P_sigma_p);
  wealth_o ~ normal(P_mu_o, P_sigma_o);
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (wealth_f1 - D_f1 * P_mu_p) / P_sigma_p;
  Dev_m = (wealth_m1 - D_m1 * P_mu_p) / P_sigma_p;
  Dev_o = (wealth_o1 - P_mu_o) / P_sigma_o;
}
