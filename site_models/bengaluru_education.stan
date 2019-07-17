
data {
  int<lower=0> no;
  int<lower=0> np;
  int<lower=0> n1;
  vector[no] ed_o;
  vector[np] ed_p;
  vector[n1] ed_o1;
  vector[n1] ed_f1;
  vector[n1] ed_m1;
  matrix[no,5] D_mu_o;
  matrix[np,2] D_mu_p;
  matrix[n1,5] D_mu_o1;
  matrix[n1,2] D_mu_f1;
  matrix[n1,2] D_mu_m1;
}
parameters {
  vector[5] P_mu_o;
  real<lower=0> sigma_o;
  vector[2] P_mu_p;
  real<lower=0> sigma_p;
}
model {
  ed_o ~ normal(D_mu_o * P_mu_o, sigma_o);
  ed_p ~ normal(D_mu_p * P_mu_p, sigma_p);
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (ed_f1 - (D_mu_f1 * P_mu_p)) / sigma_p;
  Dev_m = (ed_m1 - (D_mu_m1 * P_mu_p)) / sigma_p;
  Dev_o = (ed_o1 - (D_mu_o1 * P_mu_o)) / sigma_o;
}
