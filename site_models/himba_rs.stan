data {
  int<lower=0> nm;
  int<lower=0> no;
  int<lower=0> n1;
  vector[nm] rs_m;
  vector[no] rs_o;
  vector[n1] rs_m1;
  vector[n1] rs_o1;
  matrix[nm,2] D_m;
  matrix[no,6] D_o;
  matrix[n1,2] D_m1;
  matrix[n1,6] D_o1;
}
parameters {
  vector[2] P_mu_p;
  real<lower=0> P_sigma_p;
  vector[6] P_mu_o;
  vector[6] P_sigma_o;
}
model {
  P_mu_p ~ normal(0, 20);
  P_sigma_p ~ exponential(0.01);
  P_mu_o ~ normal(0, 20);
  P_sigma_o ~ normal(0, 20);
  rs_m ~ normal(D_m*P_mu_p, P_sigma_p);
  rs_o ~ normal(exp(D_o*P_mu_o), exp(D_o*P_sigma_o));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_m;
  Dev_m = (rs_m1 - D_m1*P_mu_p) / P_sigma_p;
  Dev_o = (rs_o1 - exp(D_o1*P_mu_o)) ./ exp(D_o1*P_sigma_o);
}
