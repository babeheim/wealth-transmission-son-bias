data {
  int<lower=0> nf;
  int<lower=0> no;
  int<lower=0> n1;
  vector[nf] livestock_f;
  vector[no] livestock_o;
  vector[n1] livestock_f1;
  vector[n1] livestock_o1;
  matrix[nf,3] D_f;
  matrix[no,3] D_o;
  matrix[n1,3] D_f1;
  matrix[n1,3] D_o1;
}
parameters {
  vector[3] P_mu_p;
  real<lower=0> P_sigma_p;
  vector[3] P_mu_o;
  real<lower=0> P_sigma_o;
}
model {
  P_mu_p ~ normal(0, 50);
  P_sigma_p ~ exponential(.05);
  P_mu_o ~ normal(0, 50);
  P_sigma_o ~ exponential(.05);
  livestock_f ~ normal(D_f*P_mu_p, P_sigma_p);
  livestock_o ~ normal(D_o*P_mu_o, P_sigma_o);
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  Dev_f = (livestock_f1 - D_f1*P_mu_p) / P_sigma_p;
  Dev_o = (livestock_o1 - D_o1*P_mu_o) / P_sigma_o;
}
