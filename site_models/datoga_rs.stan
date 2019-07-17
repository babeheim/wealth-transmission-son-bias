data {
  int<lower=0> nf;
  int<lower=0> no;
  int<lower=0> n1;
  vector[nf] rs_f;
  vector[no] rs_o;
  vector[n1] rs_f1;
  vector[n1] rs_o1;
  matrix[nf,2] D_f;
  matrix[no,4] D_o;
  matrix[n1,2] D_f1;
  matrix[n1,4] D_o1;
  vector[no] male_o;
  vector[n1] male_o1;
}
parameters {
  vector[2] P_mu_p;
  vector[2] P_sigma_p;
  vector[4] P_mu_o;
  vector[4] P_sigma_o;
  real<lower=0> max_rs_s;
  real<lower=0> max_rs_d;
}
model {
  P_mu_p ~ normal(0, 50);
  P_sigma_p ~ normal(0, 50);
  P_mu_o ~ normal(0, 50);
  P_sigma_o ~ normal(0, 50);
  max_rs_s ~ exponential(.05);
  max_rs_d ~ exponential(.05);
  rs_f ~ normal(exp(D_f * P_mu_p), exp(D_f * P_sigma_p));
  rs_o ~ normal((male_o*max_rs_s + (1-male_o)*max_rs_d) ./ (1 + exp(D_o*P_mu_o)), exp(D_o * P_sigma_o));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  Dev_f = (rs_f1 - exp(D_f1 * P_mu_p)) ./ exp(D_f1 * P_sigma_p);
  Dev_o = (rs_o1 - ((male_o1*max_rs_s + (1-male_o1)*max_rs_d) ./ (1 + exp(D_o1*P_mu_o)))) ./ (exp(D_o1 * P_sigma_o));
}
