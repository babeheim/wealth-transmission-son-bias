data {
  int<lower=0> np;
  int<lower=0> no;
  int<lower=0> n1;
  vector[np] rsp;
  vector[no] rso;
  vector[n1] rso1;
  vector[n1] rsm;
  matrix[no,4] D_o;
  matrix[n1,4] D_o1;
  vector[no] male_o;
  vector[n1] male_o1;
}
parameters {
  real<lower=0> P_mu_p;
  real<lower=0> P_sigma_p;
  vector[4] P_mu_o;
  vector[4] P_sigma_o;
  real<lower=0> max_RS_m;
  real<lower=0> max_RS_f;
  real<lower=0> max_sig_m;
  real<lower=0> max_sig_f;
}
model {
  P_mu_p ~ exponential(0.01);
  P_sigma_p ~ exponential(0.01);
  P_mu_o ~ normal(0, 100);
  P_sigma_o ~ normal(0, 100);
  max_RS_m ~ exponential(0.01);
  max_RS_f ~ exponential(0.01);
  max_sig_m ~ exponential(0.01);
  max_sig_f ~ exponential(0.01);
  rsp ~ normal(P_mu_p, P_sigma_p);
  rso ~ normal(((1-male_o)*max_RS_f + male_o*max_RS_m) ./ (1 + exp(D_o * P_mu_o)), ((1-male_o)*max_sig_f + male_o*max_sig_m) ./ (1 + exp(D_o * P_sigma_o)));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_m;
  Dev_m = (rsm - P_mu_p) / P_sigma_p;
  Dev_o = (rso1 - (((1-male_o1)*max_RS_f + male_o1*max_RS_m) ./ (1 + exp(D_o1 * P_mu_o)))) ./ (((1-male_o1)*max_sig_f + male_o1*max_sig_m) ./ (1 + exp(D_o1 * P_sigma_o)));
}
