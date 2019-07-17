data {
  int<lower=0> no;
  int<lower=0> np;
  int<lower=0> n;
  vector[np] rs_p;
  vector[no] rs_o;
  vector[n] rs_f;
  vector[n] rs_m;
  vector[n] rs_o2;
  matrix[np,2] D_p;
  matrix[no,4] D_mu_o;
  matrix[no,6] D_sig_o;
  vector[no] male_o;
  matrix[n,2] D_f;
  matrix[n,2] D_m;
  matrix[n,4] D_mu_o2;
  matrix[n,6] D_sig_o2;
  vector[n] male_o2;
}
parameters {
  vector[2] P_mu_p;
  real<lower=0> P_sigma_p;
  vector[4] P_mu_o;
  vector[6] P_sigma_o;
  real<lower=0> max_rs;
  real<lower=0> max_rs_m;
}
model {
  vector[no] E_sig_o;
  P_mu_p ~ normal(0, 20);
  P_sigma_p ~ exponential(.01);
  P_mu_o ~ normal(0, 20);
  P_sigma_o ~ normal(0, 20);
  max_rs ~ exponential(0.1);
  max_rs_m ~ exponential(0.1);
  rs_p ~ normal(D_p * P_mu_p, P_sigma_p);
  E_sig_o = exp(D_sig_o * P_sigma_o);
  rs_o ~ normal(((1-male_o)*max_rs + male_o*max_rs_m) ./ (1 + exp(-D_mu_o * P_mu_o)), E_sig_o );
}
generated quantities {
  vector[n] Dev_o;
  vector[n] Dev_f;
  vector[n] Dev_m;
  Dev_f = (rs_f - (D_f * P_mu_p)) / P_sigma_p;
  Dev_m = (rs_m - (D_m * P_mu_p)) / P_sigma_p;
  Dev_o = (rs_o2 - (((1-male_o2)*max_rs + male_o2*max_rs_m) ./ (1 + exp(-D_mu_o2 * P_mu_o)))) ./ exp(D_sig_o2 * P_sigma_o);
}
