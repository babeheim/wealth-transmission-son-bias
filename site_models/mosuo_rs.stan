data {
  int<lower=0> no;
  int<lower=0> np;
  int<lower=0> n;
  vector[np] rs_p;
  vector[no] rs_o;
  vector[n] rs_f;
  vector[n] rs_m;
  vector[n] rs_o2;
  matrix[np,4] D_mu_p;
  matrix[np,2] D_sig_p;
  matrix[no,4] D_mu_o;
  matrix[no,3] D_sig_o;
  matrix[n,4] D_mu_f;
  matrix[n,2] D_sig_f;
  matrix[n,4] D_mu_m;
  matrix[n,2] D_sig_m;
  matrix[n,4] D_mu_o2;
  matrix[n,3] D_sig_o2;
  vector[no] male_o;
  vector[n] male_o2;
  vector[no] ones;
}
parameters {
  vector[4] P_mu_p;
  vector[2] P_sigma_p;
  vector[4] P_mu_o;
  vector[3] P_sigma_o;
  real<lower=0> max_rs;
  real<lower=0> max_rs_m;
}
model {
  vector[np] E_sig_p;
  vector[no] E_sig_o;
  P_mu_p ~ normal(0, 20);
  P_sigma_p ~ normal(0, 10);
  P_mu_o ~ normal(0, 20);
  P_sigma_o ~ normal(0, 20);
  max_rs ~ exponential(0.1);
  max_rs_m ~ exponential(0.1);
  E_sig_p = exp(D_sig_p * P_sigma_p);
  rs_p ~ normal(exp(D_mu_p * P_mu_p), E_sig_p);
  E_sig_o = exp(D_sig_o * P_sigma_o);
  rs_o ~ normal(((1-male_o)*max_rs + male_o*max_rs_m) .* (ones ./ (1 + exp(-D_mu_o * P_mu_o))), E_sig_o );
}
generated quantities {
  vector[n] Dev_o;
  vector[n] Dev_f;
  vector[n] Dev_m;
  Dev_f = (rs_f - exp(D_mu_f * P_mu_p)) ./ exp(D_sig_f * P_sigma_p);
  Dev_m = (rs_m - exp(D_mu_m * P_mu_p)) ./ exp(D_sig_m * P_sigma_p);
  Dev_o = (rs_o2 - (((1-male_o2)*max_rs + male_o2*max_rs_m) .* (ones ./ (1 + exp(-D_mu_o2 * P_mu_o))))) ./ exp(D_sig_o2 * P_sigma_o);
}
