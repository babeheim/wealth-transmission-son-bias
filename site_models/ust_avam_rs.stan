data {
  int<lower=0> np;
  int<lower=0> no;
  int<lower=0> n1;
  vector[np] rs_p;
  vector[no] rs_o;
  vector[n1] rs_o1;
  vector[n1] rs_f1;
  vector[n1] rs_m1;
  matrix[np,2] D_p;
  matrix[no,4] D_mu_o;
  matrix[no,2] D_sig_o;
  matrix[n1,2] D_f1;
  matrix[n1,2] D_m1;
  matrix[n1,4] D_mu_o1;
  matrix[n1,2] D_sig_o1;
  vector[no] male_o;
  vector[n1] male_o1;
  vector[np] ones_p;
  vector[no] ones_o;
  vector[n1] ones_1;
}
parameters {
  vector[2] P_mu_p;
  vector[2] P_sigma_p;
  vector[4] P_mu_o;
  vector[2] P_sigma_o;
  real<lower=0> max_rs_p;
  real<lower=0> max_sig_p;
  real<lower=0> max_rsm_o;
  real<lower=0> max_rsf_o;
  real<lower=0> max_sig_o;
}
model {
  P_mu_p ~ normal(0, 20);
  P_sigma_p ~ normal(0, 20);
  P_mu_o ~ normal(0, 20);
  P_sigma_o ~ normal(0, 20);
  max_rs_p ~ exponential(0.1);
  max_sig_p ~ exponential(0.1);
  max_rsm_o ~ exponential(0.1);
  max_rsf_o ~ exponential(0.1);
  max_sig_o ~ exponential(0.1);
  rs_p ~ normal(
    (max_rs_p * ones_p) ./ (1 + exp(D_p * P_mu_p)),
    (max_sig_p * ones_p) ./ (1 + exp(D_p * P_sigma_p)));
  rs_o ~ normal(
    (male_o * max_rsm_o + (1 - male_o) * max_rsf_o) ./ (1 + exp(D_mu_o * P_mu_o)),
    (max_sig_o * ones_o) ./ (1 + exp(D_sig_o * P_sigma_o)));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (rs_f1 - ((max_rs_p * ones_1) ./ (1 + exp(D_f1 *
    P_mu_p)))) ./ ((max_sig_p * ones_1) ./ (1 + exp(D_f1 * P_sigma_p)));
  Dev_m = (rs_m1 - ((max_rs_p * ones_1) ./ (1 + exp(D_m1 *
    P_mu_p)))) ./ ((max_sig_p * ones_1) ./ (1 + exp(D_m1 * P_sigma_p)));
  Dev_o = (rs_o1 - ((male_o1 * max_rsm_o + (1-male_o1) *
    max_rsf_o) ./ (1 + exp(D_mu_o1 * P_mu_o)))) ./ ((max_sig_o * ones_1) ./ (1 + exp(D_sig_o1 * P_sigma_o)));
}