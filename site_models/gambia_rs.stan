data {
  int<lower=0> n;
  int<lower=0> n1;
  vector[n] rs;
  vector[n1] rs_o;
  vector[n1] rs_f;
  vector[n1] rs_m;
  matrix[n,4] D;
  matrix[n1,4] D_o;
  matrix[n1,4] D_f;
  matrix[n1,4] D_m;
  vector[n] male;
  vector[n1] male_o;
  vector[n1] ones_1;
}
parameters {
  vector[4] P_mu;
  vector[4] P_sigma;
  real<lower=0> max_rs;
  real<lower=0> max_rs_m;
  real<lower=0> max_sig;
  real<lower=0> max_sig_m;
}
model {
  P_mu ~ normal(0, 50);
  P_sigma ~ normal(0, 50);
  max_rs ~ exponential(0.1);
  max_rs_m ~ exponential(0.1);
  max_sig ~ exponential(0.1);
  max_sig_m ~ exponential(0.1);
  rs ~ normal((male*max_rs_m + (1-male)*max_rs) ./ (1+ exp(D * P_mu)) , (male*max_sig_m + (1-male)*max_sig) ./ (1+ exp(D * P_sigma)));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (rs_f - ((ones_1*max_rs_m) ./ (1+ exp(D_f * P_mu)))) ./ ((ones_1*max_sig_m) ./ (1+ exp(D_f * P_sigma)));
  Dev_m = (rs_m - ((ones_1*max_rs) ./ (1+ exp(D_m * P_mu)))) ./ ((ones_1*max_sig) ./ (1+ exp(D_m * P_sigma)));
  Dev_o = (rs_o - ((male_o*max_rs_m + (1-male_o)*max_rs) ./ (1+ exp(D_o * P_mu)))) ./ ((male_o*max_sig_m + (1-male_o)*max_sig) ./ (1+ exp(D_o * P_sigma)));
}
