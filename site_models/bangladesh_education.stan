data {
  int<lower=0> n;
  vector[n] ed;
  vector[n] male;
  matrix[n,4] D_mu;
  matrix[n,6] D_sig;
  int<lower=0> n1;
  vector[n1] ed_o;
  vector[n1] ed_f;
  vector[n1] ed_m;
  vector[n1] male_o;
  vector[n1] ones1;
  matrix[n1,4] D_mu_o;
  matrix[n1,6] D_sig_o;
  matrix[n1,4] D_mu_f;
  matrix[n1,6] D_sig_f;
  matrix[n1,4] D_mu_m;
  matrix[n1,6] D_sig_m;
}
parameters {
  vector[4] P_mu;
  vector[6] P_sig;
  real<lower=0> high;
  real<lower=0> low;
  real<lower=0> high_m;
  real<lower=0> low_m;
}
model {
  ed ~ normal((1-male)*low + male*low_m + ((1-male)*high + male*high_m) ./ (1+exp(D_mu * P_mu)), exp(D_sig * P_sig));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (ed_f - (low_m + ones1*high_m ./ (1+exp(D_mu_f * P_mu)))) ./ exp(D_sig_f * P_sig);
  Dev_m = (ed_m - (low + ones1*high ./ (1+exp(D_mu_m * P_mu)))) ./ exp(D_sig_m * P_sig);
  Dev_o = (ed_o - ((1-male_o)*low + male_o*low_m + ((1-male_o)*high + male_o*high_m) ./ (1+exp(D_mu_o * P_mu)))) ./ exp(D_sig_o * P_sig);
}
