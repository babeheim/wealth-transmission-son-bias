data {
  int<lower=0> no;
  int<lower=0> np;
  vector[no] RS_o;
  vector[no] RS_f;
  vector[no] RS_m;
  vector[np] RS_p;
  matrix[no,2] D_o;
}
parameters {
  vector[2] P_mu_o;
  vector[2] P_sig_o;
  real<lower=0> max_RS_o;
  real<lower=0> max_sig_o;
  real mu_p;
  real<lower=0> sigma_p;
}
model {
  RS_o ~ normal(max_RS_o*exp(D_o * P_mu_o) ./ (1 + exp(D_o * P_mu_o)), max_sig_o*exp(D_o * P_sig_o) ./ (1 + exp(D_o * P_sig_o)));
  RS_p ~ normal(mu_p, sigma_p);
}
generated quantities {
  vector[no] Dev_o;
  vector[no] Dev_f;
  vector[no] Dev_m;
  Dev_f = (RS_f - mu_p) / sigma_p;
  Dev_m = (RS_m - mu_p) / sigma_p;
  Dev_o = (RS_o - (max_RS_o * exp(D_o * P_mu_o) ./ (1 + exp(D_o * P_mu_o)))) ./ (max_sig_o*exp(D_o * P_sig_o) ./ (1 + exp(D_o * P_sig_o)));
}
