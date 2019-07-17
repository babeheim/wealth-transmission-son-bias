data {
  int<lower=0> n;
  int<lower=0> no;
  vector[n] ed;
  vector[no] ed_o;
  vector[no] ed_f;
  vector[no] ed_m;
  matrix[n,5] D_mu;
  matrix[n,4] D_sig;
  matrix[no,5] D_mu_o;
  matrix[no,4] D_sig_o;
  matrix[no,5] D_mu_f;
  matrix[no,4] D_sig_f;
  matrix[no,5] D_mu_m;
  matrix[no,4] D_sig_m;
}
parameters {
  vector[5] P_mu;
  vector[4] P_sig;
}
model {
  ed ~ normal(exp(D_mu * P_mu), exp(D_sig * P_sig));
}
generated quantities {
  vector[no] Dev_o;
  vector[no] Dev_f;
  vector[no] Dev_m;
  Dev_f = (ed_f - exp(D_mu_f * P_mu)) ./ exp(D_sig_f * P_sig);
  Dev_m = (ed_m - exp(D_mu_m * P_mu)) ./ exp(D_sig_m * P_sig);
  Dev_o = (ed_o - exp(D_mu_o * P_mu)) ./ exp(D_sig_o * P_sig);
}
