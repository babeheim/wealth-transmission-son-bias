data {
  int<lower=0> n;
  int<lower=0> n1;
  vector[n] skill;
  vector[n1] skill_f;
  vector[n1] skill_m;
  vector[n1] skill_o;
  matrix[n,2] D_all_sig;
  matrix[n1,2] D_f_sig;
  matrix[n1,2] D_m_sig;
  matrix[n1,2] D_o_sig;
}
parameters {
  real P_mu;
  vector[2] P_sigma;
}
model {
  vector[n] E_sig;
  P_mu ~ normal(0, 20);
  P_sigma ~ normal(0, 20);
  E_sig = exp(D_all_sig * P_sigma);
  skill ~ normal(P_mu, E_sig);
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (skill_f - P_mu) ./ exp(D_f_sig * P_sigma);
  Dev_m = (skill_m - P_mu) ./ exp(D_m_sig * P_sigma);
  Dev_o = (skill_o - P_mu) ./ exp(D_o_sig * P_sigma);
}
