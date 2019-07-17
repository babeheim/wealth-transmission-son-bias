data {
  int<lower=0> np;
  int<lower=0> no;
  int<lower=0> n1;
  vector[np] rs_p;
  vector[no] rs_o;
  vector[n1] rs_f_1;
  vector[n1] rs_m_1;
  vector[n1] rs_o_1;
  matrix[np,2] D_mu_p;
  matrix[np,3] D_sig_p;
  matrix[no,4] D_mu_o;
  matrix[no,2] D_sig_o;
  matrix[n1,2] D_mu_f_1;
  matrix[n1,3] D_sig_f_1;
  matrix[n1,2] D_mu_m_1;
  matrix[n1,3] D_sig_m_1;
  matrix[n1,4] D_mu_o_1;
  matrix[n1,2] D_sig_o_1;
  vector[no] male_o;
  vector[n1] male_1;
  vector[np] ones_p;
  vector[no] ones_o;
  vector[n1] ones_1;
}
parameters {
  vector[2] P_mu_p;
  vector[3] P_sigma_p;
  real max_RS;
  real max_RS_o_m;
  real max_RS_o_f;
  real max_sig;
  vector[4] P_mu_o;
  vector[2] P_sigma_o;
}
model {
  P_mu_p ~ normal(0, 20);
  P_sigma_p ~ normal(0, 20);
  max_RS ~ exponential(0.01);
  max_RS_o_m ~ exponential(0.01);
  max_RS_o_f ~ exponential(0.01);
  max_sig ~ exponential(0.01);
  P_mu_o ~ normal(0, 20);
  P_sigma_o ~ normal(0, 20);
  rs_p ~ normal(max_RS*ones_p ./ (1+exp(-D_mu_p*P_mu_p)), exp(D_sig_p*P_sigma_p));
  rs_o ~ normal(((1-male_o)*max_RS_o_f + male_o*max_RS_o_m).*ones_o ./ (1+exp(-D_mu_o*P_mu_o)), max_sig*ones_o ./ (1+exp(-D_sig_o*P_sigma_o)));
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  Dev_f = (rs_f_1 - max_RS*ones_1 ./ (1+exp(-D_mu_f_1*P_mu_p))) ./ exp(D_sig_f_1*P_sigma_p);
  Dev_m = (rs_m_1 - max_RS*ones_1 ./ (1+exp(-D_mu_m_1*P_mu_p))) ./ exp(D_sig_m_1*P_sigma_p);
  Dev_o = (rs_o_1 - ((1-male_1)*max_RS_o_f + male_1*max_RS_o_m) .* ones_1 ./ (1+exp(-D_mu_o_1*P_mu_o))) ./ (max_sig*ones_1 ./ (1+exp(-D_sig_o_1*P_sigma_o)));
}
