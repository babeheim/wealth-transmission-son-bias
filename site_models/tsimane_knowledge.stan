data {
  int<lower=0> n;
  int<lower=0> n1;
  vector[n] skill;
  vector[n] age;
  vector[n] male;
  int com[n];
  vector[n1] skill_f;
  vector[n1] age_f;
  int com_f[n1];
  vector[n1] skill_m;
  vector[n1] age_m;
  int com_m[n1];
  vector[n1] skill_o;
  vector[n1] age_o;
  vector[n1] male_o;
  int com_o[n1];
}
parameters {
  real B_age;
  real B_male;
  real<lower=0> sigma;
  vector[8] com_effect;
  real mu_com;
  real<lower=0> sig_com;
}
model {
  B_age ~ normal(0, 20);
  B_male ~ normal(0, 20);
  sigma ~ exponential(0.01);
  mu_com ~ normal(0, 20);
  sig_com ~ exponential(0.01);
  com_effect ~ normal(mu_com,sig_com);
  for (i in 1:n) {
    skill[i] ~ normal(com_effect[com[i]] + B_age*age[i] + B_male*male[i], sigma);
  }
}
generated quantities {
  vector[n1] Dev_o;
  vector[n1] Dev_f;
  vector[n1] Dev_m;
  for (i in 1:n1) {
    Dev_f[i] = (skill_f[i] - (com_effect[com_f[i]] + B_age*age_f[i] + B_male)) / sigma;
    Dev_m[i] = (skill_m[i] - (com_effect[com_m[i]] + B_age*age_m[i])) / sigma;
    Dev_o[i] = (skill_o[i] - (com_effect[com_o[i]] + B_age*age_o[i] + B_male*male_o[i])) / sigma;
  }
}
