
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "gambia_weight"

d0 <- read.csv("./data/gambia_weight.csv", stringsAsFactors = FALSE)

#d1: anyone who has weight and age data. 
# the purpose of this is to get a good age model.
drop <- which(is.na(d0$age) | is.na(d0$weight) | d0$age < 18)
d1 <- d0[-drop, ]
names(d1)[which(names(d1) == "sex")] <- "male"
d1$male[which(d1$male == 2)] <- 0

age_d1 <- d1$age/100

D_mu_all <- cbind(
  rep(1, nrow(d1)),
  d1$male,
  age_d1,
  age_d1^2,
  age_d1^3,
  d1$male * age_d1,
  d1$male * age_d1^2,
  d1$male * age_d1^3
)
D_sig_all <- cbind(rep(1, nrow(d1)), d1$male, age_d1, age_d1 * d1$male)

#now d2: for just the sample with parents.
drop <- which(
  is.na(d0$age) |
  is.na(d0$weightf) |
  is.na(d0$weight) |
  is.na(d0$weightm) |
  d0$age < 18
)
d2 <- d0[-drop, ]
names(d2)[which(names(d2) == "sex")] <- "male"
d2$male[which(d2$male == 2)] <- 0
d2$age <- d2$age/100
d2$agef <- d2$agef/100
d2$agem <- d2$agem/100
D_mu <- cbind(rep(1, nrow(d2)), d2$male, d2$age, d2$age^2,
  d2$age^3, d2$male * d2$age, d2$male * d2$age^2, d2$male * d2$age^3)
D_mu_f <- cbind(rep(1, nrow(d2)), rep(1, nrow(d2)), d2$agef,
  d2$agef^2, d2$agef^3, d2$agef, d2$agef^2, d2$agef^3)
D_mu_m <- cbind(rep(1, nrow(d2)), rep(0, nrow(d2)), d2$agem,
  d2$agem^2, d2$agem^3, 0, 0, 0)
D_sig <- cbind(rep(1, nrow(d2)), d2$male, d2$age, d2$male * d2$age)
D_sig_f <- cbind(rep(1, nrow(d2)), rep(1, nrow(d2)), d2$agef, d2$agef)
D_sig_m <- cbind(rep(1, nrow(d2)), rep(0, nrow(d2)), d2$agem, 0)
weight <- d2$weight
weight_f <- d2$weightf
weight_m <- d2$weightm
male <- d2$male

(mean(d2$weight[which(d2$male == 1)])-mean(d2$weight[which(d2$male == 0)]))/sd(d2$weight)
mean(d2$age * 100)

#make list for stan
data_list <- list(
  n1 = nrow(d1), 
  n2 = nrow(d2), 
  weight_d1 = d1$weight, 
  D_mu_all = D_mu_all, 
  D_sig_all = D_sig_all, 
  weight = weight, 
  weight_f = weight_f, 
  weight_m = weight_m, 
  D_mu = D_mu, 
  D_mu_f = D_mu_f, 
  D_mu_m = D_mu_m, 
  D_sig = D_sig, 
  D_sig_f = D_sig_f, 
  D_sig_m = D_sig_m, 
  male = male
)

expect_true(data_list$n1 == 2355)
expect_true(data_list$n2 == 817)



tic(paste("fit", measure, "stan model"))

if (!file.exists(file.path("temp", paste0(measure, ".stanfit")))){

  m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
    iter = n_iter, warmup = n_warmup, chains = 1, init = init)

  save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

}
toc(log = TRUE)



load(file.path("temp", paste0(measure, ".stanfit")))

#sometimes this goes through very quickly without converging at all. Just try again.
post <- extract(m1, permute=TRUE)

n_sim <- min(10000, length(post$lp__))

beta_s <- rep(NA, n_sim)
beta_d <- rep(NA, n_sim)

for(i in 1:n_sim)
{
  #find parent correlation
  dparents <- data.frame(post$Dev_f[i, ], post$Dev_m[i, ])
  dparents <- dparents[-which(duplicated(data.frame(d2$fid, d2$mid))), ]
  names(dparents) <- c("Dev_f", "Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents), 1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+r))
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(d2$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d2$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.366
beta_s_se  <- 0.042
beta_d_est <- 0.443
beta_d_se  <- 0.04

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.078)
beta_diff_se <- 0.056

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_sons = sum(d2$male == 1),
  n_daughters = sum(d2$male == 0),
  beta_son_est = sprintf("%1.3f", mean(beta_s)),
  beta_son_se  = sprintf("%1.3f", sd(beta_s)),
  beta_dau_est = sprintf("%1.3f", mean(beta_d)),
  beta_dau_se  = sprintf("%1.3f", sd(beta_d)),
  beta_diff_est = sprintf("%1.3f", mean(beta_diff)),
  beta_diff_se = sprintf("%1.3f", sd(beta_diff))
)

write_json(report, path = paste0("./temp/", measure, ".json"), pretty = TRUE)
