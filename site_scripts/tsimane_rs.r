
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "tsimane_rs"

d0 <- read.csv("./data/tsimane_rs.csv", stringsAsFactors = FALSE)

#parents only
dp <- data.frame(
  c(d0$fid, d0$mid),
  c(d0$rsf, d0$rsm),
  c(d0$agef, d0$agem)/10,
  c(rep(1, nrow(d0)),
    rep(0, nrow(d0)))
)
names(dp) <- c("id", "rs", "age", "male")
drop <- which(duplicated(dp$id) | is.na(dp$rs * dp$age * dp$male) | dp$age < 1.8)
dp <- dp[-drop, ]
#706 parents
D_mu_p <- cbind(rep(1, nrow(dp)), dp$age)
D_sig_p <- cbind(rep(1, nrow(dp)), dp$age, dp$age^2)
#706 par

#now kids
drop <-  which(
  duplicated(d0$id) |
  is.na(d0$rs * d0$age * d0$male) |
  d0$age < 18 |
  d0$id %in% d0$fid |
  d0$id %in% d0$mid
)
do <- d0[-drop, ]
do$age <- do$age/10
D_mu_o <- cbind(1 - do$male, (1 - do$male) * do$age, do$male, do$male * do$age)
D_sig_o <- cbind(rep(1, nrow(do)), do$age)
#935 kids

#folks with parents.
drop <-  which(
  duplicated(d0$id) |
  is.na(d0$rs * d0$age * d0$male * d0$rsf * d0$rsm * d0$agef * d0$agem) |
  d0$age < 18 |
  d0$id %in% d0$fid |
  d0$id %in% d0$mid
)
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
#706 people
D_mu_f_1 <- cbind(rep(1, nrow(d1)), d1$agef)
D_sig_f_1 <- cbind(rep(1, nrow(d1)), d1$agef, d1$agef^2)
D_mu_m_1 <- cbind(rep(1, nrow(d1)), d1$agem)
D_sig_m_1 <- cbind(rep(1, nrow(d1)), d1$agem, d1$agem^2)
D_mu_o_1 <- cbind(1-d1$male, (1-d1$male) * d1$age, d1$male, d1$male * d1$age)
D_sig_o_1 <- cbind(rep(1, nrow(d1)), d1$age)

(mean(d1$rs[which(d1$male == 1)])-mean(d1$rs[which(d1$male == 0)]))/sd(d1$rs)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  np = nrow(dp),
  no = nrow(do),
  n1 = nrow(d1),
  rs_p = dp$rs,
  rs_o = do$rs,
  rs_f_1 = d1$rsf,
  rs_m_1 = d1$rsm,
  rs_o_1 = d1$rs,
  D_mu_p = D_mu_p,
  D_sig_p = D_sig_p,
  D_mu_o = D_mu_o,
  D_sig_o = D_sig_o,
  D_mu_f_1 = D_mu_f_1,
  D_sig_f_1 = D_sig_f_1,
  D_mu_m_1 = D_mu_m_1,
  D_sig_m_1 = D_sig_m_1,
  D_mu_o_1 = D_mu_o_1,
  D_sig_o_1 = D_sig_o_1,
  male_o = do$male,
  male_1 = d1$male,
  ones_p = rep(1, nrow(dp)),
  ones_o = rep(1, nrow(do)),
  ones_1 = rep(1, nrow(d1))
)

expect_true(data_list$np == 706)
expect_true(data_list$no == 935)
expect_true(data_list$n1 == 841)



tic(paste("fit", measure, "stan model"))

if (!file.exists(file.path("temp", paste0(measure, ".stanfit")))){

  m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
    iter = n_iter, warmup = n_warmup, chains = 1, init = init)

  save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

}
toc(log = TRUE)



load(file.path("temp", paste0(measure, ".stanfit")))

post <- extract(m1, permute=TRUE)

n_sim <- min(10000, length(post$lp__))

beta_s <- rep(NA, n_sim)
beta_d <- rep(NA, n_sim)

for (i in 1:n_sim) {

  #find parent correlation
  dparents <- data.frame(post$Dev_f[i, ], post$Dev_m[i, ])
  dparents <- dparents[-which(duplicated(data.frame(d1$fid, d1$mid))), ]
  names(dparents) <- c("Dev_f", "Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents), 1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+r))
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]

}


# reported values here
beta_s_est <- 0.052
beta_s_se  <- 0.049
beta_d_est <- 0.088
beta_d_se  <- 0.047

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.036)
beta_diff_se <- 0.068

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_sons = sum(d1$male == 1),
  n_daughters = sum(d1$male == 0),
  beta_son_est = sprintf("%1.3f", mean(beta_s)),
  beta_son_se  = sprintf("%1.3f", sd(beta_s)),
  beta_dau_est = sprintf("%1.3f", mean(beta_d)),
  beta_dau_se  = sprintf("%1.3f", sd(beta_d)),
  beta_diff_est = sprintf("%1.3f", mean(beta_diff)),
  beta_diff_se = sprintf("%1.3f", sd(beta_diff))
)

write_json(report, path = paste0("./temp/", measure, ".json"), pretty = TRUE)
