
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "lamalera_boat_shares"

d0 <- read.csv("./data/lamalera_boat_shares.csv", stringsAsFactors = FALSE)

id <- d0$id
male <- d0$male
age <- d0$age
shares <- d0$shares06
for (i in 1:nrow(d0)) {
  if (
    !is.na(d0$fid[i]) & !(d0$fid[i] %in% id) &
    !is.na(d0$shares06p[i]) & !is.na(d0$agef[i])
  ) {
    id <- c(id, d0$fid[i])
    male <- c(male, 1)
    age <- c(age, d0$agef[i])
    shares <- c(shares, d0$shares06p[i])
  }
  if (
    !is.na(d0$mid[i]) & !(d0$mid[i] %in% id) &
    !is.na(d0$shares06p[i]) & !is.na(d0$agem[i])
  ) {
    id <- c(id, d0$mid[i])
    male <- c(male, 0)
    age <- c(age, d0$agem[i])
    shares <- c(shares, d0$shares06p[i])
  }
}

age <- age/10
d <- data.frame(id, male, age, shares)
D_mu_1 <- cbind(rep(1, nrow(d)), d$age, d$age^2)
D_sig_1 <- cbind(rep(1, nrow(d)), d$age) 

#now just those with parents
drop <- which(is.na(d0$shares06p) | is.na(d0$agem) | is.na(d0$agef))
d2 <- d0[-drop, ]
d2$age <- d2$age/10
d2$agef <- d2$agef/10
d2$agem <- d2$agem/10
D_mu_f <- cbind(rep(1, nrow(d2)), d2$agef, d2$agef^2)
D_sig_f <- cbind(rep(1, nrow(d2)), d2$agef)
D_mu_m <- cbind(rep(1, nrow(d2)), d2$agem, d2$agem^2)
D_sig_m <- cbind(rep(1, nrow(d2)), d2$agem)
D_mu_o <- cbind(rep(1, nrow(d2)), d2$age, d2$age^2)
D_sig_o <- cbind(rep(1, nrow(d2)), d2$age)

(mean(d2$shares06[which(d2$male == 1)])-mean(d2$shares06[which(d2$male == 0)]))/sd(d2$shares06)
mean(d2$age * 10)

#make list for stan
data_list <- list(
  n = nrow(d), 
  n2 = nrow(d2), 
  shares = d$shares, 
  shares_o = d2$shares06, 
  shares_p = d2$shares06p, 
  D_mu_1 = D_mu_1, 
  D_sig_1 = D_sig_1, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D_mu_m = D_mu_m, 
  D_sig_m = D_sig_m, 
  male = d2$male
)

expect_true(data_list$n == 560)
expect_true(data_list$n2 == 119)



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
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])
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
beta_s_est <- 0.152
beta_s_se  <- 0.147
beta_d_est <- 0.121
beta_d_se  <- 0.14

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.031
beta_diff_se <- 0.201

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
