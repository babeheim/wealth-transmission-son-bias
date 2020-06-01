
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "lamalera_house_quality"

d0 <- read.csv("./data/lamalera_house_quality.csv", stringsAsFactors = FALSE)

id <- d0$id
male <- d0$male
age <- d0$age
wealth <- d0$wealth
for (i in 1:nrow(d0)) {
  if (
    !is.na(d0$fid[i]) &
    !(d0$fid[i] %in% id) &
    !is.na(d0$wealthp[i]) &
    !is.na(d0$agetodayf[i]) &
    is.na(d0$ageatdeathf[i])
  ) {
    id <- c(id, d0$fid[i])
    male <- c(male, 1)
    age <- c(age, d0$agetodayf[i])
    wealth <- c(wealth, d0$wealthp[i])
  }
  if (
    !is.na(d0$mid[i]) &
    !(d0$mid[i] %in% id) &
    !is.na(d0$wealthp[i]) &
    !is.na(d0$agetodaym[i]) &
    is.na(d0$ageatdeathm[i])
  ) {
    id <- c(id, d0$mid[i])
    male <- c(male, 0)
    age <- c(age, d0$agetodaym[i])
    wealth <- c(wealth, d0$wealthp[i])
  }
}
age <- age/10
d <- data.frame(id, male, age, wealth)
D_mu_1 <- cbind(rep(1, nrow(d)), d$age, d$age^2)

has_parents <- which(
  ((!is.na(d0$agetodaym) & is.na(d0$ageatdeathm)) | (d0$agetodaym - d0$ageatdeathm) < 5) &
  ((!is.na(d0$agetodayf) & is.na(d0$ageatdeathf)) | (d0$agetodayf - d0$ageatdeathf) < 5)
)
dp <- d0[has_parents, ]
dp$agef <- dp$agetodayf
dp$agem <- dp$agetodaym
dp$age <- dp$age/10
dp$agef <- dp$agef/10
dp$agem <- dp$agem/10
D_mu_o <- cbind(rep(1, nrow(dp)), dp$age, dp$age^2)
D_mu_m <- cbind(rep(1, nrow(dp)), dp$agem, dp$agem^2)
D_mu_f <- cbind(rep(1, nrow(dp)), dp$agef, dp$agef^2)

(mean(dp$wealth[which(dp$male == 1)])-mean(dp$wealth[which(dp$male == 0)]))/sd(dp$wealth)
mean(dp$age * 10)

#make list for stan
data_list <- list(
  n = nrow(d), 
  np = nrow(dp), 
  wealth = d$wealth, 
  wealth_o = dp$wealth, 
  wealth_p = dp$wealthp, 
  D_mu_1 = D_mu_1, 
  D_mu_o = D_mu_o, 
  D_mu_f = D_mu_f, 
  D_mu_m = D_mu_m, 
  male = dp$male
)

expect_true(data_list$n == 536)
expect_true(data_list$np == 93)



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

for(i in 1:n_sim)
{
  #parent correlation essentially 1
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(dp$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(dp$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.143
beta_s_se  <- 0.139
beta_d_est <- 0.115
beta_d_se  <- 0.162

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.028
beta_diff_se <- 0.213

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_sons = sum(dp$male == 1),
  n_daughters = sum(dp$male == 0),
  beta_son_est = sprintf("%1.3f", mean(beta_s)),
  beta_son_se  = sprintf("%1.3f", sd(beta_s)),
  beta_dau_est = sprintf("%1.3f", mean(beta_d)),
  beta_dau_se  = sprintf("%1.3f", sd(beta_d)),
  beta_diff_est = sprintf("%1.3f", mean(beta_diff)),
  beta_diff_se = sprintf("%1.3f", sd(beta_diff))
)

write_json(report, path = paste0("./temp/", measure, ".json"), pretty = TRUE)
