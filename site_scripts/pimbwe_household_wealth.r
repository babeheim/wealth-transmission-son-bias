
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "pimbwe_household_wealth"

d0 <- read.csv("./data/pimbwe_household_wealth.csv", stringsAsFactors = FALSE)

# looks like parents tend to have more wealth, even when controlling for age.
# So I'm gonna do separate models.
# parent wealth VERY strongly correlated. In most cases,
# wealthm=wealthf. Very simple imputation ought to work.

d0$sex[which(d0$sex == "m")] <- 1
d0$sex[which(d0$sex == "f")] <- 0
d0$sex <- as.numeric(d0$sex)
d0$male <- d0$sex

#parents only
dp <- data.frame(
  c(d0$fid, d0$mid),
  log(c(d0$wealthf, d0$wealthm)),
  c(d0$agef, d0$agem),
  c(rep(1, nrow(d0)), rep(0, nrow(d0)))
)
names(dp) <- c("id", "wealth", "age", "male")
drop <- which(is.na(dp$age * dp$wealth))
dp <- dp[-drop, ]
dp <- dp[-which(duplicated(dp$id)), ]
dp$age <- dp$age/10
D_mu_p <- cbind(rep(1, nrow(dp)), dp$age)


#now just the kids
drop <- which(is.na(d0$age * d0$wealth) | d0$age > 70)
# dropping anyone over 70 because thats really the range
# that matters for people who actually have parents
do <- d0[-drop, ]
do <- do[-which(do$id %in% c(do$fid, do$mid)), ]
do$wealth <- log(do$wealth)
do$age <- do$age/10


#now folks with parents
drop <- which(is.na(d0$age * d0$wealth * d0$agem * d0$wealthm * d0$agef * d0$wealthf))
d1 <- d0[-drop, ]
d1$wealth <- log(d1$wealth)
d1$wealthf <- log(d1$wealthf)
d1$wealthm <- log(d1$wealthm)
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
D_mu_f <- cbind(rep(1, nrow(d1)), d1$agef)
D_mu_m <- cbind(rep(1, nrow(d1)), d1$agem)

(mean(d1$wealth[which(d1$male == 1)])-mean(d1$wealth[which(d1$male == 0)]))/sd(d1$wealth)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  np = nrow(dp),
  no = nrow(do),
  n1 = nrow(d1),
  wealth_p = dp$wealth,
  wealth_o = do$wealth,
  D_mu_p = D_mu_p,
  wealth_f = d1$wealthf,
  wealth_m = d1$wealthm,
  wealth_o1 = d1$wealth,
  D_mu_f = D_mu_f,
  D_mu_m = D_mu_m
)

expect_true(data_list$np == 173)
expect_true(data_list$no == 437)
expect_true(data_list$n1 == 169)



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

for(i in 1:n_sim) {

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
beta_s_est <- 0.332
beta_s_se  <- 0.101
beta_d_est <- (-0.081)
beta_d_se  <- 0.105

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.413
beta_diff_se <- 0.147

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
