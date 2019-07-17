
rm(list = ls())

source("./project_support.r")

measure <- "krummhorn_land"

d0 <- read.csv("./data/krummhorn_land.csv", stringsAsFactors = FALSE)

d0$kuniony <- as.numeric(d0$kuniony)
d0$male <- d0$sex
d0$male[which(d0$male == 2)] <- 0

#dp: parents only 
drop <- which(is.na(d0$kuniony))
dp <- d0[-drop, ]
dp$yearf <- dp$kuniony
dp$yearf <- (dp$yearf-1800)/100
D_mu_f <- cbind(rep(1, nrow(dp)), dp$yearf)
D_sig_f <- cbind(rep(1, nrow(dp)), dp$yearf, dp$yearf^3)


#d: offspring only
d <- d0    #no drops needed
d$year <- d$uniony
d$year <- (d$year-1800)/100
D <- cbind(rep(1, nrow(d)), d$year, d$year^2)  #same for mu and sig

#d1: those with parents
drop <- which(is.na(d0$kuniony) | is.na(d0$landf) | is.na(d0$land) | is.na(d0$uniony))
d1 <- d0[-drop, ]
d1$year <- d1$uniony
d1$yearf <- d1$kuniony
d1$yearf <- (d1$yearf-1800)/100
d1$year <- (d1$year-1800)/100
D_mu_f1 <- cbind(rep(1, nrow(d1)), d1$yearf)
D_sig_f1 <- cbind(rep(1, nrow(d1)), d1$yearf, d1$yearf^3)
D1 <- cbind(rep(1, nrow(d1)), d1$year, d1$year^2)  #same for mu and sig

(mean(d1$land[which(d1$male == 1)])-mean(d1$land[which(d1$male == 0)]))/sd(d1$land)

#make list for stan
data_list <- list(
  n = nrow(d), 
  np = nrow(dp), 
  n1 = nrow(d1), 
  land = d$land/10, 
  land_f = dp$landf/10, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D = D, 
  land_o1 = d1$land/10, 
  land_f1 = d1$landf/10, 
  D_mu_f1 = D_mu_f1, 
  D_sig_f1 = D_sig_f1, 
  D1 = D1, 
  male = d1$male
)

expect_true(data_list$n == 1602)
expect_true(data_list$np == 1482)
expect_true(data_list$n1 == 1482)



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

  #do the regression
  par_temp <- post$Dev_f[i, ]
  dtemp <- data.frame(post$Dev_o[i, ], par_temp)
  names(dtemp) <- c("Dev_o", "Dev_f")

  #sons
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]

  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]

}

# reported values here
beta_s_est <- 0.6
beta_s_se  <- 0.022
beta_d_est <- 0.521
beta_d_se  <- 0.025

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.079
beta_diff_se <- 0.032

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
