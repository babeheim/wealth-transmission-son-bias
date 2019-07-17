
rm(list = ls())

source("./project_support.r")

measure <- "kipsigis_livestock"

d0 <- read.csv("./data/kipsigis_livestock.csv", stringsAsFactors = FALSE)

#d: fit general age model
drop <- which(is.na(d0$age) | is.na(d0$livestock) | d0$age<18)
d <- d0[-drop, ]
#note that the "offspring" data includes the parents
d$age <- d$age/10
D_mu <- cbind((1-d$male), d$male, (1-d$male) * d$age, d$male * d$age, d$male * d$age^2)
D_sig <- cbind((1-d$male), d$male, d$male * d$age)

#dp: those with father data
drop <- which(is.na(d0$age) | is.na(d0$livestock) | d0$age<18 | is.na(d0$agef) | is.na(d0$livestockf))
dp <- d0[-drop, ]
dp$age <- dp$age/10
dp$agef <- dp$agef/10
D_mu_f <- cbind(rep(0, nrow(dp)), rep(1, nrow(dp)), rep(0, nrow(dp)), dp$agef, dp$agef^2)
D_sig_f <- cbind(rep(0, nrow(dp)), rep(1, nrow(dp)), dp$agef)
D_mu_o <- cbind((1-dp$male), dp$male, (1-dp$male) * dp$age, dp$male * dp$age, dp$male * dp$age^2)
D_sig_o <- cbind((1-dp$male), dp$male, dp$male * dp$age)
(mean(dp$livestock[which(dp$male == 1)])-mean(dp$livestock[which(dp$male == 0)]))/sd(dp$livestock)
mean(dp$age * 10)

#make list for stan
data_list <- list(
  n1 = nrow(d), 
  n2 = nrow(dp), 
  livestock = d$livestock, 
  D_mu = D_mu, 
  D_sig = D_sig, 
  livestock_o = dp$livestock, 
  livestock_f = dp$livestockf, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  male = dp$male
)

expect_true(data_list$n1 == 400)
expect_true(data_list$n2 == 270)



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
  #do the regression
  par_temp <- post$Dev_f[i, ]
  dtemp <- data.frame(post$Dev_o[i, ], par_temp)
  names(dtemp) <- c("Dev_o", "Dev_f")
  
  #sons
  dtemps <- dtemp[which(dp$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(dp$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.462
beta_s_se  <- 0.067
beta_d_est <- 0.093
beta_d_se  <- 0.096

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.369
beta_diff_se <- 0.116

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
