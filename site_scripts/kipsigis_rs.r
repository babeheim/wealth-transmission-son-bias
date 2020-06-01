
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "kipsigis_rs"

d0 <- read.csv("./data/kipsigis_rs.csv", stringsAsFactors = FALSE)

#d: all folks (fathers are in there)
drop <- which(is.na(d0$RS * d0$age) | d0$RS == 0)
d <- d0[-drop, ]
d <- data.frame(d$age, d$male, d$RS)
names(d) <- c("age", "male", "rs")
d$age <- d$age/10

D_mu <- cbind(1-d$male, (1-d$male) * d$age, d$male, d$male * d$age)
D_sig <- cbind(1-d$male, d$male, d$male * d$age) 


#d1: people with father data
d1 <- d0[which(!is.na(d0$RSf) & !is.na(d0$RS) & d0$RS>0), ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
D_1_mu <- cbind(1-d1$male, (1-d1$male) * d1$age, d1$male, d1$male * d1$age)
D_1_sig <- cbind(1-d1$male, d1$male, d1$male * d1$age) 
D_1_mu_f <- cbind(rep(0, nrow(d1)), rep(0, nrow(d1)), rep(1, nrow(d1)), d1$agef)
D_1_sig_f <- cbind(rep(0, nrow(d1)), rep(1, nrow(d1)), d1$agef) 
(mean(d1$RS[which(d1$male == 1)])-mean(d1$RS[which(d1$male == 0)]))/sd(d1$RS)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  n = nrow(d), 
  n1 = nrow(d1), 
  rs = d$rs, 
  rs1 = d1$RS, 
  rs1f = d1$RSf, 
  D_mu = D_mu, 
  D_sig = D_sig, 
  D_1_mu = D_1_mu, 
  D_1_sig = D_1_sig, 
  D_1_mu_f = D_1_mu_f, 
  D_1_sig_f = D_1_sig_f
)

expect_true(data_list$n == 398)
expect_true(data_list$n1 == 268)



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
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.062
beta_s_se  <- 0.109
beta_d_est <- 0.031
beta_d_se  <- 0.127

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.031
beta_diff_se <- 0.166

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
