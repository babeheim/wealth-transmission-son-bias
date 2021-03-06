
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "chewa_education"

d0 <- read.csv("./data/chewa_education.csv", stringsAsFactors = FALSE)

d0$male[which(d0$male == 2)] <- 0

#d1: everyone with age and edu data
#playing around with the data shows that every parent who has an ID and
# data is already in the "offspring" column, so we don't have to grab parents manually

drop <- which(is.na(d0$ed) | d0$age<18)
d1 <- d0[-drop, ]
d1$age <- d1$age/10
D_mu_1 <- cbind(1-d1$male, d1$male, (1-d1$male) * d1$age, d1$male * d1$age, d1$male * d1$age^2)
D_sig_1 <- cbind(1-d1$male, d1$male, (1-d1$male) * d1$age, (1-d1$male) * d1$age^2, d1$male * d1$age)

#d2: those with complete parent data

drop <- which(is.na(d0$age) | is.na(d0$ed) | is.na(d0$edm) | is.na(d0$edf) | d0$age<18 | d0$agef<18 | is.na(d0$agef) | is.na(d0$agem) | d0$agem<18)
d2 <- d0[-drop, ]
d2$age <- d2$age/10
d2$agef <- d2$agef/10
d2$agem <- d2$agem/10
D_mu_o <- cbind(1-d2$male, d2$male, (1-d2$male) * d2$age, d2$male * d2$age, d2$male * d2$age^2)
D_sig_o <- cbind(1-d2$male, d2$male, (1-d2$male) * d2$age, (1-d2$male) * d2$age^2, d2$male * d2$age)
D_mu_f <- cbind(rep(0, nrow(d2)), rep(1, nrow(d2)), rep(0, nrow(d2)), d2$agef, d2$agef^2)
D_sig_f <- cbind(rep(0, nrow(d2)), rep(1, nrow(d2)), rep(0, nrow(d2)), rep(0, nrow(d2)), d2$agef)
D_mu_m <- cbind(rep(1, nrow(d2)), rep(0, nrow(d2)), d2$agem, rep(0, nrow(d2)), rep(0, nrow(d2)))
D_sig_m <- cbind(rep(1, nrow(d2)), rep(0, nrow(d2)), d2$agem, d2$agem^2, rep(0, nrow(d2)))
ed_o <- d2$ed
ed_f <- d2$edf
ed_m <- d2$edm
male <- d2$male
(mean(d2$ed[which(d2$male == 1)])-mean(d2$ed[which(d2$male == 0)]))/sd(d2$ed)
mean(d2$age * 10)


#make list for stan
data_list <- list(
  n1 = nrow(D_mu_1), 
  n2 = nrow(D_mu_o), 
  D_mu_1 = D_mu_1, 
  D_sig_1 = D_sig_1, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D_mu_m = D_mu_m, 
  D_sig_m = D_sig_m, 
  ed = d1$ed, 
  ed_o = ed_o, 
  ed_f = ed_f, 
  ed_m = ed_m, 
  male = male
)

expect_true(data_list$n1 == 2529)
expect_true(data_list$n2 == 245)



tic(paste("fit", measure, "stan model"))

if (!file.exists(file.path("temp", paste0(measure, ".stanfit")))){

  m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
    iter = n_iter, warmup = n_warmup, chains = 1, init = init)

  save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

}
toc(log = TRUE)



load(file.path("temp", paste0(measure, ".stanfit")))

post <-extract(m1, permute=TRUE)

n_sim <- min(10000, length(post$lp__))

beta_s <- rep(NA, n_sim)
beta_d <- rep(NA, n_sim)

for (i in 1:n_sim) {

  #find parent correlation
  dparents <- data.frame(post$Dev_f[i, ], post$Dev_m[i, ])
  dparents <- dparents[-which(duplicated(data.frame(d2$idf, d2$idm))), ]
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
beta_s_est <- 0.279
beta_s_se  <- 0.066
beta_d_est <- 0.235
beta_d_se  <- 0.085

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.044
beta_diff_se <- 0.106

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
