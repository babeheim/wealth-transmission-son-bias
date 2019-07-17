
rm(list = ls())

source("./project_support.r")

measure <- "poland_education"

d0 <- read.csv("./data/poland_education.csv", stringsAsFactors = FALSE)

d0$age <- 2010-d0$yrbirth
d0$age[1:1640] <- 2009-d0$yrbirth[1:1640]

#one model for everyone
ed1 <- rep(NA, nrow(d0))
for (i in 1:nrow(d0)) {
  if(is.na(d0$ed[i])) ed1[i] <- NA
  else if(d0$ed[i] == 1) ed1[i] <- 2
  else if(d0$ed[i] == 2) ed1[i] <- 6
  else if(d0$ed[i] == 3) ed1[i] <- 8
  else if(d0$ed[i] == 4) ed1[i] <- 11
  else if(d0$ed[i] == 5) ed1[i] <- 14
}
d <- data.frame(d0, ed1)
drop <- which(d$age<25 | is.na(d$age * d$ed1 * d$male))
d <- d[-drop, -(14:21)]
d$age <- d$age/10
D_mu <- cbind((1-d$male), (1-d$male) * d$age, (1-d$male) * d$age^2,
  d$male, d$male * d$age, d$male * d$age^2)
D_sig <- cbind(rep(1, nrow(d)), d$age)


#Now need to arrange data so that father/mother data is
# in same row as offspring...

edf1 <- rep(NA, nrow(d))
edm1 <- rep(NA, nrow(d))
for (i in 1:nrow(d)) {
  if(is.na(d$edf[i])) edf1[i] <- NA
  else if(d$edf[i] == 1) edf1[i] <- 2
  else if(d$edf[i] == 2) edf1[i] <- 6
  else if(d$edf[i] == 3) edf1[i] <- 8
  else if(d$edf[i] == 4) edf1[i] <- 11
  else if(d$edf[i] == 5) edf1[i] <- 14
}
for (i in 1:nrow(d)) {
  if(is.na(d$edm[i])) edm1[i] <- NA
  else if(d$edm[i] == 1) edm1[i] <- 2
  else if(d$edm[i] == 2) edm1[i] <- 6
  else if(d$edm[i] == 3) edm1[i] <- 8
  else if(d$edm[i] == 4) edm1[i] <- 11
  else if(d$edm[i] == 5) edm1[i] <- 14
}
agef <- NULL
agem <- NULL
for (i in 1:nrow(d)) {
  agef <- c(agef, (d0$age[which(d0$id == d$idf[i])][1])/10)
  agem <- c(agem, (d0$age[which(d0$id == d$idm[i])][1])/10)
}
dp <- data.frame(d, edf1, edm1, agef, agem)
drop <- which(is.na(dp$edf1 * dp$edm1 * dp$agef * dp$agem))
dp <- dp[-drop, ]

#make list for stan
data_list <- list(
  n = nrow(d), 
  ed = d$ed1, 
  D_mu = D_mu, 
  D_sig = D_sig
)

expect_true(data_list$n == 20318)


# stan was not used!

# tic(paste("fit", measure, "stan model"))

# m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
#   iter = n_iter, warmup = n_warmup, chains = 1, init = init)

# save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

# toc(log = TRUE)



# originally n_iter = 1000 and n_warmup = 500
#post <- extract(m1, permute=TRUE)

#Using point estimates from mle2
D_mu_o <- cbind((1-dp$male), (1-dp$male) * dp$age, (1-dp$male) * dp$age^2, dp$male, dp$male * dp$age, dp$male * dp$age^2)
D_mu_f <- cbind(rep(0, nrow(dp)), rep(0, nrow(dp)), rep(0, nrow(dp)), rep(1, nrow(dp)), dp$agef, dp$agef^2)
D_mu_m <- cbind(rep(1, nrow(dp)), dp$agem, dp$agem^2, rep(0, nrow(dp)), rep(0, nrow(dp)), rep(0, nrow(dp)))
par_mu <- c(2.5021, -.0072, -.0092, 2.2984, 0.0040, -.0069)
dev_o <- (dp$ed1 - exp(D_mu_o %*% par_mu))/exp(.9519-.0211 * dp$age)
dev_f <- (dp$edf1 - exp(D_mu_f %*% par_mu))/exp(.9519-.0211 * dp$agef)
dev_m <- (dp$edm1 - exp(D_mu_m %*% par_mu))/exp(.9519-.0211 * dp$agem)


n_sim <- n_iter

#find parent correlation
dparents <- data.frame(dev_f, dev_m)
dparents <- dparents[-which(duplicated(data.frame(dp$idf, dp$idm))), ]
names(dparents) <- c("Dev_f", "Dev_m")
r <- coef(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents))[[1]]
#now do the regression
midpar_temp <- 0.5 * (dev_f+dev_m)/sqrt(1/2 * (1+r))
dtemp <- data.frame(dev_o, midpar_temp)
names(dtemp) <- c("Dev_o", "midpar_temp")
#sons
dtemps <- dtemp[which(dp$male == 1), ]
m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
beta_s <- sample.naive.posterior(m_temp, n_sim)[[1]]  
#daughters
dtempd <- dtemp[which(dp$male == 0), ]
m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
beta_d <- sample.naive.posterior(m_temp, n_sim)[[1]]

# reported values here
beta_s_est <- 0.303
beta_s_se  <- 0.012
beta_d_est <- 0.281
beta_d_se  <- 0.01

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.023
beta_diff_se <- 0.016

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_sons = sum(d$male == 1),
  n_daughters = sum(d$male == 0),
  beta_son_est = sprintf("%1.3f", mean(beta_s)),
  beta_son_se  = sprintf("%1.3f", sd(beta_s)),
  beta_dau_est = sprintf("%1.3f", mean(beta_d)),
  beta_dau_se  = sprintf("%1.3f", sd(beta_d)),
  beta_diff_est = sprintf("%1.3f", mean(beta_diff)),
  beta_diff_se = sprintf("%1.3f", sd(beta_diff))
)

write_json(report, path = paste0("./temp/", measure, ".json"), pretty = TRUE)
