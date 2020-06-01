
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "maya_weight"

d0 <- read.csv("./data/maya_weight.csv", stringsAsFactors = FALSE)

names(d0)[which(names(d0) == "sex")] <- "male"
d0$male[which(d0$male == 2)] <- 0

# using log weight! all folks:
d <- data.frame(d0$male, as.numeric(paste(d0$age))/10, log(d0$weight))
# notice log
names(d) <- c("male", "age", "weight")
drop <- which(is.na(d$age * d$male * d$weight))
d <- d[-drop, ]
D_mu <- cbind(1 - d$male, (1 - d$male) * d$age, (1 - d$male) * d$age^2,
  d$male, d$male * d$age, d$male * d$age^2)
D_sig <- cbind(rep(1, nrow(d)), d$age)

# now folks with both parents (could save a lot by imputing...)
# construct parent weights to get sample sizes
weightf <- rep(NA, nrow(d0))
weightm <- rep(NA, nrow(d0))
agef <- rep(NA, nrow(d0))
agem <- rep(NA, nrow(d0))
for (i in 1:nrow(d0)) {
  if (length(which(d0$id == d0$fid[i])) > 0) {
    weightf[i] <- d0$weight[which(d0$id == d0$fid[i])]
    agef[i] <- as.numeric(paste(d0$age[which(d0$id == d0$fid[i])]))
  }
  if (length(which(d0$id == d0$mid[i])) > 0) {
    weightm[i] <- d0$weight[which(d0$id == d0$mid[i])]
    agem[i] <- as.numeric(paste(d0$age[which(d0$id == d0$mid[i])]))
  }
}
do <- data.frame(d0, log(weightf), log(weightm), agef/10, agem/10)
do$age <- as.numeric(paste(do$age))/10
do$weight <- log(do$weight)
names(do) <- c(names(d0), "weightf", "weightm", "agef", "agem")
drop <- which(is.na(do$male * do$weight * do$weightf *
  do$weightm * do$age * do$agef * do$agem))
do <- do[-drop, ]
D_mu_o <- cbind(1 - do$male, (1 - do$male) * do$age, (1 - do$male) * do$age^2,
  do$male, do$male * do$age, do$male * do$age^2)
D_sig_o <- cbind(rep(1, nrow(do)), do$age)
D_mu_f <- cbind(rep(0, nrow(do)), rep(0, nrow(do)), rep(0, nrow(do)),
  rep(1, nrow(do)), do$agef, do$agef^2)
D_sig_f <- cbind(rep(1, nrow(do)), do$agef)
D_mu_m <- cbind(rep(1, nrow(do)), do$agem, do$agem^2, rep(0, nrow(do)),
  rep(0, nrow(do)), rep(0, nrow(do)))
D_sig_m <- cbind(rep(1, nrow(do)), do$agem)


#make list for stan
data_list <- list(
  n = nrow(d),
  no = nrow(do),
  weight = d$weight,
  weight_o = do$weight,
  weight_f = do$weightf,
  weight_m = do$weightm,
  D_mu = D_mu,
  D_sig = D_sig,
  D_mu_o = D_mu_o,
  D_sig_o = D_sig_o,
  D_mu_f = D_mu_f,
  D_sig_f = D_sig_f,
  D_mu_m = D_mu_m,
  D_sig_m = D_sig_m
)

expect_true(data_list$n == 252)
expect_true(data_list$no == 135)



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
  #find parent correlation
  dparents <- data.frame(post$Dev_f[i, ], post$Dev_m[i, ])
  dparents <- dparents[-which(duplicated(data.frame(do$fid, do$mid))), ]
  names(dparents) <- c("Dev_f", "Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents), 1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+r))
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(do$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(do$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.251
beta_s_se  <- 0.125
beta_d_est <- 0.231
beta_d_se  <- 0.109

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.02
beta_diff_se <- 0.161

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_sons = sum(do$male == 1),
  n_daughters = sum(do$male == 0),
  beta_son_est = sprintf("%1.3f", mean(beta_s)),
  beta_son_se  = sprintf("%1.3f", sd(beta_s)),
  beta_dau_est = sprintf("%1.3f", mean(beta_d)),
  beta_dau_se  = sprintf("%1.3f", sd(beta_d)),
  beta_diff_est = sprintf("%1.3f", mean(beta_diff)),
  beta_diff_se = sprintf("%1.3f", sd(beta_diff))
)

write_json(report, path = paste0("./temp/", measure, ".json"), pretty = TRUE)
