
rm(list = ls())

source("./project_support.r")

measure <- "pimbwe_rs"

d0 <- read.csv("./data/pimbwe_rs.csv", stringsAsFactors = FALSE)

# Problem: we have a lot of cases where the recorded parent's
# age is less than the offsprings, or only a few years greater.
# Normally we might say "drop anyhing where age difference is
# less than 15", since this is maybe the youngest someone
# would reproduce. But since this sampling takes place over
# a ten year interval, it's possible that parents would be
# sampled at an earlier time than offspring, and so make
# parents appear to be younger.
# So I'm going to say: if age difference between offspring
# and either parent is <5, drop them. Dropping isn't great
# - but the problem is we don't know which data point, 
# parent or offspring, is erroneous!

drop <- which(d0$daddiff<5 | d0$mumdiff<5)
d0 <- d0[-drop, ]
d0$sex[which(d0$sex == "m")] <- 1
d0$sex[which(d0$sex == "f")] <- 0
d0$sex <- as.numeric(d0$sex)
names(d0)[which(names(d0) == "sex")] <- "male"


#extract parents only
dp <- data.frame(
  c(d0$fid, d0$mid),
  c(d0$rsf, d0$rsm),
  c(d0$agef, d0$agem),
  c(rep(1, nrow(d0)),
    rep(0, nrow(d0)))
)
names(dp) <- c("id", "rs", "age", "male")
dp <- dp[-which(duplicated(dp$id)), ]
dp$age <- dp$age/10
drop <- which(is.na(dp$rs))
dp <- dp[-drop, ]
D_p <- cbind(1 - dp$male, dp$male) 

#now the kids
do <- data.frame(d0$id, d0$rs, d0$age, d0$male)
names(do) <- c("id", "rs", "age", "male")
do <- do[-which(do$id %in% dp$id), ]
do$age <- do$age/10
drop <- which(is.na(do$rs))
do <- do[-drop, ]
#removing an 80 year old outlier guy with no kids.
do <- do[-which(do$age>8), ]
D_mu_o <- cbind(1-do$male, (1-do$male) * do$age, do$male, do$male * do$age)
D_sig_o <- cbind(1-do$male, (1-do$male) * do$age, (1-do$male) * do$age^2,
  do$male, do$male * do$age, do$male * do$age^2)

#now folks with parents
drop <- which(is.na(d0$rs * d0$rsf * d0$rsm * d0$age * d0$agef * d0$agem))
d <- d0[-drop, ]
d <- data.frame(d$age/10, d$rs, d$male, d$agem/10, d$rsm, d$agem/10, d$rsm, d$fid, d$mid)
names(d) <- c("age", "rs", "male", "agem", "rsm", "agef", "rsf", "fid", "mid")
D_mu_o2 <- cbind(1-d$male, (1-d$male) * d$age, d$male, d$male * d$age)
D_sig_o2 <- cbind(1-d$male, (1-d$male) * d$age, (1-d$male) * d$age^2,
  d$male, d$male * d$age, d$male * d$age^2)
D_f <- cbind(rep(0, nrow(d)), rep(1, nrow(d)))
D_m <- cbind(rep(1, nrow(d)), rep(0, nrow(d)))

(mean(d$rs[which(d$male == 1)])-mean(d$rs[which(d$male == 0)]))/sd(d$rs)
mean(d$age * 10)

#make list for stan
data_list <- list(
  np = nrow(dp),
  no = nrow(do),
  n = nrow(d),
  rs_p = dp$rs,
  rs_o = do$rs,
  rs_f = d$rsf,
  rs_m = d$rsm,
  rs_o2 = d$rs,
  D_p = D_p,
  D_mu_o = D_mu_o,
  D_sig_o = D_sig_o,
  D_f = D_f,
  D_m = D_m,
  D_mu_o2 = D_mu_o2,
  D_sig_o2 = D_sig_o2,
  male_o = do$male,
  male_o2 = d$male
)

expect_true(data_list$n == 441)
expect_true(data_list$np == 311)
expect_true(data_list$no == 629)



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

  #parent correlation is exactly 1.
 
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(d$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]

  if (i %% 100 == 0) print(i)

}

# reported values here
beta_s_est <- (-0.037)
beta_s_se  <- 0.057
beta_d_est <- 0.035
beta_d_se  <- 0.074

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.072)
beta_diff_se <- 0.092

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
