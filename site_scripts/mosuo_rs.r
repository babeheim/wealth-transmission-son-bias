
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "mosuo_rs"

d0 <- read.csv("./data/mosuo_rs.csv", stringsAsFactors = FALSE)

d0fd <- d0[which(d0$code == "fd"), ]
d0md <- d0[which(d0$code == "md"), ]
d0fs <- d0[which(d0$code == "fs"), ]
d0ms <- d0[which(d0$code == "ms"), ]

expect_true(nrow(d0fd) == 175)
expect_true(nrow(d0md) == 362)
expect_true(nrow(d0fs) == 159)
expect_true(nrow(d0ms) == 333)

# rolling construction of d
# pull son data from fathers
sonbothpar <- as.numeric(d0fs$id %in% d0ms$id)
include <- which(sonbothpar == 1)
ds <- d0fs[include, ]
# pull son data from mothers
for (i in 1:nrow(ds)) {
  tar <- which(d0ms$id == ds$id[i])
  ds$agem[i] <- d0ms$agem[tar]
  ds$rsm[i] <- d0ms$rsm[tar]
  ds$idm[i] <- d0ms$idm[tar]
}

expect_true(nrow(ds) == 141)


# rolling construction of d
# pull daughter data from fathers
daubothpar <- as.numeric(d0fd$id %in% d0md$id)
include <- which(daubothpar == 1)
dd <- d0fd[include, ]
# pull daughter data from mothers
for (i in 1:nrow(dd)) {
  tar <- which(d0md$id == dd$id[i])
  dd$agem[i] <- d0md$agem[tar]
  dd$rsm[i] <- d0md$rsm[tar]
  dd$idm[i] <- d0md$idm[tar]
}
names(dd)[which(names(dd) == "fert")] <- "rs"

expect_true(nrow(dd) == 156)

d0 <- rbind(dd, ds)


#parents first
dp <- data.frame(
  c(d0$idf, d0$idm),
  c(d0$rsf, d0$rsm),
  c(d0$agef, d0$agem),
  c(rep(1, nrow(d0)), rep(0, nrow(d0)))
)
names(dp) <- c("id", "rs", "age", "male")
dp <- dp[-which(duplicated(dp$id)), ]
dp$age <- dp$age/10
D_mu_p <- cbind(1-dp$male, (1-dp$male) * dp$age, dp$male, dp$male * dp$age)
D_sig_p <- cbind(rep(1, nrow(dp)), dp$age)


#now offspring only
#parents first
drop <- which(d0$id %in% d0$idf | d0$id %in% d0$idm | d0$age < 18)
do <- d0[-drop, ]
do$age <- do$age/10
D_mu_o <- cbind(1 - do$male, (1 - do$male) * do$age, do$male, do$male * do$age)
D_sig_o <- cbind(rep(1, nrow(do)), do$age, do$age^2)
male_o <- do$male

#now everyone with parents
drop <- which(d0$id %in% d0$idf | d0$id %in% d0$idm | d0$age < 18)
d <- d0[-drop, ]
d$age <- d$age/10
d$agef <- d$agef/10
d$agem <- d$agem/10
D_mu_f <- cbind(rep(0, nrow(d)), rep(0, nrow(d)), rep(1, nrow(d)), d$agef)
D_mu_m <- cbind(rep(1, nrow(d)), d$agem, rep(0, nrow(d)), rep(0, nrow(d)))
D_sig_f <- cbind(rep(1, nrow(d)), d$agef)
D_sig_m <- cbind(rep(1, nrow(d)), d$agem)
D_mu_o2 <- cbind(1-d$male, (1-d$male) * d$age, d$male, d$male * d$age)
D_sig_o2 <- cbind(rep(1, nrow(d)), d$age, d$age^2)
male_o2 <- d$male

test_stat <- (mean(d$rs[which(d$male == 1)])-mean(d$rs[which(d$male == 0)]))/sd(d$rs)
expect_true(abs(test_stat - (-0.32)) < 0.01)

expect_true(abs(mean(d$age * 10) - 28.32) < 0.01)


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
  D_mu_p = D_mu_p, 
  D_sig_p = D_sig_p, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D_mu_m = D_mu_m, 
  D_sig_m = D_sig_m, 
  D_mu_o2 = D_mu_o2, 
  D_sig_o2 = D_sig_o2, 
  male_o = male_o, 
  male_o2 = male_o2, 
  ones = rep(1, nrow(do))
)

expect_true(data_list$n == 280)
expect_true(data_list$no == 280)
expect_true(data_list$np == 226)



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
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+.88)) #.88 = r
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
}

# reported values here
beta_s_est <- (-0.133)
beta_s_se  <- 0.107
beta_d_est <- 0.063
beta_d_se  <- 0.07

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.196)
beta_diff_se <- 0.129

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
