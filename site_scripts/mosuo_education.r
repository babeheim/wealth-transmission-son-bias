

rm(list = ls())

source("./project_support.r")

measure <- "mosuo_education"

d0 <- read.csv("./data/mosuo_education.csv", stringsAsFactors = FALSE)

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
  ds$edm[i] <- d0ms$edm[tar]
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
  dd$edm[i] <- d0md$edm[tar]
  dd$idm[i] <- d0md$idm[tar]
}
names(dd)[which(names(dd) == "fert")] <- "rs"

expect_true(nrow(dd) == 156)

d0 <- rbind(dd, ds)

# fit model on all folks, assuming no difference between parents and nonparents
d <- data.frame(
  c(d0$id, d0$idf, d0$idm),
  c(d0$ed, d0$edf, d0$edm),
  c(d0$age, d0$agef, d0$agem),
  c(d0$male, rep(1, nrow(d0)), rep(0, nrow(d0)))
)
names(d) <- c("id", "ed", "age", "male")
d <- d[-which(duplicated(d$id)), ]

d$age <- d$age/10
D <- cbind(1 - d$male, d$male, (1-d$male) * d$age,
  d$male * d$age, (1-d$male) * d$age^2)


#dp: only folks with parents.
dp <- d0 #no drops needed
dp$age <- dp$age/10
dp$agef <- dp$agef/10
dp$agem <- dp$agem/10
D_o <- cbind(1-dp$male, dp$male, (1-dp$male) * dp$age,
  dp$male * dp$age, (1-dp$male) * dp$age^2)
D_f <- cbind(rep(0, nrow(dp)), rep(1, nrow(dp)),
  rep(0, nrow(dp)), dp$agef, rep(0, nrow(dp)))
D_m <- cbind(rep(1, nrow(dp)), rep(0, nrow(dp)),
  dp$agem, rep(0, nrow(dp)), dp$agem^2)

expect_true(abs(mean(dp$age * 10) - 28.09) < 0.01)

test_stat <- (mean(dp$ed[which(dp$male == 1)]) - mean(dp$ed[which(dp$male == 0)]))/sd(dp$ed)
expect_true(abs(test_stat - 0.11) < 0.01)

# make list for stan
data_list <- list(
  n = nrow(d),
  np = nrow(dp),
  ed = d$ed,
  ed_o = dp$ed,
  ed_f = dp$edf,
  ed_m = dp$edm,
  D = D,
  D_o = D_o,
  D_f = D_f,
  D_m = D_m,
  male = dp$male
)

expect_true(data_list$n == 518)
expect_true(data_list$np == 297)



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
  dparents <- dparents[-which(duplicated(data.frame(dp$idf, dp$idm))), ]
  names(dparents) <- c("Dev_f", "Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents), 1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+r))
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
beta_s_est <- 0.259
beta_s_se  <- 0.085
beta_d_est <- 0.304
beta_d_se  <- 0.073

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.045)
beta_diff_se <- 0.109

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
