
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "lamalera_rs"

d0 <- read.csv("./data/lamalera_rs.csv", stringsAsFactors = FALSE)

# The plan: Use as parents all those folk who are counted as such.
# Use as offspring all people counted as such, minus those who are
# also counted as parents. Use separate models for each group.

# parents
dp <- data.frame(
  c(d0$fid, d0$mid),
  c(rep(1, nrow(d0)),
    rep(0, nrow(d0))),
  c(d0$agef, d0$agem),
  c(d0$rsf, d0$rsm)
)
names(dp) <- c("id", "male", "age", "rs")
drop <- which(is.na(dp$id) | is.na(dp$age) | is.na(dp$rs))
dp <- dp[-drop, ]
drop <- which(duplicated(dp$id))
dp <- dp[-drop, ]
dp$age <- dp$age/10
D_mu_p <- cbind(rep(1, nrow(dp)), dp$age)
# sigma is fixed; don't need data for it.

# "nonparents"
drop <- which((d0$id %in% d0$fid) | (d0$id %in% d0$mid))
d1 <- d0[-drop, ]
# names(d1)[3] <- "male" # uh oh
d1$age <- d1$age/10
D_1 <- cbind(rep(1, nrow(d1)), d1$age, d1$age^2)
# same for mu and sigma

# folks with parents
drop <- which(is.na(d0$rsf) | is.na(d0$rsm))
do <- d0[-drop, ]
# names(do)[3] <- "male" # uh oh
do$age <- do$age/10
do$agef <- do$agef/10
do$agem <- do$agem/10
D_mu_f <- cbind(rep(1, nrow(do)), do$agef)
D_mu_m <- cbind(rep(1, nrow(do)), do$agem)
D_o <- cbind(rep(1, nrow(do)), do$age, do$age^2)

(mean(do$rs[which(do$male == 1)])-mean(do$rs[which(do$male == 0)]))/sd(do$rs)
mean(do$age * 10)

# make list for stan
data_list <- list(
  np = nrow(dp), 
  n1 = nrow(d1), 
  no = nrow(do), 
  rs_p = dp$rs, 
  rs_1 = d1$rs, 
  rs_f = do$rsf, 
  rs_m = do$rsm, 
  rs_o = do$rs, 
  D_mu_p = D_mu_p, 
  D_1 = D_1, 
  D_mu_f = D_mu_f, 
  D_mu_m = D_mu_m, 
  D_o = D_o, 
  male = do$male
)

expect_true(data_list$np == 150)
expect_true(data_list$n1 == 410)
expect_true(data_list$no == 119)



# something is wrong here!

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
beta_s_est <- (-0.138)
beta_s_se  <- 0.134
beta_d_est <- 0.422
beta_d_se  <- 0.108

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.56)
beta_diff_se <- 0.175

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
