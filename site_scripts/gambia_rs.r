
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "gambia_rs"

d0 <- read.csv("./data/gambia_rs.csv", stringsAsFactors = FALSE)

d0$male <- as.numeric(d0$sex == "male")

# Executive decision: since I have no way to tell how these people are
# sampled, I'm just going to restrict my attention to people who have
# reproduced at least once.

#all folks who have reproduced
drop <- which(d0$rs == 0 | is.na(d0$age * d0$rs * d0$male) | d0$age < 18)
d <- d0[-drop, ]
d$age <- d$age/10
D <- cbind(d$male, 1 - d$male, d$male * d$age, (1 - d$male) * d$age)

#folks with parents
drop <- which(d0$rs == 0 |
  is.na(d0$age * d0$rs * d0$male * d0$rsf * d0$rsm * d0$agem * d0$agef) |
  d0$age < 18)
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
D_o <- cbind(d1$male, 1 - d1$male, d1$male * d1$age, (1 - d1$male) * d1$age)
D_f <- cbind(rep(1, nrow(d1)), rep(0, nrow(d1)), d1$agef, rep(1, nrow(d1)))
D_m <- cbind(rep(0, nrow(d1)), rep(1, nrow(d1)), d1$agem, rep(0, nrow(d1)))

(mean(d1$rs[which(d1$male == 1)])-mean(d1$rs[which(d1$male == 0)]))/sd(d1$rs)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  n = nrow(d), 
  n1 = nrow(d1), 
  rs = d$rs, 
  rs_o = d1$rs, 
  rs_f = d1$rsf, 
  rs_m = d1$rsm, 
  D = D, 
  D_o = D_o, 
  D_f = D_f, 
  D_m = D_m, 
  male = d$male, 
  male_o = d1$male, 
  ones_1 = rep(1, nrow(d1))
)

expect_true(data_list$n == 946)
expect_true(data_list$n1 == 307)



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
  dparents <- dparents[-which(duplicated(data.frame(d1$fid, d1$mid))), ]
  names(dparents) <- c("Dev_f", "Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents), 1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+r))
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.149
beta_s_se  <- 0.073
beta_d_est <- (-0.034)
beta_d_se  <- 0.045

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.184
beta_diff_se <- 0.089

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
