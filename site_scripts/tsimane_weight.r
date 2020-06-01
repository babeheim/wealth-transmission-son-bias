
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "tsimane_weight"

d0 <- read.csv("./data/tsimane_weight.csv", stringsAsFactors = FALSE)

d <- data.frame(
  c(d0$id, d0$fid, d0$mid),
  c(d0$weight, d0$weightf, d0$weightm),
  c(d0$age, d0$agef, d0$agem)/10,
  c(1-d0$female, rep(1, nrow(d0)), rep(0, nrow(d0)))
)
names(d) <- c("id", "weight", "age", "male")
d$id <- as.character(d$id)
drop <- which(
  duplicated(d$id) |
  is.na(d$id) |
  is.na(d$weight * d$age) |
  d$age < 1.8
)

expect_true(sum(is.na(d$id)) == 3506)
expect_true(sum(duplicated(d$id)) == 7570)

d <- d[-drop, ]

expect_true(nrow(d) == 1039)

drop <- which(d$weight < 32 | d$weight > 90)
d <- d[-drop, ]
D_mu <- cbind(
  1 - d$male,
  (1 - d$male) * d$age,
  (1 - d$male) * d$age^2,
  d$male,
  d$male * d$age,
  d$male * d$age^2
)
D_sig <- cbind(1 - d$male, d$male)

drop <- which(
  duplicated(d0$id) |
  is.na(d0$weight * d0$age *
    d0$weightf * d0$agef *
    d0$weightm * d0$agem) |
  d0$age < 18 |
  d0$weight < 32
)
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
d1$male <- 1 - d1$female

D_mu_o <- cbind(
  1 - d1$male,
  (1 - d1$male) * d1$age,
  (1 - d1$male) * d1$age^2,
  d1$male,
  d1$male * d1$age,
  d1$male * d1$age^2
)
D_sig_o <- cbind(1 - d1$male, d1$male)
D_mu_f <- cbind(
  rep(0, nrow(d1)),
  rep(0, nrow(d1)),
  rep(0, nrow(d1)),
  rep(1, nrow(d1)),
  d1$agef,
  d1$agef^2
)
D_sig_f <- cbind(rep(0, nrow(d1)), rep(1, nrow(d1)))
D_mu_m <- cbind(
  rep(1, nrow(d1)),
  d1$agem,
  d1$agem^2,
  rep(0, nrow(d1)),
  rep(0, nrow(d1)),
  rep(0, nrow(d1))
)
D_sig_m <- cbind(rep(1, nrow(d1)), rep(0, nrow(d1)))

(mean(d1$weight[which(d1$male == 1)])-mean(d1$weight[which(d1$male == 0)]))/sd(d1$weight)
mean(d1$age * 10)

data_list <- list(
  n = nrow(d), 
  n1 = nrow(d1), 
  weight = d$weight, 
  weight_o = d1$weight, 
  weight_f = d1$weightf, 
  weight_m = d1$weightm, 
  D_mu = D_mu, 
  D_sig = D_sig, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D_mu_m = D_mu_m, 
  D_sig_m = D_sig_m
)

expect_true(abs(mean(data_list$weight) - 57.69) < 0.01)
expect_true(abs(mean(data_list$weight_f) - 61.29) < 0.01)
expect_true(abs(mean(data_list$weight_m) - 51.75) < 0.01)

expect_true(data_list$n == 1028)
expect_true(data_list$n1 == 217)



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

  if (i %% 100 == 0) print(i)

}

expect_true(abs(mean(beta_s) - 0.354) < beta_tol)
expect_true(abs(sd(beta_s) - 0.087) < beta_tol)
expect_true(abs(mean(beta_d) - 0.224) < beta_tol)
expect_true(abs(sd(beta_d) - 0.101) < beta_tol)

beta_diff_est <- 0.131
beta_diff_se <- 0.131

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
