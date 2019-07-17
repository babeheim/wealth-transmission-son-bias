
rm(list = ls())

source("./project_support.r")

measure <- "bengaluru_networks"

d0_long <- read.csv("./data/bengaluru_networks.csv", stringsAsFactors = FALSE)

d0_mom <- d0_long[d0_long$parent_is_mother == 1, ]
d0_dad <- d0_long[d0_long$parent_is_mother == 0, ]

d0 <- data.frame(id = sort(unique(d0_long$id)))
d0$sex <- d0_long$sex[match(d0$id, d0_long$id)]
d0$wealth <- d0_long$networkwealth[match(d0$id, d0_long$id)]
d0$yb <- d0_long$yb[match(d0$id, d0_long$id)]

d0$wealthf <- d0_dad$networkwealthp[match(d0$id, d0_dad$id)]
d0$ybf <- d0_dad$ybp[match(d0$id, d0_dad$id)]
d0$idf <- d0_dad$idp[match(d0$id, d0_dad$id)]

d0$wealthm <- d0_mom$networkwealthp[match(d0$id, d0_mom$id)]
d0$ybm <- d0_mom$ybp[match(d0$id, d0_mom$id)]
d0$idm <- d0_mom$idp[match(d0$id, d0_mom$id)]

d0$male <- as.numeric(d0$sex == 1)

#parents
dp <- data.frame(
  c(d0$wealthf, d0$wealthm),
  c(d0$ybf, d0$ybm),
  c(d0$idf, d0$idm),
  c(rep(1, nrow(d0)), rep(0, nrow(d0))))
names(dp) <- c("wealth", "yb", "id", "male")
drop <- which(duplicated(dp$id) | is.na(dp$yb * dp$wealth))
dp <- dp[-drop, ]
dp$wealth <- log(dp$wealth)
dp$yb <- (dp$yb - 1950)/10
D_p <- cbind(rep(1, nrow(dp)), dp$yb)

#kids
drop <- which(duplicated(d0$id) | is.na(d0$yb * d0$wealth))
do <- d0[-drop, ]
do$wealth <- log(do$wealth)
do$yb <- (do$yb-1950)/10

#folks with parents
drop <- which(
  duplicated(d0$id) |
  is.na(d0$yb * d0$wealth * d0$wealthf * d0$wealthm * d0$ybf * d0$ybm)
)
d1 <- d0[-drop, ]
d1$wealth <- log(d1$wealth)
d1$wealthf <- log(d1$wealthf)
d1$wealthm <- log(d1$wealthm)
d1$yb <- (d1$yb-1950)/10
d1$ybf <- (d1$ybf-1950)/10
d1$ybm <- (d1$ybm-1950)/10
D_f1 <- cbind(rep(1, nrow(d1)), d1$ybf)
D_m1 <- cbind(rep(1, nrow(d1)), d1$ybm)

(mean(d1$wealth[which(d1$male == 1)])-mean(d1$wealth[which(d1$male == 0)]))/sd(d1$wealth)
mean(2002-(10 * d1$yb+1950))

#make list for stan
data_list <- list(
  np = nrow(dp), 
  no = nrow(do), 
  n1 = nrow(d1), 
  wealth_p = dp$wealth, 
  wealth_o = do$wealth, 
  wealth_o1 = d1$wealth, 
  wealth_f1 = d1$wealthf, 
  wealth_m1 = d1$wealthm, 
  D_p = D_p, 
  D_f1 = D_f1, 
  D_m1 = D_m1
)

expect_true(data_list$np == 378)
expect_true(data_list$no == 395)
expect_true(data_list$n1 == 345)



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
  dparents <- dparents[-which(duplicated(data.frame(d1$idf, d1$idm))), ]
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
beta_s_est <- 0.122
beta_s_se  <- 0.073
beta_d_est <- 0.12
beta_d_se  <- 0.073

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.002
beta_diff_se <- 0.104

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
