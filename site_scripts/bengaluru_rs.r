
rm(list = ls())

source("./project_support.r")

measure <- "bengaluru_rs"

d0 <- read.csv("./data/bengaluru_rs.csv", stringsAsFactors = FALSE)

#first: kids only.
age <- 2002-d0$child.year.birth
sonRS <- d0$son.RS
sonRS[which(is.na(sonRS))] <- 0
daughterRS <- d0$daughter.RS
daughterRS[which(is.na(daughterRS))] <- 0
RS <- sonRS + daughterRS
RS[which(RS == 0)] <- NA
male <- rep(NA, nrow(d0))
male[which(!is.na(d0$son.RS))] <- 1
male[which(!is.na(d0$daughter.RS))] <- 0
do <- data.frame(male, age, RS)
drop <- which(is.na(do$male * do$age * do$RS))
do <- do[-drop, ]
do$age <- do$age/10
do$RS <- log(do$RS)  #logs make it slightly more normal, anyway.
D_mu_o <- cbind((1-do$male), do$male, do$age)

#second: parents only.
age <- c(2002-d0$father.year.birth, 2002-d0$mother.year.birth)
RS <- c(d0$parents.RS, d0$parents.RS)
ID <- c(d0$RESPONDENT.ID, d0$RESPONDENT.ID+1000)
male <- c(rep(1, nrow(d0)), rep(0, nrow(d0)))
dp <- data.frame(ID, male, age, RS)
drop <- which(is.na(dp$male * dp$age * dp$RS))
dp <- dp[-drop, ]
drop <- which(duplicated(dp$ID))
dp <- dp[-drop, ]
dp$age <- dp$age/10
dp$RS <- log(dp$RS)  #logs make it slightly more normal, anyway.
D_mu_p <- cbind((1-dp$male), dp$male, dp$age)

#now folks with parents
age <- (2002-d0$child.year.birth)/10
sonRS <- d0$son.RS
sonRS[which(is.na(sonRS))] <- 0
daughterRS <- d0$daughter.RS
daughterRS[which(is.na(daughterRS))] <- 0
RS <- sonRS + daughterRS
RS[which(RS == 0)] <- NA
RS <- log(RS)
male <- rep(NA, nrow(d0))
male[which(!is.na(d0$son.RS))] <- 1
male[which(!is.na(d0$daughter.RS))] <- 0
agef <- (2002-d0$father.year.birth)/10
agem <- (2002-d0$mother.year.birth)/10
RSf <- log(d0$parents.RS)
RSm <- log(d0$parents.RS)
IDf <- d0$RESPONDENT.ID+1000
IDm <- d0$RESPONDENT.ID
d1 <- data.frame(male, age, RS, IDf, agef, RSf, IDm, agem, RSm)
drop <- which(is.na(d1$male * d1$age * d1$RS * d1$agef * d1$agem * d1$RSf * d1$RSm))
d1 <- d1[-drop, ]
D_mu_o1 <- cbind((1-d1$male), d1$male, d1$age)
D_mu_f1 <- cbind(rep(0, nrow(d1)), rep(1, nrow(d1)), d1$agef)
D_mu_m1 <- cbind(rep(1, nrow(d1)), rep(0, nrow(d1)), d1$agem)

#make list for stan
data_list <- list(
  n_o = nrow(do), 
  n_p = nrow(dp), 
  n_1 = nrow(d1), 
  RS_o = do$RS, 
  RS_p = dp$RS, 
  RS_o1 = d1$RS, 
  RS_f1 = d1$RSf, 
  RS_m1 = d1$RSm, 
  D_mu_o = D_mu_o, 
  D_mu_p = D_mu_p, 
  D_mu_o1 = D_mu_o1, 
  D_mu_f1 = D_mu_f1, 
  D_mu_m1 = D_mu_m1
)

expect_true(data_list$n_o == 680)
expect_true(data_list$n_p == 783)
expect_true(data_list$n_1 == 646)



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

  #parents are just about perfectly correlated
  #so I'm only going to use fathers

  #now do the regression
  dtemp <- data.frame(post$Dev_o[i, ], post$Dev_f[i, ])
  names(dtemp) <- c("Dev_o", "par_temp")

  #sons
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * par_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]

  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * par_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]

}

# reported values here
beta_s_est <- 0.214
beta_s_se  <- 0.06
beta_d_est <- 0.115
beta_d_se  <- 0.058

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.099
beta_diff_se <- 0.083

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_fathers = 0,
  n_mothers = 0,
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
