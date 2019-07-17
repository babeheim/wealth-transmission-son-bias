
rm(list = ls())

source("./project_support.r")

measure <- "bangladesh_education"

d0 <- read.csv("./data/bangladesh_education.csv", stringsAsFactors = FALSE)

#make dataframe with all folks (and no repeats), for fitting age+sex function
#all the latest children
ed <- as.numeric(paste(d0$CLD_EDU[which(d0$CLD_DUMMY == 1)]))
age <- as.numeric(paste(d0$CLD_AGE[which(d0$CLD_DUMMY == 1)]))
sex <- as.numeric(paste(d0$CLD_SEX[which(d0$CLD_DUMMY == 1)]))
#and their parents
#fathers
ed_temp <- as.numeric(paste(d0$FA_EDU[which(d0$CLD_DUMMY == 1)]))
ed_temp <- ed_temp[-which(duplicated(d0$RESP_RID[which(d0$CLD_DUMMY == 1)]))]
age_temp <- as.numeric(paste(d0$FA_AGE[which(d0$CLD_DUMMY == 1)]))
age_temp <- age_temp[-which(duplicated(d0$RESP_RID[which(d0$CLD_DUMMY == 1)]))]
ed <- c(ed, ed_temp)
age <- c(age, age_temp)
sex <- c(sex, rep(1, length(age_temp)))
#mothers
ed_temp <- as.numeric(paste(d0$MO_EDU[which(d0$CLD_DUMMY == 1)]))
ed_temp <- ed_temp[-which(duplicated(d0$RESP_RID[which(d0$CLD_DUMMY == 1)]))]
age_temp <- as.numeric(paste(d0$MO_AGE[which(d0$CLD_DUMMY == 1)]))
age_temp <- age_temp[-which(duplicated(d0$RESP_RID[which(d0$CLD_DUMMY == 1)]))]
ed <- c(ed, ed_temp)
age <- c(age, age_temp)
sex <- c(sex, rep(2, length(age_temp)))
#and now THEIR parents
#gfather
ed_temp <- as.numeric(paste(d0$FA_EDU[which(d0$SIB_DUMMY == 1)]))
ed_temp <- ed_temp[-which(duplicated(d0$RESP_RID[which(d0$SIB_DUMMY == 1)]))]
age_temp <- as.numeric(paste(d0$FA_AGE[which(d0$SIB_DUMMY == 1)]))
age_temp <- age_temp[-which(duplicated(d0$RESP_RID[which(d0$SIB_DUMMY == 1)]))]
ed <- c(ed, ed_temp)
age <- c(age, age_temp)
sex <- c(sex, rep(1, length(age_temp)))
#gmother
ed_temp <- as.numeric(paste(d0$MO_EDU[which(d0$SIB_DUMMY == 1)]))
ed_temp <- ed_temp[-which(duplicated(d0$RESP_RID[which(d0$SIB_DUMMY == 1)]))]
age_temp <- as.numeric(paste(d0$MO_AGE[which(d0$SIB_DUMMY == 1)]))
age_temp <- age_temp[-which(duplicated(d0$RESP_RID[which(d0$SIB_DUMMY == 1)]))]
ed <- c(ed, ed_temp)
age <- c(age, age_temp)
sex <- c(sex, rep(2, length(age_temp)))

d <- data.frame(age, ed, sex)
drop <- which(is.na(age * ed * sex) | age < 18)
d <- d[-drop, ]
d$sex[which(d$sex == 2)] <- 0
names(d) <- c("age", "ed", "male")
d$age <- d$age/10
D_mu <- cbind((1-d$male), (1-d$male) * d$age, d$male, d$male * d$age)
D_sig <- cbind((1-d$male), (1-d$male) * d$age, (1-d$male) * d$age^2, d$male, d$male * d$age, d$male * d$age^2)


#Now the dataset with parent data 
#start with youngest gen
#need to record family ID, for later
ed <- as.numeric(paste(d0$CLD_EDU[which(d0$CLD_DUMMY == 1)]))
age <- as.numeric(paste(d0$CLD_AGE[which(d0$CLD_DUMMY == 1)]))
sex <- as.numeric(paste(d0$CLD_SEX[which(d0$CLD_DUMMY == 1)]))
fam_id <- paste(d0$RESP_RID[which(d0$CLD_DUMMY == 1)])
edf <- as.numeric(paste(d0$FA_EDU[which(d0$CLD_DUMMY == 1)]))
agef <- as.numeric(paste(d0$FA_AGE[which(d0$CLD_DUMMY == 1)]))
edm <- as.numeric(paste(d0$MO_EDU[which(d0$CLD_DUMMY == 1)]))
agem <- as.numeric(paste(d0$MO_AGE[which(d0$CLD_DUMMY == 1)]))
#and the next highest gen
ed <- c(ed, as.numeric(paste(d0$CLD_EDU[which(d0$SIB_DUMMY == 1)])))
age <- c(age, as.numeric(paste(d0$CLD_AGE[which(d0$SIB_DUMMY == 1)])))
sex <- c(sex, as.numeric(paste(d0$CLD_SEX[which(d0$SIB_DUMMY == 1)])))
fam_id <- c(fam_id, paste(d0$RESP_RID[which(d0$SIB_DUMMY == 1)]))
edf <- c(edf, as.numeric(paste(d0$FA_EDU[which(d0$SIB_DUMMY == 1)])))
agef <- c(agef, as.numeric(paste(d0$FA_AGE[which(d0$SIB_DUMMY == 1)])))
edm <- c(edm, as.numeric(paste(d0$MO_EDU[which(d0$SIB_DUMMY == 1)])))
agem <- c(agem, as.numeric(paste(d0$MO_AGE[which(d0$SIB_DUMMY == 1)])))

# offspring table
d1 <- data.frame(age, ed, sex, agef, edf, agem, edm, fam_id)
drop <- which(is.na(age * ed * sex * agef * agem * edf * edm/100) | age < 18)
d1 <- d1[-drop, ]
d1$sex[which(d1$sex == 2)] <- 0
names(d1) <- c("age", "ed", "male", "agef", "edf", "agem", "edm", "fam_id")
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
D_mu_o <- cbind((1-d1$male), (1-d1$male) * d1$age, d1$male, d1$male * d1$age)
D_sig_o <- cbind((1-d1$male), (1-d1$male) * d1$age, (1-d1$male) * d1$age^2, d1$male, d1$male * d1$age, d1$male * d1$age^2)
D_mu_f <- cbind(rep(0, nrow(d1)), rep(0, nrow(d1)), rep(1, nrow(d1)), d1$agef)
D_sig_f <- cbind(rep(0, nrow(d1)), rep(0, nrow(d1)), rep(0, nrow(d1)), rep(1, nrow(d1)), d1$agef, d1$agef^2)
D_mu_m <- cbind(rep(1, nrow(d1)), d1$agem, rep(0, nrow(d1)), rep(0, nrow(d1)))
D_sig_m <- cbind(rep(1, nrow(d1)), d1$agem, d1$agem^2, rep(0, nrow(d1)), rep(0, nrow(d1)), rep(0, nrow(d1)))


#make list for stan
data_list <- list(
  n = nrow(d), 
  ed = d$ed, 
  male = d$male, 
  n1 = nrow(d1), 
  ed_o = d1$ed, 
  ed_f = d1$edf, 
  ed_m = d1$edm, 
  male_o = d1$male, 
  ones1 = rep(1, nrow(d1)),
  D_mu = D_mu, 
  D_sig = D_sig, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D_mu_m = D_mu_m, 
  D_sig_m = D_sig_m
)

expect_true(data_list$n == 4317)
expect_true(sum(data_list$male) == 1935)

expect_true(data_list$n1 == 1112)
expect_true(sum(data_list$male_o == 1) == 596)
expect_true(sum(data_list$male_o == 0) == 516)




tic(paste("fit", measure, "stan model"))

init <- list(
  list(
    P_mu = c(-4.62, 1.2, -7.25, 1.53), 
    P_sig = c(.73, .27, -0.04, .75, .27, -0.02), 
    high = 8, 
    low = .7, 
    high_m = 4.2, 
    low_m = 2.08
  )
)
#these are MLE estimates using bbmle2.


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
  dparents <- dparents[-which(duplicated(d1$fam_id)), ]
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
beta_s_est <- 0.468
beta_s_se  <- 0.033
beta_d_est <- 0.538
beta_d_se  <- 0.033

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.07)
beta_diff_se <- 0.05

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
