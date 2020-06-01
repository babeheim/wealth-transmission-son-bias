
#########

rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "bangladesh_income"

d0 <- read.csv("./data/bangladesh_income.csv", stringsAsFactors = FALSE)

#all folks, for fitting model

#kids
income <- c(d0$CF_INC_LOG, d0$CF_INC_LOG)
age <- c(d0$RESP_AGE, d0$HUS_AGE)/10
male <- c(rep(0, nrow(d0)), rep(1, nrow(d0)))
d <- data.frame(income, age, male) # 1886
drop <- which(is.na(d$income * d$age))
d <- d[-drop, ] # 1370
#now parents
income <- c(d0$RP_INC_LOG, d0$RP_INC_LOG, d0$HP_INC_LOG, d0$HP_INC_LOG)
age <- c(d0$RF_AGE, d0$RM_AGE, d0$HF_AGE, d0$HM_AGE)/10
age[which(is.na(age))] <- 0
age_dead <- c(d0$RF_AGE_DEATH+d0$RF_YRS_SINCE_DEATH, 
  d0$RM_AGE_DEATH+d0$RM_YRS_SINCE_DEATH, 
  d0$HF_AGE_DEATH+d0$HF_YRS_SINCE_DEATH, 
  d0$HM_AGE_DEATH+d0$HM_YRS_SINCE_DEATH)/10
age_dead[which(is.na(age_dead))] <- 0
age <- age + age_dead
age[which(age == 0)] <- NA
male <- c(rep(1, nrow(d0)), rep(0, nrow(d0)), rep(1, nrow(d0)), rep(0, nrow(d0)))
dp <- data.frame(income, age, male) # 3772
drop <- which(is.na(dp$income * dp$age))
dp <- dp[-drop, ] # 1971
#combine
dt <- rbind(d, dp)
D_mu <- cbind(dt$male, 1-dt$male, dt$age, dt$age^2)
dim(D_mu) # 3314

#now folks with parents.
income <- c(d0$CF_INC_LOG, d0$CF_INC_LOG)
age <- c(d0$RESP_AGE, d0$HUS_AGE)/10
male <- c(rep(0, nrow(d0)), rep(1, nrow(d0)))
incomef <- c(d0$RP_INC_LOG, d0$HP_INC_LOG)
#mother and father income is identical, so mother and father deviations will be very nearly identical.
# So let's just use the fathers.
agef <- c(d0$RF_AGE, d0$HF_AGE)/10
agef[which(is.na(agef))] <- 0
age_dead <- c(d0$RF_AGE_DEATH+d0$RF_YRS_SINCE_DEATH, 
  d0$HF_AGE_DEATH+d0$HF_YRS_SINCE_DEATH)/10
age_dead[which(is.na(age_dead))] <- 0
agef <- agef + age_dead
agef[which(agef == 0)] <- NA
do <- data.frame(income, age, male, incomef, agef)
drop <- which(is.na(do$income * do$age * do$incomef * do$agef))
do <- do[-drop, ]
D_mu_o <- cbind(do$male, 1-do$male, do$age, do$age^2)
D_mu_f <- cbind(rep(1, nrow(do)), rep(0, nrow(do)), do$agef, do$agef^2)

#make list for stan
data_list <- list(
  n = nrow(dt), 
  no = nrow(do), 
  inc = dt$income, 
  inc_o = do$income, 
  inc_f = do$incomef, 
  D_mu = D_mu, 
  D_mu_o = D_mu_o, 
  D_mu_f = D_mu_f
  )

expect_true(data_list$n == 3341)
expect_true(data_list$no == 882)



tic(paste("fit", measure, "stan model"))

init <- list(
  list(
    P_mu = c(12.4, 12.1, -.3, 0), 
    sigma = 1
  )
)

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
  dtemp <- data.frame(post$Dev_o[i, ], post$Dev_f[i, ])
  names(dtemp) <- c("Dev_o", "par_temp")

  #sons
  dtemps <- dtemp[which(do$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * par_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]

  #daughters
  dtempd <- dtemp[which(do$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * par_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]

}

# reported values here
beta_s_est <- 0.289
beta_s_se  <- 0.053
beta_d_est <- 0.222
beta_d_se  <- 0.042

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.067
beta_diff_se <- 0.067

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
