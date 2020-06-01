
###########

rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "bangladesh_land"

d0 <- read.csv("./data/bangladesh_land.csv", stringsAsFactors = FALSE)

d0$RF_AGE_DEATH[which(d0$RF_AGE_DEATH == 999)] <- NA
d0$HF_AGE_DEATH[which(d0$HF_AGE_DEATH == 999)] <- NA
d0$RF_YRS_SINCE_DEATH[which(d0$RF_YRS_SINCE_DEATH == 999)] <- NA
d0$HF_YRS_SINCE_DEATH[which(d0$HF_YRS_SINCE_DEATH == 999)] <- NA
d0$HF_AGE[which(d0$HF_AGE == 999)] <- NA

# get all folks who have inherited land (throw in parents) and do not have land = 0
land <- c(d0$CF_TOTAL_DEC, d0$CF_TOTAL_DEC, d0$RP_TOTAL_DEC, d0$HP_TOTAL_DEC)
male <- c(rep(0, nrow(d0)), rep(1, nrow(d0)), rep(1, nrow(d0)), rep(0, nrow(d0)))
age <- c(d0$RESP_AGE_CALC, d0$HUS_AGE, d0$RF_AGE, d0$HF_AGE)
age[which(is.na(age))] <- 0
age_dead <- c(
  rep(0, nrow(d0)),
  d0$HUS_AGE_DEATH + d0$HUS_YRS_SINCE_DEATH,
  d0$RF_AGE_DEATH + d0$RF_YRS_SINCE_DEATH,
  d0$HF_AGE_DEATH + d0$HF_YRS_SINCE_DEATH
)
age_dead[which(is.na(age_dead))] <- 0
age <- age + age_dead
age[which(age == 0)] <- NA
age <- age/10
d <- data.frame(male, land, age)
drop <- which(d0$CF_TOTAL_DEC == d0$RP_TOTAL_DEC | d0$CF_TOTAL_DEC == d0$HP_TOTAL_DEC)
#attempting to drop anyone who "lives at home" with one set of parents
drop <- c(drop, drop + nrow(d0)) #and their husbands...
d <- d[-drop, ]
drop <- which(is.na(d$land * d$age) | d$land == 0)
d <- d[-drop, ]
d$land <- log(d$land)
D_mu <- cbind(rep(1, nrow(d)), d$age, d$age^2)


#now folks with parents
land <- c(d0$CF_TOTAL_DEC, d0$CF_TOTAL_DEC)
landf <- c(d0$RP_TOTAL_DEC, d0$HP_TOTAL_DEC)
male <- c(rep(0, nrow(d0)), rep(1, nrow(d0)))
age <- c(d0$RESP_AGE_CALC, d0$HUS_AGE)
age[which(is.na(age))] <- 0
age_dead <- c(rep(0, nrow(d0)), d0$HUS_AGE_DEATH + d0$HUS_YRS_SINCE_DEATH)
age_dead[which(is.na(age_dead))] <- 0
age <- age + age_dead
age[which(age == 0)] <- NA
age <- age/10
agef <- c(d0$RF_AGE, d0$HF_AGE)
agef[which(is.na(agef))] <- 0
agef_dead <- c(d0$RF_AGE_DEATH + d0$RF_YRS_SINCE_DEATH,
  d0$HF_AGE_DEATH + d0$HF_YRS_SINCE_DEATH)
agef_dead[which(is.na(agef_dead))] <- 0
agef <- agef + agef_dead
agef[which(agef == 0)] <- NA
agef <- agef/10
d1 <- data.frame(male, land, age, landf, agef)
drop <- which(d0$CF_TOTAL_DEC == d0$RP_TOTAL_DEC | d0$CF_TOTAL_DEC == d0$HP_TOTAL_DEC)
#attempting to drop anyone who "lives at home" with one set of parents
drop <- c(drop, drop + nrow(d0)) #and their husbands...
d1 <- d1[-drop, ]
drop <- which(is.na(d1$land * d1$age * d1$landf * d1$agef) | d1$land == 0 | d1$landf == 0)
d1 <- d1[-drop, ]
d1$land <- log(d1$land)
d1$landf <- log(d1$landf)
D_mu_o <- cbind(rep(1, nrow(d1)), d1$age, d1$age^2)
D_mu_f <- cbind(rep(1, nrow(d1)), d1$agef, d1$agef^2)

#make list for stan
data_list <- list(
  n = nrow(d), 
  n1 = nrow(d1), 
  land = d$land, 
  land_o = d1$land, 
  land_f = d1$landf, 
  D_mu = D_mu, 
  D_mu_o = D_mu_o, 
  D_mu_f = D_mu_f
)

expect_true(data_list$n == 2181)
expect_true(data_list$n1 == 719)



tic(paste("fit", measure, "stan model"))

init <- list(
  list(
    P_mu = c(1.8, .3, 0), 
    sigma = 1.5
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
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * par_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * par_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]

}

# reported values here
beta_s_est <- 0.516
beta_s_se  <- 0.04
beta_d_est <- 0.372
beta_d_se  <- 0.041

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.144
beta_diff_se <- 0.056

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)







# png("./temp/figure1.png", res = 300, height = 5, width = 5, units = "in")

# D_o <- colMeans(post$Dev_o)
# D_f <- colMeans(post$Dev_f)

# par(mar=c(4, 3, 2.5, 0.25))
# par(pty="s")
# plot(D_o ~ D_f, col= ifelse(d1$male == 1, "blue", "pink"), main = "Bangladesh Land", xlab = expression(D[mp]), ylab = expression(D[s] * "  and  " * D[d]) , xlim = c(-3, 3), ylim = c(-3, 3))
# for(i in (1:1000) * 50)
# {
#   curve((beta_s[i]) * x, add=TRUE, col=rgb(0, 0, 255, alpha=10, maxColorValue=255))
#   curve((beta_d[i]) * x, add=TRUE, col=rgb(255, 192, 203, alpha=10, maxColorValue=255))
# }
# curve(mean(beta_s) * x, add=TRUE, col="blue", lwd=2)
# curve(mean(beta_d) * x, add=TRUE, col="pink", lwd=2)
# legend(x="bottomright", legend=c(expression(rho[s] * "=" * 0.52, rho[d] * "=" * 0.37)), col=c("blue", "pink"), lty=1, lwd=2)
# dev.off()


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
