
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "ust_avam_education"

d0 <- read.csv("./data/ust_avam_education.csv", stringsAsFactors = FALSE)

for (i in c("ed", "edf", "edm", "age", "agef", "agem")) {
  d0[which(d0[,i] == 9999),i] <- NA
}
d0[which(d0[,"sex"] == ""),"sex"] <- NA
# d0$sex <- as.numeric(d0$sex)-2
# names(d0)[which(names(d0)=="sex")] <- "male"
d0$male <- as.numeric(d0$sex == "M")

#all folks
d <- data.frame(
  c(d0$id,d0$idf,d0$idm),
  c(d0$ed,d0$edf,d0$edm),
  c(d0$age,d0$agef,d0$agem)/10,
  c(d0$male,rep(1,nrow(d0)),rep(0,nrow(d0))))
names(d) <- c("id","ed","age","male")
drop <- which(duplicated(d$id)| is.na(d$ed*d$age*d$male))
d <- d[-drop,]
D1 <- cbind((1-d$male), (1-d$male)*d$age, (1-d$male)*d$age^2,
  d$male, d$male*d$age, d$male*d$age^2)

#now folks with full parent data
drop <- which(duplicated(d0$id)| is.na(d0$ed*d0$age*d0$male*d0$edf*d0$agef*d0$edm*d0$agem))
d1 <- d0[-drop,]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
D_o <- cbind((1-d1$male),(1-d1$male)*d1$age,(1-d1$male)*d1$age^2,d1$male,d1$male*d1$age,d1$male*d1$age^2)
D_f <- cbind(rep(0,nrow(d1)),rep(0,nrow(d1)),rep(0,nrow(d1)),rep(1,nrow(d1)),d1$agef,d1$agef^2)
D_m <- cbind(rep(1,nrow(d1)),d1$agem,d1$agem^2,rep(0,nrow(d1)),rep(0,nrow(d1)),rep(0,nrow(d1)))

(mean(d1$ed[which(d1$male==1)])-mean(d1$ed[which(d1$male==0)]))/sd(d1$ed)

mean(d1$age*10)

#make list for stan
data_list <- list(
  n = nrow(d),
  n1 = nrow(d1),
  ed = d$ed,
  ed_o = d1$ed,
  ed_f = d1$edf,
  ed_m = d1$edm,
  D1 = D1,
  D_o = D_o,
  D_f = D_f,
  D_m = D_m
)

expect_true(data_list$n == 471)
expect_true(data_list$n1 == 100)



tic(paste("fit", measure, "stan model"))

init <- list(
  list(
    P_mu = c(6.67, 2.41, -0.38, 5.96, 2.20, -0.30), 
    P_sigma = c(1.27, -0.70, 0.11, 1.35, -0.55, 0.08)
  )
)

if (!file.exists(file.path("temp", paste0(measure, ".stanfit")))){

  m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
    iter = n_iter, warmup = n_warmup, chains = 1, init = init)

  save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

}
toc(log = TRUE)



load(file.path("temp", paste0(measure, ".stanfit")))

post <- extract(m1, permute = TRUE)

n_sim <- min(10000, length(post$lp__))

beta_s <- rep(NA, n_sim)
beta_d <- rep(NA, n_sim)

for (i in 1:n_sim) {

  #find parent correlation
  dparents <- data.frame(post$Dev_f[i,],post$Dev_m[i,])
  dparents <- dparents[-which(duplicated(data.frame(d1$idf,d1$idm))),]
  names(dparents) <- c("Dev_f","Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r*Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents),1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5*(post$Dev_f[i,]+post$Dev_m[i,])/sqrt(1/2*(1+r))
  dtemp <- data.frame(post$Dev_o[i,], midpar_temp)
  names(dtemp) <- c("Dev_o","midpar_temp")
  
  #sons
  dtemps <- dtemp[which(d1$male==1),]
  m_temp <- mle2(Dev_o ~ dnorm(beta*midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp,1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male==0),]
  m_temp <- mle2(Dev_o ~ dnorm(beta*midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp,1)[[1]]

}

# reported values here
beta_s_est <- 0.53
beta_s_se  <- 0.117
beta_d_est <- (-0.151)
beta_d_se  <- 0.181

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < 0.05) # ?
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.681
beta_diff_se <- 0.216

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
